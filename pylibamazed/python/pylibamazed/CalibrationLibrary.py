# ============================================================================
#
# This file is part of: AMAZED
#
# Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
#
# https://www.lam.fr/
#
# This software is a computer program whose purpose is to estimate the
# spectrocopic redshift of astronomical sources (galaxy/quasar/star)
# from there 1D spectrum.
#
# This software is governed by the CeCILL-C license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL-C
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL-C license and that you accept its terms.
# ============================================================================

import glob
import json
import os

import h5py
import numpy as np
import pandas as pd
from astropy.io import ascii, fits
from pylibamazed.Exception import APIException, exception_decorator
from pylibamazed.OutputSpecifications import ResultsSpecifications
from pylibamazed.Parameters import Parameters
from pylibamazed.redshift import (
    CalzettiCorrection,
    CFlagWarning,
    CLineCatalog,
    CLineCatalogsTplRatio,
    CLineRatioCatalog,
    CLog,
    CPhotBandCatalog,
    CPhotometricBand,
    CSpectrumFluxAxis_withSpectrum,
    CSpectrumFluxCorrectionCalzetti,
    CSpectrumFluxCorrectionMeiksin,
    CSpectrumSpectralAxis,
    CTemplate,
    CTemplateCatalog,
    ErrorCode,
    MeiksinCorrection,
    TAsymParams,
    VecMeiksinCorrection,
    VecTFloat64List,
    undefStr,
)

zflag = CFlagWarning.GetInstance()
zlog = CLog.GetInstance()


def _strid(waveLength, name, ltype):
    wl = round(waveLength, 2)
    return name + "_" + str(wl) + "_" + ltype


def _get_linecatalog_strid(lineCatalog_df):
    return [
        _strid(w, n, t)
        for w, n, t in zip(lineCatalog_df.WaveLength, lineCatalog_df.Name, lineCatalog_df.Type)
    ]


def load_reliability_model(model_path, parameters: Parameters, object_type):
    zlog.LogInfo(f"Loading reliability neural network for {object_type}")
    try:
        # to avoid annoying messages about gpu/cuda availability
        os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
        import keras
        from tensorflow.keras import models
    except ImportError:
        raise APIException(
            ErrorCode.RELIABILITY_NEEDS_TENSORFLOW, "Tensorflow is required to compute the reliability"
        ) from None
    ret = dict()
    model_ha = h5py.File(model_path).attrs
    keras_model_version = model_ha["keras_version"].split(".")
    keras_system_version = keras.__version__.split(".")
    if keras_model_version[0] > keras_system_version[0]:
        raise APIException(
            ErrorCode.RELIABILITY_NEEDS_TENSORFLOW,
            f"Tensorflow major version >= {keras_model_version[0]} required",
        )
    redshift_range = [model_ha["zrange_min"], model_ha["zrange_max"]]
    redshift_range_step = model_ha["zrange_step"]
    s_redshift_range = parameters.get_redshiftrange(object_type)
    s_redshift_range_step = parameters.get_redshiftstep(object_type)
    if s_redshift_range != redshift_range:
        raise APIException(
            ErrorCode.INCOMPATIBLE_PDF_MODELSHAPES,
            "redshift range of reliability model must be identical to solver one : "
            f"{redshift_range} != {s_redshift_range}",
        )
    if s_redshift_range_step != redshift_range_step:
        raise APIException(
            ErrorCode.INCOMPATIBLE_PDF_MODELSHAPES,
            "redshift step of reliability model must be identical to solver one : "
            f"{redshift_range_step} != {s_redshift_range_step}",
        )
    model = models.load_model(model_path)
    ret["model"] = model
    ret["parameters"] = dict()
    ret["parameters"]["zgrid_end"] = model_ha["zgrid_end"]
    ret["parameters"]["zrange_step"] = model_ha["zrange_step"]
    if "classes" in model_ha:
        ret["parameters"]["classes"] = json.loads(model_ha["classes"])
    else:
        ret["parameters"]["classes"] = ["failure", "success"]
    # TODO add classes here
    return ret


class CalibrationLibrary:
    """
    Class containing all calibration data necessary for the run, according to parameters.
    It includes templates, line catalogs, lambda_offsets (for line catalogs), and template ratios (i.e.
    linecatalogs with amplitude ratio for each line)

    :param parameters: parameters
    :type parameters: dict
    :param calibration_dir: path of the directory containing all calibration data
    :type calibration_dir: path
    """

    @exception_decorator
    def __init__(self, parameters: Parameters, calibration_dir):
        self.parameters = parameters
        self.calibration_dir = os.path.expanduser(calibration_dir)
        self.templates_catalogs = dict()
        self.templates_ratios = dict()
        self.line_catalogs = dict()
        self.line_catalogs_df = dict()
        for object_type in self.parameters.get_spectrum_models():
            self.line_catalogs[object_type] = dict()
            self.line_catalogs_df[object_type] = dict()

        self.line_ratio_catalog_lists = dict()
        self.lr_catalog_param = dict()
        self.lambda_offsets = dict()
        self.lsf = dict()

        self.photometric_bands = CPhotBandCatalog()

        self.calzetti = None
        self.meiksin = None
        self.reliability_models = {}
        self.reliability_parameters = dict()

    def _load_templates(self, object_type, path):
        """
        Load a template catalog from a directory
        Notes
        -----
        The template catalog must contain category directories (at least
        one category). Template files are located in category directory
        and format as two columns text file.

        :param object_type: object_type (in values in [galaxy,star,qso,linemeas])
        :type object_type: str
        :param path: path of the directory containing templates
        :type path: str
        """

        full_path = os.path.join(self.calibration_dir, path)
        # Template directory contains category directories
        categories = os.listdir(full_path)
        if not categories:
            raise APIException(
                ErrorCode.INVALID_DIRECTORY, "No template category directory found in {}".format(full_path)
            )
        for category in categories:
            category_path = os.path.join(full_path, category)
            # Category directory contains template files
            template_file_list = os.listdir(category_path)

            if not template_file_list:
                raise APIException(
                    ErrorCode.INVALID_FILEPATH, "No template file found in {}".format(category_path)
                )
            for template_filename in template_file_list:
                file_path = os.path.join(category_path, template_filename)
                # Template file is a two columns text file
                try:
                    data = np.loadtxt(file_path, unpack=True)
                except Exception:
                    raise APIException(
                        ErrorCode.INVALID_FILEPATH, "Unable to read template file {}".format(file_path)
                    ) from None
                wavelength = data[0]
                flux = data[1]

                spectralaxis = CSpectrumSpectralAxis(wavelength)
                signal = CSpectrumFluxAxis_withSpectrum(flux)
                template = CTemplate(template_filename, category, spectralaxis, signal)
                self.templates_catalogs[object_type].Add(template)

    def load_templates_catalog(self, object_type):
        self.templates_catalogs[object_type] = CTemplateCatalog()
        template_dir = self.parameters.get_template_dir(object_type)

        zlog.LogInfo(f"Loading {object_type} templates: {template_dir}")

        self._load_templates(object_type, template_dir)
        # Temporary hack before handling process flow in api
        if "all" not in self.templates_catalogs:
            self.templates_catalogs["all"] = CTemplateCatalog()
        self._load_templates("all", template_dir)

    def load_linecatalog(self, object_type, solve_method):
        line_catalog_file = os.path.join(
            self.calibration_dir, self.parameters.get_linecatalog(object_type, solve_method)
        )
        if not os.path.exists(line_catalog_file):
            raise APIException(ErrorCode.INVALID_FILEPATH, "{} cannot be found".format(line_catalog_file))

        zlog.LogInfo(f"Loading {object_type} linecatalog: {line_catalog_file}")

        nsigmasupport = self.parameters.get_nsigmasupport(object_type, solve_method)
        self.line_catalogs[object_type][solve_method] = CLineCatalog(nsigmasupport)
        try:
            line_catalog = pd.read_csv(
                line_catalog_file,
                sep="\t",
                dtype={
                    "WaveLength": float,
                    "AmplitudeGroupName": str,
                    "EnableFitWaveLengthOffset": bool,
                    "Name": str,
                    "Type": str,
                    "NominalAmplitude": float,
                },
            )
        except pd.errors.ParserError as e:
            raise APIException(
                ErrorCode.BAD_FILEFORMAT, f"bad line catalog {line_catalog_file} cause :{e}"
            ) from None
        except Exception as e:
            raise APIException(
                ErrorCode.PYTHON_API_ERROR, f"bad line catalog {line_catalog_file} cause :{e}"
            ) from None

        # force "-1" to undefStr (for compatibility)
        line_catalog.loc[line_catalog.AmplitudeGroupName == "-1", "AmplitudeGroupName"] = undefStr
        line_catalog["strId"] = _get_linecatalog_strid(line_catalog)
        if not line_catalog.strId.is_unique:
            raise APIException(
                ErrorCode.DUPLICATED_LINES,
                "some rows in the linecatalog are duplicating the\
                same line (name position, type)",
            )
        self.line_catalogs_df[object_type][solve_method] = line_catalog

        for index, row in line_catalog.iterrows():
            if row.Profile == "ASYM":
                asymParams = TAsymParams(1.0, 4.5, 0.0)
            elif row.Profile == "ASYMFIT":
                asymParams = TAsymParams(2.0, 2.0, 0.0)
            elif row.Profile == "ASYMFIXED":
                raise APIException(ErrorCode.PYTHON_API_ERROR, "Profile in linecatalog cannot be asymFixed")
            else:
                asymParams = TAsymParams(0, 0, 0)
            self.line_catalogs[object_type][solve_method].AddLineFromParams(
                row.Name,
                row.WaveLength,
                row.Type,
                row.Force,
                row.Profile,
                asymParams,
                row.AmplitudeGroupName,
                row.AmplitudeGroupValue,
                row.DispersionVelocityGroupName,
                row.WaveLengthOffset,
                row.EnableFitWaveLengthOffset,
                index,
                row.strId,
                self.meiksin,
            )

        enableIGM = self.parameters.get_solve_method_igm_fit(object_type, solve_method)
        if enableIGM:
            self.line_catalogs[object_type][solve_method].convertLineProfiles2SYMIGM(self.meiksin)

    def load_line_ratio_catalog_list(self, object_type):
        tplratio_catalog = self.parameters.get_linemodel_tplratio_catalog(object_type)
        line_ratio_catalog_list_path = os.path.join(self.calibration_dir, tplratio_catalog, "*.tsv")
        line_ratio_catalog_list = glob.glob(line_ratio_catalog_list_path)
        if not line_ratio_catalog_list:
            raise APIException(
                ErrorCode.INVALID_FILEPATH, f"Empty template ratio catalog:  {line_ratio_catalog_list_path}"
            )
        line_ratio_catalog_list.sort()

        zlog.LogInfo(f"Loading {object_type} line ratio catalogs: {tplratio_catalog}")

        self.line_ratio_catalog_lists[object_type] = CLineCatalogsTplRatio()
        n_ebmv_coeffs = 1
        if self.parameters.get_linemodel_tplratio_ismfit(object_type):
            n_ebmv_coeffs = self.parameters.get_ebmv_count()
        prior = 1.0 / (n_ebmv_coeffs * len(line_ratio_catalog_list))

        line_catalog_df = self.line_catalogs_df[object_type]["lineModelSolve"]
        line_ids = line_catalog_df["strId"].reset_index()
        line_ids["index"] = line_ids["index"].astype(pd.Int64Dtype())  # enable int missing value
        self.lr_catalog_param[object_type] = dict()
        for f in line_ratio_catalog_list:
            lr_catalog_df = pd.read_csv(
                f, sep="\t", dtype={"WaveLength": float, "Name": str, "Type": str, "NominalAmplitude": float}
            )
            name = f.split(os.sep)[-1][:-4]
            with open(
                os.path.join(
                    self.calibration_dir,
                    self.parameters.get_linemodel_tplratio_catalog(object_type),
                    name + ".json",
                )
            ) as f:
                line_ratio_catalog_parameter = json.load(f)

            self.lr_catalog_param[object_type][name] = line_ratio_catalog_parameter
            lr_catalog_df["strId"] = _get_linecatalog_strid(lr_catalog_df)
            if not lr_catalog_df.strId.is_unique:
                raise APIException(
                    ErrorCode.DUPLICATED_LINES,
                    "some rows in the lineratio catalog are\
                     duplicating the same line (name position, type)",
                )

            # set corresponding line ids from the main catalog
            lr_catalog_df = lr_catalog_df.merge(
                line_ids, on="strId", how="left", validate="1:1", indicator=True
            )
            missing_ids = lr_catalog_df["index"].isnull()
            if np.any(missing_ids):
                wrong_lines = " ; ".join([lineid for lineid in lr_catalog_df.strId[missing_ids]])
                raise APIException(
                    ErrorCode.LINE_RATIO_UNKNOWN_LINE,
                    f"unknown line ratio line in\
                     line ratio catalog : {wrong_lines}",
                )

            lr_catalog_df = (lr_catalog_df.loc[~missing_ids]).set_index("index")
            lr_catalog = CLineRatioCatalog(name, self.line_catalogs[object_type]["lineModelSolve"])
            for index, row in lr_catalog_df.iterrows():
                lr_catalog.setLineAmplitude(int(index), row.NominalAmplitude)
            lr_catalog.addVelocity("em_vel", line_ratio_catalog_parameter["velocities"]["em_vel"])
            lr_catalog.addVelocity("abs_vel", line_ratio_catalog_parameter["velocities"]["abs_vel"])
            lr_catalog.setAsymProfileAndParams(
                line_ratio_catalog_parameter["asym_params"]["profile"],
                TAsymParams(
                    line_ratio_catalog_parameter["asym_params"]["sigma"],
                    line_ratio_catalog_parameter["asym_params"]["alpha"],
                    line_ratio_catalog_parameter["asym_params"]["delta"],
                ),
            )
            lr_catalog.setPrior(prior)

            for k in range(n_ebmv_coeffs):
                lr_catalog.setIsmIndex(k)
                self.line_ratio_catalog_lists[object_type].addLineRatioCatalog(lr_catalog)

    def load_empty_line_catalog(self, object_type, method):
        self.line_catalogs[object_type][method] = CLineCatalog()

    def load_empty_line_ratio_catalog_list(self, object_type):
        self.line_ratio_catalog_lists[object_type] = CLineCatalogsTplRatio()

    def load_lsf(self):
        if self.parameters.get_lsf_type() == "gaussianVariableWidth":
            file = os.path.join(self.calibration_dir, self.parameters.get_lsf_width_file_name())
            if not os.path.isfile(file):
                raise APIException(ErrorCode.INVALID_FILEPATH, f"wrong LSF file {file}")

            zlog.LogInfo(f"Loading wavelength variable LSF: {file}")
            with fits.open(file) as hdulist:
                try:
                    lsf_hdu = hdulist[1]
                except IndexError:
                    raise APIException(
                        ErrorCode.BAD_FILEFORMAT, f"Cannot access hdu 1 for LSF in {file}"
                    ) from None
                self.lsf["wave"] = lsf_hdu.data.field(0)
                self.lsf["width"] = lsf_hdu.data.field(1)

    def load_photometric_bands(self):
        path = os.path.join(self.calibration_dir, self.parameters.get_photometry_transmission_dir())
        if not os.path.isdir(path):
            raise APIException(ErrorCode.INVALID_FILEPATH, f"Photometric transmission dir not found: {path}")
        bands = self.parameters.get_photometry_bands()
        band_paths = glob.glob(os.path.join(path, "*"))
        if not band_paths:
            raise APIException(ErrorCode.INVALID_FILEPATH, "Photometric transmission dir empty")
        for f in band_paths:
            zlog.LogInfo(f"Loading photometric transmission: {f}")
            df = pd.read_csv(f, comment="#")
            band = df.columns[1]
            if band in bands:
                if band in self.photometric_bands.GetNameList():
                    raise APIException(
                        ErrorCode.BAD_FILEFORMAT, f"Photometric transmission band duplicated: {band}"
                    )
                self.photometric_bands.Add(
                    band, CPhotometricBand(df[band].to_numpy(), df["lambda"].to_numpy())
                )
        # check all requested bands are loaded
        if len(bands) != self.photometric_bands.size():
            unknown_bands = [band for band in bands if band not in self.photometric_bands]
            raise APIException(
                ErrorCode.BAD_PARAMETER_VALUE,
                f"some bands of parameter photometryBand are unknown: {unknown_bands}",
            )

    def load_calzetti(self):
        path = os.path.join(self.calibration_dir, "ism", "SB_calzetti.dl1.txt")
        if not os.path.isfile(path):
            raise APIException(ErrorCode.INVALID_FILEPATH, f"ISM extinction file not found: {path}")
        zlog.LogInfo(f"Loading Calzetti ism extinction: {path}")
        df = ascii.read(path)
        _calzetti = CalzettiCorrection(df["lambda"], df["flux"])
        self.calzetti = CSpectrumFluxCorrectionCalzetti(
            _calzetti,
            self.parameters.get_ebmv_start(),
            self.parameters.get_ebmv_step(),
            self.parameters.get_ebmv_count(),
        )

    # Important: igm curves should be loaded in the increasing order of their extinction per bin of z,
    # i.e., from the least extinction curve to the highest extinction curve
    def load_Meiksin(self):
        zbins = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]
        columns = ["restlambda", "flux0", "flux1", "flux2", "flux3", "flux4", "flux5", "flux6"]
        meiksinCorrectionCurves = VecMeiksinCorrection()
        meiksin_root_path = os.path.join(self.calibration_dir, "igm", "IGM_variation_curves_meiksin_v3.1")
        filename_pattern = os.path.join(meiksin_root_path, "Meiksin_Var_curves_[2-7].[0,5].txt")
        zlog.LogInfo(f"Loading Meiksin igm extinction: {filename_pattern}")
        filenames = glob.glob(filename_pattern)
        if not filenames:
            raise APIException(
                ErrorCode.INVALID_FILEPATH, f"IGM extinction files not found in: {meiksin_root_path}"
            )
        filenames.sort()
        zbins = [
            float(os.path.basename(name).replace("Meiksin_Var_curves_", "").replace(".txt", ""))
            for name in filenames
        ]
        zbins.insert(0, 1.5)
        for path in filenames:
            meiksin_df = ascii.read(path)
            columns = columns[: len(meiksin_df.columns)]
            meiksin_df.rename_columns(tuple(meiksin_df.columns), tuple(columns))
            fluxcorr = VecTFloat64List([meiksin_df[col] for col in columns[1:]])
            meiksinCorrectionCurves.append(MeiksinCorrection(meiksin_df["restlambda"], fluxcorr))
        self.meiksin = CSpectrumFluxCorrectionMeiksin(meiksinCorrectionCurves, zbins)

    @exception_decorator
    def load_all(self, calibs="all"):
        """Load templates, line catalogs and template ratios for every object_type, according to parameters
        content
        """
        meiksin = calibs == "all" or "meiksin" in calibs
        calzetti = calibs == "all" or "calzetti" in calibs
        templates = calibs == "all" or "templates" in calibs
        linecatalogs = calibs == "all" or "lineCatalogs" in calibs
        lineratios = calibs == "all" or "lineratios" in calibs
        reliability = calibs == "all" or "reliability" in calibs
        try:
            if meiksin:
                self.load_Meiksin()
            if calzetti:
                self.load_calzetti()
            for object_type in self.parameters.get_spectrum_models():
                if templates:
                    self.load_templates_catalog(object_type)
                # load linecatalog for linemodelsolve

                solve_method = self.parameters.get_redshift_solver_method(object_type)
                if solve_method == "lineModelSolve":
                    if linecatalogs:
                        self.load_linecatalog(object_type, solve_method)

                    if self.parameters.is_tplratio_catalog_needed(object_type):
                        if lineratios:
                            self.load_line_ratio_catalog_list(object_type)
                            # load linecatalog for linemeassolve
                linemeas_method = self.parameters.get_linemeas_method(object_type)
                if linemeas_method == "lineMeasSolve":
                    if linecatalogs:
                        self.load_linecatalog(object_type, linemeas_method)
                # Load the reliability model
                if self.parameters.get_reliability_enabled(object_type) and reliability:
                    model_path = os.path.join(
                        self.calibration_dir, self.parameters.get_reliability_model(object_type)
                    )
                    mp = load_reliability_model(model_path, self.parameters, object_type)
                    self.reliability_models[object_type] = mp["model"]
                    self.reliability_parameters[object_type] = mp["parameters"]

            if self.parameters.get_lsf_type() != "fromSpectrumData":
                self.load_lsf()

            if self.parameters.get_photometry_transmission_dir() is not None:
                self.load_photometric_bands()
        except FileNotFoundError as e:
            raise APIException(ErrorCode.INVALID_FILEPATH, str(e)) from None

    def init(self):
        """Initialize templates (init continuum removal, init ism/igm and lsf if lsf is not spectrum dependent"""
        raise NotImplementedError("TODO after reviewing CProcessFlowContext::Init")

    def get_sub_type(self, object_type, line_ratio_catalog):
        try:
            tpl_ratio_conf = self.lr_catalog_param[object_type][line_ratio_catalog]
            return tpl_ratio_conf["sub_type"]
        except KeyError:
            raise APIException(
                ErrorCode.PYTHON_API_ERROR, f"Could not find {line_ratio_catalog} in tpl ratio catalog"
            ) from None

    @exception_decorator
    def get_lines_ids(self, attributes):
        lines_ids = dict()
        lines = None
        specs = ResultsSpecifications()
        for attr in attributes:
            attr_parts = attr.split(".")
            if len(attr_parts) == 1 or attr_parts[0] == "error" or "WarningFlags" in attr_parts[-1]:
                continue

            # Identifying attribute by dataset
            attr_name = attr_parts[-1]
            if attr_name.isnumeric():
                attr_name = attr_parts[-2]

            attr_entry = specs.get_df_by_name(attr_name)
            try:
                dataset = attr_entry["dataset"].values[0]
            except IndexError:
                continue

            if dataset == "linemeas":
                linemeas_object = attr_parts[0]
                l_method = "lineMeasSolve"
            elif dataset == "fitted_lines":
                linemeas_object = attr_parts[0]
                l_method = "lineModelSolve"
            else:
                continue
            if (
                linemeas_object in self.line_catalogs_df
                and l_method in self.line_catalogs_df[linemeas_object]
            ):
                lines = self.line_catalogs_df[linemeas_object][l_method]
            if lines is not None:
                try:
                    line_name = attr_parts[1]
                    line_id = lines[lines["Name"] == line_name].index[0]
                    lines_ids[line_name] = line_id
                except Exception:
                    raise APIException(
                        ErrorCode.PYTHON_API_ERROR, f"Could not find {line_name} in catalog"
                    ) from None
        return lines_ids
