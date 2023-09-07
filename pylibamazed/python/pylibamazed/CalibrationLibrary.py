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
import logging
import os

import h5py
import numpy as np
import pandas as pd
from astropy.io import ascii, fits
from pylibamazed.Exception import (AmazedError, AmazedErrorFromGlobalException,
                                   APIException)
from pylibamazed.Parameters import Parameters
from pylibamazed.redshift import (CalzettiCorrection, CFlagWarning,
                                  CLineCatalog, CLineCatalogsTplRatio,
                                  CLineRatioCatalog, CPhotBandCatalog,
                                  CPhotometricBand,
                                  CSpectrumFluxAxis_withSpectrum,
                                  CSpectrumFluxCorrectionCalzetti,
                                  CSpectrumFluxCorrectionMeiksin,
                                  CSpectrumSpectralAxis, CTemplate,
                                  CTemplateCatalog, ErrorCode, GlobalException,
                                  MeiksinCorrection, TAsymParams,
                                  VecMeiksinCorrection, VecTFloat64List,
                                  undefStr)

zflag = CFlagWarning.GetInstance()


def _get_linecatalog_id(row):
    wl = round(row.WaveLength, 2)
    return row.Name + "_" + str(wl) + "_" + row.Type


def load_reliability_model(model_path, parameters: Parameters, object_type):
    try:
        # to avoid annoying messages about gpu/cuda availability
        os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
        import keras
        from tensorflow.keras import models
    except ImportError:
        raise APIException(ErrorCode.RELIABILITY_NEEDS_TENSORFLOW,
                           "Tensorflow is required to compute the reliability")
    ret = dict()
    model_ha = h5py.File(model_path).attrs
    keras_model_version = model_ha["keras_version"].split(".")
    keras_system_version = keras.__version__.split(".")
    if keras_model_version[0] > keras_system_version[0]:
        raise APIException(ErrorCode.RELIABILITY_NEEDS_TENSORFLOW,
                           f"Tensorflow major version >= {keras_model_version[0]} required")
    redshift_range = [model_ha["zrange_min"],
                      model_ha["zrange_max"]]
    redshift_range_step = model_ha["zrange_step"]
    s_redshift_range = parameters.get_redshiftrange(object_type)
    s_redshift_range_step = parameters.get_redshiftstep(object_type)
    if s_redshift_range != redshift_range:
        raise APIException(
            ErrorCode.INCOMPATIBLE_PDF_MODELSHAPES,
            "redshift range of reliability model must be identical to solver one : "
            f"{redshift_range} != {s_redshift_range}"
        )
    if s_redshift_range_step != redshift_range_step:
        raise APIException(ErrorCode.INCOMPATIBLE_PDF_MODELSHAPES,
                           "redshift step of reliability model must be identical to solver one : "
                           f"{redshift_range_step} != {s_redshift_range_step}")
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

    def __init__(self, parameters: Parameters, calibration_dir):
        self.parameters = parameters
        self.calibration_dir = os.path.expanduser(calibration_dir)
        self.templates_catalogs = dict()
        self.templates_ratios = dict()
        self.line_catalogs = dict()
        self.line_catalogs_df = dict()
        for object_type in self.parameters.get_objects():
            self.line_catalogs[object_type] = dict()
            self.line_catalogs_df[object_type] = dict()

        self.line_ratio_catalog_lists = dict()
        self.lambda_offsets = dict()
        self.lsf = dict()
        try:
            self.photometric_bands = CPhotBandCatalog()
        except GlobalException as e:
            raise AmazedErrorFromGlobalException(e)

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
            raise APIException(ErrorCode.INVALID_DIRECTORY,
                               "No template category directory found in {}".format(full_path))
        for category in categories:
            category_path = os.path.join(full_path, category)
            # Category directory contains template files
            template_file_list = os.listdir(category_path)

            if not template_file_list:
                raise APIException(ErrorCode.INVALID_FILEPATH,
                                   "No template file found in {}".format(category_path))
            for template_filename in template_file_list:
                file_path = os.path.join(category_path, template_filename)
                # Template file is a two columns text file
                try:
                    data = np.loadtxt(file_path, unpack=True)
                except Exception:
                    raise APIException(ErrorCode.INVALID_FILEPATH,
                                       "Unable to read template file {}".format(file_path))
                wavelength = data[0]
                flux = data[1]

                spectralaxis = CSpectrumSpectralAxis(wavelength)
                signal = CSpectrumFluxAxis_withSpectrum(flux)
                template = CTemplate(template_filename, category,
                                     spectralaxis, signal)
                self.templates_catalogs[object_type].Add(template)

    def load_templates_catalog(self, object_type):
        self.templates_catalogs[object_type] = CTemplateCatalog()
        template_dir = self.parameters.get_template_dir(object_type)
        self._load_templates(object_type, template_dir)
        # Temporary hack before handling process flow in api
        if "all" not in self.templates_catalogs:
            self.templates_catalogs["all"] = CTemplateCatalog()
        self._load_templates("all", template_dir)

    def load_linecatalog(self, object_type, solve_method):
        logger = logging.getLogger("calibration_api")
        line_catalog_file = os.path.join(
            self.calibration_dir,
            self.parameters.get_linecatalog(object_type, solve_method)
        )
        if not os.path.exists(line_catalog_file):
            raise APIException(ErrorCode.INVALID_FILEPATH, "{} cannot be found".format(line_catalog_file))

        logger.info("Loading {} linecatalog: {}".format(object_type, line_catalog_file))

        nsigmasupport = self.parameters.get_nsigmasupport(object_type, solve_method)
        self.line_catalogs[object_type][solve_method] = CLineCatalog(nsigmasupport)
        try:
            line_catalog = pd.read_csv(
                line_catalog_file,
                sep='\t',
                dtype={"WaveLength": float, "AmplitudeGroupName": str, "EnableFitWaveLengthOffset": bool}
            )
        except pd.errors.ParserError as e:
            raise AmazedError(ErrorCode.BAD_FILEFORMAT,
                              "bad line catalog {0} cause :{1}".format(line_catalog_file, e))
        except Exception as e:
            raise Exception("bad line catalog " + line_catalog_file + " cause :" + "{}".format(e))

        # force "-1" to undefStr (for compatibility)
        line_catalog.loc[line_catalog.AmplitudeGroupName == "-1", "AmplitudeGroupName"] = undefStr

        enableIGM = self.parameters.get_solve_method_igm_fit(object_type, solve_method)

        # here should go the change of profiles if igm is applied
        self.line_catalogs_df[object_type][solve_method] = line_catalog
        for index, row in line_catalog.iterrows():
            if row.Profile == "ASYM":
                asymParams = TAsymParams(1., 4.5, 0.)
            elif row.Profile == "ASYMFIT":
                asymParams = TAsymParams(2., 2., 0.)
            elif row.Profile == "ASYMFIXED":
                raise Exception("Profile in linecatalog cannot be asymFixed")
            else:
                asymParams = TAsymParams(0, 0, 0)
            self.line_catalogs[object_type][solve_method].AddLineFromParams(row.Name,
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
                                                                            _get_linecatalog_id(row),
                                                                            self.meiksin)
            if enableIGM:
                self.line_catalogs[object_type][solve_method].convertLineProfiles2SYMIGM(self.meiksin)

    def load_line_ratio_catalog_list(self, object_type):
        logger = logging.getLogger("calibration_api")
        tplratio_catalog = self.parameters.get_lineModelSolve_tplratio_catalog(object_type)
        line_ratio_catalog_list = os.path.join(self.calibration_dir,
                                               tplratio_catalog,
                                               "*.tsv")
        line_ratio_catalog_list = glob.glob(line_ratio_catalog_list)
        if not line_ratio_catalog_list:
            raise APIException(ErrorCode.INVALID_FILEPATH, "Template ratio catalog empty")
        line_ratio_catalog_list.sort()
        logger.info("Loading {} line ratio catalogs: {}".format(object_type, tplratio_catalog))

        self.line_ratio_catalog_lists[object_type] = CLineCatalogsTplRatio()
        n_ebmv_coeffs = 1
        if self.parameters.get_lineModelSolve_tplratio_ismfit(object_type):
            n_ebmv_coeffs = self.parameters.get_ebmv_count()
        prior = 1. / (n_ebmv_coeffs * len(line_ratio_catalog_list))

        enableIGM = self.parameters.get_lineModelSolve_igmfit(object_type)

        for f in line_ratio_catalog_list:
            lr_catalog_df = pd.read_csv(f, sep='\t')
            name = f.split(os.sep)[-1][:-4]
            with open(os.path.join(self.calibration_dir,
                                   self.parameters.get_lineModelSolve_tplratio_catalog(object_type),
                                   name + ".json")) as f:
                line_ratio_catalog_parameter = json.load(f)
            for k in range(n_ebmv_coeffs):
                lr_catalog = CLineRatioCatalog(name, self.line_catalogs[object_type]["LineModelSolve"])
                for index, row in lr_catalog_df.iterrows():
                    if row.Name in list(self.line_catalogs_df[object_type]["LineModelSolve"].Name):
                        lr_catalog.setLineAmplitude(_get_linecatalog_id(row), row.NominalAmplitude)
                lr_catalog.addVelocity("em_vel", line_ratio_catalog_parameter["velocities"]["em_vel"])
                lr_catalog.addVelocity("abs_vel", line_ratio_catalog_parameter["velocities"]["abs_vel"])
                # here also we should change the profile type
                lr_catalog.setAsymProfileAndParams(
                    line_ratio_catalog_parameter["asym_params"]["profile"],
                    TAsymParams(line_ratio_catalog_parameter["asym_params"]["sigma"],
                                line_ratio_catalog_parameter["asym_params"]["alpha"],
                                line_ratio_catalog_parameter["asym_params"]["delta"],
                                )
                )
                if enableIGM:
                    lr_catalog.convertLineProfiles2SYMIGM(self.meiksin)

                lr_catalog.setIsmIndex(k)
                lr_catalog.setPrior(prior)
                self.line_ratio_catalog_lists[object_type].addLineRatioCatalog(lr_catalog)

    def load_empty_line_catalog(self, object_type, method):
        self.line_catalogs[object_type][method] = CLineCatalog()

    def load_empty_line_ratio_catalog_list(self, object_type):
        self.line_ratio_catalog_lists[object_type] = CLineCatalogsTplRatio()

    def load_lsf(self):
        if self.parameters.get_lsf_type() == "GaussianVariableWidth":
            file = os.path.join(self.calibration_dir,
                                self.parameters.get_lsf_width_file_name())
            # TODO check hdul here
            with fits.open(file) as hdul:
                self.lsf["wave"] = hdul[1].data.field(0)
                self.lsf["width"] = hdul[1].data.field(1)

    def load_photometric_bands(self):
        paths = os.path.join(self.calibration_dir,
                             self.parameters.get_photometry_transmission_dir(),
                             "*")

        bands = self.parameters.get_photometry_bands()
        for f in glob.glob(paths):
            df = pd.read_csv(f, comment='#')
            band = df.columns[1]
            if band in bands:
                self.photometric_bands.Add(band, CPhotometricBand(np.array(df[band]),
                                                                  np.array(df["lambda"])))

    def load_calzetti(self):
        df = ascii.read(os.path.join(self.calibration_dir, "ism", "SB_calzetti.dl1.txt"))
        _calzetti = CalzettiCorrection(df['lambda'], df['flux'])
        self.calzetti = CSpectrumFluxCorrectionCalzetti(_calzetti,
                                                        self.parameters.get_ebmv_start(),
                                                        self.parameters.get_ebmv_step(),
                                                        self.parameters.get_ebmv_count())

    # Important: igm curves should be loaded in the increasing order of their extinction per bin of z,
    # i.e., from the least extinction curve to the highest extinction curve
    def load_Meiksin(self):
        zbins = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]
        columns = ['restlambda', 'flux0', 'flux1', 'flux2', 'flux3', 'flux4', 'flux5', 'flux6']
        meiksinCorrectionCurves = VecMeiksinCorrection()
        for z in zbins[1:]:
            filename = f"Meiksin_Var_curves_{z}.txt"
            file_path = os.path.join(self.calibration_dir, "igm",
                                     "IGM_variation_curves_meiksin_v3.1", filename)
            if not os.path.isfile(file_path):  # mainly for unit tests
                continue
            meiksin_df = ascii.read(os.path.join(self.calibration_dir, "igm",
                                    "IGM_variation_curves_meiksin_v3.1", filename))
            columns = columns[:len(meiksin_df.columns)]
            meiksin_df.rename_columns(tuple(meiksin_df.columns), tuple(columns))
            fluxcorr = VecTFloat64List([meiksin_df[col] for col in columns[1:]])
            meiksinCorrectionCurves.append(MeiksinCorrection(meiksin_df['restlambda'], fluxcorr))
        self.meiksin = CSpectrumFluxCorrectionMeiksin(meiksinCorrectionCurves, zbins)

    def load_all(self, calibs="all"):
        """Load templates, line catalogs and template ratios for every object_type, according to parameters
        content
        """
        meiksin = calibs == "all" or "meiksin" in calibs
        calzetti = calibs == "all" or "calzetti" in calibs
        templates = calibs == "all" or "templates" in calibs
        linecatalogs = calibs == "all" or "linecatalogs" in calibs
        lineratios = calibs == "all" or "lineratios" in calibs
        reliability = calibs == "all" or "reliability" in calibs
        try:
            if meiksin:
                self.load_Meiksin()
            if calzetti:
                self.load_calzetti()
            for object_type in self.parameters.get_objects():
                if templates:
                    self.load_templates_catalog(object_type)
                # load linecatalog for linemodelsolve

                solve_method = self.parameters.get_solve_method(object_type)
                if solve_method == "LineModelSolve":
                    if linecatalogs:
                        self.load_linecatalog(object_type, solve_method)

                    if self.parameters.is_tplratio_catalog_needed(object_type):
                        if lineratios:
                            self.load_line_ratio_catalog_list(object_type)
                            # load linecatalog for linemeassolve
                linemeas_method = self.parameters.get_linemeas_method(object_type)
                if linemeas_method == "LineMeasSolve":
                    if linecatalogs:
                        self.load_linecatalog(object_type, linemeas_method)

                # Load the reliability model
                if self.parameters.get_reliability_enabled(object_type) and reliability:
                    model_path = os.path.join(self.calibration_dir,
                                              self.parameters.get_reliability_model(object_type))
                    mp = load_reliability_model(model_path,
                                                self.parameters,
                                                object_type)
                    self.reliability_models[object_type] = mp["model"]
                    self.reliability_parameters[object_type] = mp["parameters"]

            if self.parameters.get_lsf_type() != "FROMSPECTRUMDATA":
                self.load_lsf()

            if self.parameters.get_photometry_transmission_dir() is not None:
                self.load_photometric_bands()
        except GlobalException as e:
            raise AmazedErrorFromGlobalException(e)
        except FileNotFoundError as e:
            raise AmazedError(ErrorCode.INVALID_FILEPATH, str(e))
        except APIException as e:
            raise AmazedError(e.errCode, e.message)
        except Exception as e:
            raise AmazedError(ErrorCode.PYTHON_API_ERROR, str(e))

    def init(self):
        """Initialize templates (init continuum removal, init ism/igm and lsf if lsf is not spectrum dependent
        """
        raise NotImplementedError("TODO after reviewing CProcessFlowContext::Init")

    def get_sub_type(self, object_type, line_ratio_catalog):
        try:
            with open(os.path.join(
                self.calibration_dir,
                self.parameters.get_lineModelSolve_tplratio_catalog(object_type),
                line_ratio_catalog + ".json"
            )) as f:
                tpl_ratio_conf = json.load(f)
                return tpl_ratio_conf["sub_type"]
        except FileNotFoundError:
            return ""

    def get_lines_ids(self, attributes):
        lines_ids = dict()
        for attr in attributes:
            attr_parts = attr.split(".")
            try:
                if attr_parts[1] == "linemeas":
                    linemeas_object = attr_parts[0]
                    l_method = "LineMeasSolve"
                elif attr_parts[1] == "fitted_lines":
                    linemeas_object = attr_parts[0]
                    l_method = "LineModelSolve"
                else:
                    continue
                lines = self.line_catalogs_df[linemeas_object][l_method]
            except Exception:
                continue
            try:
                line_name = attr_parts[2]
                line_id = lines[lines["Name"] == line_name].index[0]
                lines_ids[line_name] = line_id
            except Exception:
                raise Exception(f"Could not find {line_name} in catalog")
        return lines_ids
