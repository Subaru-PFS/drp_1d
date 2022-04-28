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

import os
import logging
import json
import pandas as pd
from pylibamazed.redshift import (CSpectrumSpectralAxis,
                                  CSpectrumFluxAxis_withSpectrum,
                                  CTemplate, CTemplateCatalog,
                                  CLineCatalog,CLineCatalogsTplShape,
                                  CLineRatioCatalog,
                                  CPhotBandCatalog, CPhotometricBand,
                                  TAsymParams,
                                  CFlagWarning,MeiksinCorrection,
                                  CSpectrumFluxCorrectionMeiksin,
                                  VecTFloat64List, VecMeiksinCorrection,
                                  CSpectrumFluxCorrectionCalzetti, CalzettiCorrection,
                                  GlobalException, ErrorCode)
import numpy as np
from astropy.io import fits, ascii
import glob
from pylibamazed.Exception import AmazedError,AmazedErrorFromGlobalException
zflag = CFlagWarning.GetInstance()


def _get_linecatalog_id(row):
    wl = round(row.WaveLength,2)
    return row.Name + "_" + str(wl) + "_" + row.Type


class CalibrationLibrary:
    """
    Class containing all calibration data necessary for the run, according to parameters.
    It includes templates, line catalogs, lambda_offsets (for line catalogs), and template ratios (i.e. linecatalogs
    with amplitude ratio for each line

    :param parameters: parameters
    :type parameters: dict
    :param calibration_dir: path of the directory containing all calibration data
    :type calibration_dir: path
    """
    def __init__(self, parameters, calibration_dir):
        self.parameters = parameters
        self.calibration_dir = os.path.expanduser(calibration_dir)
        self.templates_catalogs = dict()
        self.templates_ratios = dict()
        self.line_catalogs = dict()
        self.line_catalogs_df = dict()
        self.line_ratio_catalog_lists = dict()
        self.lambda_offsets = dict()
        self.lsf = dict()
        self.photometric_bands = CPhotBandCatalog()
        self.calzetti = None
        self.meiksin = None
        self.reliability_models = {}

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

        logger = logging.getLogger("calibration_api")
        full_path = os.path.join(self.calibration_dir, path)
        # Template directory contains category directories
        categories = os.listdir(full_path)
        if not categories:
            raise AmazedError(ErrorCode.INVALID_DIRECTORY, "No template category directory found in {}".format(full_path))
        for category in categories:
            category_path = os.path.join(full_path, category)
            # Category directory contains template files
            template_file_list = os.listdir(category_path)

            if not template_file_list:
                raise AmazedError(ErrorCode.INVALID_FILEPATH,"No template file found in {}".format(category_path))
            for template_filename in template_file_list:
                file_path = os.path.join(category_path, template_filename)
                # Template file is a two columns text file
                try:
                    data = np.loadtxt(file_path, unpack=True)
                except Exception as e:
                    #logger.error("Unable to read template file {}".format(file_path))#todo: check if we keep the logger
                    raise AmazedError(ErrorCode.INVALID_FILEPATH,"Unable to read template file {}".format(file_path))
                wavelength = data[0]
                flux = data[1]

                spectralaxis = CSpectrumSpectralAxis(wavelength)
                signal = CSpectrumFluxAxis_withSpectrum(flux)
                template = CTemplate(template_filename, category,
                                     spectralaxis, signal)
                self.templates_catalogs[object_type].Add(template)

    def load_templates_catalog(self, object_type):
        if "template_dir" not in self.parameters[object_type]:
            # TODO create a dedicated class for setup exceptions
            raise AmazedError(ErrorCode.MISSING_PARAMETER, "Incomplete parameter file, template_dir entry mandatory")
        self.templates_catalogs[object_type] = CTemplateCatalog()
        self._load_templates(object_type, self.parameters[object_type]["template_dir"])
        # Temporary hack before handling process flow in api
        if "all" not in self.templates_catalogs:
            self.templates_catalogs["all"] = CTemplateCatalog()
        self._load_templates("all", self.parameters[object_type]["template_dir"])

    def load_linecatalog(self, object_type, method):
        logger = logging.getLogger("calibration_api")
        linemodel_params = self.parameters[object_type][method]["linemodel"]
        if "linecatalog" not in linemodel_params:
            raise AmazedError(ErrorCode.MISSING_PARAMETER, "Incomplete parameter file, {}.linemodel.linecatalog entry"
                                                          " mandatory".format(method))
        line_catalog_file = os.path.join(self.calibration_dir,linemodel_params["linecatalog"])
        if not os.path.exists(line_catalog_file):
            raise AmazedError(ErrorCode.INVALID_FILEPATH,"{} cannot be found".format(line_catalog_file))

        logger.info("Loading {} linecatalog: {}".format(object_type, line_catalog_file))

        nsigmasupport = self.parameters[object_type][method]["linemodel"]["nsigmasupport"]
        self.line_catalogs[object_type][method] = CLineCatalog(nsigmasupport)
        try:
            line_catalog = pd.read_csv( line_catalog_file, sep='\t',dtype={"WaveLength":float,
                                                                           "EnableFitWaveLengthOffset":bool})
        except pd.errors.ParserError:
            raise AmazedError(ErrorCode.BAD_FILEFORMAT, "bad line catalog {0} cause :{1}".format(line_catalog_file, e))
        except Exception as e:
            raise AmazedError(ErrorCode.BAD_FILEFORMAT,"bad line catalog {0} cause :{1}".format(line_catalog_file, e))

        self.line_catalogs_df[object_type] = line_catalog
        for index, row in line_catalog.iterrows():
            if row.Profile == "ASYM":
                asymParams = TAsymParams(1., 4.5, 0.)
            elif row.Profile == "ASYMFIT":
                asymParams = TAsymParams(2., 2., 0.)
            else:
                asymParams = TAsymParams(0, 0, 0)
            self.line_catalogs[object_type][method].AddLineFromParams(row.Name,
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
                                                             self._get_linecatalog_id(row))

    def load_line_ratio_catalog_list(self, object_type, method):
        logger = logging.getLogger("calibration_api")
        linemodel_params = self.parameters[object_type][method]["linemodel"]
        if "tplratio_catalog" not in linemodel_params:
            raise AmazedError(ErrorCode.MISSING_PARAMETER,"Missing mandatory entry: {}.linemodel.tplratio_catalog ".format(method))

        line_ratio_catalog_list = os.path.join(self.calibration_dir,
                                               linemodel_params["tplratio_catalog"],
                                               "*.tsv")
        line_ratio_catalog_list = glob.glob(line_ratio_catalog_list)
        line_ratio_catalog_list.sort()
        logger.info("Loading {} line ratio catalogs: {}".format(object_type, linemodel_params["tplratio_catalog"]))

        self.line_ratio_catalog_lists[object_type] = CLineCatalogsTplShape()
        n_ebmv_coeffs = 1
        if self.parameters[object_type][method]["linemodel"]["tplratio_ismfit"]:
            n_ebmv_coeffs = self.parameters["ebmv"]["count"]
        prior = 1./(n_ebmv_coeffs * len(line_ratio_catalog_list))

        for f in line_ratio_catalog_list:
            lr_catalog_df = pd.read_csv(f,sep='\t')
            name = f.split(os.sep)[-1][:-4]
            with open(os.path.join(self.calibration_dir,
                                   linemodel_params["tplratio_catalog"],
                                   name + ".json")) as f:
                line_ratio_catalog_parameter = json.load(f)
            for k in range(n_ebmv_coeffs):
                lr_catalog = CLineRatioCatalog(name, self.line_catalogs[object_type][method])
                for index,row in lr_catalog_df.iterrows():
                    if row.Name in list(self.line_catalogs_df[object_type].Name):
                        lr_catalog.setLineAmplitude(_get_linecatalog_id(row), row.NominalAmplitude)
                lr_catalog.addVelocity("em_vel", line_ratio_catalog_parameter["velocities"]["em_vel"])
                lr_catalog.addVelocity("abs_vel", line_ratio_catalog_parameter["velocities"]["abs_vel"])
                lr_catalog.setAsymProfileAndParams(line_ratio_catalog_parameter["asym_params"]["profile"],
                                                   TAsymParams(line_ratio_catalog_parameter["asym_params"]["sigma"],
                                                               line_ratio_catalog_parameter["asym_params"]["alpha"],
                                                               line_ratio_catalog_parameter["asym_params"]["delta"])
                                                   )
                lr_catalog.setIsmIndex(k)
                lr_catalog.setPrior(prior)
                self.line_ratio_catalog_lists[object_type].addLineRatioCatalog(lr_catalog)

    def load_empty_line_catalog(self, object_type, method):
        self.line_catalogs[object_type][method] = CLineCatalog()

    def load_empty_line_ratio_catalog_list(self, object_type):
        self.line_ratio_catalog_lists[object_type] = CLineCatalogsTplShape()

    def load_lsf(self):
        if self.parameters["LSF"]["LSFType"] == "GaussianVariableWidth":
            file = os.path.join(self.calibration_dir,
                                self.parameters["LSF"]["GaussianVariablewidthFileName"])
            # TODO check hdul here
            with fits.open(file) as hdul:
                self.lsf["wave"] = hdul[1].data.field(0)
                self.lsf["width"] = hdul[1].data.field(1)

    def load_photometric_bands(self):
        paths = os.path.join(self.calibration_dir,
                             self.parameters["photometryTransmissionDir"],
                             "*")
        if not "photometryBand" in self.parameters:
            raise AmazedError(ErrorCode.MISSING_PARAMETER,"photometryBand parameter required")
        bands = self.parameters["photometryBand"]
        if len(bands) == 0:
            raise AmazedError(ErrorCode.INVALID_PARAMETER, "photometryBand parameter is empty")
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
                                                        self.parameters["ebmv"]["start"],
                                                        self.parameters["ebmv"]["step"],
                                                        self.parameters["ebmv"]["count"])

    # Important: igm curves should be loaded in the increasing order of their extinction per bin of z,
    # i.e., from the least extinction curve to the highest extinction curve 
    def load_Meiksin(self):
        zbins = [2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0]
        columns =['restlambda', 'flux0', 'flux1', 'flux2','flux3', 'flux4', 'flux5', 'flux6']
        #we'd better create a map between meiksin curves and their zbin
        meiksinCorrectionCurves = VecMeiksinCorrection()
        for z in zbins: 
            filename = f"Meiksin_Var_curves_{z}.txt"
            meiksin_df = ascii.read(os.path.join(self.calibration_dir, "igm", "IGM_variation_curves_meiksin", filename), names=columns)
            fluxcorr = VecTFloat64List([meiksin_df[col] for col in columns[1:]])
            meiksinCorrectionCurves.append(MeiksinCorrection(meiksin_df['restlambda'], fluxcorr))
        self.meiksin = CSpectrumFluxCorrectionMeiksin(meiksinCorrectionCurves)

    def load_all(self):
        """Load templates, line catalogs and template ratios for every object_type, according to parameters content

        """
        self.load_Meiksin()
        self.load_calzetti()
        for object_type in self.parameters["objects"]:
            self.load_templates_catalog(object_type)
            #load linecatalog for linemodelsolve
            self.line_catalogs[object_type] = dict()
            method = self.parameters[object_type]["method"]            
            if method == "LineModelSolve":                
                self.load_linecatalog(object_type,method)
                if self.parameters[object_type][method]["linemodel"]["rigidity"] == "tplshape":
                    self.load_line_ratio_catalog_list(object_type, method)

            # Load the reliability model
            if self.parameters[object_type].get("enable_reliability"):
                try:
                    # to avoid annoying messages about gpu/cuda availability
                    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
                    from tensorflow.keras import models
                except ImportError:
                    zflag.warning(zflag.RELIABILITY_NEEDS_TENSORFLOW,"Tensorflow is required to compute the reliability")
                else:
                    model_path = os.path.join(self.calibration_dir,
                                              self.parameters[object_type]["reliability_model"])
                    model = models.load_model(model_path)
                    self.reliability_models[object_type] = model

        if self.parameters["LSF"]["LSFType"] != "FROMSPECTRUMDATA":
            self.load_lsf()

            if "photometryTransmissionDir" in self.parameters:
                self.load_photometric_bands()
        except AmazedError as e:
            raise e
        except GlobalException as e:
            raise AmazedErrorFromGlobalException(e)
        except FileNotFoundError as e:
            raise AmazedError(ErrorCode.INVALID_FILEPATH, str(e))
        except Exception as e:
            raise AmazedError(ErrorCode.PYTHON_API_ERROR, str(e))

    def init(self):
        """Initialize templates (init continuum removal, init ism/igm and lsf if lsf is not spectrum dependent

        """
        raise NotImplementedError("TODO after reviewing CProcessFlowContext::Init")
