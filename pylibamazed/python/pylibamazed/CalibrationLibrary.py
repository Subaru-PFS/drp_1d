import os
import logging

import pandas as pd
from pylibamazed.redshift import (CSpectrumSpectralAxis,
                                  CSpectrumFluxAxis_withSpectrum,
                                  CTemplate, CTemplateCatalog,
                                  CRayCatalog,
                                  CPhotBandCatalog, CPhotometricBand)
import numpy as np
from astropy.io import fits
import glob
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


class CalibrationLibrary:
    """
    Class containing all calibration data necessary for the run, according to parameters.
    It includes templates, line catalogs, lambda_offsets (for line catalogs), and template ratios (i.e. linecatalogs
    with amplitude ratio for each ray

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
        self.lambda_offsets = dict()
        self.lsf = dict()
        self.photometric_bands = CPhotBandCatalog()

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
            raise ValueError("No template category directory found in {}".format(full_path))

        for category in categories:
            category_path = os.path.join(full_path, category)
            # Category directory contains template files
            template_file_list = os.listdir(category_path)

            if not template_file_list:
                raise ValueError("No template file found in {}".format(category_path))

            for template_filename in template_file_list:
                file_path = os.path.join(category_path, template_filename)
                # Template file is a two columns text file
                try:
                    data = np.loadtxt(file_path, unpack=True)
                except Exception as e:
                    logger.error("Unable to read template file {}".format(file_path))
                    raise e
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
            raise Exception("Incomplete parameter file, template_dir entry mandatory")
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
            raise Exception("Incomplete parameter file, "+method+".linemodel.linecatalog entry mandatory")
        line_catalog_file = linemodel_params["linecatalog"]
        logger.info("Loading {} linecatalog: {}".format(object_type, line_catalog_file))

        self.line_catalogs[object_type] = CRayCatalog()
        self.line_catalogs[object_type].Load(os.path.join(self.calibration_dir, line_catalog_file))

    def load_empty_line_catalog(self, object_type):
        self.line_catalogs[object_type] = CRayCatalog()

    def load_lsf(self):
        if self.parameters["LSF"]["LSFType"] == "GaussianVariableWidth":
            with fits.open(self.parameters["LSF"]["GaussianVariablewidthFileName"]) as hdul:
                self.lsf["wave"] = hdul[1].data.field(0)
                self.lsf["width"] = hdul[1].data.field(1)

    def load_photometric_bands(self):
        paths = os.path.join(self.calibration_dir,
                             self.parameters["photometryTransmissionDir"],
                             "*")
        bands = self.parameters["photometryBand"]
        for f in glob.glob(paths):
            df = pd.read_csv(f, comment='#')
            band = df.columns[1]
            if band in bands:
                self.photometric_bands.Add(band, CPhotometricBand(np.array(df[band]),
                                                                  np.array(df["lambda"])))

    def load_all(self):
        """Load templates, line catalogs and template ratios for every object_type, according to parameters content

        """
        for object_type in ["star", "galaxy", "qso", "linemeas"]:
            linecatalog_loaded = False
            if self.parameters.get("enable" + object_type + "solve") is True:
                self.load_templates_catalog(object_type)
                method = self.parameters[object_type]["method"]
                if method == "linemodelsolve" or method == "linemeassolve":
                    self.load_linecatalog(object_type,method)
                    linecatalog_loaded = True
            if not linecatalog_loaded:
                self.load_empty_line_catalog(object_type)
        if "photometryTransmissionDir" in self.parameters:
            self.load_photometric_bands()

    def init(self):
        """Initialize templates (init continuum removal, init ism/igm and lsf if lsf is not spectrum dependent

        """
        raise NotImplementedError("TODO after reviewing CProcessFlowContext::Init")
