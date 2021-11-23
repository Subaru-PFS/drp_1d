import os
import logging
from pylibamazed.redshift import (CSpectrumSpectralAxis,
                                  CSpectrumFluxAxis_withSpectrum,
                                  CTemplate, CTemplateCatalog,
                                  CRayCatalog)
import numpy as np
from astropy.io import fits


class CalibrationLibrary:
    def __init__(self, parameters, calibration_dir):
        self.parameters = parameters
        self.calibration_dir = os.path.expanduser(calibration_dir)
        self.templates_catalogs = dict()
        self.templates_ratios = dict()
        self.line_catalogs = dict()
        self.lambda_offsets = dict()
        self.lsf = dict()

    def _load_templates(self, object_type, path):
        """Load a template catalog from a directory

        Parameters
        ----------
        full_path : str
            Path to template directory

        Notes
        -----
        The template catalog must contain category direcrtories (at least
        one category). Template files are located in category directory
        and format as two columns text file.

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
        #Temporary hack before handling process flow in api
        if "all" not in self.templates_catalogs:
            self.templates_catalogs["all"] = CTemplateCatalog()
        self._load_templates("all", self.parameters[object_type]["template_dir"])

    def load_linecatalog(self, object_type):
        logger = logging.getLogger("calibration_api")
        linemodel_params = self.parameters[object_type]["linemodelsolve"]["linemodel"]
        if "linecatalog" not in linemodel_params:
            raise Exception("Incomplete parameter file, linemodelsolve.linemodel.linecatalog entry mandatory")
        line_catalog_file = linemodel_params["linecatalog"]
        logger.info("Loading {} linecatalog: {}".format(object_type, line_catalog_file))

        self.line_catalogs[object_type] = CRayCatalog()
        self.line_catalogs[object_type].Load(os.path.join(self.calibration_dir, line_catalog_file))

    def load_empty_line_catalog(self,object_type):
        self.line_catalogs[object_type] = CRayCatalog()

    def load_lsf(self):
        if self.parameters["LSF"]["LSFType"] == "GaussianVariableWidth":
            with fits.open(self.parameters["LSF"]["GaussianVariablewidthFileName"]) as hdul:
                self.lsf["wave"] = hdul[1].data.field(0)
                self.lsf["width"] = hdul[1].data.field(1)

    def load(self, object_type):
        self.load_templates_catalog(object_type)
        if self.parameters[object_type]["method"] == "linemodelsolve":
            self.load_linecatalog(object_type)

    def init(self):
        '''
        Initialize templates (init continuum removal, init ism/igm and lsf if lsf is not spectrum dependent
        :return:
        '''
        pass