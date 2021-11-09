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
import os
import pandas as pd
from astropy.io import fits
import json
import numpy as np

class InputManager:
    def __init__(self, amazed_config, amazed_parameters, title =''):

        self.config = amazed_config
        self.title = title

        if not all(key in self.config for key in ["spectrum_dir", "calibration_dir"]):
            raise Exception("Given configuration is incomplete")

        if not os.access(os.path.expanduser(self.config["spectrum_dir"]), os.R_OK):
            raise Exception("Spectrum directory " + self.config["spectrum_dir"] + " does not exist")

        if not os.access(os.path.expanduser(self.config["calibration_dir"]), os.R_OK):
            raise Exception("Calibration directory " + self.config["calibration_dir"] + " does not exist")

        if amazed_parameters["enablegalaxysolve"] and "linemodelsolve" in amazed_parameters["galaxy"] \
                and "tplratio_catalog" in amazed_parameters["galaxy"]["linemodelsolve"]["linemodel"]:
            self.tpl_ratio_catalog = amazed_parameters["galaxy"]["linemodelsolve"]["linemodel"]["tplratio_catalog"]
        else:
            self.tpl_ratio_catalog = None

        self.spectra = pd.read_csv(self.config["input_file"], sep=r'\s+', header=None, dtype={2: object})
        self.spectra.set_index(self.spectra.columns[2], inplace=True)
        self.spectra.rename(columns={0: 'flux', 1: 'error'}, inplace=True)
        self.catalog_data = None
        self.amazed_parameters = amazed_parameters

    def get_spectra_dir(self):
        return os.path.expanduser(self.config["spectrum_dir"])

    def get_template_dir(self,object_type):
        return os.path.expanduser(os.path.join(self.config["calibration_dir"],
            self.amazed_parameters[object_type]["template_dir"]))

    def get_linecatalog_path(self):
        return os.path.expanduser(os.path.join(self.config["calibration_dir"],
            self.amazed_parameters["linecatalog"]))

    def get_spectrum_path(self, spectrum):
        spectrum = str(spectrum)
        return os.path.join(self.get_spectra_dir(), self.spectra.at[spectrum, 'flux'])

    def get_noise_path(self, spectrum):
        spectrum = str(spectrum)
        return os.path.join(self.get_spectra_dir(), self.spectra.at[spectrum, 'error'])

    def get_template_ratio_path(self, tplratio_file):
        return os.path.join(self.config["calibration_dir"], self.tpl_ratio_catalog, tplratio_file)

    def get_template_path(self, object_type,tpl_file):
        return os.path.join(self.get_template_dir(object_type), object_type, tpl_file)

    def get_full_catalog_path(self):
        return self.get_linecatalog_path()

    def get_mismatches_catalog_path(self):
        module_root_dir = os.path.split(__file__)[0]
        return os.path.join(module_root_dir, "resources", "mismatches.csv")

    # ************************* Data retrieval functions ************************* #

    def get_spectrum_data(self, spectrum=None):
        hdul = fits.open(self.get_spectrum_path( spectrum))
        res = dict()
        if len(hdul) > 1:
            res['lambda'] = hdul[1].data['wave'].tolist()
            res['flux'] = hdul[1].data['flux'].tolist()
        else:
            crpix = hdul[0].header["CRPIX1"]
            crval1 = hdul[0].header["CRVAL1"]
            cdelt1 = hdul[0].header["CDELT1"]
            lambda_min = crval1 - cdelt1*(crpix-1)
            if type(hdul[0].data[0]) == np.ndarray:
                nb_step = len(hdul[0].data[0])
                res['flux'] = hdul[0].data[0].tolist()
            else:
                nb_step = len(hdul[0].data)
                res['flux'] = hdul[0].data.tolist()
            res['lambda'] = [lambda_min+ i*cdelt1 for i in range(nb_step)]                
        return res

    def get_noise_data(self, spectrum=None):
        hdul = fits.open(self.get_noise_path(spectrum))
        res = dict()
        if len(hdul) > 1:
            res['lambda'] = hdul[1].data['wave'].tolist()
            res['variance'] = hdul[1].data[hdul[1].columns[1].name].tolist()
        else:
            crpix = hdul[0].header["CRPIX1"]
            crval1 = hdul[0].header["CRVAL1"]
            cdelt1 = hdul[0].header["CDELT1"]
            lambda_min = crval1 - cdelt1 * (crpix - 1)
            if type(hdul[0].data[0]) == np.ndarray:
                nb_step = len(hdul[0].data[0])
                res['variance'] = hdul[0].data[0].tolist()
            else:
                nb_step = len(hdul[0].data)
                res['variance'] = hdul[0].data.tolist()
            res['lambda'] = [lambda_min + i * cdelt1 for i in range(nb_step)]

        return res

    def get_catalog_data_df(self, catalog_path, filters = []):
        catalog_data = pd.read_csv(catalog_path, sep='\t', header=0, skiprows=1,
                                   float_precision='round_trip', index_col=None)
        for criterion in filters:
            name = criterion["field"]
            value = criterion["value"]
            c_type = criterion["type"]
            if c_type == "=" and value != "All":
                catalog_data = catalog_data[catalog_data[name] == value]

        ret = dict()
        if  "linecatalog_convert" in self.amazed_parameters and self.amazed_parameters["linecatalog_convert"]:
            s = 0.8
            catalog_data["#lambda"] = catalog_data["#lambda"] / (
                            1.0 + 8.34254 * 1e-5 + (2.406147 * 1e-2) / (130 - s*s)
                            + (1.5998 * 1e-4) / (38.9 - s*s))
        catalog_data = catalog_data.rename(columns={"Name":"name"})
        # TODO line catalogs should have EXACTLY the same headers (this one is named "Name" or "name")

        return catalog_data

    def get_catalog_data(self, catalog_path, filters=[]):
        catalog_data = self.get_catalog_data_df(catalog_path,filters)
        ret["lambda"] = catalog_data["#lambda"].to_list()
        ret["name"] = catalog_data["name"].to_list()
        ret["force"] = catalog_data["force"].to_list()
        ret["nominal_ampl"] = catalog_data["nominal_ampl"].to_list()
        ret["type"] = catalog_data["type"].to_list()
        return ret

    def get_template_ratio_data(self, tplratio_file, filters):
        return self.get_catalog_data_df(self.get_template_ratio_path(tplratio_file), filters)

    def get_full_catalog_data(self,filters):
        return self.get_catalog_data(self.get_full_catalog_path(),filters)

    def get_full_catalog_df(self,filters):
        catalog_data = pd.read_csv(self.get_full_catalog_path(), sep='\t', header=3,
                           float_precision='round_trip', index_col="id")

        for criterion in filters:
            name = criterion["field"]
            value = criterion["value"]
            c_type = criterion["type"]
            if c_type == "=" and value != "All":
                catalog_data = catalog_data[catalog_data[name] == value]

        if "linecatalog_convert" in self.amazed_parameters and self.amazed_parameters["linecatalog_convert"]:
            s = 0.8
            catalog_data["#lambda"] = catalog_data["#lambda"] / (
                    1.0 + 8.34254 * 1e-5 + (2.406147 * 1e-2) / (130 - s * s)
                    + (1.5998 * 1e-4) / (38.9 - s * s))

        return catalog_data

    def get_shifted_full_catalog_data(self, redshift,filters):
        ret = self.get_full_catalog_df(filters)
        ret = ret.rename(columns={"#lambda":"FittedRaysLambdaRest"})
        ret["FittedRaysLambda"] = (1+redshift)*ret["FittedRaysLambdaRest"]
        return json.loads(ret.to_json(orient="records"))

    def get_shifted_template_ratio_data(self, template_ratio_file, redshift,
                                        spectrum_lambda_min, spectrum_lambda_max,
                                        filters = [], full_catalog=False):
        if full_catalog:
            tsd = self.get_full_catalog_df(filters)
            tsd["#lambda"] = (1+redshift)*tsd["#lambda"]
        else:
            tsd = self.get_template_ratio_data(template_ratio_file,filters)

        l_ref_min = spectrum_lambda_min/(1+redshift)
        l_ref_max = spectrum_lambda_max/(1+redshift)
#        tsd.append(Series({"#lambda":l_ref_min,"type":"A"}), ignore_index=True)
#        tsd.append(Series({"#lambda": l_ref_max}), ignore_index=True)
        res = dict()
        tsd = tsd[tsd["#lambda"]>=l_ref_min]
        tsd = tsd[tsd["#lambda"] <= l_ref_max]
        tsd_abs = tsd[tsd["type"]== "A"]
        tsd_emi = tsd[tsd["type"] == "E"]
        res["lambda_range"] =[l_ref_min,l_ref_max]
        res["absorption"] = {"lambda": tsd_abs["#lambda"].to_list(),
                             "name": tsd_abs["name"].to_list(),
                             "force":tsd_abs["force"].to_list(),
                             "ampl": tsd_abs["nominal_ampl"].to_list()}
        res["emission"] = {"lambda": tsd_emi["#lambda"].to_list(),
                             "name": tsd_emi["name"].to_list(),
                             "force": tsd_emi["force"].to_list(),
                             "ampl": tsd_emi["nominal_ampl"].to_list()}
        return res

    def export_input_spectrum_list(self, filtered_redshifts):
        merged = self.spectra.merge(filtered_redshifts,
                                    left_index=True,
                                    right_index=True)
        merged["ProcessingID"] = merged.index

        merged.to_csv('~/input_' + self.title + '.spectrumlist',
                      columns=["flux", "error", "ProcessingID"],
                      sep='\t',
                      header=False,
                      index=False)

        return '~/input_' + self.title + '.spectrumlist'
