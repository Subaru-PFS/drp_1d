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
import json
from typing import List

import pandas as pd
from pylibamazed.Exception import APIException
from pylibamazed.redshift import ErrorCode


class Parameters:
    def __init__(self, parameters: dict):
        self.parameters = parameters
        self.check_params()

    def get_solve_methods(self, object_type) -> dict:
        method = self.get_solve_method(object_type)
        linemeas_method = self.get_linemeas_method(object_type)
        methods = []
        if method:
            methods.append(method)
        if linemeas_method:
            methods.append(linemeas_method)
        return methods

    def get_redshift_sampling(self, object_type):
        return self.get_object_params(object_type).get("redshiftsampling")

    def get_linemodel_methods(self, object_type):
        methods = []
        linemeas_method = self.get_linemeas_method(object_type)
        solve_method = self.get_solve_method(object_type)
        if linemeas_method:
            methods.append(linemeas_method)
        if solve_method == "LineModelSolve":
            methods.append(solve_method)
        return methods

    def check_lmskipsecondpass(self, object_type):
        solve_method = self.get_solve_method(object_type)
        if solve_method:
            if solve_method != "LineModelSolve":
                return False
            else:
                return self.get_linemodel_params(object_type, solve_method).get("skipsecondpass")
        return False

    def get_objects_solve_methods(self):
        ret = dict()
        for object_type in self.get_objects():
            if self.get_solve_method(object_type):
                ret[object_type] = self.get_solve_method(object_type)
        return ret

    def get_objects_linemeas_methods(self):
        ret = dict()
        for object_type in self.get_objects():
            if self.get_linemeas_method(object_type):
                ret[object_type] = self.get_linemeas_method(object_type)
        return ret

    def get_solve_method(self, object_type):
        return self.get_object_params(object_type).get("method")

    def get_linemeas_method(self, object_type):
        return self.get_object_params(object_type).get("linemeas_method")

    def get_solve_method_params(self, object_type, solve_method):
        return self.get_object_params(object_type).get(solve_method)

    def get_linemodel_params(self, object_type, solve_method) -> dict:
        return self.get_solve_method_params(object_type, solve_method).get("linemodel")

    def get_linecatalog(self, object_type, solve_method, throw=True):
        linecatalog = self.get_linemodel_params(object_type, solve_method).get("linecatalog")
        if linecatalog is None and throw:
            raise APIException(
                ErrorCode.MISSING_PARAMETER,
                "Incomplete parameter file, {}.linemodel.linecatalog entry"
                " mandatory".format(solve_method)
            )
        return linecatalog

    def get_linemodel_nsigmasupport(self, object_type, solve_method):
        return self.get_linemodel_params(object_type, solve_method).get("nsigmasupport")

    def get_linemodel_igmfit(self, object_type, solve_method):
        return self.get_linemodel_params(object_type, solve_method).get("igmfit")

    def get_tplratio_ismfit(self, object_type, solve_method):
        return self.get_linemodel_params(object_type, solve_method).get("tplratio_ismfit")

    def get_tplratio_catalog(self, object_type, solve_method="LineModelSolve", throw=True):
        catalog = self.get_linemodel_params(object_type, solve_method).get("tplratio_catalog")
        if catalog is None and throw:
            raise APIException(
                ErrorCode.MISSING_PARAMETER,
                "Missing mandatory entry: LineModelSolve.linemodel.tplratio_catalog")
        return catalog

    def get_objects(self):
        return self.parameters.get("objects")

    def get_template_dir(self, object_type: str, throw=True):
        template_dir = self.get_object_params(object_type).get("template_dir")
        if template_dir is None and throw:
            raise APIException(
                ErrorCode.MISSING_PARAMETER,
                "Incomplete parameter file, template_dir entry mandatory"
            )
        return template_dir

    def get_ebmv(self) -> dict:
        return self.parameters.get("ebmv")

    def get_ebmv_count(self):
        return self.get_ebmv().get("count")

    def get_ebmv_step(self):
        return self.get_ebmv().get("step")

    def get_ebmv_start(self):
        return self.get_ebmv().get("start")

    def set_redshiftref(self, object_type, redshift_ref) -> None:
        self.get_object_params(object_type)["redshiftref"] = redshift_ref

    def set_velocity_absorption(self, object_type: str, solve_method, velocity_abs) -> None:
        self.get_linemodel_params(object_type, solve_method)["velocityabsorption"] = velocity_abs

    def set_velocity_emission(self, object_type: str, solve_method, velocity_em) -> None:
        self.get_linemodel_params(object_type, solve_method)["velocityemission"] = velocity_em

    def load_linemeas_parameters_from_catalog(self, source_id, config):
        for object_type in config["linemeascatalog"].keys():
            lm = pd.read_csv(config["linemeascatalog"][object_type], sep='\t', dtype={'ProcessingID': object})
            lm = lm[lm.ProcessingID == source_id]
            if lm.empty:
                raise APIException(
                    ErrorCode.INVALID_PARAMETER,
                    f"Uncomplete linemeas catalog, {source_id} missing"
                )

            columns = config["linemeas_catalog_columns"][object_type]
            redshift_ref = float(lm[columns["Redshift"]].iloc[0])
            velocity_abs = float(lm[columns["VelocityAbsorption"]].iloc[0])
            velocity_em = float(lm[columns["VelocityEmission"]].iloc[0])

            self.set_redshiftref(object_type, redshift_ref)
            self.set_velocity_absorption(object_type, "LineMeasSolve", velocity_abs)
            self.set_velocity_emission(object_type, "LineMeasSolve", velocity_em)

    def load_linemeas_parameters_from_result_store(self, output, object_type):

        redshift = output.get_attribute_from_source(object_type,
                                                    self.get_solve_method(object_type),
                                                    "model_parameters",
                                                    "Redshift",
                                                    0)
        self.parameters[object_type]["redshiftref"] = redshift
        velocity_abs = output.get_attribute_from_source(object_type,
                                                        self.get_solve_method(object_type),
                                                        "model_parameters",
                                                        "VelocityAbsorption",
                                                        0)
        velocity_em = output.get_attribute_from_source(object_type,
                                                       self.get_solve_method(object_type),
                                                       "model_parameters",
                                                       "VelocityEmission",
                                                       0)
        self.set_velocity_absorption(object_type, "LineMeasSolve", velocity_abs)
        self.set_velocity_emission(object_type, "LineMeasSolve", velocity_em)

    def get_json(self):
        return json.dumps(self.parameters)

    def reliability_enabled(self, object_type):
        return self.get_object_params(object_type).get("enable_reliability")

    def get_lineratio_type(self, object_type, solve_method) -> str:
        return self.get_linemodel_params(object_type, solve_method).get("lineRatioType")

    def is_tplratio_catalog_needed(self, object_type) -> bool:
        solve_method = self.get_solve_method(object_type)
        if solve_method == "LineModelSolve":
            return self.get_lineratio_type(object_type, solve_method) in ["tplratio", "tplcorr"]
        else:
            return False

    def stage_enabled(self, object_type, stage):
        if stage == "redshift_solver":
            return self.get_solve_method(object_type) is not None
        elif stage == "linemeas_solver":
            return self.get_linemeas_method(object_type) is not None
        elif stage == "linemeas_catalog_load":
            return self.get_linemeas_method(object_type) is not None \
                and self.get_solve_method(object_type) is None
        elif stage == "reliability_solver":
            return self.reliability_enabled(object_type)
        elif stage == "sub_classif_solver":
            return self.is_tplratio_catalog_needed(object_type)
        else:
            raise Exception("Unknown stage {stage}")

    def get_filters(self):
        return self.parameters.get("filters")

    def get_additional_cols(self) -> List[str]:
        return self.parameters.get("additional_cols")

    def get_photometry_band(self, throw=True):
        band = self.parameters.get("photometryBand")
        if band is None and throw:
            raise APIException(ErrorCode.MISSING_PARAMETER, "photometryBand parameter required")
        if len(band) == 0 and throw:
            raise APIException(ErrorCode.INVALID_PARAMETER, "photometryBand parameter is empty")
        return band

    def get_lsf(self) -> dict:
        return self.parameters.get("LSF")

    def get_lsf_type(self):
        return self.get_lsf().get("LSFType")

    def get_multiobs_method(self):
        registered_methods = [None, "", "merge", "full"]
        method = self.parameters.get("multiobsmethod")
        if method not in registered_methods:
            raise APIException(
                ErrorCode.INVALID_PARAMETER,
                f"multiobsmethod must be one of {registered_methods}"
            )
        return method

    def get_airvacuum_method(self):
        return self.parameters.get("airvacuum_method", "")

    def get_gaussian_variable_width_file_name(self):
        return self.get_lsf().get("GaussianVariablewidthFileName")

    def get_photometry_transmission_dir(self):
        return self.parameters.get("photometryTransmissionDir")

    def get_reliability_model(self, object_type):
        return self.get_object_params(object_type).get("reliability_model")

    def get_object_params(self, object_type) -> dict:
        return self.parameters.get(object_type)

    def get_redshift_range(self, object_type):
        return self.get_object_params(object_type).get("redshiftrange")

    def get_redshift_step(self, object_type):
        return self.get_object_params(object_type).get("redshiftstep")

    def get_lambda_range(self):
        return self.parameters.get("lambdarange")

    def set_lsf_type(self, lsf_type):
        self.parameters["LSF"]["LSFType"] = lsf_type

    def set_lsf_param(self, param_name, data):
        self.parameters["LSF"][param_name] = data

    def to_json(self):
        return json.dumps(self.parameters)

    def check_lineameas_validity(self):
        for object_type in self.get_objects():
            method = self.get_solve_method(object_type)
            if method == "LineModelSolve":
                if self.get_linemeas_method(object_type):
                    raise APIException(
                        ErrorCode.INCOHERENT_CONFIG_OPTION,
                        "Cannot run LineMeasSolve from catalog when sequencial processing is selected"
                        "simultaneously."
                    )

    def get_photometry_bands(self, object_type, method):
        if not self.is_photometry_activated(object_type, method):
            return [None]
        param_name = "photometryBand"
        if param_name not in self.parameters:
            return [None]
        return self.parameters[param_name]

    # photometry is only activated for TFSolve and LinemodelSolve+TF.
    # photometry is not yet implemented for TplcombinationSolve
    # photometry has no sens for LineMeasSolve, since "nocontinuum"
    def is_photometry_activated(self, object_type, method):
        # TODO this should be checked elsewhere, else we duplicate code -> us try except, there are other
        # places in this class with same problem
        #        if method == "TplcombinationSolve" or self.isfftprocessingActive(object_type):
        #            return False
        if method == "LineModelSolve":
            return self.parameters[object_type][method]["linemodel"]["enablephotometry"]
        try:
            return self.parameters[object_type][method]["enablephotometry"]
        except Exception:
            return False

    # a simple check
    def isfftprocessingActive(self, object_type):
        method = self.parameters[object_type]["method"]
        if method == "TplcombinationSolve":
            return False
        if method == "LineModelSolve":
            return self.parameters[object_type][method]["linemodel"]["continuumfit"]["fftprocessing"]

        return self.parameters[object_type][method]["fftprocessing"]

    def get_linemodel_continuumfit_params(self, object_type, solve_method):
        return self.get_linemodel_params(object_type, solve_method).get("continuumfit")

    # a list of checks to do for linemodel
    def check_linemodel_params(self, object_type, solve_method):
        # nothing to check here
        if solve_method != "LineModelSolve" or solve_method != "LineMeasSolve":
            return True
        # check1
        use_loglambdasampling = self.get_linemodel_params(
            object_type, solve_method).get("useloglambdasampling")
        fftprocessing = self.get_linemodel_continuumfit_params(object_type, solve_method).get("fftprocessing")
        if use_loglambdasampling and not fftprocessing:
            raise APIException(ErrorCode.INVALID_PARAMETER,
                               "useloglambdasampling cannot be activated if fftprocessing is deactivated")

    def check_params(self):
        # iterate over objects
        for object_type in self.get_objects():
            methods = self.get_solve_methods(object_type)
            for solve_method in methods:
                self.check_linemodel_params(object_type, solve_method)
