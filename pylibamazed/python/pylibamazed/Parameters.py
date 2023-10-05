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

import pandas as pd
from pylibamazed.Exception import APIException
from pylibamazed.ParametersAccessor import ParametersAccessor
from pylibamazed.ParametersChecker import ParametersChecker
from pylibamazed.redshift import ErrorCode


class Parameters(ParametersAccessor):
    def __init__(self, parameters_dict: dict, Checker=ParametersChecker, make_checks=True):
        self.parameters = parameters_dict
        if make_checks:
            accessor = ParametersAccessor(self.parameters)
            Checker(accessor).check()

    def get_solve_methods(self, object_type) -> dict:
        method = self.get_solve_method(object_type)
        linemeas_method = self.get_linemeas_method(object_type)
        methods = []
        if method:
            methods.append(method)
        if linemeas_method:
            methods.append(linemeas_method)
        return methods

    def get_linemodel_methods(self, object_type):
        methods = []
        linemeas_method = self.get_linemeas_method(object_type)
        solve_method = self.get_solve_method(object_type)
        if linemeas_method:
            methods.append(linemeas_method)
        if solve_method == "LineModelSolve":
            methods.append(solve_method)
        return methods

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

    def is_tplratio_catalog_needed(self, object_type) -> bool:
        solve_method = self.get_solve_method(object_type)
        if solve_method == "LineModelSolve":
            return self.get_lineModelSolve_lineRatioType(object_type) in ["tplratio", "tplcorr"]
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
            return self.get_reliability_enabled(object_type)
        elif stage == "sub_classif_solver":
            return self.is_tplratio_catalog_needed(object_type)
        else:
            raise Exception("Unknown stage {stage}")

    def set_lsf_type(self, lsf_type):
        self.parameters["LSF"]["LSFType"] = lsf_type

    def set_lsf_param(self, param_name, data):
        self.parameters["LSF"][param_name] = data

    def set_redshiftref(self, object_type, redshift_ref) -> None:
        self.get_object_section(object_type)["redshiftref"] = redshift_ref

    def set_velocity_absorption(self, object_type: str, solve_method, velocity_abs) -> None:
        self.get_linemodel_section(object_type, solve_method)["velocityabsorption"] = velocity_abs

    def set_velocity_emission(self, object_type: str, solve_method, velocity_em) -> None:
        self.get_linemodel_section(object_type, solve_method)["velocityemission"] = velocity_em

    def to_json(self):
        return json.dumps(self.parameters)
  
    def is_a_redshift_solver_used(self) -> bool:
        z_solver_found: bool = False
        for object_type in self.get_objects():
            if self.get_solve_method(object_type) is not None:
                z_solver_found = True
                break
        print("z_solver_found", z_solver_found)
        return z_solver_found

    def check_linemeas_validity(self):
        for object_type in self.get_objects():
            method = self.get_solve_method(object_type)
            if method == "LineModelSolve":
                if self.get_linemeas_method(object_type):
                    raise APIException(
                        ErrorCode.INCOHERENT_CONFIG_OPTION,
                        "Cannot run LineMeasSolve from catalog when sequencial processing is selected"
                        "simultaneously."
                    )
