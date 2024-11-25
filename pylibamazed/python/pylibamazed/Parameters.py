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

from enum import Enum
import pandas as pd
from pylibamazed.Exception import APIException, exception_decorator
from pylibamazed.ParametersAccessor import ParametersAccessor
from pylibamazed.ParametersChecker import ParametersChecker
from pylibamazed.ParametersConverter import ParametersConverterSelector
from pylibamazed.ParametersExtender import ParametersExtender
from pylibamazed.redshift import ErrorCode


class Parameters(ParametersAccessor):
    defined_stages = ["redshiftSolver", "lineMeasSolver", "reliabilitySolver"]

    @exception_decorator(logging=True)
    def __init__(
        self,
        raw_params: dict,
        make_checks=True,
        accepts_v1=False,
        Checker=ParametersChecker,
        ConverterSelector=ParametersConverterSelector,
        Extender=ParametersExtender,
    ):
        version = self.get_json_schema_version(raw_params)

        if make_checks:
            Checker(raw_params, version).json_schema_check()

        converter = ConverterSelector(accepts_v1).get_converter(version)
        converted_parameters = converter().convert(raw_params)

        if make_checks:
            Checker(converted_parameters, version).custom_check()

        extended_parameters = Extender(version).extend(converted_parameters)
        self.parameters = extended_parameters
        self.remove_unused_solvers()

    def get_json_schema_version(self, raw_parameters: dict):
        version = raw_parameters.get("version")
        if version is None:
            version = 1
        if type(version) is not int:
            raise APIException(ErrorCode.INVALID_PARAMETER_FILE, "Parameter version must be an integer")
        return version

    def get_solve_methods(self, spectrum_model) -> dict:
        method = self.get_redshift_solver_method(spectrum_model)
        linemeas_method = self.get_linemeas_method(spectrum_model)
        methods = []
        if method:
            methods.append(method)
        if linemeas_method:
            methods.append(linemeas_method)
        return methods

    def get_stage_from_method(self, method: str) -> str:
        if method == "lineMeasSolve":
            return "lineMeasSolver"
        else:
            return "redshiftSolver"

    def remove_unused_solvers(self) -> None:
        for spectrum_model in self.get_spectrum_models([]):
            for stage in self.defined_stages:
                if not stage in self.get_stages(spectrum_model):
                    if stage in self.parameters.get(spectrum_model, []):
                        del self.parameters[spectrum_model][stage]

    def get_linemodel_methods(self, spectrum_model):
        methods = []
        linemeas_method = self.get_linemeas_method(spectrum_model)
        solve_method = self.get_redshift_solver_method(spectrum_model)
        if linemeas_method:
            methods.append(linemeas_method)
        if solve_method == "lineModelSolve":
            methods.append(solve_method)
        return methods

    def get_objects_solve_methods(self):
        ret = dict()
        for spectrum_model in self.get_spectrum_models():
            if self.get_redshift_solver_method(spectrum_model):
                ret[spectrum_model] = self.get_redshift_solver_method(spectrum_model)
        return ret

    def get_objects_linemeas_methods(self):
        ret = dict()
        for spectrum_model in self.get_spectrum_models():
            if self.get_linemeas_method(spectrum_model):
                ret[spectrum_model] = self.get_linemeas_method(spectrum_model)
        return ret

    def load_linemeas_parameters_from_catalog(self, source_id, config):
        for spectrum_model in config["linemeascatalog"].keys():
            lm = pd.read_csv(
                config["linemeascatalog"][spectrum_model], sep="\t", dtype={"ProcessingID": object}
            )
            lm = lm[lm.ProcessingID == source_id]
            if lm.empty:
                raise APIException(
                    ErrorCode.INVALID_PARAMETER, f"Uncomplete linemeas catalog, {source_id} missing"
                )

            columns = config["linemeas_catalog_columns"][spectrum_model]
            redshift_ref = float(lm[columns["Redshift"]].iloc[0])
            velocity_abs = float(lm[columns["VelocityAbsorption"]].iloc[0])
            velocity_em = float(lm[columns["VelocityEmission"]].iloc[0])

            self.set_redshiftref(spectrum_model, redshift_ref)
            if self.get_linemodel_section(spectrum_model, "lineMeasSolve") is not None:
                self.set_velocity_absorption(spectrum_model, "lineMeasSolve", velocity_abs)
                self.set_velocity_emission(spectrum_model, "lineMeasSolve", velocity_em)

    def load_linemeas_parameters_from_result_store(self, output, spectrum_model):
        redshift = output.get_attribute_from_source(
            spectrum_model,
            "redshiftSolver",
            self.get_redshift_solver_method(spectrum_model),
            "model_parameters",
            "Redshift",
            0,
        )
        self.parameters[spectrum_model]["redshiftref"] = redshift
        velocity_abs = output.get_attribute_from_source(
            spectrum_model,
            "redshiftSolver",
            self.get_redshift_solver_method(spectrum_model),
            "model_parameters",
            "VelocityAbsorption",
            0,
        )
        velocity_em = output.get_attribute_from_source(
            spectrum_model,
            "redshiftSolver",
            self.get_redshift_solver_method(spectrum_model),
            "model_parameters",
            "VelocityEmission",
            0,
        )
        if self.get_linemodel_section(spectrum_model, "lineMeasSolve") is not None:
            self.set_velocity_absorption(spectrum_model, "lineMeasSolve", velocity_abs)
            self.set_velocity_emission(spectrum_model, "lineMeasSolve", velocity_em)

    def is_tplratio_catalog_needed(self, spectrum_model) -> bool:
        solve_method = self.get_redshift_solver_method(spectrum_model)
        if solve_method == "lineModelSolve":
            return self.get_linemodel_line_ratio_type(spectrum_model) in ["tplRatio", "tplCorr"]
        else:
            return False

    def stage_enabled(self, spectrum_model, stage):
        if stage == "redshiftSolver":
            return self.get_redshift_solver_method(spectrum_model) is not None
        elif stage == "lineMeasSolver":
            return self.get_linemeas_method(spectrum_model) is not None
        elif stage == "linemeas_catalog_load":
            return (
                self.get_linemeas_method(spectrum_model) is not None
                and self.get_redshift_solver_method(spectrum_model) is None
            )
        elif stage == "reliabilitySolver":
            return self.get_reliability_enabled(spectrum_model)
        elif stage == "subClassifSolver":
            return self.is_tplratio_catalog_needed(spectrum_model)
        else:
            raise APIException(ErrorCode.INTERNAL_ERROR, "Unknown stage {stage}")

    def set_lsf_param(self, param_name, data):
        self.parameters["lsf"][param_name] = data

    def set_redshiftref(self, spectrum_model, redshift_ref) -> None:
        self.get_spectrum_model_section(spectrum_model)["redshiftref"] = redshift_ref

    def set_velocity_absorption(self, spectrum_model: str, solve_method, velocity_abs) -> None:
        self.get_linemodel_section(spectrum_model, solve_method)["velocityAbsorption"] = velocity_abs

    def set_velocity_emission(self, spectrum_model: str, solve_method, velocity_em) -> None:
        self.get_linemodel_section(spectrum_model, solve_method)["velocityEmission"] = velocity_em

    def to_json(self):
        return json.dumps(self.parameters)

    def is_a_redshift_solver_used(self) -> bool:
        z_solver_found: bool = False
        for spectrum_model in self.get_spectrum_models([]):
            if self.get_redshift_solver_method(spectrum_model) is not None:
                z_solver_found = True
                break
        return z_solver_found

    def check_linemeas_validity(self):
        for spectrum_model in self.get_spectrum_models():
            method = self.get_redshift_solver_method(spectrum_model)
            if "redshiftSolver" in self.get_stages(spectrum_model) and method == "lineModelSolve":
                if self.get_linemeas_method(spectrum_model):
                    raise APIException(
                        ErrorCode.INCOHERENT_CONFIG_OPTIONS,
                        "Cannot run LineMeasSolve from catalog when sequencial processing is selected"
                        "simultaneously.",
                    )

    def get_lambda_range_min(self):
        if self.get_multiobs_method() != "full":
            return self.get_lambda_range()[0]
        else:
            ret = 10 ^ 9
            for obs_id in self.get_observation_ids():
                cur = self.get_lambda_range(obs_id)[0]
                if cur < ret:
                    ret = cur
            return ret

    def get_lambda_range_max(self):
        if self.get_multiobs_method() != "full":
            return self.get_lambda_range()[0]
        else:
            ret = 0
            for obs_id in self.get_observation_ids():
                cur = self.get_lambda_range(obs_id)[1]
                if cur > ret:
                    ret = cur
            return ret
