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
import pandas as pd
from typing import Optional

from pylibamazed.Parameters import Parameters
from pylibamazed.ParametersAccessor import EVelocityType, EVelocityFitParam, ParametersAccessor, ESolveMethod
from pylibamazed.ResultStoreOutput import ResultStoreOutput
from pylibamazed.redshift import ErrorCode, CLog, WarningCode, CFlagWarning
from pylibamazed.Exception import APIException

zlog = CLog.GetInstance()
zflag = CFlagWarning.GetInstance()


class LinemeasParameters:
    def __init__(self):
        self.redshift_ref = dict()
        self.velocity = {EVelocityType.Absorption: dict(), EVelocityType.Emission: dict()}

    def _get_catalog_velocity_name(self, velocity_type: EVelocityType) -> str:
        return f"Velocity{velocity_type.value}"

    def load_from_catalogs(self, source_id: str, catalogs: dict, catalog_columns: dict):
        for spectrum_model in catalogs.keys():
            lm = pd.read_csv(catalogs[spectrum_model], sep="\t", dtype={"ProcessingID": object})
            lm = lm[lm.ProcessingID == source_id]
            if lm.empty:
                raise APIException(
                    ErrorCode.INVALID_PARAMETER, f"Uncomplete linemeas catalog, {source_id} missing"
                )

            columns = catalog_columns[spectrum_model]
            self.redshift_ref[spectrum_model] = float(lm[columns["Redshift"]].iloc[0])
            for velocity_type in self.velocity:
                self.velocity[velocity_type][spectrum_model] = float(
                    lm[columns[self._get_catalog_velocity_name(velocity_type)]].iloc[0]
                )

    def load_from_result_store(
        self, parameter: Parameters, output: ResultStoreOutput, spectrum_model: str
    ) -> None:
        redshift_solver_method: Optional[ESolveMethod] = getattr(
            parameter.get_redshift_solver_method(spectrum_model), "value", None
        )
        self.redshift_ref[spectrum_model] = output.get_attribute_from_source(
            spectrum_model,
            "redshiftSolver",
            redshift_solver_method,
            "model_parameters",
            "Redshift",
            0,
        )

        for velocity_type in self.velocity:
            if redshift_solver_method == "lineModelSolve":
                self.velocity[velocity_type][spectrum_model] = self._get_velocity_from_result_store(
                    parameter, output, spectrum_model, velocity_type
                )
            else:
                zlog.LogInfo("velocities not measured, using lineMeasSolver parameter values")
                self.velocity[velocity_type][spectrum_model] = self._get_velocity_from_linemeas_parameters(
                    parameter, spectrum_model, velocity_type
                )

    def update_parameters(self, parameter: Parameters) -> None:
        for spectrum_model in self.redshift_ref.keys():
            self._replace_nan_with_parameter(parameter, spectrum_model)
            parameter.set_redshiftref(spectrum_model, self.redshift_ref[spectrum_model])
            if parameter.get_linemodel_section(spectrum_model, ESolveMethod.LINE_MEAS) is not None:
                for velocity_type in self.velocity:
                    parameter.set_velocity(
                        spectrum_model,
                        ESolveMethod.LINE_MEAS,
                        velocity_type,
                        self.velocity[velocity_type][spectrum_model],
                    )
                self._check_velocity_range(parameter, spectrum_model)

    def _get_velocity_from_result_store(
        self,
        parameter: Parameters,
        output: ResultStoreOutput,
        spectrum_model: str,
        velocity_type: EVelocityType,
    ):
        redshift_solver: Optional[ESolveMethod] = getattr(
            parameter.get_redshift_solver_method(spectrum_model), "value", None
        )
        velocity_name = self._get_catalog_velocity_name(velocity_type)
        try:
            velocity = output.get_attribute_from_source(
                spectrum_model,
                "redshiftSolver",
                redshift_solver,
                "model_parameters",
                velocity_name,
                0,
            )
        except Exception:
            raise APIException(
                ErrorCode.OUTPUT_READER_ERROR,
                f"missing output: {spectrum_model}.redshiftSolver.{redshift_solver}.{velocity_name}",
            ) from None
        return velocity

    def _get_velocity_from_linemeas_parameters(
        self, parameter: Parameters, spectrum_model: str, velocity_type: EVelocityType
    ):
        solve_method = parameter.get_linemeas_method(spectrum_model)
        velocity = parameter.get_velocity(spectrum_model, solve_method, velocity_type)
        if velocity is None:
            raise APIException(
                ErrorCode.INVALID_PARAMETER_FILE,
                f"missing parameter {spectrum_model}.lineMeasSolver.{solve_method.value}.lineModel"
                f"{ParametersAccessor.get_velocity_name(velocity_type)}",
            )
        return velocity

    def _replace_nan_with_parameter(self, parameter: Parameters, spectrum_model: str):
        for velocity_type in self.velocity:
            if pd.isna(self.velocity[velocity_type][spectrum_model]):
                self.velocity[velocity_type][spectrum_model] = self._get_velocity_from_linemeas_parameters(
                    parameter, spectrum_model, velocity_type
                )
                zlog.LogInfo(
                    (
                        f"{ParametersAccessor.get_velocity_name(velocity_type)} "
                        "not measured, using lineMeasSolver parameter value: "
                        f"{self.velocity[velocity_type][spectrum_model]}"
                    )
                )

    def _check_velocity_range(self, parameter: Parameters, spectrum_model: str):
        velocityfit: bool = parameter.get_linemeas_velocity_fit(spectrum_model)
        if not velocityfit:
            return
        for velocity_type in EVelocityType:
            velocity = parameter.get_velocity(spectrum_model, ESolveMethod.LINE_MEAS, velocity_type)
            velocity_bound = parameter.get_velocity_fit_param(
                spectrum_model, ESolveMethod.LINE_MEAS, velocity_type, EVelocityFitParam.Min
            )
            if velocity < velocity_bound:
                warning_message = (
                    f"lineMeasSolve input velocity {parameter.get_velocity_name(velocity_type)}"
                    f"={velocity} for object {spectrum_model} is below "
                    f"{parameter.get_velocity_fit_param_name(velocity_type, EVelocityFitParam.Min)}"
                    f"={velocity_bound}:"
                    " extending velocity fit range"
                )
                zflag.warning(WarningCode.VELOCITY_FIT_RANGE, warning_message)
                parameter.set_velocity_fit_param(
                    spectrum_model, ESolveMethod.LINE_MEAS, velocity_type, EVelocityFitParam.Min, velocity
                )

            velocity_bound = parameter.get_velocity_fit_param(
                spectrum_model, ESolveMethod.LINE_MEAS, velocity_type, EVelocityFitParam.Max
            )
            if velocity > velocity_bound:
                warning_message = (
                    f"{ESolveMethod.LINE_MEAS.value} input velocity {parameter.get_velocity_name(velocity_type)}"
                    f"={velocity} for object {spectrum_model} is above "
                    f"{parameter.get_velocity_fit_param_name(velocity_type, EVelocityFitParam.Max)}"
                    f"={velocity_bound}:"
                    " extending velocity fit range"
                )
                zflag.warning(WarningCode.VELOCITY_FIT_RANGE, warning_message)
                parameter.set_velocity_fit_param(
                    spectrum_model, ESolveMethod.LINE_MEAS, velocity_type, EVelocityFitParam.Max, velocity
                )
