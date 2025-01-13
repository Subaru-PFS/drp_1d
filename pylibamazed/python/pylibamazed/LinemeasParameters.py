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
from pylibamazed.Parameters import Parameters
from pylibamazed.ResultStoreOutput import ResultStoreOutput
from pylibamazed.redshift import ErrorCode, CLog, WarningCode, CFlagWarning
from pylibamazed.Exception import APIException

zlog = CLog.GetInstance()
zflag = CFlagWarning.GetInstance()


class LinemeasParameters:
    def __init__(self):
        self.redshift_ref = dict()
        self.velocity_abs = dict()
        self.velocity_em = dict()

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
            self.velocity_abs[spectrum_model] = float(lm[columns["VelocityAbsorption"]].iloc[0])
            self.velocity_em[spectrum_model] = float(lm[columns["VelocityEmission"]].iloc[0])

    def load_from_result_store(
        self, parameter: Parameters, output: ResultStoreOutput, spectrum_model: str
    ) -> None:
        redshift_solver = parameter.get_redshift_solver_method(spectrum_model)
        self.redshift_ref[spectrum_model] = output.get_attribute_from_source(
            spectrum_model,
            "redshiftSolver",
            redshift_solver,
            "model_parameters",
            "Redshift",
            0,
        )

        if redshift_solver == "lineModelSolve":
            self.velocity_abs[spectrum_model] = self._get_velocity_from_result_store(
                parameter, output, spectrum_model, "VelocityAbsorption"
            )
            self.velocity_em[spectrum_model] = self._get_velocity_from_result_store(
                parameter, output, spectrum_model, "VelocityEmission"
            )
        else:
            zlog.LogInfo("velocities not measured, using lineMeasSolver parameter values")
            self.velocity_abs[spectrum_model] = self._get_velocity_from_parameters(
                parameter, spectrum_model, "VelocityAbsorption"
            )
            self.velocity_em[spectrum_model] = self._get_velocity_from_parameters(
                parameter, spectrum_model, "VelocityEmission"
            )

    def update_parameters(self, parameter: Parameters) -> None:
        for spectrum_model in self.redshift_ref.keys():
            self._replace_nan_with_parameter(parameter, spectrum_model)
            parameter.set_redshiftref(spectrum_model, self.redshift_ref[spectrum_model])
            if parameter.get_linemodel_section(spectrum_model, "lineMeasSolve") is not None:
                parameter.set_velocity_absorption(
                    spectrum_model, "lineMeasSolve", self.velocity_abs[spectrum_model]
                )
                parameter.set_velocity_emission(
                    spectrum_model, "lineMeasSolve", self.velocity_em[spectrum_model]
                )
                self._check_velocity_range(parameter, spectrum_model)

    def _get_velocity_from_result_store(
        self, parameter: Parameters, output: ResultStoreOutput, spectrum_model: str, velocity_name: str
    ):
        if velocity_name not in ["VelocityAbsorption", "VelocityEmission"]:
            raise APIException(ErrorCode.APIException, "invalid velocity_name")

        redshift_solver = parameter.get_redshift_solver_method(spectrum_model)
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
                f"missing output: {spectrum_model}.{redshift_solver}.{velocity_name}",
            ) from None
        return velocity

    def _get_velocity_from_parameters(self, parameter: Parameters, spectrum_model: str, velocity_name: str):
        linemeas_solver = parameter.get_linemeas_method(spectrum_model)
        if velocity_name == "VelocityAbsorption":
            velocity = parameter.get_velocity_absorption(spectrum_model, linemeas_solver)
        elif velocity_name == "VelocityEmission":
            velocity = parameter.get_velocity_emission(spectrum_model, linemeas_solver)
        return velocity

    def _replace_nan_with_parameter(self, parameter: Parameters, spectrum_model: str):
        if pd.isna(self.velocity_abs[spectrum_model]):
            self.velocity_abs[spectrum_model] = self._get_velocity_from_parameters(
                parameter, spectrum_model, "VelocityAbsorption"
            )
            zlog.LogInfo(
                (
                    "VelocityAbsorption not measured, using lineMeasSolver parameter value: "
                    + f"{self.velocity_abs[spectrum_model]}"
                )
            )
        if pd.isna(self.velocity_em[spectrum_model]):
            self.velocity_em[spectrum_model] = self._get_velocity_from_parameters(
                parameter, spectrum_model, "VelocityEmission"
            )
            zlog.LogInfo(
                (
                    "VelocityEmission not measured, using lineMeasSolver parameter value: "
                    + f"{self.velocity_em[spectrum_model]}"
                )
            )

    def _check_velocity_range(self, parameter: Parameters, spectrum_model: str):
        velocityfit: bool = parameter.get_linemeas_velocity_fit(spectrum_model)
        if not velocityfit:
            return
        params = [
            {"param": "velocityEmission", "fitparam": ["emVelocityFitMin", "emVelocityFitMax"]},
            {"param": "velocityAbsorption", "fitparam": ["absVelocityFitMin", "absVelocityFitMax"]},
        ]
        for param_dict in params:
            velocity = parameter.get_linemeas_velocity_fit_param(spectrum_model, param_dict["param"])

            velocity_bound = parameter.get_linemeas_velocity_fit_param(
                spectrum_model, param_dict["fitparam"][0]
            )
            if velocity < velocity_bound:
                warning_message = (
                    f"lineMeasSolve input velocity {param_dict['param']}={velocity}"
                    + f" for object {spectrum_model}"
                    + f" is below {param_dict['fitparam'][0]}={velocity_bound}:"
                    + " extending velocity fit range"
                )
                zflag.warning(WarningCode.VELOCITY_FIT_RANGE, warning_message)
                parameter.set_linemeas_velocity_fit_param(spectrum_model, param_dict["fitparam"][0], velocity)

            velocity_bound = parameter.get_linemeas_velocity_fit_param(
                spectrum_model, param_dict["fitparam"][1]
            )
            if velocity > velocity_bound:
                warning_message = (
                    f"lineMeasSolve input velocity {param_dict['param']}={velocity}"
                    + f" for object {spectrum_model}"
                    + f" is above {param_dict['fitparam'][1]}={velocity_bound}:"
                    + " extending velocity fit range"
                )
                zflag.warning(WarningCode.VELOCITY_FIT_RANGE, warning_message)
                parameter.set_linemeas_velocity_fit_param(spectrum_model, param_dict["fitparam"][1], velocity)
