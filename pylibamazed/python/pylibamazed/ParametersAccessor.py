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
from typing import List, Optional

from enum import Enum
from pylibamazed.Exception import APIException, exception_decorator
from pylibamazed.redshift import CLog, ErrorCode


zlog = CLog.GetInstance()


class ESolveMethod(Enum):
    LINE_MODEL = "lineModelSolve"
    TEMPLATE_FITTING = "templateFittingSolve"
    TEMPLATE_COMBINATION = "tplCombinationSolve"
    LINE_MEAS = "lineMeasSolve"


class EContinuumFit(Enum):
    FROM_FIRST_PASS = "fromFirstPass"
    RETRY_ALL = "retryAll"
    REFIT_FIRST_PASS = "reFitFirstPass"


class EVelocityType(Enum):
    Absorption = "Absorption"
    Emission = "Emission"


class EVelocityFitParam(Enum):
    Min = "Min"
    Max = "Max"
    Step = "Step"


class ParametersAccessor:
    velocity_fit_prefix_dict = {
        EVelocityType.Absorption: "abs",
        EVelocityType.Emission: "em",
    }

    def __init__(self, parameters: dict):
        self.parameters = parameters

    def get_lambda_range(self, obs_id=""):
        """Depending on multiobs method, lambda range is not of the same type:
        - mono obs or merge => lambdarange is a range [min, max]
        - full => lambdarange is a dict of ranges {obs_id: [min, max], ...}
        """
        multiobs = self.get_multiobs_method()
        if multiobs in ["", "merge"]:
            return self.parameters["lambdaRange"]
        else:
            return self.parameters["lambdaRange"][obs_id]

    def get_lambda_ranges(self):
        multiobs = self.get_multiobs_method()
        if multiobs in ["", "merge"]:
            return {"": self.parameters["lambdaRange"]}
        else:
            return self.parameters["lambdaRange"]

    def get_airvacuum_method(self):
        return self.parameters.get("airVacuumMethod", "")

    def get_photometry_transmission_dir(self) -> Optional[str]:
        return self.parameters.get("photometryTransmissionDir")

    def get_photometry_bands(self) -> List[str]:
        return self.parameters.get("photometryBand", [])

    @exception_decorator
    def get_multiobs_method(self) -> Optional[str]:
        return self.parameters.get("multiObsMethod")

    def get_spectrum_models(self, default=None) -> List[str]:
        return self.parameters.get("spectrumModels", default)

    def get_linemeas_runmode(self) -> Optional[str]:
        return self.parameters.get("lineMeasRunMode")

    def get_spectrum_model_section(self, spectrum_model, create=False) -> dict:
        spectrum_model_section = self.parameters.get(spectrum_model, {})
        if create and spectrum_model_section == {}:
            self.parameters[spectrum_model] = spectrum_model_section
        return spectrum_model_section

    def get_stages(self, spectrum_model: str) -> list:
        return self.get_spectrum_model_section(spectrum_model).get("stages", [])

    def get_redshift_solver_section(self, spectrum_model, create=False) -> dict:
        return self._get_or_create_section(
            self.get_spectrum_model_section, "redshiftSolver", create, spectrum_model
        )

    def get_redshift_solver_method(self, spectrum_model: str) -> Optional[ESolveMethod]:
        if "redshiftSolver" not in self.get_stages(spectrum_model):
            return None
        method_str = self._get_on_None(self.get_redshift_solver_section(spectrum_model), "method")
        # Not in the get because if method is defined as an empty string in the parameters json
        # we still want method to be "None"
        if method_str == "" or method_str is None:
            return None
        return ESolveMethod(method_str)

    def get_linemeas_solver_section(self, spectrum_model: str, create: bool = False) -> dict:
        return self._get_or_create_section(
            self.get_spectrum_model_section, "lineMeasSolver", create, spectrum_model
        )

    def get_linemeas_method(self, spectrum_model: str) -> Optional[ESolveMethod]:
        if "lineMeasSolver" not in self.get_stages(spectrum_model):
            return None
        method = self._get_on_None(self.get_linemeas_solver_section(spectrum_model), "method")
        if method is None:
            return None
        return ESolveMethod(method)

    def get_linemeas_solve_section(self, spectrum_model: str, create: bool = False) -> dict:
        return self._get_or_create_section(
            self.get_linemeas_solver_section,
            ESolveMethod.LINE_MEAS.value,
            create,
            spectrum_model,
        )

    def get_linemeas_dzhalf(self, spectrum_model: str) -> Optional[float]:
        return self.get_spectrum_model_section(spectrum_model).get("lineMeasDzHalf")

    def get_redshiftrange(self, spectrum_model: str) -> Optional[List[float]]:
        return self.get_spectrum_model_section(spectrum_model).get("redshiftRange")

    def get_redshiftstep(self, spectrum_model: str) -> Optional[float]:
        return self.get_spectrum_model_section(spectrum_model).get("redshiftStep")

    def get_linemeas_redshiftstep(self, spectrum_model: str) -> Optional[float]:
        return self.get_spectrum_model_section(spectrum_model).get("lineMeasRedshiftStep")

    def set_redshiftref(self, spectrum_model, redshift_ref) -> None:
        self.get_spectrum_model_section(spectrum_model)["redshiftref"] = redshift_ref

    def get_reliability_enabled(self, spectrum_model: str) -> bool:
        return "reliabilitySolver" in self.get_stages(spectrum_model)

    def get_reliability_section(self, spectrum_model: str, create: bool = False) -> dict:
        return self._get_or_create_section(
            self.get_spectrum_model_section, "reliabilitySolver", create, spectrum_model
        )

    def get_reliability_methods(self, spectrum_model: str, default=[]) -> Optional[List[str]]:
        if "reliabilitySolver" not in self.get_stages(spectrum_model):
            return default
        return self._get_on_None(self.get_reliability_section(spectrum_model), "method", default=default)

    def get_deep_learning_solver_section(self, spectrum_model: str, create: bool = False) -> str:
        return self._get_or_create_section(
            self.get_reliability_section, "deepLearningSolver", create, spectrum_model
        )

    def get_reliability_model(self, spectrum_model: str) -> str:
        return self._get_on_None(self.get_deep_learning_solver_section(spectrum_model), "reliabilityModel")

    def get_sk_learn_classifier_solver_section(self, spectrum_model: str, create: bool = False) -> str:
        return self._get_or_create_section(
            self.get_reliability_section, "skLearnClassifier", create, spectrum_model
        )

    def get_sk_learn_classifier(self, spectrum_model: str) -> str:
        return self._get_on_None(
            self.get_sk_learn_classifier_solver_section(spectrum_model), "skLearnClassifier"
        )

    def get_sk_learn_classifier_file(self, spectrum_model: str) -> str:
        return self._get_on_None(
            self.get_sk_learn_classifier_solver_section(spectrum_model), "classifierFile"
        )

    def get_template_dir(self, spectrum_model: str) -> Optional[str]:
        return self.get_spectrum_model_section(spectrum_model).get("templateDir")

    def get_template_catalog_section(self, create: bool = False) -> Optional[str]:
        template_catalog_section = self.parameters.get("templateCatalog")
        if create and template_catalog_section is None:
            template_catalog_section = {}
            self.parameters["templateCatalog"] = template_catalog_section
        return template_catalog_section

    def get_redshift_solver_method_section(self, spectrum_model: str) -> Optional[dict]:
        method = self.get_redshift_solver_method(spectrum_model)
        if not method:
            return None
        return self._get_on_None(self.get_redshift_solver_section(spectrum_model), method.value)

    def photometry_is_enabled(self):
        for spectrum_model in self.get_spectrum_models([]):
            method = self.get_redshift_solver_method(spectrum_model)
            method_section = self.get_redshift_solver_method_section(spectrum_model)
            if method == ESolveMethod.LINE_MODEL:
                method_section = self._get_on_None(method_section, "lineModel")
            if self._get_on_None(method_section, "enablePhotometry", False):
                return True
        return False

    def get_additional_cols(self, default=None) -> List[str]:
        return self.parameters.get("additionalCols") or default

    def get_filters(self, default=[], obs_id=""):
        if not obs_id:
            return self.parameters.get("filters", default)
        else:
            if self.parameters.get("filters"):
                try:
                    return self.parameters.get("filters").get(obs_id, default)
                except AttributeError:
                    return self.parameters.get("filters")
            else:
                return default

    def get_lsf(self) -> Optional[dict]:
        return self.parameters.get("lsf")

    def get_lsf_type(self):
        return self._get_on_None(self.get_lsf(), "lsfType")

    def get_lsf_width(self):
        return self._get_on_None(self.get_lsf(), "width")

    def get_lsf_resolution(self):
        return self._get_on_None(self.get_lsf(), "resolution")

    def get_lsf_sourcesize(self):
        return self._get_on_None(self.get_lsf(), "sourceSize")

    def get_lsf_width_file_name(self):
        return self._get_on_None(self.get_lsf(), "gaussianVariableWidthFileName")

    def set_lsf_param(self, param_name, data):
        self.parameters["lsf"][param_name] = data

    def get_continuum_removal_section(self, from_template_catalog: bool = False, create: bool = False):
        if from_template_catalog:
            continuum_removal = self.get_template_catalog_continuum_removal_section(create)
        else:
            continuum_removal = self.get_root_continuum_removal_section(create)
        return continuum_removal

    def get_root_continuum_removal_section(self, create: bool = False):
        continuum_removal_section = self.parameters.get("continuumRemoval")
        if create and continuum_removal_section is None:
            continuum_removal_section = {}
            self.parameters["continuumRemoval"] = continuum_removal_section
        return continuum_removal_section

    def get_template_catalog_continuum_removal_section(self, create: bool = False):
        return self._get_or_create_section(self.get_template_catalog_section, "continuumRemoval", create)

    def get_continuum_removal_method(self, nesting: bool = False):
        return self._get_on_None(self.get_continuum_removal_section(nesting), "method")

    def get_continuum_removal_median_kernel_width(self, nesting: bool = False):
        return self._get_on_None(self.get_continuum_removal_section(nesting), "medianKernelWidth")

    def get_continuum_median_kernel_reflection(self, nesting: bool = False):
        return self._get_on_None(self.get_continuum_removal_section(nesting), "medianEvenReflection")

    def _get_on_None(self, dict, key, default=None):
        if dict is None:
            return default
        else:
            return dict.get(key, default)

    def get_template_fitting_section(self, spectrum_model: str, create: bool = False) -> dict:
        return self._get_or_create_section(
            self.get_redshift_solver_section,
            ESolveMethod.TEMPLATE_FITTING.value,
            create,
            spectrum_model,
        )

    def get_template_fitting_fft(self, spectrum_model: str) -> dict:
        return self._get_on_None(self.get_template_fitting_section(spectrum_model), "fftProcessing")

    def get_template_fitting_ism(self, spectrum_model: str) -> bool:
        return self._get_on_None(self.get_template_fitting_section(spectrum_model), "ismFit")

    def get_template_fitting_igm(self, spectrum_model: str) -> bool:
        return self._get_on_None(self.get_template_fitting_section(spectrum_model), "igmFit")

    def get_template_fitting_photometry_enabled(self, spectrum_model: str) -> bool:
        return self._get_on_None(self.get_template_fitting_section(spectrum_model), "enablePhotometry")

    def get_template_fitting_photometry_section(self, spectrum_model: str, create: bool = False) -> dict:
        return self._get_or_create_section(
            self.get_template_fitting_section, "photometry", create, spectrum_model
        )

    def get_template_fitting_photometry_weight(self, spectrum_model: str) -> dict:
        return self._get_on_None(self.get_template_fitting_photometry_section(spectrum_model), "weight")

    def get_template_fitting_spectrum_component(self, spectrum_model: str) -> str:
        return self._get_on_None(
            self._get_on_None(self.get_template_fitting_section(spectrum_model), "spectrum"), "component"
        )

    def get_template_combination_section(self, spectrum_model: str, create: bool = False) -> dict:
        return self._get_or_create_section(
            self.get_redshift_solver_section,
            ESolveMethod.TEMPLATE_COMBINATION.value,
            create,
            spectrum_model,
        )

    def get_template_combination_ism(self, spectrum_model: str) -> bool:
        return self._get_on_None(self.get_template_combination_section(spectrum_model), "ismFit")

    def get_template_combination_spectrum_component(self, spectrum_model: str) -> str:
        return self._get_on_None(
            self._get_on_None(self.get_template_combination_section(spectrum_model), "spectrum"), "component"
        )

    def get_ebmv_section(self) -> Optional[dict]:
        return self.parameters.get("ebmv")

    def get_ebmv_count(self):
        return self._get_on_None(self.get_ebmv_section(), "count")

    def get_ebmv_step(self):
        return self._get_on_None(self.get_ebmv_section(), "step")

    def get_ebmv_start(self):
        return self._get_on_None(self.get_ebmv_section(), "start")

    def get_linemodel_solve_section(self, spectrum_model: str, create: bool = False) -> dict:
        return self._get_or_create_section(
            self.get_redshift_solver_section,
            ESolveMethod.LINE_MODEL.value,
            create,
            spectrum_model,
        )

    def get_linemodel_solve_linemodel_section(self, spectrum_model: str, create: bool = False) -> dict:
        return self._get_or_create_section(
            self.get_linemodel_solve_section, "lineModel", create, spectrum_model
        )

    def get_linemodel_line_width_type(self, spectrum_model: str):
        return self._get_on_None(self.get_linemodel_solve_linemodel_section(spectrum_model), "lineWidthType")

    def get_linemodel_em_velocity_fit_min(self, spectrum_model: str):
        return self.get_velocity_fit_param(
            spectrum_model, ESolveMethod.LINE_MODEL, EVelocityType.Emission, EVelocityFitParam.Min
        )

    def get_line_model_photometry(self, spectrum_model: str) -> bool:
        return self._get_on_None(
            self.get_linemodel_solve_linemodel_section(spectrum_model), "enablePhotometry"
        )

    def get_solve_method_igm_fit(self, spectrum_model: str, solve_method: ESolveMethod) -> bool:
        igmfit = None
        if solve_method == ESolveMethod.LINE_MODEL:
            igmfit = self.get_linemodel_lya_profile(spectrum_model) == "igm"
        elif solve_method == ESolveMethod.TEMPLATE_FITTING:
            igmfit = self.get_template_fitting_igmfit(spectrum_model)
        elif solve_method == ESolveMethod.TEMPLATE_COMBINATION:
            igmfit = self.get_template_combination_igmfit(spectrum_model)
        elif solve_method == ESolveMethod.LINE_MEAS:
            igmfit = self.get_linemeas_lya_profile(spectrum_model) == "igm"
        return igmfit

    def get_linemodel_lya_section(self, spectrum_model: str) -> dict:
        return self._get_on_None(self.get_linemodel_solve_linemodel_section(spectrum_model), "lya")

    def get_linemodel_lya_profile(self, spectrum_model: str) -> str:
        return self._get_on_None(self.get_linemodel_lya_section(spectrum_model), "profile")

    def get_linemodel_lya_asym_section(self, spectrum_model: str) -> str:
        return self._get_on_None(self.get_linemodel_lya_section(spectrum_model), "asymProfile")

    def get_linemeas_line_width_type(self, spectrum_model: str):
        return self._get_on_None(self.get_linemeas_linemodel_section(spectrum_model), "lineWidthType")

    def get_linemeas_em_velocity_fit_min(self, spectrum_model: str):
        return self.get_velocity_fit_param(
            spectrum_model, ESolveMethod.LINE_MEAS, EVelocityType.Emission, EVelocityFitParam.Min
        )

    def get_linemeas_lya_section(self, spectrum_model: str) -> dict:
        return self._get_on_None(self.get_linemeas_linemodel_section(spectrum_model), "lya")

    def get_linemeas_lya_profile(self, spectrum_model: str) -> str:
        return self._get_on_None(self.get_linemeas_lya_section(spectrum_model), "profile")

    def get_linemeas_lya_asym_section(self, spectrum_model: str) -> str:
        return self._get_on_None(self.get_linemeas_lya_section(spectrum_model), "asymProfile")

    def get_template_combination_igmfit(self, spectrum_model: str) -> bool:
        return self._get_on_None(self.get_template_combination_section(spectrum_model), "igmFit")

    def get_template_fitting_igmfit(self, spectrum_model: str) -> bool:
        return self._get_on_None(self.get_template_fitting_section(spectrum_model), "igmFit")

    def get_linemodel_line_ratio_type(self, spectrum_model: str) -> str:
        return self._get_on_None(self.get_linemodel_solve_linemodel_section(spectrum_model), "lineRatioType")

    def get_linemodel_rules(self, spectrum_model: str) -> str:
        return self._get_on_None(self.get_linemodel_solve_linemodel_section(spectrum_model), "rules")

    def get_linemodel_tplratio_catalog(self, spectrum_model: str) -> str:
        return self._get_on_None(
            self.get_linemodel_solve_linemodel_section(spectrum_model), "tplRatioCatalog"
        )

    def get_linemodel_tplratio_ismfit(self, spectrum_model: str) -> str:
        return self._get_on_None(self.get_linemodel_solve_linemodel_section(spectrum_model), "tplRatioIsmFit")

    def get_linemodel_continuum_component(self, spectrum_model: str) -> str:
        return self._get_on_None(
            self.get_linemodel_solve_linemodel_section(spectrum_model), "continuumComponent"
        )

    def get_linemodel_continuumfit_section(self, spectrum_model: str) -> dict:
        return self._get_on_None(self.get_linemodel_solve_linemodel_section(spectrum_model), "continuumFit")

    def get_linemodel_continuumfit_fft(self, spectrum_model: str) -> dict:
        return self._get_on_None(self.get_linemodel_continuumfit_section(spectrum_model), "fftProcessing")

    def get_firstpass_section(self, solve_method: ESolveMethod, spectrum_model: str) -> Optional[dict]:
        solve_section: dict
        if solve_method == ESolveMethod.LINE_MODEL:
            solve_section = self.get_linemodel_solve_linemodel_section(spectrum_model)
        elif solve_method == ESolveMethod.TEMPLATE_FITTING:
            solve_section = self.get_template_fitting_section(spectrum_model)
        else:
            return None
        return self._get_on_None(solve_section, "firstPass")

    def get_secondpass_section(self, solve_method: ESolveMethod, spectrum_model: str) -> Optional[dict]:
        if solve_method == ESolveMethod.LINE_MODEL:
            solve_section = self.get_linemodel_solve_linemodel_section(spectrum_model)
        elif solve_method == ESolveMethod.TEMPLATE_FITTING:
            solve_section = self.get_template_fitting_section(spectrum_model)
        else:
            return None
        return self._get_on_None(solve_section, "secondPass")

    def get_secondpass_continuumfit(self, solve_method: ESolveMethod, spectrum_model: str) -> dict:
        secondpass_section = self.get_secondpass_section(solve_method, spectrum_model)
        return self._get_on_None(secondpass_section, "continuumFit")

    def get_linemodel_continuum_reestimation(self, spectrum_model: str) -> str:
        return self._get_on_None(
            self.get_linemodel_solve_linemodel_section(spectrum_model), "continuumReestimation"
        )

    def get_linemodel_useloglambdasampling(self, spectrum_model: str) -> bool:
        return self._get_on_None(
            self.get_linemodel_solve_linemodel_section(spectrum_model), "useLogLambdaSampling"
        )

    def get_linemodel_continuumfit_ismfit(self, spectrum_model: str) -> bool:
        return self._get_on_None(self.get_linemodel_continuumfit_section(spectrum_model), "ismFit")

    def get_linemodel_continuumfit_igmfit(self, spectrum_model: str) -> bool:
        return self._get_on_None(self.get_linemodel_continuumfit_section(spectrum_model), "igmFit")

    def get_linemodel_continuumfit_fftprocessing(self, spectrum_model: str) -> bool:
        return self._get_on_None(self.get_linemodel_continuumfit_section(spectrum_model), "fftProcessing")

    def get_linemodel_firstpass_tplratio_ismfit(self, spectrum_model: str) -> bool:
        return self._get_on_None(
            self.get_firstpass_section(ESolveMethod.LINE_MODEL, spectrum_model), "tplRatioIsmFit"
        )

    def get_linemodel_firstpass_extremacount(self, spectrum_model: str) -> bool:
        return self._get_on_None(
            self.get_firstpass_section(ESolveMethod.LINE_MODEL, spectrum_model), "extremaCount"
        )

    def get_linemodel_fitting_method(self, spectrum_model: str) -> str:
        return self._get_on_None(self.get_linemodel_solve_linemodel_section(spectrum_model), "fittingMethod")

    def get_skipsecondpass(
        self, solve_method: ESolveMethod, spectrum_model: str, default: Optional[bool] = None
    ) -> bool:
        section = None
        if solve_method == ESolveMethod.LINE_MODEL:
            section = self.get_linemodel_solve_linemodel_section(spectrum_model)
        elif solve_method == ESolveMethod.TEMPLATE_FITTING:
            section = self.get_template_fitting_section(spectrum_model)
        return self._get_on_None(section, "skipSecondPass", default)

    def get_template_fitting_single_pass(self, spectrum_model: str) -> bool:
        section = self.get_template_fitting_section(spectrum_model)
        return self._get_on_None(section, "singlePass")

    def get_linemodel_nsigmasupport(self, spectrum_model: str) -> float:
        return self._get_on_None(self.get_linemodel_solve_linemodel_section(spectrum_model), "nSigmaSupport")

    def get_linemodel_improve_balmer_fit(self, spectrum_model: str) -> float:
        return self._get_on_None(
            self.get_linemodel_solve_linemodel_section(spectrum_model), "improveBalmerFit"
        )

    def get_linemodel_extremacount(self, spectrum_model: str):
        return self._get_on_None(self.get_linemodel_solve_linemodel_section(spectrum_model), "extremaCount")

    def get_linemodel_velocity_fit(self, spectrum_model: str) -> bool:
        return self._get_on_None(self.get_linemodel_solve_linemodel_section(spectrum_model), "velocityFit")

    @staticmethod
    def get_velocity_name(velocity_type: EVelocityType) -> str:
        return f"velocity{velocity_type.value}"

    def get_velocity(
        self, spectrum_model: str, solve_method: ESolveMethod, velocity_type: EVelocityType
    ) -> float:
        return self._get_on_None(
            self.get_linemodel_section(spectrum_model, solve_method),
            self.get_velocity_name(velocity_type),
        )

    @classmethod
    def get_velocity_fit_param_name(cls, velocity_type: EVelocityType, param: EVelocityFitParam) -> str:
        return f"{cls.velocity_fit_prefix_dict[velocity_type]}VelocityFit{param.value}"

    def get_velocity_fit_param(
        self,
        spectrum_model: str,
        solve_method: ESolveMethod,
        velocity_type: EVelocityType,
        param: EVelocityFitParam,
    ) -> float:
        return self._get_on_None(
            self.get_linemodel_section(spectrum_model, solve_method),
            self.get_velocity_fit_param_name(velocity_type, param),
        )

    def set_velocity(
        self, spectrum_model: str, solve_method: ESolveMethod, velocity_type: EVelocityType, value: float
    ) -> None:
        linemodel_section = self.get_linemodel_section(spectrum_model, solve_method)
        if linemodel_section is None:
            raise APIException(
                ErrorCode.INTERNAL_ERROR,
                f"Missing line model section for spectrum model {spectrum_model} and solve method {solve_method}",
            )
        linemodel_section[self.get_velocity_name(velocity_type)] = value

    def set_velocity_fit_param(
        self,
        spectrum_model: str,
        solve_method: ESolveMethod,
        velocity_type: EVelocityType,
        param: EVelocityFitParam,
        value: float,
    ) -> None:
        linemodel_section = self.get_linemodel_section(spectrum_model, solve_method)
        if linemodel_section is None:
            raise APIException(
                ErrorCode.INTERNAL_ERROR,
                f"Missing line model section for spectrum model {spectrum_model} and solve method {solve_method}",
            )
        linemodel_section[self.get_velocity_fit_param_name(velocity_type, param)] = value

    def get_linemeas_linemodel_section(self, spectrum_model: str) -> dict:
        return self._get_on_None(self.get_linemeas_solve_section(spectrum_model), "lineModel")

    def get_linemeas_line_ratio_type(self, spectrum_model: str) -> str:
        return self._get_on_None(self.get_linemeas_linemodel_section(spectrum_model), "lineRatioType")

    def get_linemeas_rules(self, spectrum_model: str) -> str:
        return self._get_on_None(self.get_linemeas_linemodel_section(spectrum_model), "rules")

    def get_linemeas_fitting_method(self, spectrum_model: str) -> str:
        return self._get_on_None(self.get_linemeas_linemodel_section(spectrum_model), "fittingMethod")

    def get_linemeas_velocity_fit(self, spectrum_model: str) -> bool:
        return self._get_on_None(self.get_linemeas_linemodel_section(spectrum_model), "velocityFit")

    def get_linemeas_nsigmasupport(self, spectrum_model: str) -> float:
        return self._get_on_None(self.get_linemeas_linemodel_section(spectrum_model), "nSigmaSupport")

    def get_nsigmasupport(self, spectrum_model: str, method: ESolveMethod) -> Optional[float]:
        nsigmasupport = None
        if method == ESolveMethod.LINE_MODEL:
            nsigmasupport = self.get_linemodel_nsigmasupport(spectrum_model)
        elif method == ESolveMethod.LINE_MEAS:
            nsigmasupport = self.get_linemeas_nsigmasupport(spectrum_model)
        return nsigmasupport

    def get_linemodel_linecatalog(self, spectrum_model: str) -> str:
        return self._get_on_None(self.get_linemodel_solve_linemodel_section(spectrum_model), "lineCatalog")

    def get_linemeas_linecatalog(self, spectrum_model: str) -> str:
        return self._get_on_None(self.get_linemeas_linemodel_section(spectrum_model), "lineCatalog")

    def get_linecatalog(self, spectrum_model: str, method: ESolveMethod) -> str:
        linecatalog = ""
        if method == ESolveMethod.LINE_MODEL:
            linecatalog = self.get_linemodel_linecatalog(spectrum_model)
        elif method == ESolveMethod.LINE_MEAS:
            linecatalog = self.get_linemeas_linecatalog(spectrum_model)
        return linecatalog

    def get_linemodel_section(self, spectrum_model: str, method: ESolveMethod) -> Optional[dict]:
        linemodel = None
        if method == ESolveMethod.LINE_MODEL:
            linemodel = self.get_linemodel_solve_linemodel_section(spectrum_model)
        elif method == ESolveMethod.LINE_MEAS:
            linemodel = self.get_linemeas_linemodel_section(spectrum_model)
        return linemodel

    def get_redshift_sampling(self, spectrum_model):
        return self.get_spectrum_model_section(spectrum_model).get("redshiftSampling")

    def get_observation_ids(self):
        try:
            return list(self.parameters["lambdaRange"].keys())
        except Exception:
            return [""]

    def get_nb_samples_min(self):
        return self.parameters["nbSamplesMin"]

    def _get_or_create_section(self, parent_section_getter, child_section_name: str, create: bool, *args):
        parent_section = parent_section_getter(*args, create)
        child_section = self._get_on_None(parent_section, child_section_name)
        if create and child_section is None:
            if parent_section is None:
                raise APIException(
                    ErrorCode.INTERNAL_ERROR,
                    f"{parent_section_getter.__name__} must not return None when called with create = true",
                )
            child_section = {}
            parent_section[child_section_name] = child_section
        return child_section
