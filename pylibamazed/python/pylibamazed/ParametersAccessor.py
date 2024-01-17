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

from typing import List


class ParametersAccessor:
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

    def get_airvacuum_method(self):
        return self.parameters.get("airVacuumMethod", "")

    def get_photometry_transmission_dir(self):
        return self.parameters.get("photometryTransmissionDir")

    def get_photometry_bands(self) -> List[str]:
        return self.parameters.get("photometryBand", [])

    def get_multiobs_method(self):
        return self.parameters.get("multiObsMethod")

    def get_spectrum_models(self, default=None):
        return self.parameters.get("spectrumModels", default)

    def get_spectrum_model_section(self, spectrum_model) -> dict:
        return self.parameters.get(spectrum_model)

    def get_stages(self, spectrum_model: str) -> list:
        return self.get_spectrum_model_section(spectrum_model).get("stages")

    def get_redshift_solver_section(self, spectrum_model) -> dict:
        return self._get_on_None(self.get_spectrum_model_section(spectrum_model), "redshiftSolver")

    def get_redshift_solver_method(self, spectrum_model: str) -> str:
        method = self._get_on_None(self.get_redshift_solver_section(spectrum_model), "method")
        # Not in the get because if method is defined as an empty string in the parameters json
        # we still want method to be "None"
        if method == "":
            method = None
        return method

    def get_linemeas_solver_section(self, spectrum_model) -> dict:
        return self._get_on_None(self.get_spectrum_model_section(spectrum_model), "lineMeasSolver")

    def get_linemeas_method(self, spectrum_model: str) -> str:
        if "lineMeasSolver" not in self.get_stages(spectrum_model):
            return None
        return self._get_on_None(self.get_linemeas_solver_section(spectrum_model), "method")

    def get_linemeas_solve_section(self, spectrum_model: str) -> dict:
        return self._get_on_None(self.get_linemeas_solver_section(spectrum_model), "lineMeasSolve")

    def get_linemeas_dzhalf(self, spectrum_model: str) -> float:
        return self.get_spectrum_model_section(spectrum_model).get("lineMeasDzHalf")

    def get_redshiftrange(self, spectrum_model: str) -> List[float]:
        return self.get_spectrum_model_section(spectrum_model).get("redshiftRange")

    def get_redshiftstep(self, spectrum_model: str) -> float:
        return self.get_spectrum_model_section(spectrum_model).get("redshiftStep")

    def get_linemeas_redshiftstep(self, spectrum_model: str) -> float:
        return self.get_spectrum_model_section(spectrum_model).get("lineMeasRedshiftStep")

    def get_reliability_enabled(self, spectrum_model: str) -> bool:
        return "reliabilitySolver" in self.get_stages(spectrum_model)

    def get_reliability_section(self, spectrum_model: str) -> dict:
        return self.get_spectrum_model_section(spectrum_model).get("reliabilitySolver")

    def get_reliability_method(self, spectrum_model: str) -> str:
        return self._get_on_None(self.get_reliability_section(spectrum_model), "method")

    def get_deep_learning_solver_section(self, spectrum_model: str) -> str:
        return self._get_on_None(self.get_reliability_section(spectrum_model), "deepLearningSolver")

    def get_reliability_model(self, spectrum_model: str) -> str:
        return self._get_on_None(self.get_deep_learning_solver_section(spectrum_model), "reliabilityModel")

    def get_template_dir(self, spectrum_model: str) -> str:
        return self.get_spectrum_model_section(spectrum_model).get("templateDir")

    def get_redshift_solver_method_section(self, spectrum_model: str) -> dict:
        method = self.get_redshift_solver_method(spectrum_model)

        return self._get_on_None(self.get_redshift_solver_section(spectrum_model), method)

    def photometry_is_enabled(self):
        for spectrum_model in self.get_spectrum_models([]):
            method = self.get_redshift_solver_method(spectrum_model)
            method_section = self.get_redshift_solver_method_section(spectrum_model)
            if method == "lineModelSolve":
                method_section = self._get_on_None(method_section, "lineModel")
            if self._get_on_None(method_section, "enablePhotometry", False):
                return True
        return False

    def get_additional_cols(self, default=None) -> List[str]:
        return self.parameters.get("additionalCols") or default

    def get_filters(self, default=None):
        return self.parameters.get("filters", default)

    def get_lsf(self) -> dict:
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

    def get_continuum_removal(self, nesting: str = None):
        dict_to_search_on = self.parameters
        if nesting:
            dict_to_search_on = dict_to_search_on.get(nesting)
        continuum_removal = self._get_on_None(dict_to_search_on, "continuumRemoval")
        return continuum_removal

    def get_continuum_removal_method(self, nesting: str = None):
        return self._get_on_None(self.get_continuum_removal(nesting), "method")

    def get_continuum_removal_median_kernel_width(self, nesting: str = None):
        return self._get_on_None(self.get_continuum_removal(nesting), "medianKernelWidth")

    def get_continuum_median_kernel_reflection(self, nesting: str = None):
        return self._get_on_None(self.get_continuum_removal(nesting), "medianEvenReflection")

    def _get_on_None(self, dict, key, default=None):
        if dict is None:
            return default
        else:
            return dict.get(key, default)

    def get_template_fitting_section(self, spectrum_model: str) -> dict:
        return self._get_on_None(self.get_redshift_solver_section(spectrum_model), "templateFittingSolve")

    def get_template_fitting_ism(self, spectrum_model: str) -> bool:
        return self._get_on_None(self.get_template_fitting_section(spectrum_model), "ismFit")

    def get_template_fitting_photometry_enabled(self, spectrum_model: str) -> bool:
        return self._get_on_None(self.get_template_fitting_section(spectrum_model), "enablePhotometry")

    def get_template_fitting_photometry_section(self, spectrum_model: str) -> dict:
        return self._get_on_None(self.get_template_fitting_section(spectrum_model), "photometry")

    def get_template_fitting_photometry_weight(self, spectrum_model: str) -> dict:
        return self._get_on_None(self.get_template_fitting_photometry_section(spectrum_model), "weight")

    def get_template_combination_section(self, spectrum_model: str) -> dict:
        return self._get_on_None(self.get_redshift_solver_section(spectrum_model), "tplCombinationSolve")

    def get_template_combination_ism(self, spectrum_model: str) -> bool:
        return self._get_on_None(self.get_template_combination_section(spectrum_model), "ismFit")

    def get_ebmv_section(self) -> dict:
        return self.parameters.get("ebmv")

    def get_ebmv_count(self):
        return self._get_on_None(self.get_ebmv_section(), "count")

    def get_ebmv_step(self):
        return self._get_on_None(self.get_ebmv_section(), "step")

    def get_ebmv_start(self):
        return self._get_on_None(self.get_ebmv_section(), "start")

    def get_linemodel_solver_section(self, spectrum_model: str) -> dict:
        return self._get_on_None(self.get_redshift_solver_section(spectrum_model), "lineModelSolve")

    def get_linemodel_solver_linemodel_section(self, spectrum_model: str) -> str:
        return self._get_on_None(self.get_linemodel_solver_section(spectrum_model), "lineModel")

    def get_solve_method_igm_fit(self, spectrum_model: str, solve_method: str) -> bool:
        igmfit = None
        if solve_method == "lineModelSolve":
            igmfit = self.get_linemodel_igmfit(spectrum_model)
        elif solve_method == "templateFittingSolve":
            igmfit = self.get_template_fitting_igmfit(spectrum_model)
        elif solve_method == "tplCombinationSolve":
            igmfit = self.get_template_combination_igmfit(spectrum_model)
        return igmfit

    def get_linemodel_igmfit(self, spectrum_model: str) -> bool:
        return self._get_on_None(self.get_linemodel_solver_linemodel_section(spectrum_model), "igmFit")

    def get_template_combination_igmfit(self, spectrum_model: str) -> bool:
        return self._get_on_None(self.get_template_combination_section(spectrum_model), "igmFit")

    def get_template_fitting_igmfit(self, spectrum_model: str) -> bool:
        return self._get_on_None(self.get_template_fitting_section(spectrum_model), "igmFit")

    def get_linemodel_line_ratio_type(self, spectrum_model: str) -> str:
        return self._get_on_None(self.get_linemodel_solver_linemodel_section(spectrum_model), "lineRatioType")

    def get_linemodel_rules(self, spectrum_model: str) -> str:
        return self._get_on_None(self.get_linemodel_solver_linemodel_section(spectrum_model), "rules")

    def get_linemodel_tplratio_catalog(self, spectrum_model: str) -> str:
        return self._get_on_None(
            self.get_linemodel_solver_linemodel_section(spectrum_model),
            "tplRatioCatalog")

    def get_linemodel_tplratio_ismfit(self, spectrum_model: str) -> str:
        return self._get_on_None(
            self.get_linemodel_solver_linemodel_section(spectrum_model),
            "tplRatioIsmFit")

    def get_linemodel_continuum_component(self, spectrum_model: str) -> str:
        return self._get_on_None(
            self.get_linemodel_solver_linemodel_section(spectrum_model),
            "continuumComponent")

    def get_linemodel_continuumfit_section(self, spectrum_model: str) -> dict:
        return self._get_on_None(self.get_linemodel_solver_linemodel_section(spectrum_model), "continuumFit")

    def get_linemodel_secondpass_section(self, spectrum_model: str) -> dict:
        return self._get_on_None(self.get_linemodel_solver_linemodel_section(spectrum_model), "secondPass")

    def get_linemodel_secondpass_continuumfit(self, spectrum_model: str) -> str:
        return self._get_on_None(self.get_linemodel_secondpass_section(spectrum_model), "continuumFit")

    def get_linemodel_continuum_reestimation(self, spectrum_model: str) -> str:
        return self._get_on_None(
            self.get_linemodel_solver_linemodel_section(spectrum_model),
            "continuumReestimation"
        )

    def get_linemodel_useloglambdasampling(self, spectrum_model: str) -> bool:
        return self._get_on_None(
            self.get_linemodel_solver_linemodel_section(spectrum_model),
            "useLogLambdaSampling"
        )

    def get_linemodel_continuumfit_ismfit(self, spectrum_model: str) -> bool:
        return self._get_on_None(self.get_linemodel_continuumfit_section(spectrum_model), "ismFit")

    def get_linemodel_continuumfit_fftprocessing(self, spectrum_model: str) -> bool:
        return self._get_on_None(self.get_linemodel_continuumfit_section(spectrum_model), "fftProcessing")

    def get_linemodel_firstpass_section(self, spectrum_model: str) -> bool:
        return self._get_on_None(self.get_linemodel_solver_linemodel_section(spectrum_model), "firstPass")

    def get_linemodel_firstpass_tplratio_ismfit(self, spectrum_model: str) -> bool:
        return self._get_on_None(self.get_linemodel_firstpass_section(spectrum_model), "tplRatioIsmFit")

    def get_linemodel_firstpass_extremacount(self, spectrum_model: str) -> bool:
        return self._get_on_None(self.get_linemodel_firstpass_section(spectrum_model), "extremaCount")

    def get_linemodel_skipsecondpass(self, spectrum_model: str) -> bool:
        return self._get_on_None(
            self.get_linemodel_solver_linemodel_section(spectrum_model), "skipSecondPass")

    def get_linemodel_nsigmasupport(self, spectrum_model: str) -> float:
        return self._get_on_None(self.get_linemodel_solver_linemodel_section(spectrum_model), "nSigmaSupport")

    def get_linemodel_improve_balmer_fit(self, spectrum_model: str) -> float:
        return self._get_on_None(
            self.get_linemodel_solver_linemodel_section(spectrum_model),
            "improveBalmerFit")

    def get_linemodel_extremacount(self, spectrum_model: str):
        return self._get_on_None(self.get_linemodel_solver_linemodel_section(spectrum_model), "extremaCount")

    def get_linemeas_linemodel_section(self, spectrum_model: str) -> dict:
        return self._get_on_None(self.get_linemeas_solve_section(spectrum_model), "lineModel")

    def get_linemeas_line_ratio_type(self, spectrum_model: str) -> str:
        return self._get_on_None(self.get_linemeas_linemodel_section(spectrum_model), "lineRatioType")

    def get_linemeas_rules(self, spectrum_model: str) -> str:
        return self._get_on_None(self.get_linemeas_linemodel_section(spectrum_model), "rules")

    def get_linemeas_fitting_method(self, spectrum_model: str) -> str:
        return self._get_on_None(self.get_linemeas_linemodel_section(spectrum_model), "fittingMethod")

    def get_linemeas_velocity_fit(self, spectrum_model: str) -> str:
        return self._get_on_None(self.get_linemeas_linemodel_section(spectrum_model), "velocityFit")

    def get_linemeas_velocity_fit_param(self, spectrum_model: str, param: str) -> float:
        return self._get_on_None(self.get_linemeas_linemodel_section(spectrum_model), param)

    def get_linemeas_nsigmasupport(self, spectrum_model: str) -> float:
        return self._get_on_None(self.get_linemeas_linemodel_section(spectrum_model), "nSigmaSupport")

    def get_nsigmasupport(self, spectrum_model: str, method: str) -> float:
        nsigmasupport = None
        if method == "lineModelSolve":
            nsigmasupport = self.get_linemodel_nsigmasupport(spectrum_model)
        elif method == "lineMeasSolve":
            nsigmasupport = self.get_linemeas_nsigmasupport(spectrum_model)
        return nsigmasupport

    def get_linemodel_linecatalog(self, spectrum_model: str) -> str:
        # print("parameters", json.dumps(self.get_spectrum_model_section(spectrum_model), indent=4))
        return self._get_on_None(self.get_linemodel_solver_linemodel_section(spectrum_model), "lineCatalog")

    def get_linemeas_linecatalog(self, spectrum_model: str) -> str:
        return self._get_on_None(self.get_linemeas_linemodel_section(spectrum_model), "lineCatalog")

    def get_linecatalog(self, spectrum_model: str, method: str) -> str:
        linecatalog = None
        # print("\n\nmethod", method)
        if method == "lineModelSolve":
            linecatalog = self.get_linemodel_linecatalog(spectrum_model)
        elif method == "lineMeasSolve":
            linecatalog = self.get_linemeas_linecatalog(spectrum_model)
        return linecatalog

    def get_linemodel_section(self, spectrum_model, method) -> dict:
        linemodel = None
        if method == "lineModelSolve":
            linemodel = self.get_linemodel_solver_linemodel_section(spectrum_model)
        elif method == "lineMeasSolve":
            linemodel = self.get_linemeas_linemodel_section(spectrum_model)
        return linemodel

    def get_redshift_sampling(self, spectrum_model):
        return self.get_spectrum_model_section(spectrum_model).get("redshiftSampling")

    def get_observation_ids(self):
        try:
            return list(self.parameters["lambdaRange"].keys())
        except Exception:
            return [""]
