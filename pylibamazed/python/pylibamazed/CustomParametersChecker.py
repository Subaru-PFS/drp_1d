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

from typing import Optional

from pylibamazed.Exception import APIException
from pylibamazed.FilterLoader import ParamJsonFilterLoader
from pylibamazed.ParametersAccessor import ParametersAccessor, EVelocityFitParam, EVelocityType
from pylibamazed.ParametersChecker import ParametersChecker
from pylibamazed.redshift import ErrorCode, WarningCode, CFlagWarning
from pylibamazed.ParametersAccessor import ESolveMethod
from pylibamazed.ParametersAccessor import EContinuumFit

zflag = CFlagWarning.GetInstance()


class CustomParametersChecker(ParametersChecker):
    def __init__(
        self,
        parameters: dict,
        FilterLoader=ParamJsonFilterLoader,
        Accessor=ParametersAccessor,
    ):
        self.filter_loader = FilterLoader()
        self.accessor = Accessor(parameters)

    def check(self):
        self._check_general()
        self._check_lsf()
        self._check_continuum_removal()
        self._check_templateCatalog_continuum_removal()
        self._check_linemeas_runmode()
        for object in self.accessor.get_spectrum_models([]):
            self._check_object(object)

    def _check_general(self):
        self._check_photometry_transmission_dir()
        self._check_photometry_band()
        for obs_id in self.accessor.get_observation_ids():
            self._check_filters(obs_id)

    def _check_photometry_transmission_dir(self):
        parameter_name = "photometryTransmissionDir"
        self._check_dependant_condition(
            self.accessor.photometry_is_enabled(),
            self.accessor.get_photometry_transmission_dir() is not None,
            parameter_name,
            parameter_name,
        )

    def _check_photometry_band(self):
        parameter_name = "photometryBand"
        self._check_dependant_condition(
            self.accessor.photometry_is_enabled(),
            len(self.accessor.get_photometry_bands()) > 0,
            parameter_name,
            parameter_name,
        )

    def _check_filters(self, obs_id: str):
        filters = self.accessor.get_filters(default=[], obs_id=obs_id)
        self._check_filters_format(filters)

        DEFAULT_COLUMN_NAMES = ["waves", "fluxes", "errors"]
        if not filters:
            return
        filter_keys = [filt["key"] for filt in filters]
        authorized_cols_names = DEFAULT_COLUMN_NAMES + self.accessor.get_additional_cols(default=[])

        for filter_name in filter_keys:
            if filter_name not in authorized_cols_names:
                raise APIException(ErrorCode.INVALID_PARAMETER_FILE, f"Unknown filter key {filter_name}")

    def _check_filters_format(self, json: list) -> None:
        if type(json) is not list:
            raise APIException(ErrorCode.INVALID_PARAMETER_FILE, "Input filters json must be a list")

        for filt in json:
            json_keys = self.filter_loader.keys
            different_keys = set(filt.keys()) != set(json_keys)
            different_length = len(filt.keys()) != len(json_keys)
            if different_keys or different_length:
                raise APIException(
                    ErrorCode.INVALID_PARAMETER_FILE,
                    f"Filters: each dictionary in json list must have exactly the following keys {json_keys}",
                )

    def _check_lsf(self) -> None:
        self._check_lsf_section()
        self._check_GaussianConstantWidth_lsf_type()
        self._check_GaussianNISPSIM_lsf_type()
        self._check_GaussianConstantResolution_lsf_type()
        self._check_GaussianVariablewidth_FileName()

    def _check_lsf_section(self) -> None:
        lsf_needed = False
        for spectrum_model in self.accessor.get_spectrum_models([]):
            if self.accessor.get_redshift_solver_method(spectrum_model) == ESolveMethod.LINE_MODEL or (
                "lineMeasSolver" in self.accessor.get_stages(spectrum_model)
            ):
                lsf_needed = True
        self._check_dependant_condition(lsf_needed, self.accessor.get_lsf() is not None, "lsf", "lsf")

    def _check_GaussianNISPSIM_lsf_type(self):
        self._check_dependant_condition(
            self.accessor.get_lsf_type() == "GaussianNISPSIM201707",
            self.accessor.get_lsf_sourcesize() is not None,
            "lsf sourceSize",
            "lsf sourceSize",
        )

    def _check_GaussianConstantWidth_lsf_type(self):
        self._check_dependant_condition(
            self.accessor.get_lsf_type() == "gaussianConstantWidth",
            self.accessor.get_lsf_width() is not None,
            "lsf width",
            "lsf width",
        )

    def _check_GaussianConstantResolution_lsf_type(self):
        self._check_dependant_condition(
            self.accessor.get_lsf_type() == "gaussianConstantResolution",
            self.accessor.get_lsf_resolution() is not None,
            "lsf resolution",
            "lsf resolution",
        )

    def _check_GaussianVariablewidth_FileName(self):
        self._check_dependant_condition(
            self.accessor.get_lsf_type() == "gaussianVariableWidth",
            self.accessor.get_lsf_width_file_name() is not None,
            "lsf gaussianVariableWidthFileName",
            "lsf gaussianVariableWidthFileName",
        )

    def _check_continuum_removal(self, fromTemplateCatalog=False) -> None:
        self._check_continuum_removal_section_presence()
        self._check_IrregularSamplingMedian_kernel_width(fromTemplateCatalog)
        self._check_IrregularSamplingMedian_kernel_reflection(fromTemplateCatalog)

    def _check_templateCatalog_continuum_removal(self) -> None:
        self._check_templateCatalog_continuum_removal_section_presence()
        self._check_IrregularSamplingMedian_kernel_width("templateCatalog")
        self._check_IrregularSamplingMedian_kernel_reflection("templateCatalog")

    def _check_templateCatalog_continuum_removal_section_presence(self):
        continuum_removal_necessity = self.template_catalog_continuum_removal_presence_condition()
        self._check_dependant_condition(
            continuum_removal_necessity,
            self.accessor.get_continuum_removal_section("templateCatalog") is not None,
            "templateCatalog continuumRemoval",
            "templateCatalog continuumRemoval",
        )

    def template_catalog_continuum_removal_presence_condition(self) -> bool:
        continuum_removal_necessity = False
        for spectrum_model in self.accessor.get_spectrum_models([]):
            if self.accessor.get_template_fitting_spectrum_component(spectrum_model) not in [
                None,
                "raw",
            ] or self.accessor.get_template_combination_spectrum_component(spectrum_model) not in [
                None,
                "raw",
            ]:
                continuum_removal_necessity = True
        return continuum_removal_necessity

    def continuum_removal_presence_condition(self) -> bool:
        continuum_removal_necessity = False
        for spectrum_model in self.accessor.get_spectrum_models([]):
            template_fitting_spectrum_component = self.accessor.get_template_fitting_spectrum_component(
                spectrum_model
            )
            template_combination_spectrum_component = (
                self.accessor.get_template_combination_spectrum_component(spectrum_model)
            )
            linemodel_continuum_component = self.accessor.get_linemodel_continuum_component(spectrum_model)
            linemodel_continuum_reestimation = self.accessor.get_linemodel_continuum_reestimation(
                spectrum_model
            )

            if (
                template_fitting_spectrum_component not in [None, "raw"]
                or template_combination_spectrum_component not in [None, "raw"]
                or linemodel_continuum_component in ["fromSpectrum", "tplFitAuto", "powerLawAuto"]
                or linemodel_continuum_reestimation not in [None, "never"]
            ):
                continuum_removal_necessity = True
        return continuum_removal_necessity

    def _check_continuum_removal_section_presence(self):
        continuum_removal_necessity = self.continuum_removal_presence_condition()
        self._check_dependant_condition(
            continuum_removal_necessity,
            self.accessor.get_continuum_removal_section() is not None,
            "continuumRemoval",
            "continuumRemoval",
        )

    def median_kernel_presence_condition(self, fromTemplateCatalog):
        return self.accessor.get_continuum_removal_method(fromTemplateCatalog) == "irregularSamplingMedian"

    def _check_IrregularSamplingMedian_kernel_width(self, fromTemplateCatalog=False):
        self._check_dependant_condition(
            self.median_kernel_presence_condition(fromTemplateCatalog),
            self.accessor.get_continuum_removal_median_kernel_width(fromTemplateCatalog),
            "continuumRemoval medianKernelWidth",
            "continuumRemoval medianKernelWidth",
        )

    def _check_IrregularSamplingMedian_kernel_reflection(self, fromTemplateCatalog=False):
        self._check_dependant_condition(
            self.median_kernel_presence_condition(fromTemplateCatalog),
            self.accessor.get_continuum_median_kernel_reflection(fromTemplateCatalog),
            "continuumRemoval medianEvenReflection",
            "continuumRemoval medianEvenReflection",
        )

    def _check_object(self, spectrum_model: str) -> None:
        self._check_linemeassolve(spectrum_model)
        self._check_object_reliability(spectrum_model)
        self._check_reliability_method(spectrum_model)
        self._check_templateFittingSolve_section(spectrum_model)
        self._check_templateCombinationSolve_section(spectrum_model)
        self._check_lineModelSolve(spectrum_model)
        self._check_template_dir(spectrum_model)

    def _check_linemeassolve(self, spectrum_model: str) -> None:
        self._check_linemeassolver_section(spectrum_model)
        self._check_linemeassolve_dzhalf(spectrum_model)
        self._check_linemeassolve_redshiftstep(spectrum_model)
        self._check_linemeassolve_lineratiotype_rules(spectrum_model)
        self._check_linemeassolve_fittingmethod_lbfgsb_velocityfit(spectrum_model)
        self._check_linemeassolve_velocity_fit_params(spectrum_model)
        self._check_linemeassolve_lya_fit(spectrum_model)

    def _check_linemeassolver_section(self, spectrum_model: str) -> None:
        self._check_dependant_condition(
            "lineMeasSolver" in self.accessor.get_stages(spectrum_model),
            self.accessor.get_linemeas_solver_section(spectrum_model) is not None,
            f"{spectrum_model} lineMeasSolver",
            f"{spectrum_model} lineMeasSolver",
        )

    def _check_template_dir(self, spectrum_model: str) -> None:
        expected_by_method = self.accessor.get_redshift_solver_method(spectrum_model) in [
            ESolveMethod.TEMPLATE_FITTING,
            ESolveMethod.TEMPLATE_COMBINATION,
        ]
        expected_by_continuum = self.accessor.get_redshift_solver_method(
            spectrum_model
        ) == ESolveMethod.LINE_MODEL and self.accessor.get_linemodel_continuum_component(spectrum_model) in [
            "tplFit",
            "tplFitAuto",
        ]

        template_dir_presence_condition = expected_by_method or expected_by_continuum

        self._check_dependant_condition(
            template_dir_presence_condition,
            self.accessor.get_template_dir(spectrum_model) is not None,
            f"{spectrum_model} templateDir",
            f"{spectrum_model} templateDir",
        )

    def _check_linemeassolve_dzhalf(self, spectrum_model: str) -> None:
        self._check_dependant_condition(
            self.accessor.get_linemeas_method(spectrum_model) is not None,
            self.accessor.get_linemeas_dzhalf(spectrum_model) is not None,
            f"lineMeasDzHalf for object {spectrum_model}",
            f"object {spectrum_model} lineMeasDzHalf",
        )

    def _check_linemeassolve_redshiftstep(self, spectrum_model: str) -> None:
        self._check_dependant_condition(
            self.accessor.get_linemeas_method(spectrum_model) is not None,
            self.accessor.get_linemeas_redshiftstep(spectrum_model) is not None,
            f"lineMeasRedshiftStep for object {spectrum_model}",
            f"object {spectrum_model} lineMeasRedshiftStep",
        )

    def _check_object_reliability(self, spectrum_model: str) -> None:
        self._check_dependant_condition(
            self.accessor.get_reliability_enabled(spectrum_model),
            self.accessor.get_reliability_section(spectrum_model) is not None,
            f"{spectrum_model} reliabilitySolver",
            f"{spectrum_model} reliabilitySolver",
        )

    def _check_reliability_method(self, spectrum_model: str) -> None:
        self._check_dependant_condition(
            "deepLearningSolver" in self.accessor.get_reliability_methods(spectrum_model),
            self.accessor.get_deep_learning_solver_section(spectrum_model) is not None,
            f"{spectrum_model} deepLearningSolver",
            f"{spectrum_model} reliabilitySolver deepLearningSolver",
        )
        self._check_dependant_condition(
            "skLearnClassifier" in self.accessor.get_reliability_methods(spectrum_model),
            self.accessor.get_sk_learn_classifier_solver_section(spectrum_model) is not None,
            f"{spectrum_model} skLearnClassifier",
            f"{spectrum_model} reliabilitySolver skLearnClassifier",
        )

    def _check_templateFittingSolve_section(self, spectrum_model: str) -> None:
        self._check_dependant_condition(
            self.accessor.get_redshift_solver_method(spectrum_model) == ESolveMethod.TEMPLATE_FITTING,
            self.accessor.get_template_fitting_section(spectrum_model) is not None,
            error_message=f"{ESolveMethod.TEMPLATE_FITTING.value} for object {spectrum_model}",
            warning_message=f"object {spectrum_model} {ESolveMethod.TEMPLATE_FITTING.value}",
        )
        self._check_templateFittingSolve_ism(spectrum_model)
        self._check_templateFittingSolve_photometry_weight(spectrum_model)
        self._check_templateFittingSolve_exclusive_fft_photometry(spectrum_model)
        self._check_templateFittingSolve_exclusive_twopass_photometry(spectrum_model)
        self._check_templateFittingSolve_exclusive_twopass_fft(spectrum_model)
        self._check_templateFittingSolve_pass_presence(spectrum_model)
        self._check_templateFittingSolve_secondpass_continuumfit(spectrum_model)
        self._check_templateFittingSolve_spectrum_consistency(spectrum_model)

    def _check_templateFittingSolve_ism(self, spectrum_model: str) -> None:
        self._check_dependant_condition(
            self.accessor.get_template_fitting_ism(spectrum_model),
            self.accessor.get_ebmv_section() is not None,
            "ebmv",
        )

    def _check_templateFittingSolve_photometry_weight(self, spectrum_model: str) -> None:
        self._check_dependant_condition(
            self.accessor.get_template_fitting_photometry_enabled(spectrum_model),
            self.accessor.get_template_fitting_photometry_weight(spectrum_model) is not None,
            f"object {spectrum_model} TemplateFittingSolve photometry weight",
            f"object {spectrum_model} TemplateFittingSolve photometry weight",
        )

    def _check_templateFittingSolve_exclusive_fft_photometry(self, spectrum_model: str) -> None:
        activateFft = self.accessor.get_template_fitting_fft(spectrum_model)
        activatePhotometry = self.accessor.get_template_fitting_photometry_enabled(spectrum_model)
        if activateFft and activatePhotometry:
            raise APIException(
                ErrorCode.INVALID_PARAMETER_FILE,
                "Template fitting: cannot activate both fft and photometry. Please deactivate "
                "templateFittingSolve.fftProcessing or templateFittingSolve.enablePhotometry "
                f"on spectrum model {spectrum_model}",
            )

    def _check_templateFittingSolve_exclusive_twopass_photometry(self, spectrum_model: str) -> None:
        activatePhotometry = self.accessor.get_template_fitting_photometry_enabled(spectrum_model)
        singlePass = self.accessor.get_template_fitting_single_pass(spectrum_model)
        if activatePhotometry and not singlePass:
            raise APIException(
                ErrorCode.INVALID_PARAMETER_FILE,
                "Template fitting: cannot activate both two pass and photometry. Please activate "
                "templateFittingSolve.singlePass or deactivate templateFittingSolve.enablePhotometry "
                f"on spectrum model {spectrum_model}",
            )

    def _check_templateFittingSolve_exclusive_twopass_fft(self, spectrum_model: str) -> None:
        activateFft = self.accessor.get_template_fitting_fft(spectrum_model)
        singlePass = self.accessor.get_template_fitting_single_pass(spectrum_model)
        if activateFft and not singlePass:
            raise APIException(
                ErrorCode.INVALID_PARAMETER_FILE,
                "Template fitting: cannot activate both two pass and fft. Please activate "
                "templateFittingSolve.singlePass or deactivate templateFittingSolve.fftProcessing "
                f"on spectrum model {spectrum_model}",
            )

    def _check_templateFittingSolve_pass_presence(self, spectrum_model: str) -> None:
        self._check_dependant_condition(
            self.accessor.get_template_fitting_single_pass(spectrum_model) is False,
            self.accessor.get_firstpass_section(ESolveMethod.TEMPLATE_FITTING, spectrum_model) is not None,
            f"object {spectrum_model} TemplateFittingSolve first pass section",
            f"object {spectrum_model} TemplateFittingSolve first pass section",
        )
        self._check_dependant_condition(
            self.accessor.get_template_fitting_single_pass(spectrum_model) is False,
            self.accessor.get_secondpass_section(ESolveMethod.TEMPLATE_FITTING, spectrum_model) is not None,
            f"object {spectrum_model} TemplateFittingSolve second pass section",
            f"object {spectrum_model} TemplateFittingSolve second pass section",
        )

    def _check_templateFittingSolve_secondpass_continuumfit(self, spectrum_model: str) -> None:
        continuum_fit = self.accessor.get_secondpass_continuumfit(
            ESolveMethod.TEMPLATE_FITTING, spectrum_model
        )
        if continuum_fit not in [
            None,
            EContinuumFit.REFIT_FIRST_PASS.value,
        ]:
            raise APIException(
                ErrorCode.INVALID_PARAMETER_FILE,
                "Second pass continuum fit for templatFittingSolve must be set to "
                f"{EContinuumFit.REFIT_FIRST_PASS.value}"
                f"(spectrum model {spectrum_model})",
            )

    def _check_templateFittingSolve_spectrum_consistency(self, spectrum_model: str):
        spectrum_component = self.accessor.get_template_fitting_spectrum_component(spectrum_model)
        self._check_dependant_condition(
            spectrum_component == "noContinuum",
            self.accessor.get_template_fitting_ism(spectrum_model) is False,
            error_message=f"noContinuum requires deactivating ism and igm, check {spectrum_model}",
            custom_error_message=True,
        )
        self._check_dependant_condition(
            spectrum_component == "noContinuum",
            self.accessor.get_template_fitting_igm(spectrum_model) is False,
            error_message=f"noContinuum requires deactivating ism and igm, check {spectrum_model}",
            custom_error_message=True,
        )

    def _check_templateCombinationSolve_section(self, spectrum_model: str) -> None:
        self._check_dependant_condition(
            self.accessor.get_redshift_solver_method(spectrum_model) == ESolveMethod.TEMPLATE_COMBINATION,
            self.accessor.get_template_combination_section(spectrum_model) is not None,
            error_message=f"{ESolveMethod.TEMPLATE_COMBINATION.value} for object {spectrum_model}",
            warning_message=f"object {spectrum_model} {ESolveMethod.TEMPLATE_COMBINATION.value}",
        )
        self._check_templateCombinationSolve_ism(spectrum_model)

    def _check_templateCombinationSolve_ism(self, spectrum_model: str) -> None:
        ism = self.accessor.get_template_combination_ism(spectrum_model)
        ebmv = self.accessor.get_ebmv_section()

        if ism and ebmv is None:
            self._raise_missing_param_error("ebmv")

    def _check_lineModelSolve(self, spectrum_model: str) -> None:
        self._check_linemodelsolve_section(spectrum_model)

        self._check_linemodelsolve_improveBalmerFit(spectrum_model)
        self._check_linemodelsolve_lineratiotype_rules(spectrum_model)

        self._check_lineratiotype_tplratio_catalog(spectrum_model)
        self._check_lineratiotype_tplratio_ismfit(spectrum_model)

        self._check_linemodelsolve_continuumreestimation(spectrum_model)

        self._check_lineModelSolve_exclusive_fft_photometry(spectrum_model)
        self._check_linemodelsolve_continuumfit_section(spectrum_model)
        self._check_linemodelsolve_continuumfit_ism(spectrum_model)
        self._check_linemodelsolve_secondpass_section(spectrum_model)
        self._check_linemodelsolve_secondpass_continuumfit(spectrum_model)
        self._check_linemodelsolve_firstpass_tplratio_ismfit(spectrum_model)
        self._check_linemodelsolve_firstpass_extremacount(spectrum_model)
        self._check_linemodelsolve_lya_fit(spectrum_model)
        self._check_linemodelsolve_useloglambdasampling(spectrum_model)
        self._check_linemodelsolve_velocity_fit_params(spectrum_model)

    def _check_linemodelsolve_section(self, spectrum_model: str):
        self._check_dependant_condition(
            self.accessor.get_redshift_solver_method(spectrum_model) == ESolveMethod.LINE_MODEL,
            self.accessor.get_linemodel_solve_section(spectrum_model) is not None,
            error_message=f"{ESolveMethod.LINE_MODEL.value} for object {spectrum_model}",
            warning_message=f"object {spectrum_model} {ESolveMethod.LINE_MODEL.value}",
        )

    def _check_linemodelsolve_lineratiotype_rules(self, spectrum_model: str):
        self._check_dependant_condition(
            self.accessor.get_linemodel_line_ratio_type(spectrum_model) == "rules",
            self.accessor.get_linemodel_rules(spectrum_model) is not None,
            error_message=f"lineModelSolve rules for object {spectrum_model}",
            warning_message=f"object {spectrum_model} lineModelSolve rules",
        )

    def _check_linemodelsolve_improveBalmerFit(self, spectrum_model: str):
        improveBalmerFit = self.accessor.get_linemodel_improve_balmer_fit(spectrum_model)
        lineRatioType = self.accessor.get_linemodel_line_ratio_type(spectrum_model)
        if improveBalmerFit and lineRatioType != "rules":
            zflag.warning(
                WarningCode.UNUSED_PARAMETER,
                f"object {spectrum_model} lineModelSolve lineRatioType must be rules to "
                "activate improveBalmerFit",
            )

    def _check_lineratiotype_tplratio_catalog(self, spectrum_model):
        self._check_dependant_condition(
            self.accessor.get_linemodel_line_ratio_type(spectrum_model) in ["tplRatio", "tplCorr"],
            self.accessor.get_linemodel_tplratio_catalog(spectrum_model) is not None,
            error_message=f"lineModelSolve tplRatioCatalog for object {spectrum_model}",
            warning_message=f"object {spectrum_model} lineModelSolve tplRatioCatalog",
        )

    def _check_lineratiotype_tplratio_ismfit(self, spectrum_model):
        self._check_dependant_condition(
            self.accessor.get_linemodel_line_ratio_type(spectrum_model) in ["tplRatio", "tplCorr"],
            self.accessor.get_linemodel_tplratio_ismfit(spectrum_model) is not None,
            error_message=f"lineModelSolve tplRatioIsmFit for object {spectrum_model}",
            warning_message=f"object {spectrum_model} lineModelSolve tplRatioIsmFit",
        )

    def linemodelsolve_continuumreestimation_presence_condition(self, spectrum_model):
        return self.accessor.get_linemodel_fitting_method(spectrum_model) == "hybrid"

    def _check_linemodelsolve_continuumreestimation(self, spectrum_model):
        self._check_dependant_condition(
            self.accessor.get_linemodel_continuum_component(spectrum_model) == "fromSpectrum",
            self.accessor.get_linemodel_continuum_reestimation(spectrum_model) is not None,
            error_message=f"object {spectrum_model} lineModelSolve continuumReestimation",
            warning_message=f"object {spectrum_model} lineModelSolve continuumReestimation",
        )

    def _check_lineModelSolve_exclusive_fft_photometry(self, spectrum_model: str) -> None:
        activateContinuumFft = self.accessor.get_linemodel_continuumfit_fft(spectrum_model)
        activatePhotometry = self.accessor.get_line_model_photometry(spectrum_model)
        if activateContinuumFft and activatePhotometry:
            raise APIException(
                ErrorCode.INVALID_PARAMETER_FILE,
                "Line model solve: cannot activate both fft and photometry. Please deactivate "
                "lineModelSolve.continuumFit.fftProcessing or lineModelSolve.enablePhotometry "
                f"on spectrum model {spectrum_model}",
            )

    def _check_linemodelsolve_continuumfit_section(self, spectrum_model):
        self._check_dependant_condition(
            self.accessor.get_linemodel_continuum_component(spectrum_model)
            in ["tplFit", "tplFitAuto", "powerLaw", "powerLawAuto"],
            self.accessor.get_linemodel_continuumfit_section(spectrum_model) is not None,
            error_message=f"object {spectrum_model} lineModelSolve continuumFit",
            warning_message=f"object {spectrum_model} continuumFit section",
        )

    def _check_linemodelsolve_secondpass_section(self, spectrum_model):
        self._check_dependant_condition(
            self.accessor.get_skipsecondpass(ESolveMethod.LINE_MODEL, spectrum_model) is False,
            self.accessor.get_secondpass_section(ESolveMethod.LINE_MODEL, spectrum_model) is not None,
            f"object {spectrum_model} lineModelSolve secondPass",
            f"object {spectrum_model} lineModelSolve secondPass",
        )

    def _check_linemodelsolve_secondpass_continuumfit(self, spectrum_model):
        condition = self.accessor.get_secondpass_section(
            ESolveMethod.LINE_MODEL, spectrum_model
        ) is not None and self.accessor.get_linemodel_continuum_component(spectrum_model) in [
            "tplFit",
            "tplFitAuto",
            "powerLaw",
            "powerLawAuto",
        ]
        self._check_dependant_condition(
            condition,
            self.accessor.get_secondpass_continuumfit(ESolveMethod.LINE_MODEL, spectrum_model) is not None,
            error_message=f"object {spectrum_model} lineModelSolve secondpass continuumFit",
        )

    def _check_linemodelsolve_continuumfit_ism(self, spectrum_model: str) -> None:
        self._check_dependant_condition(
            self.accessor.get_linemodel_continuumfit_ismfit(spectrum_model),
            self.accessor.get_ebmv_section() is not None,
            "ebmv",
        )

    def _check_linemodelsolve_firstpass_tplratio_ismfit(self, spectrum_model: str):
        self._check_dependant_condition(
            self.accessor.get_linemodel_line_ratio_type(spectrum_model) in ["tplRatio", "tplCorr"],
            self.accessor.get_linemodel_firstpass_tplratio_ismfit(spectrum_model) is not None,
            f"object {spectrum_model} lineModelSolve firstpass tplRatioIsmFit",
            f"object {spectrum_model} lineModelSolve firstpass tplRatioIsmFit",
        )

    def _check_linemodelsolve_firstpass_extremacount(self, spectrum_model: str):
        lm_extremacount = self.accessor.get_linemodel_extremacount(spectrum_model)
        fp_extremacount = self.accessor.get_linemodel_firstpass_extremacount(spectrum_model)
        if lm_extremacount and fp_extremacount:
            if fp_extremacount < lm_extremacount:
                raise APIException(
                    ErrorCode.INVALID_PARAMETER_FILE,
                    "linemodel.firstpass.extremaCount is lower than linemodel.extremaCount for object "
                    f"{spectrum_model}",
                )

    def _check_linemodelsolve_lya_fit(self, spectrum_model: str):
        if self.accessor.get_redshift_solver_method(spectrum_model) == ESolveMethod.LINE_MODEL:
            self._check_dependant_condition(
                self.accessor.get_linemodel_lya_profile(spectrum_model) == "asym",
                self.accessor.get_linemodel_lya_asym_section(spectrum_model) is not None,
                error_message=f"{ESolveMethod.LINE_MODEL.value} linemodel lya asymProfile section for object {spectrum_model}",
                warning_message=f"object {spectrum_model} {ESolveMethod.LINE_MODEL.value} lineModel lya asymProfile section",
            )

    def _check_linemeas_runmode(self):
        for spectrum_model in self.accessor.get_spectrum_models([]):
            if "lineMeasSolver" in self.accessor.get_stages(spectrum_model):
                self._check_dependant_condition(
                    True,
                    self.accessor.get_linemeas_runmode() is not None,
                    error_message=f"LineMeasSolve stage present, lineMeasRunMode parameter is mandatory",
                )
                break

    def useloglambdasampling_presence_condition(self, spectrum_model):
        return self.accessor.get_linemodel_continuum_component(spectrum_model) in [
            "tplFit",
            "tplFitAuto",
        ] and self.accessor.get_linemodel_continuumfit_fftprocessing(spectrum_model)

    def _check_linemodelsolve_useloglambdasampling(self, spectrum_model: str) -> None:
        self._check_dependant_condition(
            self.useloglambdasampling_presence_condition(spectrum_model),
            self.accessor.get_linemodel_useloglambdasampling(spectrum_model) is not None,
            f"{spectrum_model} lineModelSolve lineModel useLogLambdaSampling",
            f"{spectrum_model} lineModelSolve lineModel useLogLambdaSampling",
        )

    def _check_linemodelsolve_velocity_fit_params(self, spectrum_model: str):
        velocity_fit: bool = self.accessor.get_linemodel_velocity_fit(spectrum_model)
        for velocity_type in EVelocityType:
            for param in EVelocityFitParam:
                message = (
                    f"lineModelSolve {self.accessor.get_velocity_fit_param_name(velocity_type, param)}"
                    f" for object {spectrum_model}"
                )
                self._check_dependant_condition(
                    velocity_fit,
                    self.accessor.get_velocity_fit_param(
                        spectrum_model, ESolveMethod.LINE_MODEL, velocity_type, param
                    )
                    is not None,
                    error_message=message,
                    warning_message=message,
                )
        self._check_linemodelsolve_em_velocity_fit_min(spectrum_model)
        self._check_linemodelsolve_velocity_in_range(spectrum_model)

    def _check_linemodelsolve_em_velocity_fit_min(self, spectrum_model: str):
        velocityMin = self.accessor.get_linemodel_em_velocity_fit_min(spectrum_model)
        self._check_dependant_condition(
            self.accessor.get_linemodel_line_width_type(spectrum_model) == "velocityDriven",
            velocityMin is not None and velocityMin > 0,
            error_message=f"{spectrum_model} lineModelSolve lineModel emVelocityFit min must be > 0 when lineWidthType is velocityDriven",
            custom_error_message=True,
        )

    def _check_linemodelsolve_velocity_in_range(self, spectrum_model: str):
        velocityfit: bool = self.accessor.get_linemodel_velocity_fit(spectrum_model)
        if not velocityfit:
            return
        for velocity_type in EVelocityType:
            velocity = self.accessor.get_velocity(spectrum_model, ESolveMethod.LINE_MODEL, velocity_type)
            velocity_fit_min = self.accessor.get_velocity_fit_param(
                spectrum_model, ESolveMethod.LINE_MODEL, velocity_type, EVelocityFitParam.Min
            )
            velocity_fit_max = self.accessor.get_velocity_fit_param(
                spectrum_model, ESolveMethod.LINE_MODEL, velocity_type, EVelocityFitParam.Max
            )
            velocity_in_range = velocity >= velocity_fit_min and velocity <= velocity_fit_max
            self._check_dependant_condition(
                velocityfit,
                velocity_in_range,
                error_message=(
                    f"{ESolveMethod.LINE_MODEL.value}  {self.accessor.get_velocity_name(velocity_type)}"
                    f"={velocity} for object {spectrum_model} is outside "
                    f"[{self.accessor.get_velocity_fit_param_name(velocity_type, EVelocityFitParam.Min)}"
                    f" {self.accessor.get_velocity_fit_param_name(velocity_type, EVelocityFitParam.Max)}]"
                    f"=[{velocity_fit_min}, {velocity_fit_max}]"
                ),
                warning_message=None,
                custom_error_message=True,
            )

    def _check_linemeassolve_em_velocity_fit_min(self, spectrum_model: str):
        velocityMin = self.accessor.get_linemeas_em_velocity_fit_min(spectrum_model)
        self._check_dependant_condition(
            self.accessor.get_linemeas_line_width_type(spectrum_model) == "velocityDriven",
            velocityMin is not None and velocityMin > 0,
            error_message=f"{spectrum_model} lineMeasSolve lineModel emVelocityFit min must be > 0 when lineWidthType is velocityDriven",
            custom_error_message=True,
        )

    def _check_linemeassolve_lineratiotype_rules(self, spectrum_model: str):
        self._check_dependant_condition(
            self.accessor.get_linemeas_line_ratio_type(spectrum_model) == "rules",
            self.accessor.get_linemeas_rules(spectrum_model) is not None,
            error_message=f"lineMeasSolve rules for object {spectrum_model}",
            warning_message=f"object {spectrum_model} LineMeasSolve rules",
        )

    def _check_linemeassolve_fittingmethod_lbfgsb_velocityfit(self, spectrum_model: str):
        self._check_dependant_condition(
            self.accessor.get_linemeas_fitting_method(spectrum_model) == "lbfgsb",
            self.accessor.get_linemeas_velocity_fit(spectrum_model) is not None,
            error_message=f"lineMeasSolve velocityFit for object {spectrum_model}",
            warning_message=f"object {spectrum_model} LineMeasSolve velocityFit",
        )

    def _check_linemeassolve_velocity_fit_params(self, spectrum_model: str):
        velocity_fit: bool = self.accessor.get_linemeas_velocity_fit(spectrum_model)
        for velocity_type in EVelocityType:
            for param in [EVelocityFitParam.Min, EVelocityFitParam.Max]:
                message = (
                    f"lineMeasSolve {self.accessor.get_velocity_fit_param_name(velocity_type, param)}"
                    f" for object {spectrum_model}"
                )
                self._check_dependant_condition(
                    velocity_fit,
                    self.accessor.get_velocity_fit_param(
                        spectrum_model, ESolveMethod.LINE_MEAS, velocity_type, param
                    )
                    is not None,
                    error_message=message,
                    warning_message=message,
                )
        self._check_linemeassolve_em_velocity_fit_min(spectrum_model)
        self._check_linemeassolve_velocity_in_range(spectrum_model)
        self._check_linemeassolve_piped_linemodel_velocity(spectrum_model)

    def _check_linemeassolve_velocity_in_range(self, spectrum_model: str):
        velocity_fit: bool = self.accessor.get_linemeas_velocity_fit(spectrum_model)
        if not velocity_fit:
            return
        for velocity_type in EVelocityType:
            velocity = self.accessor.get_velocity(spectrum_model, ESolveMethod.LINE_MEAS, velocity_type)
            velocity_fit_min = self.accessor.get_velocity_fit_param(
                spectrum_model, ESolveMethod.LINE_MEAS, velocity_type, EVelocityFitParam.Min
            )
            velocity_fit_max = self.accessor.get_velocity_fit_param(
                spectrum_model, ESolveMethod.LINE_MEAS, velocity_type, EVelocityFitParam.Max
            )
            velocity_in_range = velocity >= velocity_fit_min and velocity <= velocity_fit_max
            self._check_dependant_condition(
                velocity_fit,
                velocity_in_range,
                error_message=(
                    f"{ESolveMethod.LINE_MEAS.value} {self.accessor.get_velocity_name(velocity_type)}"
                    f"={velocity} for object {spectrum_model} is outside "
                    f"[{self.accessor.get_velocity_fit_param_name(velocity_type, EVelocityFitParam.Min)},"
                    f" {self.accessor.get_velocity_fit_param_name(velocity_type, EVelocityFitParam.Max)}]="
                    f"[{velocity_fit_min}, {velocity_fit_max}]"
                ),
                warning_message=None,
                custom_error_message=True,
            )

    def _check_linemeassolve_piped_linemodel_velocity(self, spectrum_model: str):
        stages = self.accessor.get_stages(spectrum_model)
        redshift_solver_method = self.accessor.get_redshift_solver_method(spectrum_model)
        linemeas_piped_linemodel = (
            "redshiftSolver" in stages
            and "lineMeasSolver" in stages
            and redshift_solver_method == ESolveMethod.LINE_MODEL
        )
        linemeas_velocity_fit = self.accessor.get_linemeas_velocity_fit(spectrum_model)
        if not linemeas_piped_linemodel or not linemeas_velocity_fit:
            return
        linemodel_velocity_fit = self.accessor.get_linemodel_velocity_fit(spectrum_model)
        for velocity_type in EVelocityType:
            linemodel_velocity = self.accessor.get_velocity(
                spectrum_model, ESolveMethod.LINE_MODEL, velocity_type
            )
            velocity_fit_min = self.accessor.get_velocity_fit_param(
                spectrum_model, ESolveMethod.LINE_MEAS, velocity_type, EVelocityFitParam.Min
            )
            velocity_fit_max = self.accessor.get_velocity_fit_param(
                spectrum_model, ESolveMethod.LINE_MEAS, velocity_type, EVelocityFitParam.Max
            )
            velocity_in_range = (
                linemodel_velocity >= velocity_fit_min and linemodel_velocity <= velocity_fit_max
            )
            self._check_dependant_condition(
                linemeas_piped_linemodel and linemeas_velocity_fit and not linemodel_velocity_fit,
                velocity_in_range,
                error_message=(
                    f"{ESolveMethod.LINE_MODEL.value} {self.accessor.get_velocity_name(velocity_type)}"
                    f"={linemodel_velocity} for object {spectrum_model} is outside lineMeasSolve "
                    f"[{self.accessor.get_velocity_fit_param_name(velocity_type, EVelocityFitParam.Min)},"
                    f" {self.accessor.get_velocity_fit_param_name(velocity_type, EVelocityFitParam.Max)}]"
                    f"=[{velocity_fit_min}, {velocity_fit_max}]"
                ),
                warning_message=None,
                custom_error_message=True,
            )

            if not linemodel_velocity_fit:
                continue

            linemodel_velocity_fit_min = self.accessor.get_velocity_fit_param(
                spectrum_model, ESolveMethod.LINE_MODEL, velocity_type, EVelocityFitParam.Min
            )
            linemodel_velocity_fit_max = self.accessor.get_velocity_fit_param(
                spectrum_model, ESolveMethod.LINE_MODEL, velocity_type, EVelocityFitParam.Max
            )
            velocity_in_range = (
                linemodel_velocity_fit_min >= velocity_fit_min
                and linemodel_velocity_fit_max <= velocity_fit_max
            )
            self._check_dependant_condition(
                linemeas_piped_linemodel and linemeas_velocity_fit and linemodel_velocity_fit,
                velocity_in_range,
                error_message=(
                    "lineModelSolve velocity fit range "
                    f"[{self.accessor.get_velocity_fit_param_name(velocity_type, EVelocityFitParam.Min)},"
                    f" {self.accessor.get_velocity_fit_param_name(velocity_type, EVelocityFitParam.Max)}]"
                    f"=[{linemodel_velocity_fit_min}, {linemodel_velocity_fit_max}]"
                    " is not included in lineMeasSolve velocity fit range"
                    f"=[{velocity_fit_min}, {velocity_fit_max}] for object {spectrum_model}"
                ),
                warning_message=None,
                custom_error_message=True,
            )

    def _check_linemeassolve_lya_fit(self, spectrum_model: str):
        if self.accessor.get_linemeas_method(spectrum_model) == ESolveMethod.LINE_MEAS:
            self._check_dependant_condition(
                self.accessor.get_linemeas_lya_profile(spectrum_model) == "asym",
                self.accessor.get_linemeas_lya_asym_section(spectrum_model) is not None,
                error_message=f"{ESolveMethod.LINE_MEAS.value} linemodel lya asymProfile section for object {spectrum_model}",
                warning_message=f"object {spectrum_model} {ESolveMethod.LINE_MEAS.value} linemodel lya asymProfile section",
            )

    def _check_dependant_condition(
        self,
        triggering_condition: bool,
        dependant_condition: bool,
        error_message: Optional[str] = None,
        warning_message: Optional[str] = None,
        custom_error_message: bool = False,
        custom_warning_message: bool = False,
    ):
        # Checks that if triggering_condition is valid, then dependant_condition too. Otherwise error.
        # If a warning message is specified, raises an error if dependant_condition is present without its
        # triggering_condition
        if triggering_condition and not dependant_condition:
            if error_message is not None:
                if custom_error_message:
                    raise APIException(ErrorCode.INVALID_PARAMETER_FILE, error_message)
                else:
                    self._raise_missing_param_error(error_message)
        elif not triggering_condition and dependant_condition:
            if warning_message is not None:
                if custom_warning_message:
                    zflag.warning(WarningCode.UNUSED_PARAMETER, warning_message)
                else:
                    self._add_unused_parameter_warning(warning_message)

    def _add_unused_parameter_warning(self, param_name):
        zflag.warning(WarningCode.UNUSED_PARAMETER, f"Unused parameter {param_name}")

    def _raise_missing_param_error(self, param_name: str) -> None:
        raise APIException(ErrorCode.INVALID_PARAMETER_FILE, f"Missing parameter {param_name}")
