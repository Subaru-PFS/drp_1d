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
import os

from jsonschema import RefResolver, validate
from jsonschema.exceptions import ValidationError
from pylibamazed.Exception import APIException
from pylibamazed.FilterLoader import ParamJsonFilterLoader
from pylibamazed.Paths import module_root_dir
from pylibamazed.redshift import CFlagWarning, ErrorCode, WarningCode

zflag = CFlagWarning.GetInstance()


class JsonSchemaFileAccessor:
    def get_json_schema(self, version: dict) -> dict:
        with open(self.json_schema_filename(version)) as schemaFile:
            jsonSchema = json.load(schemaFile)
        return jsonSchema

    def json_schema_filename(self, version: str) -> str:
        return os.path.join(self._json_schema_path(version), "general.json")

    def json_schema_filesystem_path(self, version: str) -> str:
        return f"file:///{self._json_schema_path(version)}/"

    def _json_schema_path(self, version: str) -> str:
        return os.path.join(module_root_dir, "resources", f"jsonschema-v{version}")


class ParametersChecker:

    jsonSchemaFilename = os.path.join(module_root_dir, "resources", "jsonschema-v1", "general.json")

    def __init__(self, FilterLoader=ParamJsonFilterLoader, FileAccessor=JsonSchemaFileAccessor):
        self.filter_loader = FilterLoader()
        self.file_accesor = FileAccessor()

    def json_schema_check(self, parameters_dict: dict, version: int) -> None:
        # Parameters dict must be at "raw parameters dict"
        # How to link the different json schema files
        jsonSchema = self.file_accesor.get_json_schema(version)
        resolver = RefResolver(self.file_accesor.json_schema_filesystem_path(version), jsonSchema)

        try:
            validate(instance=parameters_dict, schema=jsonSchema, resolver=resolver)
        except ValidationError as e:
            raise APIException(
                ErrorCode.INVALID_PARAMETER_FILE,
                e.message
            )

    def custom_check(self, accessor):
        self.accessor = accessor
        self._check_general()
        self._check_lsf()
        self._check_continuum_removal()
        self._check_continuum_removal("templateCatalog")
        for object in accessor.get_spectrum_models([]):
            self._check_object(object)

    def _check_general(self):
        self._check_photometry_transmission_dir()
        self._check_photometry_band()
        self._check_filters()

    def _check_photometry_transmission_dir(self):
        parameter_name = "photometryTransmissionDir"
        self._check_dependant_parameter_presence(
            self.accessor.photometry_is_enabled(),
            self.accessor.get_photometry_transmission_dir() is not None,
            parameter_name,
            parameter_name
        )

    def _check_photometry_band(self):
        parameter_name = "photometryBand"
        self._check_dependant_parameter_presence(
            self.accessor.photometry_is_enabled(),
            len(self.accessor.get_photometry_bands()) > 0,
            parameter_name,
            parameter_name
        )

    def _check_filters(self):
        filters = self.accessor.get_filters(default=[])
        self._check_filters_format(filters)

        DEFAULT_COLUMN_NAMES = ["Wave", "Flux", "Err"]
        if not filters:
            return
        filter_keys = [filt["key"] for filt in filters]
        authorized_cols_names = DEFAULT_COLUMN_NAMES + self.accessor.get_additional_cols(default=[])

        for filter_name in filter_keys:
            if filter_name not in authorized_cols_names:
                raise APIException(
                    ErrorCode.INVALID_PARAMETER_FILE,
                    f"Unknown filter key {filter_name}"
                )

    def _check_filters_format(self, json: json) -> None:
        if type(json) is not list:
            raise APIException(
                ErrorCode.INVALID_PARAMETER_FILE,
                "Input filters json must be a list"
            )

        for filt in json:
            json_keys = self.filter_loader.keys
            different_keys = set(filt.keys()) != set(json_keys)
            different_length = len(filt.keys()) != len(json_keys)
            if different_keys or different_length:
                raise APIException(
                    ErrorCode.INVALID_PARAMETER_FILE,
                    f"Filters: each dictionary in json list must have exactly the following keys {json_keys}"
                )

    def _check_lsf(self) -> None:
        self._check_GaussianConstantWidth_lsf_type()
        self._check_GaussianNISPSIM_lsf_type()
        self._check_GaussianConstantResolution_lsf_type()
        self._check_GaussianVariablewidth_FileName()

    def _check_GaussianNISPSIM_lsf_type(self):
        self._check_dependant_parameter_presence(
            self.accessor.get_lsf_type() == "GaussianNISPSIM201707",
            self.accessor.get_lsf_sourcesize() is not None,
            "lsf sourceSize",
            "lsf sourceSize"
        )

    def _check_GaussianConstantWidth_lsf_type(self):
        self._check_dependant_parameter_presence(
            self.accessor.get_lsf_type() == "gaussianConstantWidth",
            self.accessor.get_lsf_width() is not None,
            "lsf width",
            "lsf width"
        )

    def _check_GaussianConstantResolution_lsf_type(self):
        self._check_dependant_parameter_presence(
            self.accessor.get_lsf_type() == "gaussianConstantResolution",
            self.accessor.get_lsf_resolution() is not None,
            "lsf resolution",
            "lsf resolution"
        )

    def _check_GaussianVariablewidth_FileName(self):
        self._check_dependant_parameter_presence(
            self.accessor.get_lsf_type() == "gaussianVariableWidth",
            self.accessor.get_lsf_width_file_name() is not None,
            "lsf gaussianVariableWidthFileName",
            "lsf gaussianVariableWidthFileName"
        )

    def _check_continuum_removal(self, nested=None) -> None:
        self._check_IrregularSamplingMedian_kernel_width(nested)
        self._check_IrregularSamplingMedian_kernel_reflection(nested)

    def _check_IrregularSamplingMedian_kernel_width(self, nested=None):
        self._check_dependant_parameter_presence(
            self.accessor.get_continuum_removal_method(nested) == "irregularSamplingMedian",
            self.accessor.get_continuum_removal_median_kernel_width(nested) is not None,
            "continuumRemoval medianKernelWidth",
            "continuumRemoval medianKernelWidth"
        )

    def _check_IrregularSamplingMedian_kernel_reflection(self, nested=None):
        self._check_dependant_parameter_presence(
            self.accessor.get_continuum_removal_method(nested) == "irregularSamplingMedian",
            self.accessor.get_continuum_median_kernel_reflection(nested) is not None,
            "continuumRemoval medianEvenReflection",
            "continuumRemoval medianEvenReflection"
        )

    def _check_object(self, object_type: str) -> None:
        self._check_linemeassolve_dzhalf(object_type)
        self._check_linemeassolve_redshiftstep(object_type)
        self._check_object_reliability(object_type)
        self._check_templateFittingSolve_section(object_type)
        self._check_templateCombinationSolve_section(object_type)
        self._check_lineModelSolve(object_type)
        self._check_template_dir(object_type)

    def _check_template_dir(self, object_type: str) -> None:
        template_dir_presence_condition = \
            self.accessor.get_redshift_solver_method(object_type) in ["templateFittingSolve",
                                                                      "tplCombinationSolve"] or (
                self.accessor.get_redshift_solver_method(object_type) == "lineModelSolve"
                and self.accessor.get_linemodel_continuum_component(object_type) in ["tplFit", "tplFitAuto"]
            )
        self._check_dependant_parameter_presence(
            template_dir_presence_condition,
            self.accessor.get_template_dir(object_type) is not None,
            f"{object_type} templateDir",
            f"{object_type} templateDir"
        )

    def _check_linemeassolve_dzhalf(self, object_type: str) -> None:
        self._check_dependant_parameter_presence(
            self.accessor.get_linemeas_method(object_type) is not None,
            self.accessor.get_linemeas_dzhalf(object_type) is not None,
            f"lineMeasDzHalf for object {object_type}",
            f"object {object_type} lineMeasDzHalf"
        )

    def _check_linemeassolve_redshiftstep(self, object_type: str) -> None:
        self._check_dependant_parameter_presence(
            self.accessor.get_linemeas_method(object_type) is not None,
            self.accessor.get_linemeas_redshiftstep(object_type) is not None,
            f"lineMeasRedshiftStep for object {object_type}",
            f"object {object_type} lineMeasRedshiftStep"
        )

    def _check_object_reliability(self, object_type: str) -> None:
        self._check_dependant_parameter_presence(
            self.accessor.get_reliability_enabled(object_type),
            self.accessor.get_reliability_model(object_type) is not None,
            error_message=f"reliabilityModel for object {object_type}",
            warning_message=f"object {object_type} reliability_enabled"
        )

    def _check_templateFittingSolve_section(self, object_type: str) -> None:
        self._check_dependant_parameter_presence(
            self.accessor.get_redshift_solver_method(object_type) == "templateFittingSolve",
            self.accessor.get_template_fitting_section(object_type) is not None,
            error_message=f"templateFittingSolve for object {object_type}",
            warning_message=f"object {object_type} templateFittingSolve"
        )

        self._check_templateFittingSolve_ism(object_type)
        self._check_templateFittingSolve_photometry_weight(object_type)

    def _check_templateFittingSolve_ism(self, object_type: str) -> None:
        self._check_dependant_parameter_presence(
            self.accessor.get_template_fitting_ism(object_type),
            self.accessor.get_ebmv_section() is not None,
            "ebmv"
        )

    def _check_templateFittingSolve_photometry_weight(self, object_type: str) -> None:
        self._check_dependant_parameter_presence(
            self.accessor.get_template_fitting_photometry_enabled(object_type),
            self.accessor.get_template_fitting_photometry_weight(object_type) is not None,
            f"object {object_type} TemplateFittingSolve photometry weight",
            f"object {object_type} TemplateFittingSolve photometry weight"
        )

    def _check_templateCombinationSolve_section(self, object_type: str) -> None:
        self._check_dependant_parameter_presence(
            self.accessor.get_redshift_solver_method(object_type) == "tplCombinationSolve",
            self.accessor.get_template_combination_section(object_type) is not None,
            error_message=f"tplCombinationSolve for object {object_type}",
            warning_message=f"object {object_type} tplCombinationSolve"
        )
        self._check_templateCombinationSolve_ism(object_type)

    def _check_templateCombinationSolve_ism(self, object_type: str) -> None:
        ism = self.accessor.get_template_combination_ism(object_type)
        ebmv = self.accessor.get_ebmv_section()

        if ism and ebmv is None:
            self._raise_missing_param_error("ebmv")

    def _check_lineModelSolve(self, object_type: str) -> None:

        self._check_linemodelsolve_section(object_type)

        self._check_linemodelsolve_improveBalmerFit(object_type)
        self._check_linemodelsolve_lineratiotype_rules(object_type)

        self._check_lineratiotype_tplratio_catalog(object_type)
        self._check_lineratiotype_tplratio_ismfit(object_type)

        self._check_linemodelsolve_continuumremoval(object_type)
        self._check_linemodelsolve_continuumreestimation(object_type)

        self._check_linemodelsolve_continuumfit_section(object_type)
        self._check_linemodelsolve_continuumfit_ism(object_type)
        self._check_linemodelsolve_secondpass_section(object_type)
        self._check_linemodelsolve_secondpass_continuumfit(object_type)
        self._check_linemodelsolve_firstpass_tplratio_ismfit(object_type)
        self._check_linemodelsolve_firstpass_extremacount(object_type)

        self._check_linemeassolve_lineratiotype_rules(object_type)
        self._check_linemeassolve_fittingmethod_lbfgsb_velocityfit(object_type)
        self._check_linemeassolve_velocityfit_params(object_type)

    def _check_linemodelsolve_section(self, object_type: str):
        self._check_dependant_parameter_presence(
            self.accessor.get_redshift_solver_method(object_type) == "lineModelSolve",
            self.accessor.get_linemodel_solver_section(object_type) is not None,
            error_message=f"lineModelSolve for object {object_type}",
            warning_message=f"object {object_type} lineModelSolve"
        )

    def _check_linemodelsolve_lineratiotype_rules(self, object_type: str):
        self._check_dependant_parameter_presence(
            self.accessor.get_linemodel_line_ratio_type(object_type) == "rules",
            self.accessor.get_linemodel_rules(object_type) is not None,
            error_message=f"lineModelSolve rules for object {object_type}",
            warning_message=f"object {object_type} lineModelSolve rules"
        )

    def _check_linemodelsolve_improveBalmerFit(self, object_type: str):
        improveBalmerFit = self.accessor.get_linemodel_improve_balmer_fit(object_type)
        lineRatioType = self.accessor.get_linemodel_line_ratio_type(object_type)
        if improveBalmerFit and lineRatioType != "rules":
            zflag.warning(
                WarningCode.UNUSED_PARAMETER.value,
                f"object {object_type} lineModelSolve lineRatioType must be rules to "
                "activate improveBalmerFit"
            )

    def _check_lineratiotype_tplratio_catalog(self, object_type):
        self._check_dependant_parameter_presence(
            self.accessor.get_linemodel_line_ratio_type(object_type) in ["tplRatio", "tplCorr"],
            self.accessor.get_linemodel_tplratio_catalog(object_type) is not None,
            error_message=f"lineModelSolve tplRatioCatalog for object {object_type}",
            warning_message=f"object {object_type} lineModelSolve tplRatioCatalog"
        )

    def _check_lineratiotype_tplratio_ismfit(self, object_type):
        self._check_dependant_parameter_presence(
            self.accessor.get_linemodel_line_ratio_type(object_type) in ["tplRatio", "tplCorr"],
            self.accessor.get_linemodel_tplratio_ismfit(object_type) is not None,
            error_message=f"lineModelSolve tplRatioIsmFit for object {object_type}",
            warning_message=f"object {object_type} lineModelSolve tplRatioIsmFit"
        )

    def _check_linemodelsolve_continuumremoval(self, object_type):
        self._check_dependant_parameter_presence(
            self.accessor.get_linemodel_continuum_component(object_type) in [
                "fromSpectrum", "tplFitAuto"],
            self.accessor.get_continuum_removal() is not None,
            error_message="continuumRemoval"
        )

    def _check_linemodelsolve_continuumreestimation(self, object_type):
        self._check_dependant_parameter_presence(
            self.accessor.get_linemodel_continuum_component(object_type) == "fromSpectrum",
            self.accessor.get_linemodel_continuum_reestimation(object_type) is not None,
            error_message=f"object {object_type} lineModelSolve continuumReestimation",
            warning_message=f"object {object_type} lineModelSolve continuumReestimation"
        )

    def _check_linemodelsolve_continuumfit_section(self, object_type):
        self._check_dependant_parameter_presence(
            self.accessor.get_linemodel_continuum_component(object_type) in ["tplFit", "tplFitAuto"],
            self.accessor.get_linemodel_continuumfit_section(object_type) is not None,
            error_message=f"object {object_type} lineModelSolve continuumFit"
        )

    def _check_linemodelsolve_secondpass_section(self, object_type):
        self._check_dependant_parameter_presence(
            self.accessor.get_linemodel_skipsecondpass(object_type) is False,
            self.accessor.get_linemodel_secondpass_section(object_type) is not None,
            f"object {object_type} lineModelSolve secondPass",
            f"object {object_type} lineModelSolve secondPass"
        )

    def _check_linemodelsolve_secondpass_continuumfit(self, object_type):
        self._check_dependant_parameter_presence(
            self.accessor.get_linemodel_continuum_component(object_type) in ["tplFit", "tplFitAuto"],
            self.accessor.get_linemodel_secondpass_continuumfit(object_type) is not None,
            error_message=f"object {object_type} lineModelSolve secondpass continuumFit"
        )

    def _check_linemodelsolve_continuumfit_ism(self, object_type: str) -> None:
        self._check_dependant_parameter_presence(
            self.accessor.get_linemodel_continuumfit_ismfit(object_type),
            self.accessor.get_ebmv_section() is not None,
            "ebmv"
        )

    def _check_linemodelsolve_firstpass_tplratio_ismfit(self, object_type: str):
        self._check_dependant_parameter_presence(
            self.accessor.get_linemodel_line_ratio_type(object_type) in ["tplRatio", "tplCorr"],
            self.accessor.get_linemodel_firstpass_tplratio_ismfit(object_type) is not None,
            f"object {object_type} lineModelSolve firstpass tplRatioIsmFit",
            f"object {object_type} lineModelSolve firstpass tplRatioIsmFit"
        )

    def _check_linemodelsolve_firstpass_extremacount(self, object_type: str):
        lm_extremacount = self.accessor.get_linemodel_extremacount(object_type)
        fp_extremacount = self.accessor.get_linemodel_firstpass_extremacount(object_type)
        if lm_extremacount and fp_extremacount:
            if fp_extremacount < lm_extremacount:
                raise APIException(
                    ErrorCode.INVALID_PARAMETER_FILE,
                    "linemodel.firstpass.extremaCount is lower than linemodel.extremaCount for object "
                    f"{object_type}"
                )

    def _check_linemeassolve_lineratiotype_rules(self, object_type: str):
        self._check_dependant_parameter_presence(
            self.accessor.get_linemeas_line_ratio_type(object_type) == "rules",
            self.accessor.get_linemeas_rules(object_type) is not None,
            error_message=f"lineMeasSolve rules for object {object_type}",
            warning_message=f"object {object_type} LineMeasSolve rules"
        )

    def _check_linemeassolve_fittingmethod_lbfgsb_velocityfit(self, object_type: str):
        self._check_dependant_parameter_presence(
            self.accessor.get_linemeas_fitting_method(object_type) == "lbfgsb",
            self.accessor.get_linemeas_velocity_fit(object_type) is not None,
            error_message=f"lineMeasSolve velocityFit for object {object_type}",
            warning_message=f"object {object_type} LineMeasSolve velocityFit"
        )

    def _check_linemeassolve_velocityfit_params(self, object_type: str):
        params = ["emVelocityFitMin", "emVelocityFitMax", "absVelocityFitMin", "absVelocityFitMax"]
        velocityfit: bool = self.accessor.get_linemeas_velocity_fit(object_type)
        for param in params:
            self._check_dependant_parameter_presence(
                velocityfit,
                self.accessor.get_linemeas_velocity_fit_param(object_type, param) is not None,
                error_message=f"lineMeasSolve {param} for object {object_type}",
                warning_message=f"object {object_type} LineMeasSolve {param}"
            )

    def _check_dependant_parameter_presence(self, triggering_condition: bool, dependant_condition: bool,
                                            error_message: str = None,
                                            warning_message: str = None,
                                            custom_error_message: str = False,
                                            custom_warning_message: str = False
                                            ):
        # Checks that if triggering_condition is valid, then dependant_condition too. Otherwise error.
        # If a warning message is specified, raises an error if dependant_condition is present without its
        # triggering_condition
        if triggering_condition and not dependant_condition:
            if error_message is not None:
                if custom_error_message:
                    raise APIException(
                        ErrorCode.INVALID_PARAMETER_FILE,
                        error_message
                    )
                else:
                    self._raise_missing_param_error(error_message)
        elif not triggering_condition and dependant_condition:
            if warning_message is not None:
                if custom_warning_message:
                    zflag.warning(
                        WarningCode.UNUSED_PARAMETER.value,
                        warning_message
                    )
                else:
                    self._add_unused_parameter_warning(warning_message)

    def _add_unused_parameter_warning(self, param_name):
        zflag.warning(
            WarningCode.UNUSED_PARAMETER.value,
            f"Unused parameter {param_name}"
        )

    def _raise_missing_param_error(self, param_name: str) -> None:
        raise APIException(
            ErrorCode.INVALID_PARAMETER_FILE,
            f"Missing parameter {param_name}"
        )
