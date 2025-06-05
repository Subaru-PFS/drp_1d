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
import os
import os.path
import copy
from contextlib import contextmanager, suppress
from decorator import decorator
from typing import Optional, Any

from pylibamazed.CalibrationLibrary import CalibrationLibrary
from pylibamazed.Exception import APIException, exception_decorator
from pylibamazed.Parameters import Parameters
from pylibamazed.ParametersAccessor import ESolveMethod
from pylibamazed.Spectrum import Spectrum
from datetime import datetime
import resource

# NB: DO NOT REMOVE - these libs are used in globals
from pylibamazed.redshift import CClassificationSolve  # noqa F401
from pylibamazed.redshift import CLineMatchingSolve  # noqa F401
from pylibamazed.redshift import CLineMeasSolve  # noqa F401
from pylibamazed.redshift import CLineModelSolve  # noqa F401
from pylibamazed.redshift import CTemplateFittingSolve  # noqa F401
from pylibamazed.redshift import CTplCombinationSolve  # noqa F401
from pylibamazed.redshift import AmzException, CFlagWarning, CLog, CProcessFlowContext, ErrorCode, ScopeType
from pylibamazed.AbstractReliabilitySolver import (
    get_reliability_solver_from_name,
    get_reliability_dataset_suffix,
)
import pylibamazed.DeepLearningSolve  # noqa F401
import pylibamazed.SkLearnSolve  # noqa F401
from pylibamazed.ResultStoreOutput import ResultStoreOutput
from pylibamazed.ScopeManager import get_scope_spectrum_model, get_scope_stage, push_scope
from pylibamazed.SubType import SubType
from pylibamazed.LinemeasParameters import LinemeasParameters
from enum import Enum

zflag = CFlagWarning.GetInstance()
zlog = CLog.GetInstance()


class EProcessingMode(Enum):
    TWO_OR_SINGLE_PASS = "normal"
    FIRST_PASS_ONLY = "firstPass"
    SECOND_PASS_AND_PDF = "secondPass"


class ProcessFlowException(Exception):
    pass


def _get_time():
    ret = dict()
    ret["clock"] = datetime.now()
    ret["user"] = resource.getrusage(resource.RUSAGE_SELF).ru_utime
    ret["system"] = resource.getrusage(resource.RUSAGE_SELF).ru_stime
    return ret


# TO BE USED at ScopeType.METHOD
@contextmanager
def store_flags(name="warningFlag"):
    try:
        yield
    finally:
        resultStore = CProcessFlowContext.GetInstance().GetResultStore()
        resultStore.StoreScopedFlagResult(name)
        zflag.resetFlag()


class ProcessFlow:
    @exception_decorator
    def __init__(self, config: dict[str, Any], parameters: Parameters):
        _check_config(config)
        _check_lineMeasValidity(config, parameters)
        self.parameters: Parameters = parameters
        zlog.LogInfo("Loading all needed calibration files :")
        self.calibration_library = CalibrationLibrary(parameters, config["calibration_dir"])
        self.calibration_library.load_all()
        zlog.LogInfo("Calibration files loaded.")
        self.process_flow_context = CProcessFlowContext.GetInstance()
        self.process_flow_context.reset()
        self.config = config
        if "linemeascatalog" not in self.config:
            self.config["linemeascatalog"] = {}

        self.extended_results = config["extended_results"]

        # save context warning flag to reinject at each spectrum
        resultStore = self.process_flow_context.GetResultStore()
        resultStore.StoreScopedFlagResult("context_warningFlag")
        zflag.resetFlag()
        self.context_warning_Flag = resultStore.GetFlagLogResult("", "", "", "context_warningFlag")

    @property
    def _scope_spectrum_model(self):
        return get_scope_spectrum_model()

    @property
    def _scope_stage(self):
        return get_scope_stage()

    @decorator
    def _store_exception(func, self, *args, **kwargs):
        try:
            new_func = exception_decorator(logging=True)(func)
            return new_func(self, *args, **kwargs)
        except AmzException as e:
            self.rso.store_error(e, self._scope_spectrum_model, self._scope_stage)
            raise ProcessFlowException from e

    @decorator
    def _store_perfs(func, self, *args, **kwargs):
        start = _get_time()
        try:
            return func(self, *args, **kwargs)
        finally:
            end = _get_time()
            perfs = dict()
            for k in start.keys():
                perfs[k] = end[k] - start[k]
            perfs["clock"] = perfs["clock"].total_seconds()
            self.rso.store_perfs(
                self._scope_spectrum_model,
                self._scope_stage,
                perfs,
                kwargs.get("mode", EProcessingMode("normal")).value,
            )

    @exception_decorator
    def run(self, spectrum: Spectrum) -> ResultStoreOutput:
        resultStore = self.process_flow_context.GetResultStore()
        self.rso = ResultStoreOutput(
            resultStore, self.parameters, auto_load=False, extended_results=self.extended_results
        )
        try:
            self._initialize(spectrum)
        except ProcessFlowException:
            self.rso.load_root()
            return self.rso

        processing_mode = EProcessingMode("normal")
        if self.parameters.second_pass_after_classification():
            processing_mode = EProcessingMode("firstPass")

        # loop on spectrum models (galaxy, star, qso, ...)
        for spectrum_model in self.parameters.get_spectrum_models():
            with push_scope(spectrum_model, ScopeType.SPECTRUMMODEL):
                with suppress(ProcessFlowException):
                    self._process_spectrum_model(processing_mode)

        if self.parameters.is_a_redshift_solver_used():
            with suppress(ProcessFlowException):
                classif_model = self._run_classification_solver()
                with push_scope(classif_model, ScopeType.SPECTRUMMODEL):
                    if processing_mode == EProcessingMode.FIRST_PASS_ONLY:
                        self._run_second_pass_after_classification(classif_model)
                    # Running linemeas only on classified model (if any)
                    if self.parameters.get_linemeas_runmode() == "classif":
                        self._run_linemeas_after_classification(classif_model)

        with suppress(ProcessFlowException):
            self._load_result_store()

        return self.rso

    @_store_perfs
    @store_flags(name="init_warningFlag")
    @_store_exception
    def _initialize(self, spectrum: Spectrum) -> None:
        zlog.LogInfo("Context initialization")
        self.process_flow_context.reset()

        # inject in the clean resultStore the contextFlag
        self.rso.results_store.StoreScopedGlobalResult("context_warningFlag", self.context_warning_Flag)

        for object_type in self.parameters.get_spectrum_models():
            if object_type in self.calibration_library.line_catalogs:
                for method in self.parameters.get_linemodel_methods_str(object_type):
                    self.process_flow_context.setLineCatalog(
                        object_type, method, self.calibration_library.line_catalogs[object_type][method]
                    )
            if object_type in self.calibration_library.line_ratio_catalog_lists:
                self.process_flow_context.setLineRatioCatalogCatalog(
                    object_type, self.calibration_library.line_ratio_catalog_lists[object_type]
                )
        self.process_flow_context.setTemplateCatalog(self.calibration_library.templates_catalogs["all"])
        self.process_flow_context.setPhotBandCatalog(self.calibration_library.photometric_bands)
        self.process_flow_context.setFluxCorrectionMeiksin(self.calibration_library.meiksin)
        self.process_flow_context.setFluxCorrectionCalzetti(self.calibration_library.calzetti)

        spectrum.init()
        spectrum.push_in_context()

        parameters = copy.deepcopy(self.parameters)
        if self.config.get("linemeascatalog"):
            lp = LinemeasParameters()
            lp.load_from_catalogs(
                spectrum.source_id,
                self.config["linemeascatalog"],
                self.config["linemeas_catalog_columns"],
            )
            lp.update_parameters(parameters)

        self.process_flow_context.LoadParameterStore(parameters.to_json())
        self.process_flow_context.Init()

    def _process_spectrum_model(self, mode: EProcessingMode) -> None:
        spectrum_model = self._scope_spectrum_model

        redshift_solver_method = self.parameters.get_redshift_solver_method(spectrum_model)
        linemeas_method = self.parameters.get_linemeas_method(spectrum_model)

        if redshift_solver_method:
            self._run_redshift_solver(redshift_solver_method.value, mode)
            if self.parameters.is_tplratio_catalog_needed(spectrum_model):
                with suppress(ProcessFlowException):
                    self._run_sub_classification_solver()

            if self.parameters.get_reliability_enabled(spectrum_model):
                with suppress(ProcessFlowException):
                    self._run_reliability_solver()

            if self.parameters.get_linemeas_runmode() == "all" and linemeas_method:
                self._run_load_linemeas_params()
                self._run_linemeas_solver(linemeas_method.value)

        elif linemeas_method:  # linemeas alone
            self._run_linemeas_solver(linemeas_method.value)

    def _run_linemeas_after_classification(self, classif_model: str) -> None:
        linemeas_method = getattr(self.parameters.get_linemeas_method(classif_model), "value", None)
        if linemeas_method is None:
            return
        self._run_load_linemeas_params()
        self._run_linemeas_solver(linemeas_method)

    def _run_second_pass_after_classification(self, classif_model: str) -> None:
        # for spectrum_model in self.parameters.get_spectrum_models():
        redshift_solver_method = self.parameters.get_redshift_solver_method(classif_model)
        if redshift_solver_method != ESolveMethod.LINE_MODEL:
            return
        self._process_spectrum_model(EProcessingMode("secondPass"))

    @push_scope("redshiftSolver", ScopeType.STAGE)
    @_store_perfs()
    @_store_exception
    def _run_redshift_solver(
        self, method: str, mode: Optional[EProcessingMode] = EProcessingMode("normal")
    ) -> None:
        self._run_cpp_method(method, mode)

    @push_scope("lineMeasSolver", ScopeType.STAGE)
    @_store_perfs
    @_store_exception
    def _run_linemeas_solver(self, method: str) -> None:
        self._run_cpp_method(method)

    @push_scope("linemeas_catalog_load", ScopeType.STAGE)
    @_store_exception
    def _run_load_linemeas_params(self) -> None:
        parameters = copy.deepcopy(self.parameters)
        lp = LinemeasParameters()
        lp.load_from_result_store(parameters, self.rso, self._scope_spectrum_model)
        lp.update_parameters(parameters)
        self.process_flow_context.LoadParameterStore(parameters.to_json())

    @push_scope("reliabilitySolver", ScopeType.STAGE)
    @_store_perfs
    @_store_exception  # NOTE: an exception raised in one reliability method prevent the following to be run
    def _run_reliability_solver(self) -> None:
        reliability_methods = self.parameters.get_reliability_methods(self._scope_spectrum_model)
        if reliability_methods is None:
            raise APIException(
                ErrorCode.INVALID_PARAMETER, f"{self._scope_spectrum_model}.reliabilitySolver.method empty"
            )
        for solver_name in reliability_methods:
            zlog.LogInfo(f"run reliability with {solver_name}")
            solver = get_reliability_solver_from_name(solver_name)
            solverArgs = [self._scope_spectrum_model, self.parameters, self.calibration_library]
            rel = self._run_python_method(solver, self.process_flow_context, classArgs=solverArgs)
            dataset_suffix = get_reliability_dataset_suffix(solver_name)
            # rel = solver(self._scope_spectrum_model, self.parameters, self.calibration_library)
            self.rso.object_results[self._scope_spectrum_model][f"reliability{dataset_suffix}"] = dict()
            self.rso.object_results[self._scope_spectrum_model][f"reliability{dataset_suffix}"][
                "Reliability"
            ] = rel

    @push_scope("subClassifSolver", ScopeType.STAGE)
    @_store_perfs
    @_store_exception
    def _run_sub_classification_solver(self) -> None:
        subTypeArgs = [self._scope_spectrum_model, self.parameters, self.calibration_library]
        sub_types = self._run_python_method(SubType, self.process_flow_context, classArgs=subTypeArgs)
        self.rso.object_results[self._scope_spectrum_model]["model_parameters"] = []
        for rank in range(len(sub_types)):
            self.rso.object_results[self._scope_spectrum_model]["model_parameters"].append(dict())
            self.rso.object_results[self._scope_spectrum_model]["model_parameters"][rank][
                "SubType"
            ] = sub_types[rank]

    def _all_redshift_solver_failed(self) -> bool:
        for obj in self.parameters.get_spectrum_models():
            if not self.rso.has_error(obj, "redshiftSolver"):
                return False
        return True

    @push_scope("classification", ScopeType.SPECTRUMMODEL)
    @push_scope("classification", ScopeType.STAGE)
    @_store_perfs
    @_store_exception
    def _run_classification_solver(self) -> str:
        if self._all_redshift_solver_failed():
            raise APIException(
                ErrorCode.NO_CLASSIFICATION, "Classification not run because all redshiftSolver failed"
            )
        self._run_cpp_method("classificationSolve")
        classif_model = self.rso.get_attribute_from_source("root", None, None, "classification", "Type")
        return classif_model

    @push_scope("load_result_store", ScopeType.STAGE)
    @_store_exception
    def _load_result_store(self) -> None:
        self.rso.load_all()

    def _run_cpp_method(self, method: str, mode: Optional[EProcessingMode] = None) -> None:
        method_to_solver = {
            "classificationSolve": "CClassificationSolve",
            "lineMeasSolve": "CLineMeasSolve",
            "lineModelSolve": "CLineModelSolve",
            "templateFittingSolve": "CTemplateFittingSolve",
            "tplCombinationSolve": "CTplCombinationSolve",
        }
        if method_to_solver[method] not in globals():
            raise APIException(ErrorCode.INVALID_PARAMETER, "Unknown method {}".format(method))
        solver_method = globals()[method_to_solver[method]]
        solver = solver_method()
        zlog.LogInfo(f"Running method {method} on mode {mode}")
        if mode is not None:
            if mode == EProcessingMode.FIRST_PASS_ONLY:
                solver.initForClassificationAfterFirstPass()
            elif mode == EProcessingMode.SECOND_PASS_AND_PDF:
                solver.setRunSecondPassFromResultStore()
        solver.Compute()

    def _run_python_method(
        self, methodClass: type, *args, classArgs: list = [], classKwargs: dict = {}, **kwargs
    ):
        with push_scope(methodClass.__name__, ScopeType.METHOD):
            solver = methodClass(*classArgs, **classKwargs)
            with store_flags():
                results = solver.Compute(*args, **kwargs)
        return results


def _check_config(config: dict) -> None:
    if "calibration_dir" not in config:
        raise APIException(ErrorCode.MISSING_CONFIG_OPTION, "Config must contain 'calibration_dir' key")
    if not os.path.exists(config["calibration_dir"]):
        raise APIException(
            ErrorCode.INVALID_DIRECTORY,
            "Calibration directory {} does not exist".format(config["calibration_dir"]),
        )
    if "linemeascatalog" in config:
        if "linemeas_catalog_columns" not in config:
            raise APIException(
                ErrorCode.MISSING_CONFIG_OPTION,
                "Missing linemeas_catalog_columns key in linemeascatalog config-option",
            )
        for object_type in config["linemeascatalog"].keys():
            if object_type not in config["linemeas_catalog_columns"]:
                raise APIException(
                    ErrorCode.INCOHERENT_CONFIG_OPTIONS,
                    "Missing category {} in linemeas_catalog_columns ".format(object_type),
                )

            for attr in ["Redshift", "VelocityAbsorption", "VelocityEmission"]:
                if attr not in config["linemeas_catalog_columns"][object_type]:
                    raise APIException(
                        ErrorCode.ATTRIBUTE_NOT_SUPPORTED,
                        f"Not supported Attribute {object_type} in Config['linemeas_catalog_columns'][{attr}]",
                    )


def _check_lineMeasValidity(config: dict, parameters: Parameters) -> None:
    for spectrum_model in parameters.get_spectrum_models():
        if parameters.is_linemeas_alone(spectrum_model):
            if "linemeascatalog" not in config or spectrum_model not in config["linemeascatalog"].keys():
                raise APIException(
                    ErrorCode.INCOHERENT_CONFIG_OPTIONS,
                    "Cannot run lineMeasSolver alone without linemeascatalog in config "
                    f"for model {spectrum_model}.",
                )
        elif parameters.is_linemeas_piped(spectrum_model):
            if "linemeascatalog" in config and spectrum_model in config["linemeascatalog"].keys():
                raise APIException(
                    ErrorCode.INCOHERENT_CONFIG_OPTIONS,
                    "Cannot run lineMeasSolver from catalog when redshiftSolver stage "
                    f"is enabled for model {spectrum_model}",
                )
