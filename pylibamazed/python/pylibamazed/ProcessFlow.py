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
import functools
import os
import os.path
from contextlib import contextmanager, suppress

from pylibamazed.CalibrationLibrary import CalibrationLibrary
from pylibamazed.Exception import APIException, exception_decorator
from pylibamazed.Parameters import Parameters
from pylibamazed.Spectrum import Spectrum

# NB: DO NOT REMOVE - these libs are used in globals
from pylibamazed.redshift import CClassificationSolve  # noqa F401
from pylibamazed.redshift import CLineMatchingSolve  # noqa F401
from pylibamazed.redshift import CLineMeasSolve  # noqa F401
from pylibamazed.redshift import CLineModelSolve  # noqa F401
from pylibamazed.redshift import CTemplateFittingSolve  # noqa F401
from pylibamazed.redshift import CTplCombinationSolve  # noqa F401
from pylibamazed.redshift import AmzException, CFlagWarning, CLog, CProcessFlowContext, ErrorCode, ScopeType
from pylibamazed.Reliability import ReliabilitySolve
from pylibamazed.ResultStoreOutput import ResultStoreOutput
from pylibamazed.ScopeManager import get_scope_spectrum_model, get_scope_stage, push_scope
from pylibamazed.SubType import SubType

zflag = CFlagWarning.GetInstance()
zlog = CLog.GetInstance()


class ProcessFlowException(Exception):
    pass


class ProcessFlow:
    @exception_decorator(logging=True)
    def __init__(self, config, parameters: Parameters):
        _check_config(config)
        self.parameters = parameters
        zlog.LogInfo("Loading all needed calibration files :")
        self.calibration_library = CalibrationLibrary(parameters, config["calibration_dir"])
        self.calibration_library.load_all()
        zlog.LogInfo("Calibration files loaded.")
        self.process_flow_context = CProcessFlowContext.GetInstance()
        self.process_flow_context.reset()
        self.config = config
        if "linemeascatalog" not in self.config:
            self.config["linemeascatalog"] = {}
        else:
            _check_LinemeasValidity(config, parameters)

        self.extended_results = config["extended_results"]

        self.store_flags(name="context_warningFlag")

        # save context warning flag to reinject at each spectrum
        resultStore = self.process_flow_context.GetResultStore()
        self.context_warning_Flag = resultStore.GetFlagLogResult("", "", "", "context_warningFlag")

    def store_exception_handler(func):
        @functools.wraps(func)
        def wrapper(self, rso, *args, **kwargs):
            try:
                return func(self, rso, *args, **kwargs)
            except AmzException as e:
                e.LogError()
                rso.store_error(e, get_scope_spectrum_model(), get_scope_stage())
                raise ProcessFlowException from e
            except Exception as e:
                api_exception = APIException.fromException(e)
                api_exception.LogError()
                rso.store_error(api_exception, get_scope_spectrum_model(), get_scope_stage())
                raise ProcessFlowException from e

        return wrapper

    @exception_decorator(logging=True)
    def run(self, spectrum: Spectrum):
        resultStore = self.process_flow_context.GetResultStore()
        rso = ResultStoreOutput(
            resultStore, self.parameters, auto_load=False, extended_results=self.extended_results
        )
        try:
            with self.store_flags_handler(name="init_warningFlag"):
                self.initialize(rso, spectrum)
        except ProcessFlowException:
            rso.load_root()
            return rso

        # loop on spectrum models (galaxy, star, qso, ...)
        for spectrum_model in self.parameters.get_spectrum_models():
            with push_scope(spectrum_model, ScopeType.SPECTRUMMODEL):
                with suppress(ProcessFlowException):
                    self.process_spectrum_model(rso)

        if self.parameters.is_a_redshift_solver_used():
            with suppress(ProcessFlowException):
                self.run_classification_solver(rso)
                # Running linemeas only on classified model (if any)
                self._run_linemeas_after_classification(rso)

        with suppress(ProcessFlowException):
            self.load_result_store(rso)

        return rso

    def _run_linemeas_after_classification(self, rso: ResultStoreOutput) -> None:
        if self.parameters.get_linemeas_runmode() == "classif":
            classif_model = rso.get_attribute_from_source("root", None, None, "classification", "Type")
            linemeas_method = self.parameters.get_linemeas_method(classif_model)
            if linemeas_method is None:
                return
            with push_scope(classif_model, ScopeType.SPECTRUMMODEL):
                self.run_load_linemeas_params(rso)
                self.run_linemeas_solver(rso, linemeas_method)

    def process_spectrum_model(self, rso):
        spectrum_model = self.scope_spectrum_model

        redshift_solver_method = self.parameters.get_redshift_solver_method(spectrum_model)
        linemeas_method = self.parameters.get_linemeas_method(spectrum_model)
        linemeas_alone = self.config["linemeascatalog"] and spectrum_model in self.config["linemeascatalog"]
        pipe_redshift_linemeas = linemeas_method and not linemeas_alone

        if redshift_solver_method:
            self.run_redshift_solver(rso, redshift_solver_method)

            if self.parameters.is_tplratio_catalog_needed(spectrum_model):
                with suppress(ProcessFlowException):
                    self.run_sub_classification_solver(rso)

            if (
                self.parameters.get_reliability_enabled(spectrum_model)
                and spectrum_model in self.calibration_library.reliability_models
            ):
                with suppress(ProcessFlowException):
                    self.run_reliability_solver(rso)

            if self.parameters.get_linemeas_runmode() == "all" and pipe_redshift_linemeas:
                self.run_load_linemeas_params(rso)
                self.run_linemeas_solver(rso, linemeas_method)

        elif linemeas_method and linemeas_alone:
            self.run_linemeas_solver(rso, linemeas_method)

    @store_exception_handler
    def initialize(self, rso, spectrum: Spectrum):
        zlog.LogInfo("Context initialization")
        self.process_flow_context.reset()
        
        # inject in the clean resultStore the contextFlag
        rso.results_store.StoreScopedGlobalResult("context_warningFlag", self.context_warning_Flag)

        for object_type in self.parameters.get_spectrum_models():
            if object_type in self.calibration_library.line_catalogs:
                for method in self.parameters.get_linemodel_methods(object_type):
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

        if self.config.get("linemeascatalog"):
            self.parameters.load_linemeas_parameters_from_catalog(spectrum.source_id, self.config)
        self.process_flow_context.LoadParameterStore(self.parameters.to_json())
        self.process_flow_context.Init()

    @push_scope("redshiftSolver", ScopeType.STAGE)
    @store_exception_handler
    def run_redshift_solver(self, rso, method):
        self.run_method(method)

    @push_scope("lineMeasSolver", ScopeType.STAGE)
    @store_exception_handler
    def run_linemeas_solver(self, rso, method):
        self.run_method(method)

    @push_scope("linemeas_catalog_load", ScopeType.STAGE)
    @store_exception_handler
    def run_load_linemeas_params(self, rso):
        self.parameters.load_linemeas_parameters_from_result_store(rso, self.scope_spectrum_model)
        self.process_flow_context.LoadParameterStore(self.parameters.to_json())

    @push_scope("reliabilitySolver", ScopeType.STAGE)
    @store_exception_handler
    def run_reliability_solver(self, rso):
        rel = ReliabilitySolve(self.scope_spectrum_model, self.parameters, self.calibration_library)
        rso.object_results[self.scope_spectrum_model]["reliability"] = dict()
        rso.object_results[self.scope_spectrum_model]["reliability"]["Reliability"] = rel.Compute(
            self.process_flow_context
        )

    @push_scope("subClassifSolver", ScopeType.STAGE)
    @store_exception_handler
    def run_sub_classification_solver(self, rso):
        sub_type = SubType(self.scope_spectrum_model, self.parameters, self.calibration_library)
        sub_types = sub_type.Compute(self.process_flow_context)
        rso.object_results[self.scope_spectrum_model]["model_parameters"] = []
        for rank in range(len(sub_types)):
            rso.object_results[self.scope_spectrum_model]["model_parameters"].append(dict())
            rso.object_results[self.scope_spectrum_model]["model_parameters"][rank]["SubType"] = sub_types[
                rank
            ]

    def all_redshift_solver_failed(self, rso):
        answer = True
        for obj in self.parameters.get_spectrum_models():
            if not rso.has_error(obj, "redshiftSolver"):
                answer = False
                break
        return answer

    @push_scope("classification", ScopeType.SPECTRUMMODEL)
    @push_scope("classification", ScopeType.STAGE)
    @store_exception_handler
    def run_classification_solver(self, rso):
        if self.all_redshift_solver_failed(rso):
            raise APIException(
                ErrorCode.NO_CLASSIFICATION, "Classification not run because all redshiftSolver failed"
            )
        self.run_method("classificationSolve")

    @push_scope("load_result_store", ScopeType.STAGE)
    @store_exception_handler
    def load_result_store(self, rso):
        rso.load_all()

    def run_method(self, method):
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
        solver.Compute()

    @property
    def scope_spectrum_model(self):
        return get_scope_spectrum_model()

    @property
    def scope_stage(self):
        return get_scope_stage()

    @contextmanager
    def store_flags_handler(self, name="warningFlag"):
        try:
            yield
        finally:
            self.store_flags(name)

    def store_flags(self, name="warningFlag"):
        resultStore = self.process_flow_context.GetResultStore()
        resultStore.StoreScopedFlagResult(name)
        zflag.resetFlag()


def _check_config(config):
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


def _check_LinemeasValidity(config, parameters):
    if not config["linemeascatalog"]:
        return
    parameters.check_linemeas_validity()
