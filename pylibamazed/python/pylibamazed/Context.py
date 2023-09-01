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

from pylibamazed.CalibrationLibrary import CalibrationLibrary
from pylibamazed.Exception import (AmazedError, AmazedErrorFromGlobalException,
                                   APIException)
from pylibamazed.Parameters import Parameters
# NB: DO NOT REMOVE - these libs are used in globals
from pylibamazed.redshift import CClassificationSolve  # noqa F401
from pylibamazed.redshift import CLineMatchingSolve  # noqa F401
from pylibamazed.redshift import CLineMeasSolve  # noqa F401
from pylibamazed.redshift import CLineModelSolve  # noqa F401
from pylibamazed.redshift import CReliabilitySolve  # noqa F401
from pylibamazed.redshift import CTemplateFittingSolve  # noqa F401
from pylibamazed.redshift import CTplcombinationSolve  # noqa F401
from pylibamazed.redshift import (CFlagWarning, CLog, CProcessFlowContext,
                                  ErrorCode, GlobalException)
from pylibamazed.Reliability import Reliability
from pylibamazed.ResultStoreOutput import ResultStoreOutput
from pylibamazed.SubType import SubType

zflag = CFlagWarning.GetInstance()
zlog = CLog.GetInstance()


class Context:

    def __init__(self, config, parameters: Parameters):
        try:
            _check_config(config)
            self.parameters = parameters
            self.calibration_library = CalibrationLibrary(parameters,
                                                          config["calibration_dir"])
            self.calibration_library.load_all()
            self.process_flow_context = None
            self.config = config
            if "linemeascatalog" not in self.config:
                self.config["linemeascatalog"] = {}
            else:
                _check_LinemeasValidity(config, parameters)

            self.extended_results = config["extended_results"]
        except GlobalException as e:
            raise AmazedErrorFromGlobalException(e)
        except APIException as e:
            raise AmazedError(e.errCode, e.message)
        except Exception as e:
            raise AmazedError(ErrorCode.PYTHON_API_ERROR, str(e))

    def init_context(self):
        self.process_flow_context = CProcessFlowContext.GetInstance()
        self.process_flow_context.reset()
        for object_type in self.parameters.get_objects():
            if object_type in self.calibration_library.line_catalogs:
                for method in self.parameters.get_linemodel_methods(object_type):
                    self.process_flow_context.setLineCatalog(
                        object_type,
                        method,
                        self.calibration_library.line_catalogs[object_type][method]
                    )
            if object_type in self.calibration_library.line_ratio_catalog_lists:
                self.process_flow_context.setLineRatioCatalogCatalog(
                    object_type,
                    self.calibration_library.line_ratio_catalog_lists[object_type]
                )

        self.process_flow_context.setTemplateCatalog(self.calibration_library.templates_catalogs["all"])
        self.process_flow_context.setPhotBandCatalog(self.calibration_library.photometric_bands)
        self.process_flow_context.setfluxCorrectionMeiksin(self.calibration_library.meiksin)
        self.process_flow_context.setfluxCorrectionCalzetti(self.calibration_library.calzetti)

    def run(self, spectrum_reader):
        try:
            resultStore = CProcessFlowContext.GetInstance().GetResultStore()
            rso = ResultStoreOutput(resultStore,
                                    self.parameters,
                                    auto_load=False,
                                    extended_results=self.extended_results)
        except Exception as e:
            rso.store_error(AmazedError(ErrorCode.PYTHON_API_ERROR, str(e)), None, "init")
        context_warningFlagRecorded = False
        try:
            self.init_context()
            spectrum_reader.init()
            if self.config.get("linemeascatalog"):
                self.parameters.load_linemeas_parameters_from_catalog(spectrum_reader.source_id, self.config)
            self.process_flow_context.LoadParameterStore(self.parameters.get_json())
            self.process_flow_context.Init()
            # store flag in root object
            resultStore.StoreFlagResult("context_warningFlag", zflag.getBitMask())
            context_warningFlagRecorded = True
            zflag.resetFlag()
        except GlobalException as e:
            rso.store_error(AmazedErrorFromGlobalException(e), None, "init")
            rso.load_root()
            return rso
        except APIException as e:
            if not context_warningFlagRecorded:
                resultStore.StoreFlagResult("context_warningFlag", zflag.getBitMask())
                zflag.resetFlag()
            rso.store_error(AmazedError(e.errCode, e.message), None, "init")
            rso.load_root()
            return rso
        except Exception as e:
            if not context_warningFlagRecorded:
                resultStore.StoreFlagResult("context_warningFlag", zflag.getBitMask())
                zflag.resetFlag()
            rso.store_error(AmazedError(ErrorCode.PYTHON_API_ERROR, str(e)), None, "init")
            rso.load_root()
            return rso

        for object_type in self.parameters.get_objects():
            linemeas_params_from_solver = not self.config["linemeascatalog"] \
                and object_type not in self.config["linemeascatalog"]
            method = self.parameters.get_solve_method(object_type)
            if method:
                self.run_redshift_solver(rso, object_type, "redshift_solver")

            if self.parameters.get_linemeas_method(object_type):
                if linemeas_params_from_solver and not rso.has_error(object_type, "redshift_solver"):
                    self.run_load_linemeas_params(rso, object_type, "linemeas_catalog_load")
                if not rso.has_error(object_type, "linemeas_catalog_load"):
                    self.run_linemeas_solver(rso, object_type, "linemeas_solver")

            if self.parameters.is_tplratio_catalog_needed(object_type) \
                    and not rso.has_error(object_type, "redshift_solver"):
                self.run_sub_classification(rso, object_type, "sub_classif_solver")

            if self.parameters.reliability_enabled(object_type) \
                and object_type in self.calibration_library.reliability_models \
                    and not rso.has_error(object_type, "redshift_solver"):
                self.run_reliability(rso, object_type, "reliability_solver")

        self.run_classification(rso, None, "classification")

        self.load_result_store(rso, None, "load_result_store")

        return rso

    def run_method_exception_handler(func):
        @functools.wraps(func)
        def inner_function(self, *args, **kwargs):
            try:
                rso = args[0]
                object_type = args[1]
                stage = args[2]
                func(self, *args, **kwargs)
            except GlobalException as e:
                rso.store_error(AmazedErrorFromGlobalException(e), object_type, stage)
            except APIException as e:
                rso.store_error(AmazedError(e.errCode, e.message), object_type, stage)
            except Exception as e:
                rso.store_error(AmazedError(ErrorCode.PYTHON_API_ERROR, str(e)),
                                object_type,
                                stage)
        return inner_function

    @run_method_exception_handler
    def run_redshift_solver(self, rso, object_type, stage):
        self.run_method(object_type,
                        self.parameters.get_solve_method(object_type))

    @run_method_exception_handler
    def run_linemeas_solver(self, rso, object_type, stage):
        self.run_method(object_type,
                        self.parameters.get_linemeas_method(object_type))

    @run_method_exception_handler
    def run_load_linemeas_params(self, rso, object_type, stage):
        self.parameters.load_linemeas_parameters_from_result_store(rso, object_type)
        self.process_flow_context.LoadParameterStore(self.parameters.get_json())

    @run_method_exception_handler
    def run_reliability(self, rso, object_type, stage):
        rel = Reliability(object_type, self.parameters, self.calibration_library)
        rso.object_results[object_type]['reliability'] = dict()
        rso.object_results[object_type]['reliability']['Reliability'] = rel.Compute(self.process_flow_context)

    @run_method_exception_handler
    def run_sub_classification(self, rso, object_type, stage):
        sub_type = SubType(object_type,
                           self.parameters,
                           self.calibration_library)
        sub_types = sub_type.Compute(self.process_flow_context)
        rso.object_results[object_type]['model_parameters'] = []
        for rank in range(len(sub_types)):
            rso.object_results[object_type]['model_parameters'].append(dict())
            rso.object_results[object_type]['model_parameters'][rank]['SubType'] = sub_types[rank]

    @run_method_exception_handler
    def run_classification(self, rso, object_type, stage):
        enable_classification = False
        for object_type in self.parameters.get_objects():
            if self.parameters.get_solve_method(object_type) \
                    and not rso.has_error(object_type, "redshift_solver"):
                enable_classification = True
                break

        if enable_classification:
            self.run_method("classification", "ClassificationSolve")
        else:
            raise APIException(ErrorCode.NO_CLASSIFICATION,
                               "Classification not run because all redshift_solver failed")

    @run_method_exception_handler
    def load_result_store(self, rso, object_type, stage):
        rso.load_all()

    def run_method(self, object_type, method):
        if "C" + method not in globals():
            raise APIException(ErrorCode.INVALID_PARAMETER, "Unknown method {}".format(method))
        solver_method = globals()["C" + method]
        solver = solver_method(self.process_flow_context.m_ScopeStack,
                               object_type)
        solver.Compute()


def _check_config(config):
    if "calibration_dir" not in config:
        raise APIException(ErrorCode.MISSING_CONFIG_OPTION, "Config must contain 'calibration_dir' key")
    if not os.path.exists(config["calibration_dir"]):
        raise APIException(ErrorCode.INVALID_DIRECTORY,
                           "Calibration directory {} does not exist".format(config["calibration_dir"]))
    if "linemeascatalog" in config:
        if "linemeas_catalog_columns" not in config:
            raise APIException(ErrorCode.MISSING_CONFIG_OPTION,
                               "Missing linemeas_catalog_columns key in linemeascatalog config-option")
        for object_type in config["linemeascatalog"].keys():
            if object_type not in config["linemeas_catalog_columns"]:
                raise APIException(ErrorCode.INCOHERENT_CONFIG_OPTIONS,
                                   "Missing category {} in linemeas_catalog_columns ".format(object_type))

            for attr in ["Redshift", "VelocityAbsorption", "VelocityEmission"]:
                if attr not in config["linemeas_catalog_columns"][object_type]:
                    raise APIException(
                        ErrorCode.ATTRIBUTE_NOT_SUPPORTED,
                        f"Not supported Attribute {object_type} in Config['linemeas_catalog_columns'][{attr}]"
                    )


def _check_LinemeasValidity(config, parameters):
    if not config["linemeascatalog"]:
        return
    parameters.check_lineameas_validity()
