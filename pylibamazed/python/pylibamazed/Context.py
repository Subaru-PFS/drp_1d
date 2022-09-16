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
import os.path
import sys,traceback
from pylibamazed.CalibrationLibrary import CalibrationLibrary
from pylibamazed.ResultStoreOutput import ResultStoreOutput
from pylibamazed.Reliability import Reliability
from pylibamazed.Parameters import Parameters
from pylibamazed.SubType import SubType
from pylibamazed.redshift import (CProcessFlowContext,
                                  CLineModelSolve,
                                  CTemplateFittingSolve,
                                  CTplcombinationSolve,
                                  CClassificationSolve,
                                  CReliabilitySolve,
                                  CLineMeasSolve,
                                  CLineMatchingSolve,
                                  CFlagWarning,CLog,
                                  GlobalException, ErrorCode)
import pandas as pd
import json
import os

from pylibamazed.Exception import (AmazedError,
                                   AmazedErrorFromGlobalException,
                                   APIException)
zflag = CFlagWarning.GetInstance()
zlog = CLog.GetInstance()


class Context:

    def __init__(self, config, parameters):
        try:
            _check_config(config)
            self.calibration_library = CalibrationLibrary(parameters,
                                                          config["calibration_dir"])
            self.calibration_library.load_all()
            self.process_flow_context = None
            self.config = config
            if "linemeascatalog" not in self.config:
                self.config["linemeascatalog"] = {}
            else:
                _check_LinemeasValidity(config, parameters)
            self.parameters = Parameters(parameters,config)

            self.extended_results = config["extended_results"]
        except GlobalException as e:
            raise AmazedErrorFromGlobalException(e)
        except APIException as e:
            raise AmazedError(e.errCode,e.message)
        except Exception as e:
            raise AmazedError(ErrorCode.PYTHON_API_ERROR, str(e))
            
    def init_context(self):
        self.process_flow_context = CProcessFlowContext.GetInstance()
        self.process_flow_context.reset()
        for object_type in self.parameters.get_objects():
            if object_type in self.calibration_library.line_catalogs:
                for method in self.parameters.get_linemodel_methods(object_type):
                    self.process_flow_context.setLineCatalog(object_type, method,
                                                             self.calibration_library.line_catalogs[object_type][method])
            if object_type in self.calibration_library.line_ratio_catalog_lists:
                self.process_flow_context.setLineRatioCatalogCatalog(object_type,
                                                                     self.calibration_library.line_ratio_catalog_lists[object_type])

        self.process_flow_context.setTemplateCatalog(self.calibration_library.templates_catalogs["all"])
        self.process_flow_context.setPhotBandCatalog(self.calibration_library.photometric_bands)
        self.process_flow_context.setfluxCorrectionMeiksin(self.calibration_library.meiksin)
        self.process_flow_context.setfluxCorrectionCalzetti(self.calibration_library.calzetti)

    def run(self, spectrum_reader):
        resultStore = CProcessFlowContext.GetInstance().GetResultStore()
        rso = None
        context_warningFlagRecorded=False
        try:
            self.init_context()
            spectrum_reader.init()
            if self.config.get("linemeascatalog"):
                self.parameters.load_linemeas_parameters_from_catalog(spectrum_reader.source_id)
            self.process_flow_context.LoadParameterStore(self.parameters.get_json())
            self.process_flow_context.Init()

            #store flag in root object
            
            rso = ResultStoreOutput(self.process_flow_context.GetResultStore(),
                                    self.parameters,
                                    auto_load = False,
                                    extended_results = self.extended_results)
            
            resultStore.StoreFlagResult( "context_warningFlag", zflag.getBitMask())
            zflag.resetFlag()
            context_warningFlagRecorded = True
        except GlobalException as e:
            if not context_warningFlagRecorded:
                resultStore.StoreFlagResult( "context_warningFlag", zflag.getBitMask())
                zflag.resetFlag()
            rso.store_error(AmazedErrorFromGlobalException(e),None,"init")
            return rso
        except APIException as e:
            if not context_warningFlagRecorded:
                resultStore.StoreFlagResult( "context_warningFlag", zflag.getBitMask())
                zflag.resetFlag()
            rso.store_error(AmazedError(e.errCode,e.message), None, "init")
            return rso
        except Exception as e:
            if not context_warningFlagRecorded:
                resultStore.StoreFlagResult( "context_warningFlag", zflag.getBitMask())
                zflag.resetFlag()
            rso.store_error(AmazedError(ErrorCode.PYTHON_API_ERROR, str(e)),None,"init")
            return rso

        enable_classification = False
        reliabilities = dict()
        sub_types = dict()
        for object_type in self.parameters.get_objects():
            linemeas_params_from_solver = not self.config["linemeascatalog"] and not object_type in self.config["linemeascatalog"]
            method = self.parameters.get_solve_method(object_type)
            solver_success = True
            if method:
                try:
                    self.run_method(object_type, method)
                    enable_classification = True
                except GlobalException as e:
                    rso.store_error(AmazedErrorFromGlobalException(e),object_type,"redshift_solver")
                    # if 'warningFlag' not in rso.object_results[object_type]:
                    #     rso.object_results[object_type]['warningFlag'] = dict()
                    # if not rso.object_results[object_type]['warningFlag']:
                    #     rso.object_results[object_type]['warningFlag'][method+'WarningFlags'] = 0
                    linemeas_params_from_solver = False
                    solver_success = False

            linemeas_parameters_loaded=True
            if self.parameters.get_linemeas_method(object_type):
                if linemeas_params_from_solver :
                    try:
                        output = ResultStoreOutput(resultStore,
                                                   self.parameters,
                                                   auto_load=False,
                                                   extended_results = self.extended_results)

                        self.parameters.load_linemeas_parameters_from_result_store(output,object_type)
                        self.process_flow_context.LoadParameterStore(self.parameters.get_json())
                    except APIException as e:
                        rso.store_error(AmazedError(e.errCode,e.message), object_type,"linemeas_catalog_load")
                        linemeas_parameters_loaded = False
                try:
                    if linemeas_parameters_loaded:
                        self.run_method(object_type, self.parameters.get_linemeas_method(object_type))
                except GlobalException as e:
                    rso.store_error(AmazedErrorFromGlobalException(e),object_type,"linemeas_solver")


            if self.parameters.reliability_enabled(object_type) and object_type in self.calibration_library.reliability_models and solver_success:
                try:
                    rel = Reliability(object_type, self.parameters,self.calibration_library)
                    reliabilities[object_type] = rel.Compute(self.process_flow_context)
                except APIException as e:
                    rso.store_error(AmazedError(e.errCode,e.message),object_type, "reliability_solver" )

            if self.parameters.lineratio_catalog_enabled(object_type) and solver_success:
                try:
                    sub_classif = SubType(object_type,
                                          self.parameters,
                                          self.calibration_library)
                    sub_types[object_type] = sub_classif.Compute(self.process_flow_context)
                except APIException as e:
                    rso.store_error(AmazedError(e.errCode,e.message), object_type,"sub_classif_solver")


        if enable_classification:
            try:
                self.run_method("classification", "ClassificationSolve")
            except GlobalException as e:
                rso.store_error(AmazedErrorFromGlobalException(e),None,"classification")

        try:
            rso.load_all()
            for object_type in reliabilities.keys():
                rso.object_results[object_type]['reliability'] = dict()
                rso.object_results[object_type]['reliability']['Reliability'] = reliabilities[object_type]
            for object_type in sub_types.keys():
                for rank in range(len(sub_types[object_type])):
                    rso.object_results[object_type]['model_parameters'][rank]['SubType'] = sub_types[object_type][rank]
        except APIException as e:
            rso.store_error(AmazedError(e.errCode,e.message),None,
                            "result_store_fill")
        except Exception as e:
            rso.store_error(AmazedError(ErrorCode.OUTPUT_READER_ERROR,"{}".format(e)),None,
                            "result_store_fill")

        return rso

    def run_method(self, object_type, method):

        if "C" + method not in globals():
            raise APIException(ErrorCode.INVALID_PARAMETER,"Unknown method {}".format(method))
        solver_method = globals()["C" + method]
        solver = solver_method(self.process_flow_context.m_ScopeStack,
                               object_type)
        solver.Compute()


def _check_config(config):
    if "calibration_dir" not in config:
        raise APIException(ErrorCode.MISSING_CONFIG_OPTION,"Config must contain 'calibration_dir' key")
    if not os.path.exists(config["calibration_dir"]):
        raise APIException(ErrorCode.INVALID_DIRECTORY,"Calibration directory {} does not exist".format(config["calibration_dir"]))
    if "linemeascatalog" in config:
        if "linemeas_catalog_columns" not in config:
            raise APIException(ErrorCode.MISSING_CONFIG_OPTION,"Missing linemeas_catalog_columns key in linemeascatalog config-option")
        for object_type in config["linemeascatalog"].keys():
            if object_type not in config["linemeas_catalog_columns"]:
                raise APIException(ErrorCode.INCOHERENT_CONFIG_OPTIONS,"Missing category {} in linemeas_catalog_columns ".format(object_type))

            for attr in ["Redshift", "VelocityAbsorption", "VelocityEmission"]:
                if attr not in config["linemeas_catalog_columns"][object_type]:
                    raise APIException(ErrorCode.ATTRIBUTE_NOT_SUPPORTED,"Not supported Attribute {0} in Config['linemeas_catalog_columns'][{1}]".format(object_type, attr))


def _check_LinemeasValidity(config, parameters):
    if not config["linemeascatalog"]:
        return
    for object_type in parameters["objects"]:
        method = parameters[object_type]["method"]
        if method == "LineModelSolve":
            if "linemeas_method" in parameters[object_type] and parameters[object_type]["linemeas_method"]:
                raise APIException(ErrorCode.INCOHERENT_CONFIG_OPTION, "Cannot run LineMeasSolve from catalog when sequencial processing is selected simultaneously.")
