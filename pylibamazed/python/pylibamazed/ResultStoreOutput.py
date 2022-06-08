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
from pylibamazed.AbstractOutput import AbstractOutput

import numpy as np
import pandas as pd
import resource
from collections import defaultdict
from pylibamazed.redshift import (PC_Get_Float64Array, PC_Get_Int32Array, CLog, ErrorCode)
#from pylibamazed.redshift import PC_Get_Float32Array #TODO ask Ali how to define this latter in the lib
zlog = CLog.GetInstance()
from pylibamazed.Exception import AmazedError,APIException


class ResultStoreOutput(AbstractOutput): 
    def __init__(self, result_store, parameters, auto_load=True, extended_results=True):
        AbstractOutput.__init__(self, parameters, extended_results=extended_results)
        self.results_store = result_store
        self.parameters = parameters
        
        if auto_load:
            self.load_all()

        

    def _get_attribute_from_result_store(self,object_type,method,data_spec,rank):
        operator_result = self._get_operator_result(object_type, method, data_spec,rank)
        if data_spec.dimension == "mono":
            if "[object_type]" in data_spec.OperatorResult_name:
                operator_result_name = data_spec.OperatorResult_name.replace("[object_type]","")
                return getattr(operator_result, operator_result_name)[object_type]
            else:
                return getattr(operator_result, data_spec.OperatorResult_name)
        else:
            if data_spec.hdf5_type == 'f8':
                return PC_Get_Float64Array(getattr(operator_result, data_spec.OperatorResult_name))
            if data_spec.hdf5_type == 'f4':#this should be Float32Array..to be corrected
                return PC_Get_Float64Array(getattr(operator_result, data_spec.OperatorResult_name))
            if data_spec.hdf5_type == 'i':
                return PC_Get_Int32Array(getattr(operator_result, data_spec.OperatorResult_name))

    def get_attribute_from_source(self, object_type, method, dataset, attribute,  rank=None):
        rs = self.results_specifications
        rs = rs[rs["name"] == attribute]
        rs = rs[rs["dataset"] == dataset]
        attribute_info = rs.iloc[0]
        
        return self._get_attribute_from_result_store(object_type,
                                                     method,
                                                     attribute_info,
                                                     rank)

    def has_attribute_in_source(self,object_type,method, dataset, attribute, rank=None):
        rs = self.results_specifications
        rs = rs[rs["name"] == attribute]
        rs = rs[rs["dataset"] == dataset]
        
        attribute_info = rs.iloc[0]
        
        if rank is not None:
            method = self.parameters.get_solve_method(object_type)
            if self.results_store.HasCandidateDataset(object_type,
                                                      method,
                                                      attribute_info.ResultStore_key,
                                                      attribute_info.dataset):
                operator_result = self._get_operator_result(object_type, method,attribute_info, rank)
            else:
                return False
        else:
            try:
                operator_result = self._get_operator_result(object_type, method, attribute_info, rank=None)
            except Exception as e:
                return False
        if "[object_type]" in attribute_info.OperatorResult_name:
            or_name = attribute_info.OperatorResult_name.replace("[object_type]", "")
            if hasattr(operator_result, or_name):
                return object_type in getattr(operator_result, or_name)
            else:
                return False
        else:
            return hasattr(operator_result, attribute_info.OperatorResult_name)

    def has_dataset_in_source(self, object_type, method, dataset):
        return self.results_store.HasDataset(object_type,
                                             method,
                                             dataset)

    def has_candidate_dataset_in_source(self, object_type, method, dataset):
        rs = self.results_specifications
        rs = rs[ rs.dataset == dataset]
        rs_key = rs["ResultStore_key"].unique()[0]

        return self.results_store.HasCandidateDataset(object_type,
                                                      method,
                                                      rs_key,
                                                      dataset)
    
    def get_nb_candidates_in_source(self, object_type, method):
        return self.results_store.getNbRedshiftCandidates(object_type, method)

    def _get_operator_result(self, object_type, method, attribute_info, rank=None):
        if attribute_info.level == "root":
            if attribute_info.ResultStore_key == "context_warningFlag":
                return self.results_store.GetFlagResult(attribute_info.dataset,
                                                        attribute_info.dataset,
                                                        attribute_info.ResultStore_key)
            else :
                or_type = self.results_store.GetGlobalResultType(attribute_info.dataset,
                                                                attribute_info.dataset,
                                                                attribute_info.ResultStore_key)
                if or_type == "CClassificationResult":
                    return self.results_store.GetClassificationResult(attribute_info.dataset,
                                                                    attribute_info.dataset,
                                                                    attribute_info.ResultStore_key)
                else:
                    raise APIException(ErrorCode.OutputReaderError,"Unknown OperatorResult type {}".format(str(or_type)))
        elif attribute_info.level == "method":
            or_type = self.results_store.GetGlobalResultType(object_type,
                                                             method,
                                                             attribute_info.ResultStore_key)
            if or_type == "CFlagLogResult":
                return  self.results_store.GetFlagResult(object_type,
                                                        method,
                                                        attribute_info.ResultStore_key)
            else:
                raise Exception("Unknown OperatorResult type " + or_type)                                     
        elif attribute_info.level == "object":
            or_type = self.results_store.GetGlobalResultType(object_type,
                                                             method,
                                                             attribute_info.ResultStore_key)

            if or_type == "CPdfMargZLogResult":
                return self.results_store.GetPdfMargZLogResult(object_type,
                                                               method,
                                                               attribute_info.ResultStore_key)
            elif or_type == "CLineModelSolution":
                return self.results_store.GetLineModelSolution(object_type,
                                                               method,
                                                               attribute_info.ResultStore_key)
            elif or_type == "CModelSpectrumResult":
                return self.results_store.GetModelSpectrumResult(object_type,
                                                                 method,
                                                                 attribute_info.ResultStore_key)
            else:
                raise APIException(ErrorCode.OutputReaderError,"Unknown OperatorResult type {}".format(str(or_type)))
        elif attribute_info.level == "candidate":
            or_type = self.results_store.GetCandidateResultType(object_type,
                                                                method,
                                                                attribute_info.ResultStore_key,
                                                                attribute_info.dataset)
            if or_type == "TLineModelResult":
                firstpass_result = "Firstpass" in attribute_info.name
                return self.results_store.GetLineModelResult(object_type,
                                                             method,
                                                             attribute_info.ResultStore_key,
                                                             rank, firstpass_result)
            elif or_type == "TExtremaResult":
                return self.results_store.GetExtremaResult(object_type,
                                                           method,
                                                           attribute_info.ResultStore_key,
                                                           rank)
            elif or_type == "CModelSpectrumResult":
                return self.results_store.GetModelSpectrumResult(object_type,
                                                                 method,
                                                                 attribute_info.ResultStore_key,
                                                                 rank)
            elif or_type == "CSpectraFluxResult":
                return self.results_store.GetSpectraFluxResult(object_type,
                                                               method,
                                                               attribute_info.ResultStore_key,
                                                               rank)
            elif or_type == "CModelFittingResult":
                return self.results_store.GetModelFittingResult(object_type,
                                                                method,
                                                                attribute_info.ResultStore_key,
                                                                rank)
            elif or_type == "TTplCombinationResult":
                return self.results_store.GetTplCombinationResult(object_type,
                                                                method,
                                                                attribute_info.ResultStore_key,
                                                                rank)
            elif or_type == "CLineModelSolution":
                return self.results_store.GetLineModelSolution(object_type,
                                                               method,
                                                               attribute_info.ResultStore_key,
                                                               rank)
            else:
                raise APIException(ErrorCode.OutputReaderError,"Unknown OperatorResult type {}".format(or_type))

