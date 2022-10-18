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
from pylibamazed.redshift import (PC_Get_Float64Array,
                                  PC_Get_Float32Array,
                                  PC_Get_Int32Array,
                                  CLog, ErrorCode)
from pylibamazed.Exception import AmazedError,APIException


zlog = CLog.GetInstance()


class ResultStoreOutput(AbstractOutput): 
    def __init__(self, result_store, parameters, auto_load=True, extended_results=True):
        AbstractOutput.__init__(self, parameters, extended_results=extended_results)
        self.results_store = result_store
        self.parameters = parameters
        if auto_load:
            self.load_all()
        

    def _get_attribute_from_result_store(self,object_type,method,data_spec,rank):
        operator_result = self._get_operator_result(object_type, method, data_spec,rank)
        if "[object_type]" in data_spec.OperatorResult_name:
            operator_result_name = data_spec.OperatorResult_name.replace("[object_type]","")
        else:
            operator_result_name = data_spec.OperatorResult_name
        attr = getattr(operator_result, operator_result_name)
        attr_type = type(attr).__name__
        if attr_type == "TMapFloat64":
            return attr[object_type]
        elif attr_type == "TFloat64List":
            return PC_Get_Float64Array(attr)
        elif attr_type == "TFloat32List":
            return PC_Get_Float32Array(attr)
        elif attr_type == "TInt32List":
            return PC_Get_Int32Array(attr)
        else:
            return attr

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
        if type(attribute_info.ResultStore_key) != str:
            return False
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
        if object_type:
            return self.results_store.HasDataset(object_type,
                                                 method,
                                                 dataset)
        else:
            return self.results_store.HasDataset(dataset,
                                                 dataset,
                                                 "solveResult")

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
                return self.results_store.GetFlagLogResult(attribute_info.dataset,
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
                    raise APIException(ErrorCode.OUTPUT_READER_ERROR,"Unknown OperatorResult type {}".format(str(or_type)))
        elif attribute_info.level == "object" or attribute_info.level == "method":
            or_type = self.results_store.GetGlobalResultType(object_type,
                                                             method,
                                                             attribute_info.ResultStore_key)
            getter = getattr(self.results_store,"Get"+or_type[1:])
            return getter(object_type,
                          method,
                          attribute_info.ResultStore_key)
        elif attribute_info.level == "candidate":
            or_type = self.results_store.GetCandidateResultType(object_type,
                                                                method,
                                                                attribute_info.ResultStore_key,
                                                                attribute_info.dataset)
            if or_type == "TLineModelResult":
                firstpass_result = "Firstpass" in attribute_info["name"]
                return self.results_store.GetLineModelResult(object_type,
                                                             method,
                                                             attribute_info.ResultStore_key,
                                                             rank, firstpass_result)
            else:
                getter = getattr(self.results_store,"Get"+or_type[1:])
                return getter(object_type,
                              method,
                              attribute_info.ResultStore_key,
                              rank)
        else:
            raise APIException(ErrorCode.OUTPUT_READER_ERROR,
                               "Unknown level {}".format(attribute_info.level))

