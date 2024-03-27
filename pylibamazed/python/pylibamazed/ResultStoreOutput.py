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

import numpy as np
from pylibamazed.AbstractOutput import AbstractOutput, spectrum_model_stages
from pylibamazed.Exception import APIException
from pylibamazed.redshift import CLog, ErrorCode

zlog = CLog.GetInstance()


class ResultStoreOutput(AbstractOutput):
    def __init__(self, result_store, parameters, auto_load=True, extended_results=True):
        AbstractOutput.__init__(self, parameters, extended_results=extended_results)
        self.results_store = result_store
        self.parameters = parameters
        if auto_load:
            self.load_all()

    def _get_attribute_from_result_store(self, object_type, stage, method, data_spec,
                                         rank, band_name=None, obs_id=None):
        operator_result = self._get_operator_result(object_type, stage, method, data_spec, rank)
        band_name_hook = "[band_name]"
        object_type_hook = "[object_type]"
        obs_id_hook = "[obs_id]"
        if object_type_hook in data_spec.OperatorResult_name:
            operator_result_name = data_spec.OperatorResult_name.replace(object_type_hook, "")
        elif band_name_hook in data_spec.OperatorResult_name:
            operator_result_name = data_spec.OperatorResult_name.replace(band_name_hook, "")
        elif obs_id_hook in data_spec.OperatorResult_name:
            operator_result_name = data_spec.OperatorResult_name.replace(obs_id_hook, "")
        else:
            operator_result_name = data_spec.OperatorResult_name
        if "." in operator_result_name:
            o = data_spec.OperatorResult_name.split(".")
            attr = getattr(getattr(operator_result, o[0]), o[1])
        else:
            attr = getattr(operator_result, operator_result_name)

        attr_type = type(attr).__name__
        if attr_type == "TMapFloat64":
            if band_name is not None:
                return attr[band_name]
            return attr[object_type]
        if attr_type == "TMapTFloat64List":
            if obs_id is not None:
                return attr.to_numpy(obs_id)
        elif (attr_type == "TFloat64List"
                or attr_type == "TInt32List"
                or attr_type == "TBoolList"):
            return attr.to_numpy()
        elif attr_type == "TStringList":
            return np.array(attr)
        elif attr_type == "CMask":
            return attr.getMaskList().to_numpy()
        else:
            return attr

    def get_attribute_from_source(
        self,
        object_type,
        stage,
        method,
        dataset,
        attribute,
        rank=None,
        band_name=None,
        obs_id=None
    ):
        rs = self.results_specifications.get_df_by_name(attribute)
        rs = rs[rs["dataset"] == dataset]
        attribute_info = rs.iloc[0]
        return self._get_attribute_from_result_store(
            object_type,
            stage,
            method,
            attribute_info,
            rank=rank,
            band_name=band_name,
            obs_id=obs_id
        )

    def has_attribute_in_source(self, object_type, stage, method, dataset, attribute,
                                rank=None, band_name=None, obs_id=None):
        rs = self.results_specifications.get_df_by_name(attribute)
        rs = rs[rs["dataset"] == dataset]

        attribute_info = rs.iloc[0]
        if type(attribute_info.ResultStore_key) is not str:
            return False
        if rank is not None:
            method = self.parameters.get_redshift_solver_method(object_type)
            if self.results_store.HasCandidateDataset(
                object_type,
                stage,
                method,
                attribute_info.ResultStore_key,
                attribute_info.dataset.replace("<ObsID>", "")
            ):
                operator_result = self._get_operator_result(object_type, stage, method, attribute_info, rank)
            else:
                return False
        else:
            try:
                operator_result = self._get_operator_result(
                    object_type, stage, method, attribute_info, rank=None)
            except Exception:
                return False
        if "[object_type]" in attribute_info.OperatorResult_name:
            or_name = attribute_info.OperatorResult_name.replace("[object_type]", "")
            if hasattr(operator_result, or_name):
                return object_type in getattr(operator_result, or_name)
            else:
                return False
        elif "[band_name]" in attribute_info.OperatorResult_name:
            or_name = attribute_info.OperatorResult_name.replace("[band_name]", "")
            if hasattr(operator_result, or_name):
                o = getattr(operator_result, or_name)
                if o:
                    return band_name in o
                else:
                    return False
            else:
                return False
        elif "[obs_id]" in attribute_info.OperatorResult_name:
            or_name = attribute_info.OperatorResult_name.replace("[obs_id]", "")
            if hasattr(operator_result, or_name):
                o = getattr(operator_result, or_name)
                if o:
                    return obs_id in o
                else:
                    return False
            else:
                return False
        elif "." in attribute_info.OperatorResult_name:
            o = attribute_info.OperatorResult_name.split(".")
            has_o = hasattr(operator_result, o[0])
            if has_o:
                return hasattr(getattr(operator_result, o[0]), o[1])
            else:
                return False
        else:
            return hasattr(operator_result, attribute_info.OperatorResult_name)

    def has_dataset_in_source(self, object_type, stage, method, dataset):
        if object_type:
            return self.results_store.HasDataset(object_type,
                                                 stage,
                                                 method,
                                                 dataset)
        else:
            return self.results_store.HasDataset(dataset,
                                                 dataset,
                                                 dataset,
                                                 "solveResult")

    def has_candidate_dataset_in_source(self, object_type, stage, method, dataset):
        rs = self.results_specifications.get_df_by_dataset(dataset)
        rs_key = rs["ResultStore_key"].unique()[0]

        return self.results_store.HasCandidateDataset(
            object_type,
            stage,
            method,
            rs_key,
            dataset.replace("<ObsID>", "")
        )

    def get_nb_candidates_in_source(self, object_type, stage, method):
        return self.results_store.getNbRedshiftCandidates(object_type, stage, method)

    def _get_operator_result(self, object_type, stage, method, attribute_info, rank=None):
        if attribute_info.level == "root":
            if attribute_info.ResultStore_key in ["context_warningFlag", "warningFlag"]:
                return self.results_store.GetFlagLogResult(attribute_info.dataset,
                                                           attribute_info.dataset,
                                                           attribute_info.dataset,
                                                           attribute_info.ResultStore_key)
            else:
                or_type = self.results_store.GetGlobalResultType(attribute_info.dataset,
                                                                 attribute_info.dataset,
                                                                 attribute_info.dataset,
                                                                 attribute_info.ResultStore_key)
                if or_type == "CClassificationResult":
                    return self.results_store.GetClassificationResult(attribute_info.dataset,
                                                                      attribute_info.dataset,
                                                                      attribute_info.dataset,
                                                                      attribute_info.ResultStore_key)
                else:
                    raise APIException(ErrorCode.OUTPUT_READER_ERROR,
                                       "Unknown OperatorResult type {}".format(str(or_type)))
        elif attribute_info.level == "object" or attribute_info.level == "method":
            or_type = self.results_store.GetGlobalResultType(object_type,
                                                             stage,
                                                             method,
                                                             attribute_info.ResultStore_key)
            getter = getattr(self.results_store, "Get" + or_type[1:])
            return getter(object_type,
                          stage,
                          method,
                          attribute_info.ResultStore_key)
        elif attribute_info.level == "candidate":
            or_type = self.results_store.GetCandidateResultType(object_type,
                                                                stage,
                                                                method,
                                                                attribute_info.ResultStore_key,
                                                                attribute_info.dataset.replace("<ObsID>", ""))
            if or_type == "TLineModelResult":
                firstpass_result = "Firstpass" in attribute_info["name"]
                return self.results_store.GetLineModelResult(object_type,
                                                             stage,
                                                             method,
                                                             attribute_info.ResultStore_key,
                                                             attribute_info.dataset,
                                                             rank, firstpass_result)
            else:
                getter = getattr(self.results_store, "Get" + or_type[1:])
                return getter(object_type,
                              stage,
                              method,
                              attribute_info.ResultStore_key,
                              attribute_info.dataset.replace("<ObsID>", ""),
                              rank)
        else:
            raise APIException(ErrorCode.OUTPUT_READER_ERROR,
                               "Unknown level {}".format(attribute_info.level))

    def store_error(self, amz_exception, object_type, stage):
        full_name = self.get_error_full_name(object_type, stage)
        fn_path = amz_exception.getFileName()
        if len(fn_path.split(os.sep)) > 3:
            fn_path = "/".join(fn_path.split(os.sep)[-3:])
        self.errors[full_name] = dict()
        self.errors[full_name]["code"] = ErrorCode(amz_exception.getErrorCode()).name
        self.errors[full_name]["message"] = amz_exception.getMessage()
        self.errors[full_name]["line"] = amz_exception.getLine()
        self.errors[full_name]["filename"] = fn_path
        self.errors[full_name]["method"] = amz_exception.getMethod()

        # propagate errors to dependant stages
        if not object_type and stage == "init":
            for ot in self.object_types:
                for dependant_stage in spectrum_model_stages:
                    if self.parameters.stage_enabled(ot, dependant_stage):
                        self.store_consequent_error(ot, dependant_stage, stage)
        elif object_type is not None:
            for dependant_stage in spectrum_model_stages[stage]:
                if self.parameters.stage_enabled(object_type, dependant_stage):
                    self.store_consequent_error(object_type, dependant_stage, stage)

    def store_consequent_error(self, object_type, stage, causing_stage):
        full_name = self.get_error_full_name(object_type, stage)
        self.errors[full_name] = dict()
        self.errors[full_name]["code"] = ErrorCode(ErrorCode.STAGE_NOT_RUN_BECAUSE_OF_PREVIOUS_FAILURE).name
        self.errors[full_name]["message"] = f"not run because {causing_stage} failed"
        self.errors[full_name]["line"] = -1
        self.errors[full_name]["filename"] = ""
        self.errors[full_name]["method"] = ""
