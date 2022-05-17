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
from pylibamazed.r_specifications import rspecifications
import numpy as np
import pandas as pd
import resource
from collections import defaultdict
from pylibamazed.redshift import PC_Get_Float64Array, PC_Get_Int32Array, CLog

zlog = CLog.GetInstance()


class ResultStoreOutput(AbstractOutput): 
    def __init__(self, result_store, parameters, results_specifications=rspecifications, auto_load=True, extended_results=True):
        AbstractOutput.__init__(self)
        self.results_store = result_store
        self.parameters = parameters
        self.operator_results = dict()
        self.extended_results = extended_results 
        self.results_specifications = pd.read_csv(results_specifications,
                                     sep='\t',
                                     dtype={'format': object}
                                     )
        
        self.object_types = self.parameters["objects"]
           
        for object_type in self.object_types:
            self.object_results[object_type] = dict()
            self.object_dataframes[object_type] = dict()
            self.operator_results[object_type] = dict()
        if auto_load:
            self.load_all()

    def load_root(self):
        level= "root"
        rs, root_datasets = self.filter_datasets(level)
        for ds in root_datasets:
            ds_attributes = self.filter_dataset_attributes(ds)
            self.root_results[ds]=dict()
            if ds == "context_warningFlag" :
                for index, ds_row in ds_attributes.iterrows():
                    self.root_results[ds][ds_row["hdf5_name"]] = self._get_attribute_from_result_store(ds_row)
            elif self.results_store.HasDataset(ds,ds,"solveResult"):
                for index, ds_row in ds_attributes.iterrows():
                    if "<" in ds_row["hdf5_name"]:
                        for object_type in self.parameters["objects"]:
                            attr_name = ds_row["hdf5_name"].replace("<ObjectType>", object_type)
                            self.root_results[ds][attr_name] = self._get_attribute_from_result_store(ds_row,
                                                                                                object_type=object_type)
                    else:
                        self.root_results[ds][ds_row["hdf5_name"]] = self._get_attribute_from_result_store(ds_row)

    def load_method_level(self, object_type):
        level = "method"
        rs, object_datasets = self.filter_datasets(level)
        for ds in object_datasets:
            methods = self.get_solve_methods(object_type)
            self.object_results[object_type][ds] = dict()
            for method in methods:
                if self.results_store.HasDataset(object_type,
                                                 method,
                                                 ds):
                    ds_attributes = self.filter_dataset_attributes(ds)                    
                    for index, ds_row in ds_attributes.iterrows():
                        attr_name = ds_row["hdf5_name"]
                        if "<MethodType>" in ds_row["hdf5_name"]:
                            attr_name = ds_row["hdf5_name"].replace("<MethodType>", method)
                        methodIsLineMeasSolve = (method=="LineMeasSolve")
                        attr = self._get_attribute_from_result_store(ds_row,object_type, None, methodIsLineMeasSolve)
                        self.object_results[object_type][ds][attr_name] = attr

    def load_object_level(self, object_type):
        level = "object"
        rs, object_datasets = self.filter_datasets(level)
        for ds in object_datasets:
            methods = self.get_solve_methods(object_type)
            if self.get_solve_method(object_type):
                method = self.get_solve_method(object_type)
            else:
                method = self.parameters[object_type]["linemeas_method"]
            for method in methods:
                if self.results_store.HasDataset(object_type,
                                                 method,
                                                 ds):
                    ds_attributes = self.filter_dataset_attributes(ds)
                    dimension = ds_attributes["dimension"].iat[0]
                    if dimension == "multi":
                        self.object_results[object_type][ds] = dict()
                        self.object_dataframes[object_type][ds] = pd.DataFrame()
                        for index, ds_row in ds_attributes.iterrows():
                            attr = self._get_attribute_from_result_store(ds_row,object_type)
                            self.object_results[object_type][ds][ds_row["hdf5_name"]] = attr
                            self.object_dataframes[object_type][ds][ds_row["hdf5_name"]] = attr
                    else:
                        self.object_results[object_type][ds] = dict()
                        for index, ds_row in ds_attributes.iterrows():
                            methodIsLineMeasSolve = (method=="LineMeasSolve")
                            attr = self._get_attribute_from_result_store(ds_row,object_type, None, methodIsLineMeasSolve)
                            self.object_results[object_type][ds][ds_row["hdf5_name"]] = attr
                    
    def load_candidate_level(self, object_type):
        if not self.parameters[object_type]["method"]:
            return
        level = "candidate"
        rs, candidate_datasets = self.filter_datasets(level)
        for ds in candidate_datasets:
            ds_attributes = self.filter_dataset_attributes(ds)
            rs_key = ds_attributes["ResultStore_key"].unique()[0]
            if not self.results_store.HasCandidateDataset(object_type,
                                                          self.get_solve_method(object_type),
                                                          rs_key,
                                                          ds):
                continue

            nb_candidates = self.results_store.getNbRedshiftCandidates(object_type,
                                                                       self.get_solve_method(object_type))
            candidates = []
            candidates_df = []  # useless if dataset attributes dimension is not multi
            dimension = None
            for rank in range(nb_candidates):
                candidates.append(dict())
                candidates_df.append(pd.DataFrame())
                for index, ds_row in ds_attributes.iterrows():
                    #            for attr in list(ds_attributes["hdf5_name"]):
                    if self.has_attribute_in_result_store(ds_row,object_type, rank):
                        attr = self._get_attribute_from_result_store(ds_row, object_type, rank)
                        candidates[rank][ds_row["hdf5_name"]] = attr
                        dimension = ds_row["dimension"]
                        if dimension == "multi":
                            candidates_df[rank][ds_row["hdf5_name"]] = candidates[rank][ds_row["hdf5_name"]]
            self.object_results[object_type][ds] = candidates
            if dimension == "multi":
                self.object_dataframes[object_type][ds] = candidates_df
            else:
                res = defaultdict(list)
                {res[key].append(sub[key]) for sub in candidates for key in sub}
                self.object_dataframes[object_type][ds] = pd.DataFrame(res)
                self.object_dataframes[object_type][ds]["Rank"] = range(nb_candidates)
    
    def filter_datasets(self, level):
        rs = self.results_specifications
        #filter by level
        rs = rs[rs["level"] == level]
        all_datasets = list(rs["hdf5_dataset"].unique())
        
        #filter by extended_results
        if self.extended_results:
            return rs, all_datasets

        #a dataset is considered as debug if all its elements have debug = True
        filtered_datasets = []
        for ds in all_datasets:
            ds_attributes = rs[rs["hdf5_dataset"]==ds]
            debug = all(ds_row["debug"] == True for index, ds_row in ds_attributes.iterrows())
            if not debug:
                filtered_datasets.append(ds) 
        
        return rs, filtered_datasets

    def filter_dataset_attributes(self, ds_name):  
        rs = self.results_specifications 
        ds_attributes = rs[rs["hdf5_dataset"]==ds_name]   
        #filter ds_attributes by debug column
        if self.extended_results:
            return ds_attributes
        filtered_df = ds_attributes[ds_attributes["debug"]==self.extended_results]     
        return filtered_df

    def write_hdf5_root(self, hdf5_spectrum_node):
        level = "root"
        rs, root_datasets = self.filter_datasets(level)
                
        for ds in root_datasets:
            if ds in self.root_results:
                if not self.root_results[ds]:
                    continue
                ds_attributes = self.filter_dataset_attributes(ds)
                dsg = hdf5_spectrum_node.create_group(ds)
                for index, ds_row in ds_attributes.iterrows():
                    if "<" in ds_row["hdf5_name"]:
                        for object_type in self.parameters["objects"]:
                            attr_name = ds_row["hdf5_name"].replace("<ObjectType>", object_type)
                            dsg.attrs[attr_name] = self.root_results[ds][attr_name]
                    else:
                        if ds_row["hdf5_name"] in self.root_results[ds]:
                            # we know we only have dimension = mono here
                            dsg.attrs[ds_row["hdf5_name"]] = self.root_results[ds][ds_row["hdf5_name"]]

    def write_hdf5_object_level(self, object_type, object_results_node):
        level = "object"
        rs, object_datasets = self.filter_datasets(level)
        for ds in object_datasets:
            ds_attributes = self.filter_dataset_attributes(ds)
            ds_datatype = np.dtype([(row["hdf5_name"], row["hdf5_type"]) for index, row in ds_attributes.iterrows()])
            # TODO we should add dynamic column(s) to results_specification to specify if the attribute has been loaded
            if self.has_dataset(object_type, ds):
                ds_size = self.get_dataset_size(object_type, ds)
                if ds_size > 1:
                    if ds == "firstpass_pdf":#compress firstpass pdf
                        object_results_node.create_dataset(ds,
                                                        (ds_size,),
                                                        ds_datatype,
                                                        self.object_dataframes[object_type][ds].to_records(),
                                                        compression="lzf")
                    else:
                        object_results_node.create_dataset(ds,
                                (ds_size,),
                                ds_datatype,
                                self.object_dataframes[object_type][ds].to_records())

                else:
                    object_results_node.create_group(ds)
                    for index,ds_row in ds_attributes.iterrows():
                        if self.has_attribute(object_type,ds,ds_row["hdf5_name"]):
                            object_results_node.get(ds).attrs[ds_row["hdf5_name"]] = self.object_results[object_type][ds][ds_row["hdf5_name"]]
                    
    def write_hdf5(self,hdf5_root,spectrum_id):
        obs = hdf5_root.create_group(spectrum_id)
        self.write_hdf5_root(obs)

        for object_type in self.object_types:
            object_results = obs.create_group(object_type) #h5
            self.write_hdf5_object_level(object_type, object_results)
            
            # warning flag
            rs = self.results_specifications
            rs = rs[rs["level"] == "method"]
            methods_datasets = list(rs["hdf5_dataset"].unique())
            for ds in methods_datasets:
                if self.has_dataset(object_type,ds):
                    ds_attributes = rs[rs["hdf5_dataset"] == ds] 
                    object_results.create_group(ds)
                    methods = self.get_solve_methods(object_type)
                    for method in methods:                                               
                        for index,ds_row in ds_attributes.iterrows():
                            if "<MethodType>" in ds_row["hdf5_name"]:
                                attr_name = ds_row["hdf5_name"].replace("<MethodType>", method)
                            else :
                                attr_name = ds_row["hdf5_name"]
                            if self.has_attribute(object_type,ds,attr_name):
                                object_results.get(ds).attrs[attr_name] = self.object_results[object_type][ds][attr_name]

            if self.get_solve_method(object_type):

                candidates = object_results.create_group("candidates")
                level = "candidate"
                rs, candidate_datasets = self.filter_datasets(level)
                nb_candidates = len(self.object_results[object_type]["model_parameters"])
                for rank in range(nb_candidates):
                    candidate = candidates.create_group(self.get_candidate_group_name(rank))
                    for ds in candidate_datasets:
                        if self.has_dataset(object_type,ds):
                            ds_attributes = rs[rs["hdf5_dataset"]==ds]
                            ds_dim = ds_attributes["dimension"].unique()[0]
                            if ds_dim == "mono":
                                candidate.create_group(ds)

                for ds in candidate_datasets:
                    if not ds in self.object_results[object_type]:
                        continue
                    ds_attributes = rs[rs["hdf5_dataset"]==ds]

                    # TODO change here when we will deal with 2D ranking or model ranking
                    dimension=None
                    for rank in range(nb_candidates):
                        candidate = candidates.get(self.get_candidate_group_name(rank))
                        for index,ds_row in ds_attributes.iterrows():
                            dimension = ds_row["dimension"]
                            if self.has_attribute(object_type,
                                                  ds_row.hdf5_dataset,
                                                  ds_row.hdf5_name,
                                                  rank) and dimension == "mono":
                                candidate.get(ds).attrs[ds_row["hdf5_name"]] = self.object_results[object_type][ds][rank][ds_row["hdf5_name"]]
                        if dimension == "multi":
                            ds_datatype = np.dtype([(row["hdf5_name"],row["hdf5_type"]) for index, row in ds_attributes.iterrows()])
                            ds_size = self.get_dataset_size(object_type,ds,rank)
                            candidate.create_dataset(ds,
                                                     (ds_size,),
                                                     ds_datatype,
                                                     self.object_dataframes[object_type][ds][rank].to_records())

    def _get_attribute_from_result_store(self,data_spec,object_type=None,rank=None, linemeas = None):
        operator_result = self.get_operator_result(data_spec,object_type,rank, linemeas)
        if data_spec.dimension == "mono":
            if "[object_type]" in data_spec.OperatorResult_name:
                operator_result_name = data_spec.OperatorResult_name.replace("[object_type]","")
                return getattr(operator_result, operator_result_name)[object_type]
            elif "[method_type]" in data_spec.OperatorResult_name:
                operator_result_name = data_spec.OperatorResult_name.replace("[method_type]","")
                return getattr(operator_result, operator_result_name)
            else:
                return getattr(operator_result, data_spec.OperatorResult_name)
        else:
            if data_spec.hdf5_type == 'f8':
                return PC_Get_Float64Array(getattr(operator_result, data_spec.OperatorResult_name))
            if data_spec.hdf5_type == 'f4':#this should be Float32Array..to be corrected
                return PC_Get_Float64Array(getattr(operator_result, data_spec.OperatorResult_name))
            if data_spec.hdf5_type == 'i':
                return PC_Get_Int32Array(getattr(operator_result, data_spec.OperatorResult_name))

    def get_attribute_from_result_store(self, attribute_name, object_type, rank):
        rs = self.results_specifications
        rs = rs[rs["hdf5_name"] == attribute_name]
        return self._get_attribute_from_result_store(rs.iloc[0], object_type, rank)

    def has_attribute_in_result_store(self,data_spec,object_type,rank=0):
        if rank is not None:
            if self.results_store.HasCandidateDataset(object_type,
                                                      self.get_solve_method(object_type),
                                                      data_spec.ResultStore_key,
                                                      data_spec.hdf5_dataset):
                operator_result = self.get_operator_result(data_spec, object_type, rank)
            else:
                return False
        else:
            operator_result = self.get_operator_result(data_spec, object_type, rank=None)
        return hasattr(operator_result, data_spec.OperatorResult_name)

    def get_operator_result(self, data_spec, object_type, rank = None, linemeas=None):
        if object_type is not None:
            if data_spec.hdf5_dataset in self.operator_results[object_type]:
                if rank is not None:
                    if rank not in self.operator_results[object_type][data_spec.hdf5_dataset]:
                        self.operator_results[object_type][data_spec.hdf5_dataset][rank] = self.load_operator_result(
                            data_spec,
                            object_type,
                            rank)
                    return self.operator_results[object_type][data_spec.hdf5_dataset][rank]
                elif rank is None:
                    return self.operator_results[object_type][data_spec.hdf5_dataset]
            else:
                if rank is not None:
                    self.operator_results[object_type][data_spec.hdf5_dataset] = dict()
                    self.operator_results[object_type][data_spec.hdf5_dataset][rank] = self.load_operator_result(data_spec,
                                                                                                                 object_type,
                                                                                                                 rank)
                    return self.operator_results[object_type][data_spec.hdf5_dataset][rank]
                else:
                    self.operator_results[object_type][data_spec.hdf5_dataset] = self.load_operator_result(data_spec,
                                                                                                           object_type,
                                                                                                           rank, linemeas)
                    return self.operator_results[object_type][data_spec.hdf5_dataset]
        else:
            if data_spec.hdf5_dataset in self.operator_results:
                return self.operator_results[data_spec.hdf5_dataset]
            else:
                self.operator_results[data_spec.hdf5_dataset] = self.load_operator_result(data_spec,
                                                                                          object_type,
                                                                                          rank)
                return self.operator_results[data_spec.hdf5_dataset]

    def load_operator_result(self, data_spec, object_type, rank=None, linemeas = None):
        if data_spec.level == "root":
            if data_spec.ResultStore_key == "context_warningFlag":
                return self.results_store.GetFlagResult(data_spec.hdf5_dataset,
                                                        data_spec.hdf5_dataset,
                                                        data_spec.ResultStore_key)
            else :
                or_type = self.results_store.GetGlobalResultType(data_spec.hdf5_dataset,
                                                                data_spec.hdf5_dataset,
                                                                data_spec.ResultStore_key)
                if or_type == "CClassificationResult":
                    return self.results_store.GetClassificationResult(data_spec.hdf5_dataset,
                                                                    data_spec.hdf5_dataset,
                                                                    data_spec.ResultStore_key)
                elif or_type == "CReliabilityResult":
                    return self.results_store.GetReliabilityResult(data_spec.hdf5_dataset,
                                                                    data_spec.hdf5_dataset,
                                                                    data_spec.ResultStore_key)
                else:
                    raise Exception("Unknown OperatorResult type " + or_type)
        elif data_spec.level == "method":
            if linemeas :
                method = self.parameters[object_type]["linemeas_method"]
            else: 
                method = self.get_solve_method(object_type)
            or_type = self.results_store.GetGlobalResultType(object_type,
                                                             method,
                                                             data_spec.ResultStore_key)
            if or_type == "CFlagLogResult":
                return  self.results_store.GetFlagResult(object_type,
                                                        method,
                                                        data_spec.ResultStore_key)
            else:
                raise Exception("Unknown OperatorResult type " + or_type)                                     

        elif data_spec.level == "object":
            if "linemeas" in data_spec.hdf5_dataset :
                method = self.parameters[object_type]["linemeas_method"]
            else: 
                method = self.get_solve_method(object_type)

            or_type = self.results_store.GetGlobalResultType(object_type,
                                                             method,
                                                             data_spec.ResultStore_key)

            if or_type == "CPdfMargZLogResult":
                return self.results_store.GetPdfMargZLogResult(object_type,
                                                               method,
                                                               data_spec.ResultStore_key)
            elif or_type == "CLineModelSolution":
                return self.results_store.GetLineModelSolution(object_type,
                                                               method,
                                                               data_spec.ResultStore_key)
            elif or_type == "CModelSpectrumResult":
                return self.results_store.GetModelSpectrumResult(object_type,
                                                                 method,
                                                                 data_spec.ResultStore_key)
            else:
                raise Exception("Unknown OperatorResult type " + or_type)
        elif data_spec.level == "candidate":
            method = self.get_solve_method(object_type)
            or_type = self.results_store.GetCandidateResultType(object_type,
                                                                method,
                                                                data_spec.ResultStore_key,
                                                                data_spec.hdf5_dataset)
            if or_type == "TLineModelResult":
                return self.results_store.GetLineModelResult(object_type,
                                                             method,
                                                             data_spec.ResultStore_key,
                                                             rank)
            elif or_type == "TExtremaResult":
                return self.results_store.GetExtremaResult(object_type,
                                                           method,
                                                           data_spec.ResultStore_key,
                                                           rank)
            elif or_type == "CModelSpectrumResult":
                return self.results_store.GetModelSpectrumResult(object_type,
                                                                 method,
                                                                 data_spec.ResultStore_key,
                                                                 rank)
            elif or_type == "CSpectraFluxResult":
                return self.results_store.GetSpectraFluxResult(object_type,
                                                               method,
                                                               data_spec.ResultStore_key,
                                                               rank)
            elif or_type == "CModelFittingResult":
                return self.results_store.GetModelFittingResult(object_type,
                                                                method,
                                                                data_spec.ResultStore_key,
                                                                rank)
            elif or_type == "TTplCombinationResult":
                return self.results_store.GetTplCombinationResult(object_type,
                                                                method,
                                                                data_spec.ResultStore_key,
                                                                rank)
            elif or_type == "CLineModelSolution":
                return self.results_store.GetLineModelSolution(object_type,
                                                               method,
                                                               data_spec.ResultStore_key,
                                                               rank)
            else:
                raise Exception("Unknown OperatorResult type " + or_type)

