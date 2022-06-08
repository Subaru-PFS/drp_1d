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
import pandas as pd
from pylibamazed.r_specifications import rspecifications
from pylibamazed.Parameters import Parameters
import numpy as np


class AbstractOutput:

    def __init__(self,
                 parameters,
                 results_specifications=rspecifications,
                 extended_results = True):
        self.parameters = parameters
        self.spectrum_id = ''
        self.root_results = dict()
        self.object_results = dict()
        self.extended_results = extended_results 
        self.results_specifications = pd.read_csv(results_specifications,
                                                  sep='\t'
                                                  )
        self.object_types = self.parameters.get_objects()

        for object_type in self.object_types:
            self.object_results[object_type] = dict()

        
    def load_all(self):
        self.load_root()
        for object_type in self.object_types:
            self.load_object_level(object_type)
            self.load_method_level(object_type)
            self.load_candidate_level(object_type)

    def get_attribute(self,object_type, dataset, attribute, rank = None):
        if object_type:
            if rank is None:
                return self.object_results[object_type][dataset][attribute]
            else:
                return self.object_results[object_type][dataset][rank][attribute]
        else:
            return self.root_results[dataset][attribute]
    
    def has_dataset(self, object_type, dataset):
        if object_type in self.object_results:
            return dataset in self.object_results[object_type]
        else:
            return False

    def has_attribute(self,object_type, dataset,attribute,rank=None):
        if object_type in self.object_results and dataset in self.object_results[object_type]:
            if rank is None:
                return attribute in self.object_results[object_type][dataset]
            else:
                return attribute in self.object_results[object_type][dataset][rank]
        else:
            return False

    def get_dataset_size(self, object_type, dataset, rank = None):
        if rank is None:
            if dataset in self.object_results[object_type]:
                first_attr = next(iter(self.object_results[object_type][dataset].values()))
            else:
                raise APIException("Dataset " + dataset + " does not exist")
        else:
            first_attr = next(iter(self.object_results[object_type][dataset][rank].values()))
        if type(first_attr).__module__ == np.__name__:
            return len(first_attr)
        else:
            return 1

    def get_available_datasets(self, level,object_type=None):
        if level == "root":
            return self.root_results.keys()
        elif level == "object":
            l = []
            for d in self.object_results[object_type].keys():
                if type(self.object_results[object_type][d]) == dict:
                    l.append(d)
            return l
        elif level == "candidate":
            l = []
            for d in self.object_results[object_type].keys():
                if type(self.object_results[object_type][d]) == list:
                    l.append(d)
            return l
        else:
            raise APIException("Unknown level " + level)

    def get_candidate_data(self, object_type, rank, data_name):
        mp = self.object_results[object_type]["model_parameters"][rank][data_name]
        return mp

    def get_dataset(self, dataset, object_type=None, rank=None):
        if object_type:
            if rank is not None:
                return self.object_results[object_type][dataset][rank]
            else:
                return self.object_results[object_type][dataset]
        else:
            return self.root_results[dataset]
        
    # TODO more robust version, should iterate over candidate datasets and check existence
    def get_nb_candidates(self,object_type):
        return len(self.object_results[object_type]["model"])
            
    def get_attribute_from_source(self,object_type, method, dataset, attribute ,rank=None):
        raise NotImplementedError("Implement in derived class")

    def has_attribute_in_source(self,object_type, method, dataset, attribute,rank=None):
        raise NotImplementedError("Implement in derived class")    

    def has_dataset_in_source(self, object_type, method, dataset):
        raise NotImplementedError("Implement in derived class")

    def has_candidate_dataset_in_source(self, object_type, method, dataset):
        raise NotImplementedError("Implement in derived class")    

    def get_nb_candidates_in_source(self, object_type, method):
        raise NotImplementedError("Implement in derived class")

    def get_level(self, dataset):
        rs = self.results_specifications
        rs = rs[rs.hdf5_dataset == dataset]
        return rs.level.unique()[0]
    
    def filter_datasets(self, level):
        rs = self.results_specifications
        # filter by level
        rs = rs[rs["level"] == level]
        all_datasets = list(rs["hdf5_dataset"].unique())
        
        # filter by extended_results
        if self.extended_results:
            return rs, all_datasets

        # a dataset is considered as debug if all its elements have debug = True
        filtered_datasets = []
        for ds in all_datasets:
            ds_attributes = rs[rs["hdf5_dataset"]==ds]
            extended_results = all(ds_row["extended_results"] == True for index, ds_row in ds_attributes.iterrows())
            if not extended_results:
                filtered_datasets.append(ds) 
        
        return rs, filtered_datasets

    def filter_dataset_attributes(self, ds_name):  
        rs = self.results_specifications 
        ds_attributes = rs[rs["hdf5_dataset"]==ds_name]   
        #filter ds_attributes by extended_results column
        if self.extended_results:
            return ds_attributes
        filtered_df = ds_attributes[ds_attributes["extended_results"]==self.extended_results]     
        return filtered_df

    # root is every first level data excluding self.objects
    # (currently, only classification)
    def load_root(self):
        level = "root"
        rs, root_datasets = self.filter_datasets(level)
        for ds in root_datasets:
            ds_attributes = self.filter_dataset_attributes(ds)
            self.root_results[ds] = dict()
            for index, ds_row in ds_attributes.iterrows():
                if "<" in ds_row["hdf5_name"]:
                    for object_type in self.parameters.get_objects():
                        if self.has_attribute_in_source(object_type,
                                                        None,
                                                        ds,
                                                        ds_row.hdf5_name):
                            attr = self.get_attribute_from_source(object_type,
                                                                  None,
                                                                  ds,
                                                                  ds_row.hdf5_name)
                            attr_name = ds_row["hdf5_name"].replace("<ObjectType>", object_type)
                            self.root_results[ds][attr_name] = attr
                else:
                    if self.has_attribute_in_source(object_type,
                                                    None,
                                                    ds_row.hdf5_dataset,
                                                    ds_row.hdf5_name):
                        self.root_results[ds][ds_row["hdf5_name"]] = self.get_attribute_from_source("root", None, ds_row.hdf5_dataset,ds_row.hdf5_name)

    def load_object_level(self, object_type):
        level = "object"
        rs, object_datasets = self.filter_datasets(level)
        for dataset in object_datasets:
            methods = self.parameters.get_solve_methods(object_type)
            for method in methods:
                if self.has_dataset_in_source(object_type,
                                              method,
                                              dataset):
                    self.object_results[object_type][dataset] = dict()
                    self.fill_object_dataset(object_type, method, dataset)
                    
    def fill_object_dataset(self, object_type, method, dataset):
        ds_attributes = self.filter_dataset_attributes(dataset)
        for index, ds_row in ds_attributes.iterrows():
            attr_name = ds_row.hdf5_name
            if self.has_attribute_in_source(object_type, method, dataset,attr_name):
                attr = self.get_attribute_from_source(object_type,
                                                      method,
                                                      dataset,
                                                      attr_name)
                self.object_results[object_type][dataset][attr_name] = attr

    def load_method_level(self, object_type):
        level = "method"
        rs, object_datasets = self.filter_datasets(level)
        for ds in object_datasets:
            methods = self.parameters.get_solve_methods(object_type)
            self.object_results[object_type][ds] = dict()
            for method in methods:
                if self.has_dataset_in_source(object_type,
                                              method,
                                              ds):
                    ds_attributes = self.filter_dataset_attributes(ds)                    
                    for index, ds_row in ds_attributes.iterrows():
                        attr_name = ds_row["hdf5_name"]
                        if "<MethodType>" in ds_row["hdf5_name"]:
                            attr_name = ds_row["hdf5_name"].replace("<MethodType>", method)
                        attr = self.get_attribute_from_source(object_type,
                                                              method,
                                                              ds_row.hdf5_dataset,
                                                              ds_row.hdf5_name)
                        self.object_results[object_type][ds][attr_name] = attr

    def load_candidate_level(self, object_type):
        method = self.parameters.get_solve_method(object_type)
        if not method:
            return
        level = "candidate"
        rs, candidate_datasets = self.filter_datasets(level)
        for ds in candidate_datasets:
            ds_attributes = self.filter_dataset_attributes(ds).copy()
            rs_key = ds_attributes["ResultStore_key"].unique()[0]
            if not self.has_candidate_dataset_in_source(object_type,
                                              method,
                                              ds):
                continue

            nb_candidates = self.get_nb_candidates_in_source(object_type,
                                                             method)
            candidates = []
            for rank in range(nb_candidates):
                candidates.append(dict())
                for index, ds_row in ds_attributes.iterrows():
                    attr_name = ds_row["hdf5_name"]
                    if self.has_attribute_in_source(object_type,
                                                    method,
                                                    ds,
                                                    attr_name,
                                                    rank):
                        attr = self.get_attribute_from_source(object_type,
                                                              method,
                                                              ds,
                                                              attr_name,
                                                              rank)
                        candidates[rank][attr_name] = attr
            self.object_results[object_type][ds] = candidates

                
                    
    def get_candidate_group_name(self,rank):
        return "candidate" + chr(rank+65) # 0=A, 1=B,....


