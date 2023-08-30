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
import numpy as np
import pandas as pd
from pylibamazed.Exception import APIException
from pylibamazed.Parameters import Parameters
from pylibamazed.r_specifications import rspecifications
from pylibamazed.redshift import CLog, ErrorCode

RootStages = ["init", "classification", "load_result_store"]
ObjectStages = ["redshift_solver", "linemeas_catalog_load",
                "linemeas_solver", "reliability_solver", "sub_classif_solver"]

zlog = CLog.GetInstance()


class AbstractOutput:

    def __init__(self,
                 parameters: Parameters,
                 results_specifications=rspecifications,
                 extended_results=True):
        self.parameters = parameters
        self.spectrum_id = ''
        self.root_results = dict()
        self.object_results = dict()
        self.extended_results = extended_results
        self.results_specifications = pd.read_csv(results_specifications,
                                                  sep='\t'
                                                  )
        self.object_types = self.parameters.get_objects()
        self.errors = dict()
        for object_type in self.object_types:
            self.object_results[object_type] = dict()
        self.load_errors()
        self.cache = False

    def get_attribute_from_source(self, object_type, method, dataset, attribute,
                                  rank=None, band_name=None, obs_id=None):
        raise NotImplementedError("Implement in derived class")

    def has_attribute_in_source(self, object_type, method, dataset, attribute,
                                rank=None, band_name=None, obs_id=None):
        raise NotImplementedError("Implement in derived class")

    def has_dataset_in_source(self, object_type, method, dataset):
        raise NotImplementedError("Implement in derived class")

    def has_candidate_dataset_in_source(self, object_type, method, dataset):
        raise NotImplementedError("Implement in derived class")

    def get_nb_candidates_in_source(self, object_type, method):
        raise NotImplementedError("Implement in derived class")

    def load_errors(self):
        pass

    def has_error(self, object_type, stage):
        return self.get_error_full_name(object_type, stage) in self.errors

    def get_error(self, object_type, stage):
        return self.errors[self.get_error_full_name(object_type, stage)]

    def get_error_full_name(self, object_type, stage):
        if object_type:
            return f"{stage}_{object_type}"
        else:
            return stage

    def load_all(self):
        self.load_root()
        if self.has_error(None, "init"):
            return
        for object_type in self.object_types:
            self.load_object_level(object_type)
            self.load_method_level(object_type)
            self.load_candidate_level(object_type)
        self.cache = True

    def get_attribute_short(self, attribute, lines_ids):

        attr_parts = attribute.split(".")
        root = attr_parts[0]
        data = attr_parts[-1]
        rank = None
        if root == "classification":
            return self.get_attribute(None, "classification", data, None)
        elif root == "error":
            if self.has_error(attr_parts[1], attr_parts[2]):
                return self.get_error(attr_parts[1], attr_parts[2])[attr_parts[3]]
            elif self.has_error(None, attr_parts[1]):
                return self.get_error(None, attr_parts[1])[attr_parts[2]]
            else:
                return None
        elif root == "ContextWarningFlags":
            return self.get_attribute(None, "context_warningFlag", "ContextWarningFlags")
        elif "WarningFlags" in data:
            return self.get_attribute(root, "warningFlag", data)
        else:
            object_type = root
            if len(attr_parts) == 2:
                rank = 0
            elif len(attr_parts) == 3:
                rank = int(attr_parts[1])
        if len(attr_parts) < 4:
            if self.has_attribute(object_type, "model_parameters", data, rank):
                return self.get_attribute(object_type, "model_parameters", data, rank)
            else:
                rank = None
                for dataset in ["linemeas_parameters", "reliability"]:
                    if self.has_attribute(object_type, dataset, data, rank):
                        return self.get_attribute(object_type, dataset, data, rank)
        else:
            dataset = attr_parts[1]
            line_name = attr_parts[2]
            col_name = attr_parts[3]
            object_type = root
            if line_name not in lines_ids:
                raise APIException(ErrorCode.INTERNAL_ERROR, f"Line {line_name}  not found in {lines_ids}")
            if len(attr_parts) == 4:
                rank = None
                index_col = "LinemeasLineID"
            else:
                rank = int(attr_parts[4])
                index_col = "FittedLineID"
            fitted_lines_attr = self.get_attribute(object_type, dataset, col_name, rank)
            fitted_lines_idx = self.get_attribute(object_type, dataset, index_col, rank)
            df = pd.DataFrame({"idx": fitted_lines_idx, col_name: fitted_lines_attr}).set_index("idx")
            return df.at[lines_ids[line_name], col_name]
        return None

    def get_attribute(self, object_type, dataset, attribute, rank=None):
        if not self.cache:
            method = self.get_method(object_type, dataset)
            return self.get_attribute_from_source(object_type, method, dataset, attribute, rank)
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

    def get_method(self, object_type, dataset):
        if object_type is None:
            return None
        if dataset == "linemeas":
            return self.parameters.get_linemeas_method(object_type)
        else:
            return self.parameters.get_solve_method(object_type)

    def has_attribute(self, object_type, dataset, attribute, rank=None):
        if not self.cache:
            method = self.get_method(object_type, dataset)
            return self.has_attribute_in_source(object_type, method, dataset, attribute, rank)
        if not object_type:
            if dataset in self.root_results:
                return attribute in self.root_results[dataset]
            else:
                return False
        elif object_type in self.object_results and dataset in self.object_results[object_type]:
            if rank is None:
                return attribute in self.object_results[object_type][dataset]
            else:
                if type(self.object_results[object_type][dataset]) == list:
                    if len(self.object_results[object_type][dataset]) > rank:
                        return attribute in self.object_results[object_type][dataset][rank]
                    else:
                        return False
                else:
                    return False
        else:
            return False

    def get_dataset_size(self, object_type, dataset, rank=None):
        first_attr = None
        if rank is None:
            if dataset in self.object_results[object_type]:
                if not self.object_results[object_type][dataset]:
                    return 0
                first_attr = next(iter(self.object_results[object_type][dataset].values()))
            else:
                raise APIException(ErrorCode.INTERNAL_ERROR, "Dataset " + dataset + " does not exist")
        else:
            if len(self.object_results[object_type][dataset][rank]):
                first_attr = next(iter(self.object_results[object_type][dataset][rank].values()))
        if first_attr is None:
            return 0
        if type(first_attr) == np.ndarray:
            return len(first_attr)
        else:
            return 1

    def get_available_datasets(self, level, object_type=None):
        if level == "root":
            return self.root_results.keys()
        elif level == "object":
            datasets = []
            for d in self.object_results[object_type].keys():
                if type(self.object_results[object_type][d]) == dict:
                    datasets.append(d)
            return datasets
        elif level == "candidate":
            datasets = []
            for d in self.object_results[object_type].keys():
                if type(self.object_results[object_type][d]) == list:
                    datasets.append(d)
            return datasets
        else:
            raise APIException(ErrorCode.INTERNAL_ERROR, "Unknown level " + level)

    def get_candidate_data(self, object_type, rank, data_name):
        mp = self.object_results[object_type]["model_parameters"][rank][data_name]
        return mp

    def get_dataset(self, object_type, dataset, rank=None):
        if object_type:
            if rank is not None:
                return self.object_results[object_type][dataset][rank]
            else:
                return self.object_results[object_type][dataset]
        else:
            return self.root_results[dataset]

    # TODO more robust version, should iterate over candidate datasets and check existence
    def get_nb_candidates(self, object_type):
        available_datasets = self.get_available_datasets("candidate", object_type)
        if len(available_datasets) > 0:
            return len(self.object_results[object_type][available_datasets[0]])
        else:
            return 0

    def get_level(self, dataset):
        rs = self.results_specifications
        rs = rs[rs.dataset == dataset]
        return rs.level.unique()[0]

    def filter_datasets(self, level):
        rs = self.results_specifications
        # filter by level
        rs = rs[rs["level"] == level]
        all_datasets = list(rs["dataset"].unique())

        # filter by extended_results
        if self.extended_results:
            return rs, all_datasets

        # a dataset is considered as debug if all its elements have debug = True
        filtered_datasets = []
        for ds in all_datasets:
            ds_attributes = rs[rs["dataset"] == ds]
            extended_results = all(ds_row["extended_results"] for index,
                                   ds_row in ds_attributes.iterrows())
            if not extended_results:
                filtered_datasets.append(ds)

        return rs, filtered_datasets

    def filter_dataset_attributes(self, ds_name, object_type=None, method=None):
        rs = self.results_specifications
        ds_attributes = rs[rs["dataset"] == ds_name]
        # filter ds_attributes by extended_results column
        skipsecondpass = False
        if object_type is not None:
            skipsecondpass = self.parameters.get_lineModelSolve_skipsecondpass(object_type)

        # retrieve results which are not firstpass results
        # exclude firstpass results for methods other than LineModelSolve
        notLinemodel = False
        if method:
            notLinemodel = "LineModelSolve" not in method
        if skipsecondpass or notLinemodel:
            filtered_df = ds_attributes[~ds_attributes["name"].str.contains("Firstpass", na=True)]
        else:
            filtered_df = ds_attributes

        if self.extended_results:
            return filtered_df
        filtered_df = filtered_df.loc[~ds_attributes["extended_results"]]
        return filtered_df

    # root is every first level data excluding self.objects
    # (currently, only classification)
    def load_root(self):
        level = "root"
        rs, root_datasets = self.filter_datasets(level)
        for ds in root_datasets:
            skip = not self.has_dataset_in_source(None, None, ds)
            skip = skip and "warning" not in ds
            if skip:
                zlog.LogDebug("skipping " + ds)
                continue
            ds_attributes = self.filter_dataset_attributes(ds)
            self.root_results[ds] = dict()
            for index, ds_row in ds_attributes.iterrows():
                if "<" in ds_row["name"]:
                    for object_type in self.parameters.get_objects():
                        if self.has_attribute_in_source(object_type,
                                                        None,
                                                        ds,
                                                        ds_row["name"]):
                            attr = self.get_attribute_from_source(object_type,
                                                                  None,
                                                                  ds,
                                                                  ds_row["name"])
                            attr_name = ds_row["name"].replace("<ObjectType>", object_type)
                            self.root_results[ds][attr_name] = attr
                else:
                    if self.has_attribute_in_source(None,
                                                    None,
                                                    ds_row.dataset,
                                                    ds_row["name"]):
                        self.root_results[ds][ds_row["name"]] = self.get_attribute_from_source(
                            "root", None, ds_row.dataset, ds_row["name"])

    def load_object_level(self, object_type):
        level = "object"
        rs, object_datasets = self.filter_datasets(level)
        for dataset in object_datasets:
            methods = self.parameters.get_solve_methods(object_type)
            for method in methods:
                if self.has_dataset_in_source(object_type,
                                              method,
                                              dataset.replace("<ObsID>", "")):

                    if "<ObsID>" in dataset:
                        for obs_id in self.parameters.get_observation_ids():
                            self.object_results[object_type][dataset.replace("<ObsID>", obs_id)] = dict()
                            self.fill_object_dataset(object_type, method, dataset, obs_id)
                    else:
                        self.object_results[object_type][dataset] = dict()
                        self.fill_object_dataset(object_type, method, dataset)

    def fill_object_dataset(self, object_type, method, dataset, obs_id=""):
        ds_attributes = self.filter_dataset_attributes(dataset)
        for index, ds_row in ds_attributes.iterrows():
            attr_name = ds_row["name"]
            if self.has_attribute_in_source(object_type, method, dataset, attr_name, obs_id=obs_id):
                attr = self.get_attribute_from_source(object_type,
                                                      method,
                                                      dataset,
                                                      attr_name,
                                                      obs_id=obs_id)
                self.object_results[object_type][dataset.replace("<ObsID>", obs_id)][attr_name] = attr

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
                        attr_name = ds_row["name"]
                        if "<MethodType>" in ds_row["name"]:
                            attr_name = ds_row["name"].replace("<MethodType>", method)
                        attr = self.get_attribute_from_source(object_type,
                                                              method,
                                                              ds_row.dataset,
                                                              ds_row["name"])
                        self.object_results[object_type][ds][attr_name] = attr

    def load_candidate_level(self, object_type):
        method = self.parameters.get_solve_method(object_type)
        if not method:
            return
        level = "candidate"

        rs, candidate_datasets = self.filter_datasets(level)
        for ds in candidate_datasets:
            if not self.has_candidate_dataset_in_source(object_type,
                                                        method,
                                                        ds):
                continue
            multiobs_ds = "<ObsID>" in ds
            if not multiobs_ds:
                cds = self.build_candidate_dataset(object_type,
                                                   method,
                                                   ds)
                self.object_results[object_type][ds] = cds
            else:
                for obs_id in self.parameters.get_observation_ids():
                    ocds = self.build_candidate_dataset(object_type,
                                                        method,
                                                        ds,
                                                        obs_id)
                    self.object_results[object_type][ds.replace("<ObsID>", obs_id)] = ocds

    def build_candidate_dataset(self, object_type, method, dataset, obs_id=""):
        nb_candidates = self.get_nb_candidates_in_source(object_type,
                                                         method)
        ds_attributes = self.filter_dataset_attributes(dataset,
                                                       object_type,
                                                       method).copy()
        candidates = []
        for rank in range(nb_candidates):
            candidates.append(dict())
            for index, ds_row in ds_attributes.iterrows():
                attr_name = ds_row["name"]
                if "<BandName>" in attr_name:
                    # get phot bands from params
                    bands = self.parameters.get_photometry_bands()
                    for band in bands:
                        attr = self.get_attribute_wrapper(object_type,
                                                          method,
                                                          dataset,
                                                          attr_name,
                                                          rank=rank,
                                                          band_name=band)
                        if attr is not None:
                            attr_name_ = band
                            candidates[rank][attr_name_] = attr
                if "<ObsID>" in attr_name:
                    for obsId in self.parameters.get_observation_ids():
                        attr = self.get_attribute_wrapper(object_type,
                                                          method,
                                                          dataset,
                                                          attr_name,
                                                          rank=rank,
                                                          obs_id=obs_id)
                        if attr is not None:
                            attr_name_ = attr_name.replace("<ObsID>", obs_id)
                            candidates[rank][attr_name_] = attr
                else:
                    attr = self.get_attribute_wrapper(object_type,
                                                      method,
                                                      dataset,
                                                      attr_name,
                                                      rank=rank,
                                                      obs_id=obs_id)
                    if attr is not None:
                        candidates[rank][attr_name] = attr
        return candidates

    def get_attribute_wrapper(self, object_type, method, ds, attr_name,
                              rank=None, band_name=None, obs_id=None):
        attr = None
        if self.has_attribute_in_source(object_type,
                                        method,
                                        ds,
                                        attr_name,
                                        rank=rank,
                                        band_name=band_name,
                                        obs_id=obs_id):
            attr = self.get_attribute_from_source(object_type,
                                                  method,
                                                  ds,
                                                  attr_name,
                                                  rank,
                                                  band_name=band_name,
                                                  obs_id=obs_id)
        return attr

    def get_candidate_group_name(self, rank):
        return "candidate" + chr(rank + 65)  # 0=A, 1=B,....

    def get_attributes(self, attributes, lines_ids):
        ret = dict()
        ret["ProcessingID"] = self.spectrum_id
        for attribute in attributes:
            try:
                value = self.get_attribute_short(attribute, lines_ids)
                if value is not None:
                    ret[attribute] = value
            except Exception as e:
                zlog.LogDebug(f"could not extract {attribute} : {e}")
        return ret
