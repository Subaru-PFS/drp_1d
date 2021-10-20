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
from .AbstractOutput import AbstractOutput
from .OutputSpecifications import results_specifications
import numpy as np
import pandas as pd
import abc
from collections import defaultdict
from pandas.core.dtypes.common import is_numeric_dtype, is_string_dtype


def get_attribute_from_attrs(hdf5_group, attr_name, is_string):
    return hdf5_group.attrs[attr_name]


class Hdf5Output(AbstractOutput):

    def __init__(self,input_manager,hdf5_group, parameters):
        AbstractOutput.__init__(self,input_manager)
        self.hdf5_group = hdf5_group
        self.spectrum_id = hdf5_group.name[1:]
        self.parameters = parameters
        self.add_sup_columns = True
        self.object_types = []
        if self.parameters["enablegalaxysolve"]:
            self.object_types.append("galaxy")
        if self.parameters["enablestarsolve"]:
            self.object_types.append("star")
        if self.parameters["enableqsosolve"]:
            self.object_types.append("qso")
        if self.parameters["enablelinemeassolve"]:
            self.object_types.append("linemeas")
        for object_type in self.object_types:
            self.object_results[object_type] = dict()
            self.object_dataframes[object_type] = dict()
        self.load_all()

    def load_root(self):
        rs = results_specifications
        rs = rs[rs["level"] == "root"]
        root_datasets = list(rs["hdf5_dataset"].unique())
        for ds in root_datasets:
            if ds in self.hdf5_group.keys():
                ds_attributes = rs[rs["hdf5_dataset"]==ds]
                self.root_results[ds]=dict()
                for index,ds_row in ds_attributes.iterrows():
                    if self.hdf5_group.get(ds).attrs.__contains__(ds_row["hdf5_name"]):
                        self.root_results[ds][ds_row["hdf5_name"]] = self.get_attribute_from_h5(ds_row)

    def load_object_level(self, object_type):
        rs = results_specifications
        rs = rs[rs["level"] == "object"]
        object_datasets = list(rs["hdf5_dataset"].unique())
        for ds in object_datasets:
            if ds in self.hdf5_group.get(object_type).keys():
                ds_attributes = rs[rs["hdf5_dataset"] == ds]
                # TODO handle case where attribute dimension != multi (here we only have pdf for now)
                self.object_results[object_type][ds] = dict()
                self.object_dataframes[object_type][ds] = pd.DataFrame(np.array(self.hdf5_group.get(object_type).get(ds)))
            
    def load_candidate_level(self, object_type):
        rs = results_specifications
        rs = rs[rs["level"] == "candidate"]
        candidate_datasets = list(rs["hdf5_dataset"].unique())

        for ds in candidate_datasets:
            ds_attributes = rs[rs["hdf5_dataset"] == ds]
            rs_key = ds_attributes["ResultStore_key"].unique()[0]
            if not self.has_candidate_dataset(object_type, ds):
                continue
            #                self.object_results[ds]=dict()
            # self.object_results[ds]["dataframe"]=pd.DataFrame()
            # TODO change here when we will deal with 2D ranking or model ranking
            nb_candidates = self.get_nb_candidates(object_type)
            candidates = []
            candidates_df = []  # useless if dataset attributes dimension is not multi
            dimension = ds_attributes["dimension"].unique()[0]
            for rank in range(nb_candidates):
                candidates.append(dict())
                if dimension == "multi":
                    candidates_df.append(pd.DataFrame(np.array(self.get_h5_candidate_dataset(object_type, ds, rank))))
                else:
                    for index, ds_row in ds_attributes.iterrows():
                        if self.has_candidate_attribute_in_h5(object_type, ds, rank, ds_row["hdf5_name"]):
                            candidates[rank][ds_row["hdf5_name"]] = self.get_attribute_from_h5(ds_row,object_type,rank)
            self.object_results[object_type][ds] = candidates
            if dimension == "mono":
                res = defaultdict(list)
                {res[key].append(sub[key]) for sub in candidates for key in sub}
                self.object_dataframes[object_type][ds] = pd.DataFrame(res)
                self.object_dataframes[object_type][ds]["Rank"] = range(nb_candidates)
            else:
                self.object_dataframes[object_type][ds] = candidates_df

    def get_attribute_from_h5(self, attribute_spec, object_type=None, rank=None):
        if attribute_spec["dimension"] == "mono":
            if attribute_spec["level"] == "root":
                return get_attribute_from_attrs(self.hdf5_group.get(attribute_spec["hdf5_dataset"]),
                                                attribute_spec["hdf5_name"],
                                                attribute_spec["c_type"]=="string")
            elif attribute_spec["level"] == "object":
                return get_attribute_from_attrs(self.hdf5_group.get(object_type).get(attribute_spec["hdf5_dataset"]),
                                                attribute_spec["hdf5_name"],
                                                attribute_spec["c_type"] == "string")
            elif attribute_spec["level"] == "candidate":
                return get_attribute_from_attrs(self.get_h5_candidate_dataset(object_type,attribute_spec["hdf5_dataset"],rank),
                                                attribute_spec["hdf5_name"],
                                                attribute_spec["c_type"] == "string")
            else:
                raise Exception("Unknown level + " + attribute_spec["level"])
        else:
            raise Exception("get_attribute_from_h5 does not support dimension=multi (not necessary)")

    def load_candidates_results_supplementary_data(self, object_type):
       self.object_dataframes[object_type]["model_parameters"].set_index("Rank", inplace=True, drop=False)
       if self.reference_redshift is not None:
           self.object_dataframes[object_type]["model_parameters"]["abs_deltaZ"] = abs(
               self.object_dataframes[object_type]["model_parameters"]["Redshift"] -
                self.reference_redshift)
           reference_values = dict()
           for col in self.object_dataframes[object_type]["model_parameters"].columns:
                if is_numeric_dtype(self.object_dataframes[object_type]["model_parameters"][col]):
                    reference_values[col] = np.nan
                if is_string_dtype(self.object_dataframes[object_type]["model_parameters"][col]):
                    reference_values[col] = ""
                else:
                    reference_values[col] = []
           reference_values["Redshift"] = self.reference_redshift
           reference_values["Rank"] = -1
           self.object_dataframes[object_type]["model_parameters"] = self.object_dataframes[object_type]["model_parameters"].append(
                pd.Series(reference_values),
                ignore_index=True)
       manual_values = dict()
       for col in self.object_dataframes[object_type]["model_parameters"].columns:
            if is_numeric_dtype(self.object_dataframes[object_type]["model_parameters"][col]):
                manual_values[col] = np.nan
            if is_string_dtype(self.object_dataframes[object_type]["model_parameters"][col]):
                manual_values[col] = ""
            else:
                manual_values[col] = []
       manual_values["Redshift"]=0
       manual_values["Rank"]=-2
       self.object_dataframes[object_type]["model_parameters"] =self.object_dataframes[object_type]["model_parameters"].append(pd.Series(manual_values),ignore_index=True)

    def has_candidate_dataset(self,object_type, dataset):
        if "candidateA" in self.hdf5_group.get(object_type).get("candidates").keys():
            return dataset in self.hdf5_group.get(object_type).get("candidates").get("candidateA").keys()
        else:
            return False

    def get_nb_candidates(self, object_type):
        return len(self.hdf5_group.get(object_type).get("candidates"))

    def get_h5_candidate_dataset(self,object_type, dataset, rank):
        return self.hdf5_group.get(object_type).get("candidates").get(self.get_candidate_group_name(rank)).get(dataset)

    def has_candidate_attribute_in_h5(self,object_type, dataset, rank, attribute):
        return self.get_h5_candidate_dataset(object_type,dataset,rank).attrs.__contains__(attribute)
