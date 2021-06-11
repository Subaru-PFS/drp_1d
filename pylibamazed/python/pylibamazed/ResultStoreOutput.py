from .AbstractOutput import AbstractOutput
from .OutputSpecifications import results_specifications
import numpy as np
import pandas as pd
import resource
from collections import defaultdict
from pylibamazed.redshift import PC_Get_Float64Array,PC_Get_Int32Array,CLog

zlog = CLog.GetInstance()


class ResultStoreOutput(AbstractOutput):

    def __init__(self, input_manager, result_store, parameters):
        AbstractOutput.__init__(self, input_manager)
        self.results_store = result_store
        self.parameters = parameters

        self.object_types = ["galaxy"]
        if self.parameters["enablestellarsolve"] == "yes":
            self.object_types.append("star")
        if self.parameters["enableqsosolve"] == "yes":
            self.object_types.append("qso")

        for object_type in self.object_types:
            self.object_results[object_type] = dict()
            self.object_dataframes[object_type] = dict()
        self.load_all()

    def load_root(self):
        rs = results_specifications
        rs = rs[rs["level"] == "root"]
        root_datasets = list(rs["hdf5_dataset"].unique())
        for ds in root_datasets:
            ds_attributes = rs[rs["hdf5_dataset"] == ds]
            self.root_results[ds]=dict()
            for index,ds_row in ds_attributes.iterrows():
                self.root_results[ds][ds_row["hdf5_name"]] = self.get_attribute_from_result_store(ds_row)

    def load_object_level(self, object_type):
        rs = results_specifications
        rs = rs[rs["level"] == "object"]
        object_datasets = list(rs["hdf5_dataset"].unique())
        for ds in object_datasets:
            ds_attributes = rs[rs["hdf5_dataset"] == ds]
            # TODO handle case where attribute dimension != multi (here we only have pdf for now)
            self.object_results[object_type][ds] = dict()
            self.object_dataframes[object_type][ds] = pd.DataFrame()
            for index, ds_row in ds_attributes.iterrows():
#                if self.results_store.hasAttribute(object_type, ds, attr):
                self.object_results[object_type][ds][ds_row["hdf5_name"]] = self.get_attribute_from_result_store(ds_row,object_type)
                self.object_dataframes[object_type][ds][ds_row["hdf5_name"]] = self.get_attribute_from_result_store(ds_row,object_type)

    def load_candidate_level(self, object_type):
        rs = results_specifications
        rs = rs[rs["level"] == "candidate"]
        candidate_datasets = list(rs["hdf5_dataset"].unique())
        for ds in candidate_datasets:
            ds_attributes = rs[rs["hdf5_dataset"] == ds]
            rs_key = ds_attributes["ResultStore_key"].unique()[0]
            if not self.results_store.HasCandidateDataset(object_type,
                                                          self.get_solve_method(object_type),
                                                          rs_key,
                                                          ds):
                continue
            #                self.object_results[ds]=dict()
            # self.object_results[ds]["dataframe"]=pd.DataFrame()
            # TODO change here when we will deal with 2D ranking or model ranking
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
                        candidates[rank][ds_row["hdf5_name"]] = self.get_attribute_from_result_store(ds_row,object_type,rank)
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

    def write_hdf5_root(self, hdf5_spectrum_node):
        rs = results_specifications
        rs = rs[rs["level"] == "root"]
        root_datasets = list(rs["hdf5_dataset"].unique())
                                   
        for ds in root_datasets:
            ds_attributes = rs[rs["hdf5_dataset"]==ds]
            dsg = hdf5_spectrum_node.create_group(ds)
            for index, ds_row in ds_attributes.iterrows():
                if ds_row["hdf5_name"] in self.root_results[ds]:
                    # we know we only have dimension = mono here
                    dsg.attrs[ds_row["hdf5_name"]] = self.root_results[ds][ds_row["hdf5_name"]]

    def write_hdf5_object_level(self, object_type, object_results_node):
        rs = results_specifications
        rs = rs[rs["level"] == "object"]
        object_datasets = list(rs["hdf5_dataset"].unique())
        for ds in object_datasets:
            ds_attributes = rs[rs["hdf5_dataset"] == ds]
            ds_datatype = np.dtype([(row["hdf5_name"], row["hdf5_type"]) for index, row in ds_attributes.iterrows()])
            # TODO we should add dynamic column(s) to results_specification to specify if the attribute has been loaded
            if self.has_dataset(object_type, ds):
                ds_size = self.get_dataset_size(object_type, ds)
                if ds_size > 1:
                    object_results_node.create_dataset(ds,
                                                  (ds_size,),
                                                  ds_datatype,
                                                  self.object_dataframes[object_type][ds].to_records()
                                                  )
                else:
                    raise Exception("h5write : TODO handle dataset with of size 1 in object scope")

    def write_hdf5(self,hdf5_root,spectrum_id):
        obs = hdf5_root.create_group(spectrum_id)
        obs.attrs["user_time"] = resource.getrusage(resource.RUSAGE_SELF).ru_utime
        obs.attrs["system_time"] = resource.getrusage(resource.RUSAGE_SELF).ru_stime
        obs.attrs["memory_used"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

        self.write_hdf5_root(obs)

        for object_type in self.object_types:
            object_results = obs.create_group(object_type) #h5
            self.write_hdf5_object_level(object_type, object_results)
            candidates = object_results.create_group("candidates")
            rs = results_specifications
            rs = rs[rs["level"] == "candidate"]
            candidate_datasets = list(rs["hdf5_dataset"].unique())
            nb_candidates = self.results_store.getNbRedshiftCandidates(object_type,self.get_solve_method(object_type))
            for rank in range(nb_candidates):
                candidate = candidates.create_group(self.get_candidate_group_name(rank))
                for ds in candidate_datasets:
                    if self.has_dataset(object_type,ds):
                        ds_attributes = rs[rs["hdf5_dataset"]==ds]
                        ds_dim = ds_attributes["dimension"].unique()[0]
                        if ds_dim == "mono":
                            candidate.create_group(ds)

            for ds in candidate_datasets:
                ds_attributes = rs[rs["hdf5_dataset"]==ds]

                # TODO change here when we will deal with 2D ranking or model ranking
                dimension=None
                for rank in range(nb_candidates):
                    candidate = candidates.get(self.get_candidate_group_name(rank))
                    for index,ds_row in ds_attributes.iterrows():
                        if self.has_attribute_in_result_store(ds_row, object_type, rank):
                            dimension=ds_row["dimension"]
                            if dimension == "mono":
                                candidate.get(ds).attrs[ds_row["hdf5_name"]] = self.object_results[object_type][ds][rank][ds_row["hdf5_name"]]
                    if dimension == "multi":
                        ds_datatype = np.dtype([(row["hdf5_name"],row["hdf5_type"]) for index, row in ds_attributes.iterrows()])
                        ds_size = self.get_dataset_size(object_type,ds,rank)
                        candidate.create_dataset(ds,
                                                 (ds_size,),
                                                 ds_datatype,
                                                 self.object_dataframes[object_type][ds][rank].to_records())

    def get_attribute_from_result_store(self,data_spec,object_type=None,rank=None):
        operator_result = self.get_operator_result(data_spec,object_type,rank)
        if data_spec.dimension == "mono":
            return getattr(operator_result, data_spec.OperatorResult_name)
        else:
            if data_spec.hdf5_type == 'f8':
                return PC_Get_Float64Array(getattr(operator_result, data_spec.OperatorResult_name))
            else:
                return PC_Get_Int32Array(getattr(operator_result, data_spec.OperatorResult_name))

    def has_attribute_in_result_store(self,data_spec,object_type,rank=0):
        if self.results_store.HasCandidateDataset(object_type,
                                                  self.get_solve_method(object_type),
                                                  data_spec.ResultStore_key,
                                                  data_spec.hdf5_dataset):
            operator_result = self.get_operator_result(data_spec, object_type, rank)
            return hasattr(operator_result, data_spec.OperatorResult_name)
        else:
            return False

    def get_operator_result(self, data_spec, object_type, rank = None):
        if data_spec.level == "root":
            or_type = self.results_store.GetGlobalResultType(data_spec.hdf5_dataset,
                                                             data_spec.hdf5_dataset,
                                                             data_spec.ResultStore_key)
            if or_type == "CClassificationResult":
                return self.results_store.GetClassificationResult(data_spec.hdf5_dataset,
                                                                  data_spec.hdf5_dataset,
                                                                  data_spec.ResultStore_key)
            else:
                raise Exception("Unknown OperatorResult type " + or_type)
        elif data_spec.level == "object":
            or_type = self.results_store.GetGlobalResultType(object_type,
                                                             self.get_solve_method(object_type),
                                                             data_spec.ResultStore_key)

            if or_type == "CPdfMargZLogResult":
                return self.results_store.GetPdfMargZLogResult(object_type,
                                                           self.get_solve_method(object_type),
                                                           data_spec.ResultStore_key)
            else:
                raise Exception("Unknown OperatorResult type " + or_type)
        elif data_spec.level == "candidate":
            or_type = self.results_store.GetCandidateResultType(object_type,
                                                                self.get_solve_method(object_type),
                                                                data_spec.ResultStore_key,
                                                                data_spec.hdf5_dataset)
            if or_type == "TLineModelResult":
                return self.results_store.GetLineModelResult(object_type,
                                                             self.get_solve_method(object_type),
                                                             data_spec.ResultStore_key,
                                                             rank)
            elif or_type == "TExtremaResult":
                return self.results_store.GetExtremaResult(object_type,
                                                           self.get_solve_method(object_type),
                                                           data_spec.ResultStore_key,
                                                           rank)
            elif or_type == "CModelSpectrumResult":
                return self.results_store.GetModelSpectrumResult(object_type,
                                                                 self.get_solve_method(object_type),
                                                                 data_spec.ResultStore_key,
                                                                 rank)
            elif or_type == "CSpectraFluxResult":
                return self.results_store.GetSpectraFluxResult(object_type,
                                                               self.get_solve_method(object_type),
                                                               data_spec.ResultStore_key,
                                                               rank)
            elif or_type == "CModelFittingResult":
                return self.results_store.GetModelFittingResult(object_type,
                                                                self.get_solve_method(object_type),
                                                                data_spec.ResultStore_key,
                                                                rank)
            else:
                raise Exception("Unknown OperatorResult type " + or_type)

