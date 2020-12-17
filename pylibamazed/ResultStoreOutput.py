from .AbstractOutput import AbstractOutput
from .OutputSpecifications import candidate_specifications
import numpy as np
import pandas as pd
import abc
import resource
class ResultStoreOutput(AbstractOutput):

    def __init__(self,input_manager,result_store, parameters):
        AbstractOutput.__init__(self,input_manager)
        self.result_store = result_store
        self.parameters = parameters
        self.object_types = ["galaxy"]
        if self.parameters["enablestellarsolve"] == "yes":
            self.object_types.append("star")

    def load_pdf(self, object_type):
        if object_type not in self.pdf.keys():
            pdf_proba = self.result_store.Get_Float64ArrayData(object_type, "all", "pdf_probaLog")
            pdf_zgrid = self.result_store.Get_Float64ArrayData(object_type, "all", "pdf_zgrid")
            self.pdf[object_type] = pd.DataFrame()
            self.pdf[object_type]["probaLog"] = pdf_proba
            self.pdf[object_type]["zGrid"] = pdf_zgrid

    def load_candidates_results(self, object_type):
        if object_type not in self.candidates_results.keys():
            self.candidates_results[object_type] = self.load_candidates_dataset_df("model_parameters",object_type)

            # self.candidates_results[object_type][col].apply(
            #         lambda x: x.decode('utf-8'))
            # cand_specs = candidate_specifications[(candidate_specifications["method"] == "all") |
            #                                       (candidate_specifications["method"].str.contains(solve_method))]
            # cand_specs = cand_specs[cand_specs["dataset"] == "model_parameters"]
            # 
            # comp_type = []
            # for index, row in cand_specs[cand_specs["dimension"] == "mono"].iterrows():
            #     comp_type.append((row["name"], row["hdf5_type"]))
            # for index, row in cand_specs[cand_specs["dimension"] == "multi"].iterrows():
            #     comp_type.append((row["name"], "object"))
            # 
            # nb_candidates = self.result_store.Get_Int32Data(object_type, "all", "NbCandidates")
            # candidates_ra = np.recarray((nb_candidates,),dtype=comp_type)
            # for rank in range(nb_candidates):
            #     for index, row in cand_specs.iterrows():
            #         candidates_ra[row["name"]][rank] = self.get_attribute(row,object_type,rank,solve_method)
            # self.candidates_results[object_type] = pd.DataFrame(candidates_ra)

    @abc.abstractmethod
    def load_rays_info(self):
        if self.rays_infos is None:
            self.rays_infos = np.recarray((1,),[("snrHa","f8"),("lfHa","f8"),("snrOII","f8"),("lfOII","f8")])
            self.rays_infos["snrHa"][0]= self.result_store.Get_Float64Data("all", "all", "snrHa")
            self.rays_infos["lfHa"][0]= self.result_store.Get_Float64Data("all", "all", "lfHa")
            self.rays_infos["snrOII"][0] = self.result_store.Get_Float64Data("all", "all", "snrOII")
            self.rays_infos["lfOII"][0] = self.result_store.Get_Float64Data("all", "all", "lfOII")

    @abc.abstractmethod
    def load_fitted_rays(self, object_type,rank):
        if rank not in self.fitted_rays:
            self.fitted_rays[object_type][rank] = self.load_candidate_dataset_df(object_type,
                                                                                 "fitted_rays",
                                                                                 rank)

    @abc.abstractmethod
    def load_fitted_continuum_by_rank(self, object_type, rank):
        if rank not in self.best_continuum[object_type]:
            self.best_continuum[object_type][rank] = self.load_candidate_dataset_df(object_type,
                                                                                    "continuum",
                                                                                    rank)
            
    @abc.abstractmethod
    def load_model_by_rank(self, object_type, rank):
        if rank not in self.model[object_type]:
            self.model[object_type][rank] = self.load_candidate_dataset_df(object_type, "model", rank)

    @abc.abstractmethod
    def load_classification(self):
        self.classification["Type"]=self.result_store.Get_StringData("classification", "all", "Type")
        self.classification["GalaxyProba"]=self.result_store.Get_Float64Data("classification", "all", "ProbGalaxy")
        self.classification["StarProba"]=self.result_store.Get_Float64Data("classification", "all", "ProbStar")
        self.classification["QSOProba"]=self.result_store.Get_Float64Data("classification", "all", "ProbQSO")
        self.classification["GalaxyEvidence"] = self.result_store.Get_Float64Data("classification", "all", "EvidenceGalaxy")
        self.classification["StarEvidence"] = self.result_store.Get_Float64Data("classification", "all", "EvidenceStar")
        self.classification["QSOEvidence"] = self.result_store.Get_Float64Data("classification", "all", "EvidenceQSO")

    def load_all_fitted_rays(self,object_type):
        nb_candidates = self.result_store.Get_Int32Data(object_type, "all", "NbCandidates")
        for rank in range(nb_candidates):
            self.load_fitted_rays(object_type,rank)

    def load_all_best_continuum(self,object_type):
        nb_candidates = self.result_store.Get_Int32Data(object_type, "all", "NbCandidates")
        for rank in range(nb_candidates):
            self.load_fitted_continuum_by_rank(object_type,rank)

    def load_all_models(self,object_type):
        nb_candidates = self.result_store.Get_Int32Data(object_type, "all", "NbCandidates")
        for rank in range(nb_candidates):
            self.load_model_by_rank(object_type,rank)

    def load_all(self):
        self.load_classification()
        self.load_rays_info()
        for object_type in self.object_types:
            self.load_candidates_results(object_type)
            self.load_all_models(object_type)
            self.load_all_best_continuum(object_type)
            self.load_all_fitted_rays(object_type)
            self.load_pdf(object_type)

    def get_attribute(self,spec_row, object_type, rank):
        if spec_row["method"] == "all":
            method = "all"
        else:
            method = self.get_solve_method(object_type)
        if spec_row["c_type"] == 'int':
            return self.result_store.Get_Int32CandidateData(object_type, method, rank, spec_row["name"])
        elif spec_row["c_type"] == 'string':
            return self.result_store.Get_StringCandidateData(object_type, method, rank, spec_row["name"])
        elif spec_row["c_type"] == 'double':
            return self.result_store.Get_Float64CandidateData(object_type, method, rank, spec_row["name"])
        elif spec_row["c_type"] == 'list[double]':
            return self.result_store.Get_Float64ArrayCandidateData(object_type, method, rank, spec_row["name"])
        elif spec_row["c_type"] == 'list[int]':
            return self.result_store.Get_Int32ArrayCandidateData(object_type, method, rank, spec_row["name"])
        else:
            raise Exception("unknown c type " + spec_row["c_type"])
        
    def load_candidates_dataset(self,dataset, object_type):
        solve_method = self.get_solve_method(object_type)
        cand_specs = candidate_specifications[(candidate_specifications["method"] == "all") |
                                              (candidate_specifications["method"].str.contains(solve_method))]
        cand_specs = cand_specs[cand_specs["dataset"] == dataset]

        comp_type = []
        for index, row in cand_specs[cand_specs["dimension"] == "mono"].iterrows():
            comp_type.append((row["name"], row["hdf5_type"]))
        for index, row in cand_specs[cand_specs["dimension"] == "multi"].iterrows():
            comp_type.append((row["name"], "object"))

        nb_candidates = self.result_store.Get_Int32Data(object_type, "all", "NbCandidates")
        self.nb_candidates[object_type]=nb_candidates
        ra = np.recarray((nb_candidates,), dtype=comp_type)
        for rank in range(nb_candidates):
            for index, row in cand_specs.iterrows():
                ra[row["name"]][rank] = self.get_attribute(row, object_type, rank)

        return ra

    def load_candidates_dataset_df(self,dataset, object_type):
        solve_method = self.get_solve_method(object_type)
        cand_specs = candidate_specifications[(candidate_specifications["method"] == "all") |
                                              (candidate_specifications["method"].str.contains(solve_method))]
        cand_specs = cand_specs[cand_specs["dataset"] == dataset]
        df = pd.DataFrame()
        comp_type = []

        nb_candidates = self.result_store.Get_Int32Data(object_type, "all", "NbCandidates")
        self.nb_candidates[object_type]=nb_candidates
        columns_data = dict()
        for index, row in cand_specs.iterrows():
            columns_data[row["name"]]=[]
        for rank in range(nb_candidates):
            for index, row in cand_specs.iterrows():
                columns_data[row["name"]].append(self.get_attribute(row, object_type, rank))
        return pd.DataFrame(columns_data)

    def load_candidate_dataset(self,object_type, dataset, rank):
        solve_method=self.get_solve_method(object_type)
        cand_specs = candidate_specifications[(candidate_specifications["method"] == "all") |
                                              (candidate_specifications["method"].str.contains(solve_method))]
        cand_specs = cand_specs[cand_specs["dataset"] == dataset]

        comp_type = []
        for index, row in cand_specs[cand_specs["dimension"] == "mono"].iterrows():
            comp_type.append((row["name"], row["hdf5_type"]))
        for index, row in cand_specs[cand_specs["dimension"] == "multi"].iterrows():
            comp_type.append((row["name"], "object"))

        nb_records = self.get_attribute(next(cand_specs.iterrows())[1],object_type,rank).size
        ra = np.recarray((nb_records,), dtype=comp_type)
        for index, row in cand_specs.iterrows():
            ra[row["name"]] = self.get_attribute(row, object_type, rank)
        return ra

    def load_candidate_dataset_df(self,object_type, dataset, rank):
        solve_method=self.get_solve_method(object_type)
        cand_specs = candidate_specifications[(candidate_specifications["method"] == "all") |
                                              (candidate_specifications["method"].str.contains(solve_method))]
        cand_specs = cand_specs[cand_specs["dataset"] == dataset]

        comp_type = []
        for index, row in cand_specs[cand_specs["dimension"] == "mono"].iterrows():
            comp_type.append((row["name"], row["hdf5_type"]))
        for index, row in cand_specs[cand_specs["dimension"] == "multi"].iterrows():
            comp_type.append((row["name"], "object"))

        nb_records = self.get_attribute(next(cand_specs.iterrows())[1],object_type,rank).size
        df = pd.DataFrame()
        for index, row in cand_specs.iterrows():
            df[row["name"]] = self.get_attribute(row, object_type, rank)
        return df

    def get_solve_method(self,object_type):
        if object_type == "galaxy":
            return self.parameters["method"]
        else:
            return "chisquare2solve"

    def write_candidate_dataset_to_hdf5(self, dataset,rank , datatype):
        data = self.__getattr__(dataset)[object_type][rank].to_records()
    # model = np.array()
        data_size = self.model[object_type][rank].index.size
        candidate.create_dataset(dataset,
                                 (data_size,),
                                 cand_datatype,
                                 data=model)
        
    def write_hdf5(self,hdf5_root,spectrum_id,zlog):
        zlog.LogInfo("[hdf5] writing " + str(spectrum_id))
        obs = hdf5_root.create_group(spectrum_id)
        obs.attrs["user_time"] = resource.getrusage(resource.RUSAGE_SELF).ru_utime
        obs.attrs["system_time"] = resource.getrusage(resource.RUSAGE_SELF).ru_stime
        obs.attrs["memory_used"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        self.load_all()

        classification = obs.create_dataset('classification', (1,), np.dtype([('Type', 'S1'),
                                                                              ('GalaxyProba','f8'),
                                                                              ('StarProba', 'f8'),
                                                                              ('QSOProba', 'f8'),
                                                                              ('GalaxyEvidence', 'f8'),
                                                                              ('StarEvidence', 'f8'),
                                                                              ('QSOEvidence', 'f8'),
                                                                              ]))
        classification[...] = [(self.classification["Type"],
                                self.classification["GalaxyProba"],
                                self.classification["StarProba"],
                                self.classification["QSOProba"],
                                self.classification["GalaxyEvidence"],
                                self.classification["StarEvidence"],
                                self.classification["QSOEvidence"]
                                )]

        rays_infos = obs.create_dataset('rays_info', (1,), np.dtype([('snrHa', 'f8'), ('lfHa', 'f8'),
                                                                     ('snrOII', 'f8'), ('lfOII', 'f8')]))
        rays_infos[...] = self.rays_infos

        for object_type in self.object_types:
                # zlog.LogDebug("Writing " + object_type + " results")
                object_result = obs.create_group(object_type)
                pdf = self.pdf[object_type].to_records()
                grid_size = self.pdf[object_type].index.size
                object_result.create_dataset("pdf",
                                             (grid_size,),
                                             np.dtype([('zGrid', 'f8'), ('probaLog', 'f8')]),
                                             data=pdf)

                cand_specs = candidate_specifications[candidate_specifications["object_type"].str.contains(object_type)]
                # TODO solve method retrieval should be reviewed after parameters.json restructuration
                solve_method = self.get_solve_method(object_type)
                cand_specs = cand_specs[
                    (cand_specs["method"] == "all") | (cand_specs["method"].str.contains(solve_method))]
                datasets = list(cand_specs["dataset"].unique())
                mono_specs = []
                multi_specs = []
                cand_types = dict()
                cand_multi_types = dict()

                for dataset in datasets:
                    cs = cand_specs[cand_specs["dataset"] == dataset]
                    mocs = cs[cs["dimension"] == "mono"]
                    mucs = cs[cs["dimension"] == "multi"]
                    if mucs.empty:
                        mono_specs.append(mocs)
                        cand_types[dataset]=np.dtype([(row["name"], row["hdf5_type"]) for index, row in cs.iterrows()])
                    else:
                        multi_specs.append(mucs)
                        cand_multi_types[dataset]=np.dtype([(row["name"],row["hdf5_type"]) for index, row in cs.iterrows()])

                candidates = object_result.create_group("candidates")
                nb_cands = self.nb_candidates[object_type]
                zlog.LogDebug("writing " + str(nb_cands) + "candidates")
                for rank in range(nb_cands):
                    zlog.LogDebug("writing candidate " + chr(rank+65))
                    candidate = candidates.create_group(self.get_candidate_group_name(rank))
                    model_parameters = candidate.create_dataset("model_parameters",
                                                                (1,), 
                                                                cand_types["model_parameters"])
#                                                                data=self.get_candidate_results(object_type,rank).to_records())
                    model_parameters[0]=self.get_candidate_results(object_type,rank).to_records(index=False)[0]

                    model = self.model[object_type][rank].to_records()
                    # model = np.array()
                    lambda_grid_size = self.model[object_type][rank].index.size
                    candidate.create_dataset("model",
                                             (lambda_grid_size,),
                                              cand_multi_types["model"],
                                              data=model)

                    if solve_method == "linemodel":
                        best_continuum = self.best_continuum[object_type][rank]
                        candidate.create_dataset("continuum",
                                                 (best_continuum.index.size,),
                                                 cand_multi_types["continuum"],
                                                 data=best_continuum.to_records(index=False))
                        fitted_rays = self.fitted_rays[object_type][rank]
                        candidate.create_dataset("fitted_rays",
                                                 (fitted_rays.index.size,),
                                                 cand_multi_types["fitted_rays"],
                                                 data=fitted_rays.to_records(index=False))



