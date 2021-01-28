from .AbstractOutput import AbstractOutput
from .OutputSpecifications import candidate_specifications
import numpy as np
import pandas as pd
import abc
from pandas.core.dtypes.common import is_numeric_dtype, is_string_dtype

class Hdf5Output(AbstractOutput):

    def __init__(self,input_manager,hdf5_group, parameters):
        AbstractOutput.__init__(self,input_manager)
        self.hdf5_group = hdf5_group
        self.spectrum_id = hdf5_group.name[1:]
        self.parameters = parameters
        self.add_sup_columns = True
        self.object_types = ["galaxy"]
        if self.parameters["enablestellarsolve"] == "yes":
            self.object_types.append("star")

    def load_classification(self):
        classification = np.array(self.hdf5_group.get("classification"))
        self.classification["Type"]=classification["Type"][0].decode('utf-8')
        self.classification["GalaxyProba"]=classification["GalaxyProba"][0]
        self.classification["StarProba"]=classification[ "StarProba"][0]
        self.classification["QSOProba"]=classification[ "QSOProba"][0]
        self.classification["GalaxyEvidence"] = classification[ "GalaxyEvidence"][0]
        self.classification["StarEvidence"] = classification["StarEvidence"][0]
        self.classification["QSOEvidence"] = classification["QSOEvidence"][0]

    def load_pdf(self, object_type):
        if object_type not in self.pdf.keys():
            self.pdf[object_type] = pd.DataFrame(np.array(self.hdf5_group.get(object_type + "/pdf")))

    def load_candidates_results(self, object_type):
        if object_type not in self.candidates_results.keys():
            first = True
            candidates = self.hdf5_group.get(object_type + "/candidates")
            params_names = []
            for candidate in candidates.keys():
                cand = candidates.get(candidate)
                params = np.array(cand.get("model_parameters"))
                if first:
                    for p in params.dtype.names:
                        if np.issubdtype(params.dtype.fields[p][0], bytes):
                            params_names.append(p)
                    self.candidates_results[object_type] = pd.DataFrame(params)
                    first = False
                else:
                    self.candidates_results[object_type] = self.candidates_results[object_type].append(
                        pd.DataFrame(params))

            self.candidates_results[object_type] = pd.DataFrame(self.candidates_results[object_type])
            for param in params_names:
                self.candidates_results[object_type][param] = self.candidates_results[object_type][param].apply(
                    lambda x: x.decode('utf-8'))
            self.candidates_results[object_type].set_index("Rank", inplace=True, drop=False)
            if self.add_sup_columns:
                self.load_pdf(object_type)
                proba = []
                ranks = []
                proba_df = pd.DataFrame()
                for z, rank in zip(self.candidates_results[object_type]["Redshift"],
                                   self.candidates_results[object_type]["Rank"]):
                    self.pdf[object_type]["dist_z"] = abs(self.pdf[object_type]["zGrid"] - z)
                    proba.append(self.pdf[object_type][self.pdf[object_type]["dist_z"] ==
                                                       self.pdf[object_type]["dist_z"].min()]
                                 ["probaLog"].iloc[0])
                    ranks.append(rank)
                proba_df["proba"] = proba
                proba_df["ranks"] = ranks
                proba_df.set_index("ranks", inplace=True)

                self.candidates_results[object_type] = pd.merge(self.candidates_results[object_type], proba_df,
                                                                left_index=True, right_index=True)

                # Find closest indexes for [z-3dz,z+3z] in self.pdf[object_type][
                intg_area_indexes = []
                ranks = []
                intg_area_indexes_df = pd.DataFrame()
                for z, rank, dz in zip(self.candidates_results[object_type]["Redshift"],
                                       self.candidates_results[object_type]["Rank"],
                                       self.candidates_results[object_type]["RedshiftError"]):
                    zmin = min(enumerate(self.pdf[object_type]["zGrid"]),
                               key=lambda x: abs(z - 3 * dz - x[1]))
                    zmax = min(enumerate(self.pdf[object_type]["zGrid"]),
                               key=lambda x: abs(z + 3 * dz - x[1]))
                    intg_area_indexes.append([zmin[0], zmax[0]])
                    ranks.append(rank)

                intg_area_indexes_df["Rank"] = ranks
                intg_area_indexes_df["IntgAreaIndexes"] = intg_area_indexes
                intg_area_indexes_df.set_index("Rank", inplace=True)
                self.candidates_results[object_type] = pd.merge(self.candidates_results[object_type],
                                                                intg_area_indexes_df,
                                                                left_index=True, right_index=True)

                if self.reference_redshift is not None:
                    self.candidates_results[object_type]["abs_deltaZ"] = abs(
                        self.candidates_results[object_type]["Redshift"] -
                        self.reference_redshift)
                    reference_values = dict()
                    for col in self.candidates_results[object_type].columns:
                        if is_numeric_dtype(self.candidates_results[object_type][col]):
                            reference_values[col] = np.nan
                        if is_string_dtype(self.candidates_results[object_type][col]):
                            reference_values[col] = ""
                        else:
                            reference_values[col] = []
                    reference_values["Redshift"] = self.reference_redshift
                    reference_values["Rank"] = -1
                    self.candidates_results[object_type] = self.candidates_results[object_type].append(
                        pd.Series(reference_values),
                        ignore_index=True)
                manual_values = dict()
                for col in self.candidates_results[object_type].columns:
                    if is_numeric_dtype(self.candidates_results[object_type][col]):
                        manual_values[col] = np.nan
                    if is_string_dtype(self.candidates_results[object_type][col]):
                        manual_values[col] = ""
                    else:
                        manual_values[col] = []
                manual_values["Redshift"]=0
                manual_values["Rank"]=-2
                self.candidates_results[object_type] = self.candidates_results[object_type].append(pd.Series(manual_values),ignore_index=True)

    def load_rays_info(self):
        if self.rays_infos is None:
            self.rays_infos = np.array(self.hdf5_group.get("rays_info"))

    def load_fitted_rays(self, object_type, rank):
        if rank not in self.fitted_rays[object_type].keys():
            # TODO [bug] group should be retrieved by looking in parameters.rank rather than by name
            cand = self.hdf5_group.get("galaxy/candidates/" + self.get_candidate_group_name(rank))
            fitted_rays = np.array(cand.get("fitted_rays"))
            df_fr = pd.DataFrame(fitted_rays)
            df_fr.set_index('FittedRaysID', inplace=True)
            df_catalog = self.input_manager.get_full_catalog_df([])
            df_full = pd.merge(df_fr, df_catalog, left_index=True, right_index=True)
            self.fitted_rays[object_type][rank] = df_full[df_full["FittedRaysLambda"] > 0]

    def load_fitted_continuum_by_rank(self, object_type,rank):
        if rank not in self.best_continuum[object_type]:
            cand = self.hdf5_group.get("galaxy/candidates/"  + self.get_candidate_group_name(rank))
            self.best_continuum[object_type][rank] = np.array(cand.get("continuum"))

    def load_model_by_rank(self, object_type, rank):
        if rank not in self.model[object_type]:
            cand = self.hdf5_group.get(object_type + "/candidates/"  + self.get_candidate_group_name(rank))
            self.model[object_type][rank] = np.array(cand.get("model"))

    def load_nb_candidates(self,object_type):
        self.nb_candidates[object_type] = len(self.hdf5_group.get(object_type + "/candidates"))

