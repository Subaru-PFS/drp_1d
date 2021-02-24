import abc
from .InputManager import InputManager
import pandas as pd


class AbstractOutput:

    def __init__(self,input_manager):
        self.spectrum_id = ''
        self.root_results = dict()
        self.root_dataframes = dict()
        self.object_results = dict()
        self.object_dataframes = dict()
        self.reference_redshift = None # TODO should be a map or np array with all attributes presents in the ref file
        self.manual_redshift = None
        
    def load_all(self):
        pass

    def get_solve_method(self,object_type):
        return self.parameters[object_type]["method"]

    def get_attribute(self,object_type, dataset, attribute, rank = None):
        if object_type in ["star","qso","galaxy"]:
            if rank is None:
                return self.object_results[object_type][dataset][attribute]
            else:
                return self.object_results[object_type][dataset][rank][attribute]
        else:
            return self.root_results[dataset][attribute]
    
    def get_candidate_data(self, object_type, rank, data_name):
        mp = self.object_dataframes[object_type]["model_parameters"]
        return mp[mp['Rank'] == rank][data_name].iat[0]
    
    def get_candidate_results(self, object_type, rank, columns=[]):
        cr = self.object_dataframes[object_type]["model_parameters"]
        if not columns:
            return cr[cr['Rank'] == rank]
        else:
            return cr[cr['Rank'] == rank][columns]

    def get_candidates_results(self, object_type, columns=[]):
        cr = self.object_dataframes[object_type]["model_parameters"]
        if not columns:
            return cr
        else:
            return cr[columns]

    def get_pdf(self,object_type):
        return self.object_dataframes[object_type]["pdf"]

    def get_rays_info(self,object_type="galaxy"):
        return self.object_results[object_type]["rays_info"]

    def get_rays_info_df(self,object_type="galaxy"):
        return self.object_results[object_type]["rays_info"]

    def get_classification_type(self):
        return self.root_results["classification"]["Type"]

    def get_classification(self):
        return self.root_results["classification"]
    
    def get_fitted_continuum_by_rank(self, object_type, rank):
        self.load_fitted_continuum_by_rank(object_type,rank)
        return self.best_continuum[object_type][rank]

    def get_fitted_model_by_rank(self, object_type, rank):
        return self.object_dataframes[object_type]["model"][rank]

    def get_fitted_rays_by_rank(self, object_type, rank):
        return self.object_dataframes[object_type]["fitted_rays"][rank]

    def get_candidate_group_name(self,rank):
        return "candidate" + chr(rank+65) # 0=A, 1=B,....

    def get_nb_candidates(self,object_type):
        self.load_nb_candidates(object_type)
        return self.nb_candidates[object_type]
    
    def set_manual_redshift_value(self,val):
        self.manual_redshift=val

    def get_reliability(self,object_type):
        return self.object_results[object_type]["reliability"]

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
            return self.object_dataframes[object_type][dataset].index.size
        else:
            if type(self.object_dataframes[object_type][dataset]) == pd.DataFrame:
                return 1
            else:
                return self.object_dataframes[object_type][dataset][rank].index.size
