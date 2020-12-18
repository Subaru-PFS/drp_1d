import abc
from .InputManager import InputManager
import pandas as pd

class AbstractOutput:

    def __init__(self,input_manager):
        self.spectrum_id = ''
        self.candidates_results = dict()
        self.rays_infos = None
        self.fitted_rays = dict()
        self.fitted_rays["galaxy"] = dict()
        self.fitted_rays["qso"] = dict()
        self.fitted_rays["star"] = dict() # stars do not have fitted rays, be we want an empty slot for consistency
        self.model = dict()
        self.model["galaxy"] = dict()
        self.model["qso"] = dict()
        self.model["star"] = dict()
        self.pdf = dict()
        self.input_manager = input_manager
        self.best_continuum = dict()
        self.best_continuum["galaxy"] = dict()
        self.best_continuum["qso"] = dict()
        self.best_continuum["star"] = dict()
        self.classification = dict()
        self.nb_candidates = dict()
        self.reference_redshift = None # TODO should be a map or np array with all attributes presents in the ref file

    @abc.abstractmethod
    def load_pdf(self, object_type):
        pass

    @abc.abstractmethod
    def load_candidates_results(self, object_type):
        pass

    @abc.abstractmethod
    def load_rays_info(self):
        pass

    @abc.abstractmethod
    def load_fitted_rays(self, object_type,rank):
        pass
    
    @abc.abstractmethod
    def load_fitted_continuum_by_rank(self, rank):
        pass
    
    @abc.abstractmethod
    def load_model_by_rank(self, object_type, rank):
        pass

    @abc.abstractmethod
    def load_classification(self):
        pass


    def get_candidate_data(self, object_type, rank, data_name):
        self.load_candidates_results(object_type)
        return self.candidates_results[object_type][self.candidates_results[object_type]['Rank'] == rank][data_name].at[0]
    
    def get_candidate_results(self, object_type, rank, columns=[]):
        self.load_candidates_results(object_type)
        if not columns:
            return self.candidates_results[object_type][self.candidates_results[object_type]['Rank'] == rank]
        else:
            return self.candidates_results[object_type][self.candidates_results[object_type]['Rank'] == rank][columns]

    def get_candidates_results(self, object_type, columns=[]):
        self.load_candidates_results(object_type)
        if not columns:
            return self.candidates_results[object_type]
        else:
            return self.candidates_results[object_type][columns]

    def get_pdf(self,object_type):
        self.load_pdf(object_type)
        return self.pdf[object_type]

    def get_rays_info(self):
        self.load_rays_info()
        return self.rays_infos

    def get_rays_info_df(self):
        return pd.DataFrame(self.get_rays_info())

    def get_classification_type(self):
        return self.classification["Type"]

    def get_fitted_continuum_by_rank(self, object_type, rank):
        self.load_fitted_continuum_by_rank(object_type,rank)
        return self.best_continuum[object_type][rank]

    def get_fitted_model_by_rank(self, object_type, rank):
        self.load_model_by_rank(object_type,rank)
        return self.model[object_type][rank]

    def get_fitted_rays_by_rank(self, object_type, rank):
        self.load_fitted_rays(object_type, rank)
        return self.fitted_rays[object_type][rank]
    
    def get_candidate_group_name(self,rank):
        return "candidate" + chr(rank+65) # 0=A, 1=B,....