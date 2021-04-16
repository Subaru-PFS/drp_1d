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
        self.manual_redshift = None
        
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

    @abc.abstractmethod
    def load_nb_candidates(self,object_type):
        pass

    def load_all_fitted_rays(self,object_type):
        for rank in range(self.nb_candidates[object_type]):
            self.load_fitted_rays(object_type,rank)

    def load_all_best_continuum(self,object_type):
        for rank in range(self.nb_candidates[object_type]):
            self.load_fitted_continuum_by_rank(object_type,rank)

    def load_all_models(self,object_type):
        for rank in range(self.nb_candidates[object_type]):
            self.load_model_by_rank(object_type,rank)

    def load_all(self):
        self.load_classification()
        self.load_rays_info()
        for object_type in self.object_types:
            self.load_nb_candidates(object_type)
            self.load_candidates_results(object_type)
            self.load_all_models(object_type)
            self.load_pdf(object_type)
            if self.get_solve_method(object_type) == "linemodelsolve":
                self.load_all_best_continuum(object_type)
                self.load_all_fitted_rays(object_type)


    def get_solve_method(self,object_type):
        return self.parameters[object_type]["method"]

    def get_candidate_data(self, object_type, rank, data_name):
        self.load_candidates_results(object_type)
        return self.candidates_results[object_type][self.candidates_results[object_type]['Rank'] == rank][data_name].iat[0]
    
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

    def get_classification(self):
        if not self.classification:
            self.load_classification()
        return self.classification

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

    def get_nb_candidates(self,object_type):
        self.load_nb_candidates(object_type)
        return self.nb_candidates[object_type]
    
    def set_manual_redshift_value(self,val):
        self.manual_redshift=val

    def compare_candidates_results(self,other,object_type,quiet):
        if not self.candidates_results[object_type].equals(other.candidates_results[object_type]):
            if not quiet:
                print("difference in " + object_type + "candidates results")
                print(self.candidates_results[object_type].compare(other.candidates_results[object_type]))
            return False
        else:
            return True
        
    def compare_pdf(self,other,object_type,quiet):
        pdf = self.pdf[object_type]
        opdf = other.pdf[object_type]
        if not pdf.equals(opdf):
            if not quiet:
                print("difference in " + object_type + " pdf")
            return False
        else:
            return True

    def compare_fitted_rays(self,other,object_type,rank,quiet):
        fr = self.fitted_rays[object_type][rank]
        ofr = other.fitted_rays[object_type][rank]
        if not fr.equals(ofr):
            if not quiet:
                print(" difference in " + object_type + " fitted rays for candidate n° " + str(rank))
                print(fr.compare(ofr))
            return False
        else:
            return True

    def compare_best_continuum(self,other,object_type,rank,quiet):
        bc = self.best_continuum[object_type][rank]
        obc = other.best_continuum[object_type][rank]
        if not (bc==obc).all():
            if not quiet:
                print(" difference in " + object_type + " best continuum for candidate n° " + str(rank))
                return False
        else:
            return True

    def compare_model(self,other,object_type,rank,quiet):
        m = self.model[object_type][rank]
        om = other.model[object_type][rank]
        if not (m==om).all():
            if not quiet:
                print(" difference in " + object_type + " model for candidate n° " + str(rank))
            return False
        else:
            return True

    def compare_classification(self,other,quiet):
        equal = True
        for k in self.classification.keys():
            if self.classification[k] != other.classification[k]:
                if not quiet:
                    print("classification " + k + " different : " + str(self.classification[k]) + " vs " + str(other.classification[k]))
                equal = False
        return equal
        
    def diff(self,other,object_types=["galaxy","star","qso"],quiet=False):
        equal = True
        for object_type in object_types:
            equal = self.compare_candidates_results(other,object_type,quiet)
            equal = equal and self.compare_pdf(other,object_type,quiet)                           
            nb_cand = self.nb_candidates[object_type]
            if nb_cand != other.nb_candidates[object_type]:
                equal = False
                if not quiet:
                    print(object_type + " results do not have same number of candidates : " + str(nb_cand) + " vs " + str(other.nb_candidates[object_type]) + ". Terminating diff operation")
                return equal
            for cand in range(nb_cand):
                equal = equal and self.compare_model(other,object_type,cand,quiet)
            if object_type == 'galaxy':
                for cand in range(nb_cand):
                    equal = equal and self.compare_fitted_rays(other,object_type,cand,quiet)
                for cand in range(nb_cand):
                    equal = equal and self.compare_best_continuum(other,object_type,cand,quiet)
        equal = equal and self.compare_classification(other,quiet)
        return equal

    
