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
import abc
import pandas as pd


class AbstractOutput:

    def __init__(self):
        self.spectrum_id = ''
        self.root_results = dict()
        self.root_dataframes = dict()
        self.object_results = dict()
        self.object_dataframes = dict()
        self.reference_redshift = None # TODO should be a map or np array with all attributes presents in the ref file
        self.manual_redshift = None
        
    def load_all(self):
        self.load_root()
        for object_type in self.object_types:
            self.load_object_level(object_type)
            self.load_candidate_level(object_type)

    def get_solve_methods(self,object_type):
        method = self.parameters[object_type]["method"]
        linemeas_method = self.parameters[object_type]["linemeas_method"]
        methods = []
        if method:
            methods.append(method)
        if linemeas_method:
            methods.append(linemeas_method)
        return methods

    def get_solve_method(self, object_type):
        return self.parameters[object_type]["method"]

    def get_attribute(self,object_type, dataset, attribute, rank = None):
        if object_type:
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

    def get_classification_type(self):
        return self.root_results["classification"]["Type"]

    def get_classification(self):
        return self.root_results["classification"]
    
    def get_fitted_continuum_by_rank(self, object_type, rank):
        return self.object_dataframes[object_type]["continuum"][rank]

    def get_fitted_model_by_rank(self, object_type, rank, method):
        if method == "LineMeasSolve":
            return self.object_dataframes[object_type]["linemeas_model"]
        else:
            return self.object_dataframes[object_type]["model"][rank]

    def get_fitted_lines_by_rank(self, object_type, rank, method):
        if method == "LineMeasSolve":
            return self.object_dataframes[object_type]["linemeas"]
        return self.object_dataframes[object_type]["fitted_lines"][rank]

    def get_candidate_group_name(self,rank):
        return "candidate" + chr(rank+65) # 0=A, 1=B,....

    def get_nb_candidates(self,object_type):
        return len(self.object_dataframes[object_type]["model"])
    
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

# TODO this method should not exit, dataset size should be inferred in an other way
    def get_dataset_size(self, object_type, dataset, rank = None):
        if rank is None:
            if dataset in self.object_dataframes[object_type]:
                return self.object_dataframes[object_type][dataset].index.size
            else:
                return 1
        else:
            if type(self.object_dataframes[object_type][dataset]) == pd.DataFrame:
                return 1
            else:
                return self.object_dataframes[object_type][dataset][rank].index.size
