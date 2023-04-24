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
import json
from pylibamazed.Exception import APIException
from pylibamazed.redshift import ErrorCode

class Parameters:
    def __init__(self, parameters, config=None):
        self.parameters = parameters
        self.config = config
        
    def get_solve_methods(self,object_type):
        method = self.get_solve_method(object_type)
        linemeas_method = self.get_linemeas_method(object_type)
        methods = []
        if method:
            methods.append(method)
        if linemeas_method:
            methods.append(linemeas_method)
        return methods
        
    def get_redshift_sampling(self,object_type):
        return self.parameters[object_type]["redshiftsampling"]

    def get_linemodel_methods(self, object_type):
        methods = []
        linemeas_method = self.get_linemeas_method(object_type)
        solve_method = self.get_solve_method(object_type)
        if linemeas_method:
            methods.append(linemeas_method) 
        if solve_method and solve_method == "LineModelSolve":
            methods.append(solve_method) 
        return methods
    
    def check_lmskipsecondpass(self, object_type):
        solve_method = self.get_solve_method(object_type)
        if solve_method:
            if solve_method != "LineModelSolve":
                return False 
            else:
                return self.parameters[object_type][solve_method]["linemodel"]["skipsecondpass"]
        return False

    def get_solve_method(self, object_type):
        return self.parameters[object_type]["method"]

    def get_linemeas_method(self, object_type):
        return self.parameters[object_type]["linemeas_method"]

    def get_objects(self):
        return self.parameters["objects"]

    def load_linemeas_parameters_from_catalog(self, source_id):
        for object_type in self.config["linemeascatalog"].keys():
            lm = pd.read_csv(self.config["linemeascatalog"][object_type], sep='\t', dtype={'ProcessingID': object})
            lm = lm[lm.ProcessingID == source_id]
            if lm.empty:
                raise APIException(ErrorCode.INVALID_PARAMETER,f"Uncomplete linemeas catalog, {source_id} missing")
            
            columns = self.config["linemeas_catalog_columns"][object_type]
            redshift_ref = float(lm[columns["Redshift"]].iloc[0])
            velocity_abs = float(lm[columns["VelocityAbsorption"]].iloc[0])
            velocity_em = float(lm[columns["VelocityEmission"]].iloc[0])
            self.parameters[object_type]["redshiftref"] = redshift_ref
            self.parameters[object_type]["LineMeasSolve"]["linemodel"]["velocityabsorption"] = velocity_abs
            self.parameters[object_type]["LineMeasSolve"]["linemodel"]["velocityemission"] = velocity_em

    def load_linemeas_parameters_from_result_store(self, output, object_type):

        redshift = output.get_attribute_from_source(object_type,
                                                    self.get_solve_method(object_type),
                                                    "model_parameters",
                                                    "Redshift",
                                                    0)
        self.parameters[object_type]["redshiftref"] = redshift
        vel_a = output.get_attribute_from_source(object_type,
                                                 self.get_solve_method(object_type),
                                                 "model_parameters",
                                                 "VelocityAbsorption",
                                                 0)
        vel_e = output.get_attribute_from_source(object_type,
                                                 self.get_solve_method(object_type),
                                                 "model_parameters",
                                                 "VelocityEmission",
                                                 0)
        self.parameters[object_type]["LineMeasSolve"]["linemodel"]["velocityabsorption"] = vel_a
        self.parameters[object_type]["LineMeasSolve"]["linemodel"]["velocityemission"] = vel_e
        
    def get_json(self):
        return json.dumps(self.parameters)
    
    def reliability_enabled(self, object_type):
        return self.parameters[object_type].get("enable_reliability")

    def lineratio_catalog_enabled(self, object_type):
        if self.get_solve_method(object_type) == "LineModelSolve" :
            return self.parameters[object_type]["LineModelSolve"]["linemodel"]["lineRatioType"] == "tplratio"
        else:
            return False
        
    def stage_enabled(self, object_type, stage):
        if stage == "redshift_solver":
            return self.get_solve_method(object_type) is not None
        elif stage == "linemeas_solver":
            return self.get_linemeas_method(object_type) is not None
        elif stage == "linemeas_catalog_load":
            return self.get_linemeas_method(object_type) is not None and self.get_solve_method(object_type) is None
        elif stage == "reliability_solver":
            return self.reliability_enabled(object_type)
        elif stage == "sub_classif_solver":
            return self.lineratio_catalog_enabled(object_type)
        else:
            raise Exception("Unknown stage {stage}") 
       
    def get_filters(self):
        return self.parameters.get("filters")
