import pandas as pd
import json


class Parameters():

    def __init__(self, parameters, config):
        self.parameters = parameters
        self.config = config
        self.calibration_dir = config["calibration_dir"]
        self.parameters["calibrationDir"]=config["calibration_dir"]
        
    def get_solve_methods(self,object_type):
        method = self.parameters[object_type]["method"]
        linemeas_method = self.parameters[object_type]["linemeas_method"]
        methods = []
        if method:
            methods.append(method)
        if linemeas_method:
            methods.append(linemeas_method)
        return methods

    def get_linemodel_methods(self, object_type):
        methods = []
        if self.parameters[object_type]["linemeas_method"]:
            methods.append(self.parameters[object_type]["linemeas_method"]) 
        if self.parameters[object_type]["method"] and self.parameters[object_type]["method"] == "LineModelSolve":
            methods.append(self.parameters[object_type]["method"]) 
        return methods
    
    def check_lmskipsecondpass(self, object_type):
        if self.parameters[object_type]["method"]:
            if self.parameters[object_type]["method"] != "LineModelSolve":
                return False 
            else:
                return self.parameters[object_type][self.parameters[object_type]["method"]]["linemodel"]["skipsecondpass"]
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
            columns = self.config["linemeas_catalog_columns"][object_type]
            redshift_ref = lm[columns["Redshift"]].iloc[0]
            velocity_abs = lm[columns["VelocityAbsorption"]].iloc[0]
            velocity_em = lm[columns["VelocityEmission"]].iloc[0]
            #            zlog.LogInfo("Linemeas on " + spectrum_reader.source_id + " with redshift " + str(redshift_ref))
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
        

