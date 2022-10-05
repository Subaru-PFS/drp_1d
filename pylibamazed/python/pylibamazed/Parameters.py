import pandas as pd
import json
import h5py
import os
from pylibamazed.Exception import APIException
from pylibamazed.redshift import ErrorCode

class Parameters():

    def __init__(self, parameters, config):
        self.parameters = parameters
        self.config = config
        self.calibration_dir = config["calibration_dir"]
        self.parameters["calibrationDir"]=config["calibration_dir"]
        self.check_reliability()
        
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
       
    def check_reliability(self):
        for ot in self.get_objects():
            if self.reliability_enabled(ot):
                try:
                    # to avoid annoying messages about gpu/cuda availability
                    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
                    from tensorflow.keras import models
                    import keras
                    model_path = os.path.join(self.calibration_dir,
                                              self.parameters[ot]["reliability_model"])
                    model_ha = h5py.File(model_path).attrs
                    keras_model_version = model_ha["keras_version"].split(".")
                    keras_system_version = keras.__version__.split(".")
                    if keras_model_version[0] > keras_system_version[0]:
                        raise APIException(ErrorCode.RELIABILITY_NEEDS_TENSORFLOW,
                                           f"Tensorflow major version >= {keras_model_version[0]} required")
                    redshift_range = [model_ha["zrange_min"],
                                      model_ha["zrange_max"]]
                    redshift_range_step = model_ha["zrange_step"]
                    s_redshift_range = self.parameters[ot]["redshiftrange"]
                    s_redshift_range_step = self.parameters[ot]["redshiftstep"]
                    if s_redshift_range != redshift_range:
                        raise APIException(ErrorCode.INCOMPATIBLE_PDF_MODELSHAPES,f"redshift range of reliability model must be identical to solver one : {redshift_range} != {s_redshift_range}")
                    if s_redshift_range_step != redshift_range_step:
                        raise APIException(ErrorCode.INCOMPATIBLE_PDF_MODELSHAPES,f"redshift step of reliability model must be identical to solver one : {redshift_range_step} != {s_redshift_range_step}")
                except ImportError:
                    raise APIException(ErrorCode.RELIABILITY_NEEDS_TENSORFLOW,"Tensorflow is required to compute the reliability")

