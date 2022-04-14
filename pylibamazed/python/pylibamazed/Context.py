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
from distutils.log import error
import os.path

from pylibamazed.CalibrationLibrary import CalibrationLibrary
from pylibamazed.ResultStoreOutput import ResultStoreOutput
from pylibamazed.Reliability import Reliability

from pylibamazed.redshift import (CProcessFlowContext,
                                  CLineModelSolve,
                                  CTemplateFittingSolve,
                                  CTplcombinationSolve,
                                  CClassificationSolve,
                                  CReliabilitySolve,
                                  CLineMeasSolve,
                                  CLineMatchingSolve,
                                  CFlagWarning)
import pandas as pd
import json
import os

zflag = CFlagWarning.GetInstance()
class Context:

    def __init__(self, config, parameters):
        _check_config(config)
        self.calibration_library = CalibrationLibrary(parameters, config["calibration_dir"])
        self.calibration_library.load_all()
        self.process_flow_context = None
        self.config = config
        if "linemeascatalog" not in self.config:
            self.config["linemeascatalog"] = {}
        else:
            _check_LinemeasValidity(config, parameters)
        self.parameters = parameters
        self.parameters["calibrationDir"]=config["calibration_dir"]
        self.extended_results = config["extended_results"]

    def init_context(self):
        self.process_flow_context = CProcessFlowContext.GetInstance()
        for object_type in self.parameters["objects"]:
            if object_type in self.calibration_library.line_catalogs:
                self.process_flow_context.setLineCatalog(object_type,
                                                         self.calibration_library.line_catalogs[object_type])
            if object_type in self.calibration_library.line_ratio_catalog_lists:
                self.process_flow_context.setLineRatioCatalogCatalog(object_type,
                                                                     self.calibration_library.line_ratio_catalog_lists[
                                                                         object_type])
            if "linemeas_method" in self.parameters[object_type] and self.parameters[object_type]["linemeas_method"]:
                self.process_flow_context.setLineCatalog(object_type, self.calibration_library.line_catalogs[object_type])

        self.process_flow_context.setTemplateCatalog(self.calibration_library.templates_catalogs["all"])
        self.process_flow_context.setPhotBandCatalog(self.calibration_library.photometric_bands)
        self.process_flow_context.setfluxCorrectionMeiksin(self.calibration_library.meiksin)
        self.process_flow_context.setfluxCorrectionCalzetti(self.calibration_library.calzetti)


    def run(self, spectrum_reader):
        self.init_context()
        spectrum_reader.init()
        if self.config.get("linemeascatalog"):
            for object_type in self.config["linemeascatalog"].keys():
                lm = pd.read_csv(self.config["linemeascatalog"][object_type], sep='\t', dtype={'ProcessingID': object})
                lm = lm[lm.ProcessingID == spectrum_reader.source_id]
                columns = self.config["linemeas_catalog_columns"][object_type]
                redshift_ref = lm[columns["Redshift"]].iloc[0]
                velocity_abs = lm[columns["VelocityAbsorption"]].iloc[0]
                velocity_em = lm[columns["VelocityEmission"]].iloc[0]
    #            zlog.LogInfo("Linemeas on " + spectrum_reader.source_id + " with redshift " + str(redshift_ref))
                self.parameters[object_type]["redshiftref"] = redshift_ref
                self.parameters[object_type]["LineMeasSolve"]["linemodel"]["velocityabsorption"] = velocity_abs
                self.parameters[object_type]["LineMeasSolve"]["linemodel"]["velocityemission"] = velocity_em
        self.process_flow_context.LoadParameterStore(json.dumps(self.parameters))
        self.process_flow_context.Init()

        #store flag in root object
        resultStore = self.process_flow_context.GetResultStore()
        resultStore.StoreFlagResult( "context_warningFlag", zflag.getBitMask())
        zflag.resetFlag()

        enable_classification = False
        reliabilities = dict()
        for object_type in self.parameters["objects"]:
            method = self.parameters[object_type]["method"]
            if method:
                self.run_method(object_type, method)
                enable_classification = True

            if "linemeas_method" in self.parameters[object_type] and self.parameters[object_type]["linemeas_method"]:
                if not self.config["linemeascatalog"]:
                    output = ResultStoreOutput(self.process_flow_context.GetResultStore(),
                                               self.parameters,
                                               auto_load=False,
                                               extended_results = self.extended_results)

                    self.parameters[object_type]["redshiftref"] = output.get_attribute_from_result_store("Redshift",
                                                                                                         object_type,
                                                                                                         0)
                    vel_a = output.get_attribute_from_result_store("VelocityAbsorption",
                                                                   object_type,
                                                                   0)
                    vel_e = output.get_attribute_from_result_store("VelocityEmission",
                                                                   object_type,
                                                                   0)
                    self.parameters[object_type]["LineMeasSolve"]["linemodel"]["velocityabsorption"] = vel_a
                    self.parameters[object_type]["LineMeasSolve"]["linemodel"]["velocityemission"] = vel_e
                    self.process_flow_context.LoadParameterStore(json.dumps(self.parameters))
                self.run_method(object_type, self.parameters[object_type]["linemeas_method"])

            if self.parameters[object_type].get("enable_reliability") and object_type in self.calibration_library.reliability_models:
                rel = Reliability(object_type, self.parameters,self.calibration_library, extended_results = self.extended_results)
                reliabilities[object_type] = rel.Compute(self.process_flow_context)

        if enable_classification:
            self.run_method("classification", "ClassificationSolve")

        rso = ResultStoreOutput(self.process_flow_context.GetResultStore(), self.parameters, extended_results = self.extended_results)
        for object_type in reliabilities.keys():
            rso.object_results[object_type]['reliability'] = dict()
            rso.object_results[object_type]['reliability']['Reliability'] = reliabilities[object_type]

        self.process_flow_context.reset()
        return rso

    def run_method(self, object_type, method):

        if "C" + method not in globals():
            raise Exception("Unkown method " + method)
        solver_method = globals()["C" + method]
        solver = solver_method(self.process_flow_context.m_ScopeStack,
                               object_type)
        solver.Compute()


def _check_config(config):
    if "calibration_dir" not in config:
        raise Exception("Config must contain 'calibration_dir' key")
    if not os.path.exists(config["calibration_dir"]):
        raise Exception("Calibration directory {} does not exist".format(config["calibration_dir"]))
    if "linemeascatalog" in config:
        if "linemeas_catalog_columns" not in config:
            raise Exception("With a linemeas catalog Config must contain a linemeas_catalog_columns key")
        for object_type in config["linemeascatalog"].keys():
            if object_type not in config["linemeas_catalog_columns"]:
                raise Exception("Config['linemeas_catalog_columns'] misses category {}".format(object_type))
            for attr in ["Redshift", "VelocityAbsorption", "VelocityEmission"]:
                if attr not in config["linemeas_catalog_columns"][object_type]:
                    raise Exception("Config['linemeas_catalog_columns'][{}] misses attribute".format(object_type, attr))

def _check_LinemeasValidity(config, parameters):
    if not config["linemeascatalog"]:
        return
    for object_type in parameters["objects"]:
        method = parameters[object_type]["method"]
        if method == "LineModelSolve":
            if "linemeas_method" in parameters[object_type] and parameters[object_type]["linemeas_method"]:
                raise Exception("Cannot run LineMeasSolve from catalog when sequencial processing is selected simultaneously.")
