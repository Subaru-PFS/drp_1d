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
import os.path

from pylibamazed.CalibrationLibrary import CalibrationLibrary
from pylibamazed.ResultStoreOutput import ResultStoreOutput

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
        self.parameters = parameters
        self.parameters["calibrationDir"]=config["calibration_dir"]
        self.reliability_models = {}

    def init_context(self):
        self.process_flow_context = CProcessFlowContext()
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

            if parameters[object_type].get("enable_reliability"):
                # Load model
                try:
                    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
                    from tensorflow.keras import models
                except ImportError:
                    import warnings
                    warnings.warn("Tensorflow is required to compute the reliability")
                else:
                    model_path = os.path.join(self.parameters["calibrationDir"],
                                              parameters[object_type]["reliability_model"])
                    model = models.load_model(model_path)
                    self.reliability_models[object_type] = model

        self.process_flow_context.setTemplateCatalog(self.calibration_library.templates_catalogs["all"])
        self.process_flow_context.setPhotBandCatalog(self.calibration_library.photometric_bands)

    def run(self, spectrum_reader):
        self.init_context()
        spectrum_reader.init()
        self.process_flow_context.setSpectrum(spectrum_reader.get_spectrum())
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
        enable_reliability = False
        for object_type in self.parameters["objects"]:
            method = self.parameters[object_type]["method"]
            if method:
                self.run_method(object_type, method)
                enable_classification = True
                enable_reliability = True
            if "linemeas_method" in self.parameters[object_type] and self.parameters[object_type]["linemeas_method"]:
                if not self.config["linemeascatalog"]:
                    output = ResultStoreOutput(self.process_flow_context.GetResultStore(),
                                               self.parameters,
                                               auto_load=False)

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

            if self.parameters[object_type].get("enable_reliability"):
                output = ResultStoreOutput(self.process_flow_context.GetResultStore(),
                                           self.parameters,
                                           auto_load=False)
                pdf = output.get_attribute_from_result_store("PDFProbaLog", object_type, 0)
                model = self.reliability_models[object_type]
                if pdf.shape[0] != model.input_shape[1]:
                    raise ValueError('PDF and model shapes are not compatible')
                # The model needs a PDF, not LogPDF
                reliability = model.predict(np.exp(pdf[None, :, None]))[0, 1]
                output.object_results[object_type]['reliability'] = reliability

        if enable_classification:
            self.run_method("classification", "ClassificationSolve")
        if enable_reliability:
            self.run_method("reliability", "ReliabilitySolve")

        rso = ResultStoreOutput(self.process_flow_context.GetResultStore(), self.parameters)
        del self.process_flow_context
        return rso

    def run_method(self, object_type, method):

        if "C" + method not in globals():
            raise Exception("Unkown method " + method)
        solver_method = globals()["C" + method]
        solver = solver_method(self.process_flow_context.m_ScopeStack,
                               object_type)
        solver.Compute(self.process_flow_context)


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
