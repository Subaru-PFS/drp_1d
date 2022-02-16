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

from pylibamazed.CalibrationLibrary import CalibrationLibrary
from pylibamazed.ResultStoreOutput import ResultStoreOutput

from pylibamazed.redshift import (CProcessFlowContext,
                                  CLineModelSolve,
                                  CMethodTemplateFittingSolve,
                                  CMethodTplcombinationSolve,
                                  CClassificationSolve,
                                  CReliabilitySolve,
                                  CLineMeasSolve)
import pandas as pd
import json


class Context:

    def __init__(self, config, parameters):
        self.calibration_library = CalibrationLibrary(parameters, config["calibration_dir"])
        self.calibration_library.load_all()
        self.process_flow_context = CProcessFlowContext()
        self.config = config
        self.parameters = parameters
        self.parameters["calibrationDir"]=config["calibration_dir"]
        for object_type in self.parameters["objects"]:
            if object_type in self.calibration_library.line_catalogs:
                self.process_flow_context.setLineCatalog(object_type,
                                                         self.calibration_library.line_catalogs[object_type])
            if object_type in self.calibration_library.line_ratio_catalog_lists:
                self.process_flow_context.setLineRatioCatalogCatalog(object_type,
                                                                     self.calibration_library.line_ratio_catalog_lists[
                                                                         object_type])
            if "linemeasmethod" in parameters[object_type] and parameters[object_type]["linemeasmethod"]:
                self.process_flow_context.setLineCatalog(object_type, self.calibration_library.line_catalogs["linemeas"])

        self.process_flow_context.setTemplateCatalog(self.calibration_library.templates_catalogs["all"])
        self.process_flow_context.setPhotBandCatalog(self.calibration_library.photometric_bands)

    def run(self, spectrum_reader):
        spectrum_reader.init()
        self.process_flow_context.setSpectrum(spectrum_reader.get_spectrum())
        if "linemeascatalog" in self.config:
            for object_type in self.config["linemeascatalog"].keys():
                lm = pd.read_csv(self.config["linemeascatalog"][object_type], sep='\t', dtype={'ProcessingID': object})
                redshift_ref = lm[lm.ProcessingID == spectrum_reader.source_id]["g.Redshift"].iloc[0]
                velocity_abs = lm[lm.ProcessingID == spectrum_reader.source_id]["g.VelocityAbsorption"].iloc[0]
                velocity_em = lm[lm.ProcessingID == spectrum_reader.source_id]["g.VelocityEmission"].iloc[0]
    #            zlog.LogInfo("Linemeas on " + spectrum_reader.source_id + " with redshift " + str(redshift_ref))
                self.parameters[object_type]["redshiftref"] = redshift_ref
                self.parameters[object_type]["linemeassolve"]["linemodel"]["velocityabsorption"] = velocity_abs
                self.parameters[object_type]["linemeassolve"]["linemodel"]["velocityemission"] = velocity_em
        self.process_flow_context.LoadParameterStore(json.dumps(self.parameters))
        self.process_flow_context.Init()

        enable_classification = False
        enable_reliability = False
        for object_type in self.parameters["objects"]:
            method = self.parameters[object_type]["method"]
            if method:
                self.run_method(object_type, method)
                enable_classification = True
                enable_reliability = True
            if "linemeas_method" in self.parameters[object_type] and self.parameters[object_type]["linemeas_method"]:
                if "linemeascatalog" not in self.config:
                    output = ResultStoreOutput(self.process_flow_context.GetResultStore(),
                                               self.parameters,
                                               auto_load=False)

                    self.parameters[object_type]["redshiftref"] = output.get_attribute_from_result_store("Redshift",
                                                                                                         "galaxy",
                                                                                                         0)
                    vel_a = output.get_attribute_from_result_store("VelocityAbsorption",
                                                                   "galaxy",
                                                                   0)
                    vel_e = output.get_attribute_from_result_store("VelocityEmission",
                                                                   "galaxy",
                                                                   0)
                    self.parameters[object_type]["linemeassolve"]["linemodel"]["velocityabsorption"] = vel_a
                    self.parameters[object_type]["linemeassolve"]["linemodel"]["velocityemission"] = vel_e
                    self.process_flow_context.LoadParameterStore(json.dumps(self.parameters))
                self.run_method(object_type, self.parameters[object_type]["linemeas_method"])
        if enable_classification:
            self.run_method("classification", "classification")
        if enable_reliability:
            self.run_method("reliability", "reliability")

        return ResultStoreOutput(self.process_flow_context.GetResultStore(), self.parameters)

    def run_method(self, object_type, method):

        if method == "linemodelsolve":
            solver = CLineModelSolve(self.process_flow_context.m_ScopeStack,
                                     object_type,
                                     self.config["calibration_dir"])
        elif method == "templatefittingsolve":
            solver = CMethodTemplateFittingSolve(self.process_flow_context.m_ScopeStack, object_type)
        elif method == "tplcombinationsolve":
            solver = CMethodTplcombinationSolve(self.process_flow_context.m_ScopeStack, object_type)
        elif method == "linemeassolve":
            solver = CLineMeasSolve(self.process_flow_context.m_ScopeStack,
                                    object_type,
                                    self.config["calibration_dir"])
        elif method == "classification":
            solver = CClassificationSolve(self.process_flow_context.m_ScopeStack, object_type)
        elif method == "reliability":
            solver = CReliabilitySolve(self.process_flow_context.m_ScopeStack, object_type)
        else:
            raise Exception("Unkown method " + method)
        solver.Compute(self.process_flow_context)
