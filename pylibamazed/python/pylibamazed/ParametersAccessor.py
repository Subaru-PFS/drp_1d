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

from typing import List


class ParametersAccessor:
    def __init__(self, parameters: dict):
        self.parameters = parameters

    def get_lambda_range(self, obs_id=""):
        """Depending on multiobs method, lambda range is not of the same type:
        - mono obs or merge => lambdarange is a range [min, max]
        - full => lambdarange is a dict of ranges {obs_id: [min, max], ...}
        """
        multiobs = self.get_multiobs_method()
        if multiobs in ["", "merge"]:
            return self.parameters["lambdarange"]
        else:
            return self.parameters["lambdarange"][obs_id]

    def get_airvacuum_method(self):
        return self.parameters.get("airvacuum_method", "")

    def get_photometry_transmission_dir(self):
        return self.parameters.get("photometryTransmissionDir")

    def get_photometry_bands(self) -> List[str]:
        return self.parameters.get("photometryBand",[])

    def get_multiobs_method(self):
        return self.parameters.get("multiobsmethod")

    def get_objects(self, default=None):
        return self.parameters.get("objects", default)

    def get_object_section(self, object_type) -> dict:
        return self.parameters.get(object_type)

    def get_solve_method(self, object_type: str) -> str:
        return self.get_object_section(object_type).get("method")

    def get_linemeas_method(self, object_type: str) -> str:
        return self.get_object_section(object_type).get("linemeas_method")

    def get_linemeas_dzhalf(self, object_type: str) -> float:
        return self.get_object_section(object_type).get("linemeas_dzhalf")

    def get_redshiftrange(self, obejct_type: str) -> List[float]:
        return self.get_object_section(obejct_type).get("redshiftrange")

    def get_redshiftstep(self, object_type: str) -> float:
        return self.get_object_section(object_type).get("redshiftstep")

    def get_linemeas_redshiftstep(self, object_type: str) -> float:
        return self.get_object_section(object_type).get("linemeas_redshiftstep")

    def get_reliability_enabled(self, object_type: str) -> str:
        return self.get_object_section(object_type).get("enable_reliability")

    def get_reliability_model(self, object_type: str) -> str:
        return self.get_object_section(object_type).get("reliability_model")

    def get_template_dir(self, object_type: str):
        return self.get_object_section(object_type).get("template_dir")

    def photometry_is_enabled(self):
        for obj in self.get_objects([]):
            method = self.get_solve_method(obj)
            method_params: dict = (self.parameters[obj]).get(method)
            if method == "LineModelSolve":
                method_params = self._get_on_None(method_params, "linemodel")
            if self._get_on_None(method_params, "enablephotometry", False):
                return True
        return False

    def get_additional_cols(self, default=None) -> List[str]:
        return self.parameters.get("additional_cols") or default

    def get_filters(self, default=None):
        return self.parameters.get("filters", default)

    def get_lsf(self) -> dict:
        return self.parameters.get("LSF")

    def get_lsf_type(self):
        return self._get_on_None(self.get_lsf(), "LSFType")

    def get_lsf_width(self):
        return self._get_on_None(self.get_lsf(), "width")

    def get_lsf_resolution(self):
        return self._get_on_None(self.get_lsf(), "resolution")

    def get_lsf_sourcesize(self):
        return self._get_on_None(self.get_lsf(), "sourcesize")

    def get_lsf_width_file_name(self):
        return self._get_on_None(self.get_lsf(), "GaussianVariablewidthFileName")

    def get_continuum_removal(self, nesting: str = None):
        dict_to_search_on = self.parameters
        if nesting:
            dict_to_search_on = dict_to_search_on.get(nesting)
        continuum_removal = self._get_on_None(dict_to_search_on, "continuumRemoval")
        return continuum_removal

    def get_continuum_removal_method(self, nesting: str = None):
        return self._get_on_None(self.get_continuum_removal(nesting), "method")

    def get_continuum_removal_median_kernel_width(self, nesting: str = None):
        return self._get_on_None(self.get_continuum_removal(nesting), "medianKernelWidth")

    def get_continuum_median_kernel_reflection(self, nesting: str = None):
        return self._get_on_None(self.get_continuum_removal(nesting), "medianEvenReflection")

    def _get_on_None(self, dict, key, default=None):
        if dict is None:
            return default
        else:
            return dict.get(key, default)

    def get_templateFittingSolve_section(self, object_type: str) -> dict:
        return self._get_on_None(self.get_object_section(object_type), "TemplateFittingSolve")

    def get_templateFittingSolve_ism(self, object_type: str) -> bool:
        return self._get_on_None(self.get_templateFittingSolve_section(object_type), "ismfit")

    def get_templateFittingSolve_photometry_enabled(self, object_type: str) -> bool:
        return self._get_on_None(self.get_templateFittingSolve_section(object_type), "enablephotometry")

    def get_templateFittingSolve_photometry_section(self, object_type: str) -> dict:
        return self._get_on_None(self.get_templateFittingSolve_section(object_type), "photometry")

    def get_templateFittingSolve_photometry_weight(self, object_type: str) -> dict:
        return self._get_on_None(self.get_templateFittingSolve_photometry_section(object_type), "weight")

    def get_templateCombinationSolve_section(self, object_type: str) -> dict:
        return self._get_on_None(self.get_object_section(object_type), "TplcombinationSolve")

    def get_templateCombinationSolve_ism(self, object_type: str) -> bool:
        return self._get_on_None(self.get_templateCombinationSolve_section(object_type), "ismfit")

    def get_ebmv_section(self) -> dict:
        return self.parameters.get("ebmv")

    def get_ebmv_count(self):
        return self._get_on_None(self.get_ebmv_section(), "count")

    def get_ebmv_step(self):
        return self._get_on_None(self.get_ebmv_section(), "step")

    def get_ebmv_start(self):
        return self._get_on_None(self.get_ebmv_section(), "start")

    def get_lineModelSolve_section(self, object_type: str) -> dict:
        return self._get_on_None(self.get_object_section(object_type), "LineModelSolve")

    def get_lineModelSolve_linemodel_section(self, object_type: str) -> str:
        return self._get_on_None(self.get_lineModelSolve_section(object_type), "linemodel")

    def get_solve_method_igm_fit(self, object_type: str, solve_method: str) -> bool:
        igmfit = None
        if solve_method == "LineModelSolve":
            igmfit = self.get_lineModelSolve_igmfit(object_type)
        elif solve_method == "TemplateFittingSolve":
            igmfit = self.get_templateFittingSolve_igmfit(object_type)
        elif solve_method == "TplcombinationSolve":
            igmfit = self.get_templateCombinationSolve_igmfit(object_type)
        return igmfit

    def get_lineModelSolve_igmfit(self, object_type: str) -> bool:
        return self._get_on_None(self.get_lineModelSolve_linemodel_section(object_type), "igmfit")

    def get_templateCombinationSolve_igmfit(self, object_type: str) -> bool:
        return self._get_on_None(self.get_templateCombinationSolve_section(object_type), "igmfit")

    def get_templateFittingSolve_igmfit(self, object_type: str) -> bool:
        return self._get_on_None(self.get_templateFittingSolve_section(object_type), "igmfit")

    def get_lineModelSolve_lineRatioType(self, object_type: str) -> str:
        return self._get_on_None(self.get_lineModelSolve_linemodel_section(object_type), "lineRatioType")

    def get_lineModelSolve_rules(self, object_type: str) -> str:
        return self._get_on_None(self.get_lineModelSolve_linemodel_section(object_type), "rules")

    def get_lineModelSolve_tplratio_catalog(self, object_type: str) -> str:
        return self._get_on_None(self.get_lineModelSolve_linemodel_section(object_type), "tplratio_catalog")

    def get_lineModelSolve_tplratio_ismfit(self, object_type: str) -> str:
        return self._get_on_None(self.get_lineModelSolve_linemodel_section(object_type), "tplratio_ismfit")

    def get_lineModelSolve_continuumComponent(self, object_type: str) -> str:
        return self._get_on_None(self.get_lineModelSolve_linemodel_section(object_type), "continuumcomponent")

    def get_lineModelSolve_continuumfit_section(self, object_type: str) -> dict:
        return self._get_on_None(self.get_lineModelSolve_linemodel_section(object_type), "continuumfit")

    def get_lineModelSolve_secondpass_section(self, object_type: str) -> dict:
        return self._get_on_None(self.get_lineModelSolve_linemodel_section(object_type), "secondpass")

    def get_lineModelSolve_secondpass_continuumfit(self, object_type: str) -> str:
        return self._get_on_None(self.get_lineModelSolve_secondpass_section(object_type), "continuumfit")

    def get_lineModelSolve_continuumreestimation(self, object_type: str) -> str:
        return self._get_on_None(
            self.get_lineModelSolve_linemodel_section(object_type),
            "continuumreestimation"
        )

    def get_lineModelSolve_useloglambdasampling(self, object_type: str) -> bool:
        return self._get_on_None(
            self.get_lineModelSolve_linemodel_section(object_type),
            "useloglambdasampling"
        )

    def get_lineModelSolve_continuumfit_ismfit(self, object_type: str) -> bool:
        return self._get_on_None(self.get_lineModelSolve_continuumfit_section(object_type), "ismfit")

    def get_lineModelSolve_continuumfit_fftprocessing(self, object_type: str) -> bool:
        return self._get_on_None(self.get_lineModelSolve_continuumfit_section(object_type), "fftprocessing")

    def get_lineModelSolve_firstpass_section(self, object_type: str) -> bool:
        return self._get_on_None(self.get_lineModelSolve_linemodel_section(object_type), "firstpass")

    def get_lineModelSolve_firstpass_tplratio_ismfit(self, object_type: str) -> bool:
        return self._get_on_None(self.get_lineModelSolve_firstpass_section(object_type), "tplratio_ismfit")

    def get_lineModelSolve_skipsecondpass(self, object_type: str) -> bool:
        return self._get_on_None(
            self.get_lineModelSolve_linemodel_section(object_type), "skipsecondpass")

    def get_lineModelSolve_nsigmasupport(self, object_type: str) -> float:
        return self._get_on_None(self.get_lineModelSolve_linemodel_section(object_type), "nsigmasupport")

    def get_lineModelSolve_improveBalmerFit(self, object_type: str) -> float:
        return self._get_on_None(self.get_lineModelSolve_linemodel_section(object_type), "improveBalmerFit")

    def get_lineMeasSolve_section(self, object_type: str) -> dict:
        return self._get_on_None(self.get_object_section(object_type), "LineMeasSolve")

    def get_lineMeasSolve_linemodel_section(self, object_type: str) -> dict:
        return self._get_on_None(self.get_lineMeasSolve_section(object_type), "linemodel")

    def get_lineMeasSolve_lineRatioType(self, object_type: str) -> str:
        return self._get_on_None(self.get_lineMeasSolve_linemodel_section(object_type), "lineRatioType")

    def get_lineMeasSolve_rules(self, object_type: str) -> str:
        return self._get_on_None(self.get_lineMeasSolve_linemodel_section(object_type), "rules")

    def get_lineMeasSolve_fittingmethod(self, object_type: str) -> str:
        return self._get_on_None(self.get_lineMeasSolve_linemodel_section(object_type), "fittingmethod")

    def get_lineMeasSolve_velocityfit(self, object_type: str) -> str:
        return self._get_on_None(self.get_lineMeasSolve_linemodel_section(object_type), "velocityfit")

    def get_lineMeasSolve_velocityfit_param(self, object_type: str, param: str) -> float:
        return self._get_on_None(self.get_lineMeasSolve_linemodel_section(object_type), param)

    def get_lineMeasSolve_nsigmasupport(self, object_type: str) -> float:
        return self._get_on_None(self.get_lineMeasSolve_linemodel_section(object_type), "nsigmasupport")

    def get_nsigmasupport(self, object_type: str, method: str) -> float:
        nsigmasupport = None
        if method == "LineModelSolve":
            nsigmasupport = self.get_lineModelSolve_nsigmasupport(object_type)
        elif method == "LineMeasSolve":
            nsigmasupport = self.get_lineMeasSolve_nsigmasupport(object_type)
        return nsigmasupport

    def get_lineModelSolve_linecatalog(self, object_type: str) -> str:
        return self._get_on_None(self.get_lineModelSolve_linemodel_section(object_type), "linecatalog")

    def get_lineMeasSolve_linecatalog(self, object_type: str) -> str:
        return self._get_on_None(self.get_lineMeasSolve_linemodel_section(object_type), "linecatalog")

    def get_linecatalog(self, object_type: str, method: str) -> str:
        linecatalog = None
        if method == "LineModelSolve":
            linecatalog = self.get_lineModelSolve_linecatalog(object_type)
        elif method == "LineMeasSolve":
            linecatalog = self.get_lineMeasSolve_linecatalog(object_type)
        return linecatalog

    def get_linemodel_section(self, object_type, method) -> dict:
        linemodel = None
        if method == "LineModelSolve":
            linemodel = self.get_lineModelSolve_linemodel_section(object_type)
        elif method == "LineMeasSolve":
            linemodel = self.get_lineMeasSolve_linemodel_section(object_type)
        return linemodel

    def get_redshift_sampling(self, object_type):
        return self.get_object_section(object_type).get("redshiftsampling")

    def get_observation_ids(self):
        try:
            return list(self.parameters["lambdarange"].keys())
        except Exception:
            return [""]
