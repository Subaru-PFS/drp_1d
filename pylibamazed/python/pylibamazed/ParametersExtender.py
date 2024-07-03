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

import copy

from pylibamazed.ParametersAccessor import ParametersAccessor
from pylibamazed.ParametersChecker import ParametersChecker


class ParametersExtender:
    default_continuum_removal_method = "default"
    default_median_kernel_width = -1
    default_median_even_reflection = True

    default_continuum_removal_section = {
        "method": default_continuum_removal_method,
        "medianKernelWidth": default_median_kernel_width,
        "medianEvenReflection": default_median_even_reflection,
    }

    def __init__(self, Accessor=ParametersAccessor, Checker=ParametersChecker):
        self.Accessor = Accessor
        self.Checker = Checker

    def extend(self, parameters: dict):
        self.initial_accessor = self.Accessor(parameters)
        self.extended_accessor = self.Accessor(copy.deepcopy(parameters))
        self.checker = self.Checker(parameters)

        self._extend_continuum_removal()
        self._extend_template_catalog_continuum_removal()
        self._extend_continuum_reestimation()
        self._extend_useloglambdasampling()
        self._extend_mediankernel()

        extended_params = self.extended_accessor.parameters

        self.initial_accessor = None
        self.extended_accessor = None
        self.checker = None

        return extended_params

    def _extend_continuum_removal(self):
        is_present = self.checker.continuum_removal_presence_condition()
        if not is_present:
            self.extended_accessor.parameters["continuumRemoval"] = self.default_continuum_removal_section

    def _extend_template_catalog_continuum_removal(self):
        is_present = self.checker.template_catalog_continuum_removal_presence_condition()
        if not is_present:
            self.extended_accessor.get_template_catalog_section(True)
            self.extended_accessor.get_template_catalog_section()[
                "continuumRemoval"
            ] = self.default_continuum_removal_section

    def _extend_continuum_reestimation(self):
        for spectrum_model in self.initial_accessor.get_spectrum_models([]):
            is_present = self.checker.linemodelsolve_continuumreestimation_presence_condition(spectrum_model)
            if not is_present:
                self.extended_accessor.get_linemodel_solve_linemodel_section(True)
                self.extended_accessor.get_linemodel_solve_linemodel_section(spectrum_model, True)[
                    "continuumReestimation"
                ] = "no"

    def _extend_useloglambdasampling(self):
        for spectrum_model in self.initial_accessor.get_spectrum_models([]):
            is_present = self.checker.useloglambdasampling_presence_condition(spectrum_model)
            if not is_present:
                self.extended_accessor.get_linemodel_solve_linemodel_section(True)
                self.extended_accessor.get_linemodel_solve_linemodel_section(spectrum_model, True)[
                    "useLogLambdaSampling"
                ] = False

    def _extend_mediankernel(self):
        for fromTemplateCatalog in [False, True]:
            is_present = self.checker.median_kernel_presence_condition(fromTemplateCatalog)
            if not is_present:
                self.extended_accessor.get_continuum_removal_section(fromTemplateCatalog, True)[
                    "medianKernelWidth"
                ] = self.default_median_kernel_width
                self.extended_accessor.get_continuum_removal_section(fromTemplateCatalog)[
                    "medianEvenReflection"
                ] = self.default_median_even_reflection
