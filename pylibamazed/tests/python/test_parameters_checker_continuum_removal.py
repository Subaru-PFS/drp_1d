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
import pytest
from pylibamazed.Exception import APIException
from pylibamazed.ParametersChecker import ParametersChecker
from tests.python.utils import (DictUtils, WarningUtils,
                                make_parameter_dict_at_redshift_solver_level)


class TestContinuumRemoval:

    class TestMethod:
        def test_ok_if_method_is_IrregularSamplingMedian_with_medianKernelWidth_medianEvenReflection(
                self, zflag):
            param_dict = make_parameter_dict_at_redshift_solver_level(**{
                "method": "templateFittingSolve",
                "templateFittingSolve": {
                    "spectrum": {
                        "component": "sth"
                    }
                }
            })
            param_dict["galaxy"]["templateDir"] = "sth"
            param_dict["templateCatalog]"] = {
                "continuumRemoval": {
                    "method": "irregularSamplingMedian",
                    "medianKernelWidth": 1,
                    "medianEvenReflection": 1,
                }
            }
            param_dict["continuumRemoval"] = {
                "method": "irregularSamplingMedian",
                "medianKernelWidth": 1,
                "medianEvenReflection": 1,
            }
            param_dict["templateCatalog"] = {"continuumRemoval": {
                "method": "irregularSamplingMedian",
                "medianKernelWidth": 1,
                "medianEvenReflection": 1,
            }}
            ParametersChecker(param_dict).custom_check()
            assert not WarningUtils.has_any_warning()

        @pytest.mark.parametrize('nesting', [None, "templateCatalog"])
        def test_error_if_method_is_IrregularSamplingMedian_without_medianKernelWidth(self, nesting):
            parametersDict = DictUtils.make_nested({
                "continuumRemoval": {
                    "method": "irregularSamplingMedian",
                    "medianEvenReflection": 1
                }
            }, nesting)
            with pytest.raises(APIException, match="Missing parameter continuumRemoval medianKernelWidth"):
                ParametersChecker(parametersDict).custom_check()

        @pytest.mark.parametrize('nesting', [None, "templateCatalog"])
        def test_error_if_method_is_IrregularSamplingMedian_without_medianEvenReflection(self, nesting):
            parametersDict = DictUtils.make_nested({
                "continuumRemoval": {
                    "method": "irregularSamplingMedian",
                    "medianKernelWidth": 1
                }
            }, nesting)
            with pytest.raises(APIException, match="Missing parameter continuumRemoval medianEvenReflection"):
                ParametersChecker(parametersDict).custom_check()

    class TestMedianEvenWidth:
        @pytest.mark.parametrize('nesting', [None, "templateCatalog"])
        def test_warning_present_and_method_is_not_IrregularSamplingMedian(self, zflag, nesting):
            parametersDict = DictUtils.make_nested({
                "continuumRemoval": {
                    "method": "sth",
                    "medianKernelWidth": 1,
                }
            }, nesting)

            ParametersChecker(parametersDict).custom_check()
            assert WarningUtils.has_any_warning()

    class TestMedianEvenReflection:
        @pytest.mark.parametrize('nesting', [None, "templateCatalog"])
        def test_warning_present_and_method_is_not_IrregularSamplingMedian(self, zflag, nesting):
            parametersDict = DictUtils.make_nested({
                "continuumRemoval": {
                    "method": "sth",
                    "medianEvenReflection": 1
                }
            }, nesting)

            ParametersChecker(parametersDict).custom_check()
            assert WarningUtils.has_any_warning()


class TestBaseContinuumRemoval:
    class TestSectionPresence:
        def test_error_if_absent_but_templateFittingSolve_spectrum_component_is_not_raw(self):
            kwargs = {
                "method": "templateFittingSolve",
                "templateFittingSolve": {
                    "spectrum": {
                        "component": "sth"
                    }
                }
            }
            parametersDict = make_parameter_dict_at_redshift_solver_level(**kwargs)
            parametersDict["galaxy"]["templateDir"] = "sth"
            with pytest.raises(APIException, match="Missing parameter continuumRemoval"):
                ParametersChecker(parametersDict).custom_check()

        def test_error_if_absent_but_tplCombinationSolve_spectrum_component_is_not_raw(self):
            kwargs = {
                "method": "tplCombinationSolve",
                "tplCombinationSolve": {
                    "spectrum": {
                        "component": "sth"
                    }
                }
            }
            parametersDict = make_parameter_dict_at_redshift_solver_level(**kwargs)
            parametersDict["galaxy"]["templateDir"] = "sth"
            with pytest.raises(APIException, match="Missing parameter continuumRemoval"):
                ParametersChecker(parametersDict).custom_check()

        def test_error_if_absent_but_lineModelSolve_continuumComponent_is_fromSpectrum(self):
            kwargs = {
                "method": "lineModelSolve",
                "lineModelSolve": {
                    "lineModel": {
                        "continuumComponent": "fromSpectrum"
                    }
                }
            }
            parametersDict = make_parameter_dict_at_redshift_solver_level(**kwargs)
            with pytest.raises(APIException, match="Missing parameter continuumRemoval"):
                ParametersChecker(parametersDict).custom_check()

        def test_necessity_for_any_object(self):
            # First object which does not need continuumRemoval
            parametersDict = make_parameter_dict_at_redshift_solver_level()

            # Second object whivh needs it
            parametersDict2 = make_parameter_dict_at_redshift_solver_level(None, "qso", **{
                "method": "templateFittingSolve",
                "templateFittingSolve": {
                    "spectrum": {
                        "component": "sth"
                    }
                }
            })
            parametersDict = parametersDict | parametersDict2
            parametersDict["qso"]["templateDir"] = "sth"
            with pytest.raises(APIException, match="Missing parameter continuumRemoval"):
                ParametersChecker(parametersDict).custom_check()
