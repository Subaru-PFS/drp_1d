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
from pylibamazed.ParametersAccessor import ParametersAccessor
from pylibamazed.ParametersChecker import ParametersChecker
from tests.python.utils import DictUtils, WarningUtils


@pytest.mark.parametrize('nesting', [None, "templateCatalog"])
class TestContinuumRemoval:

    class TestMethod:
        def test_ok_if_method_is_IrregularSamplingMedian_with_medianKernelWidth_medianEvenReflection(
                self, zflag, nesting):
            parametersDict = DictUtils.make_nested({
                "continuumRemoval": {
                    "method": "IrregularSamplingMedian",
                    "medianKernelWidth": 1,
                    "medianEvenReflection": 1,
                }
            }, nesting)

            accessor = ParametersAccessor(parametersDict)
            ParametersChecker(accessor).custom_check()
            assert not WarningUtils.has_warning(zflag)

        def test_error_if_method_is_IrregularSamplingMedian_without_medianKernelWidth(self, nesting):
            parametersDict = DictUtils.make_nested({
                "continuumRemoval": {
                    "method": "IrregularSamplingMedian",
                    "medianEvenReflection": 1
                }
            }, nesting)
            accessor = ParametersAccessor(parametersDict)

            with pytest.raises(APIException, match="Missing parameter continuumRemoval medianKernelWidth"):
                ParametersChecker(accessor).custom_check()

        def test_error_if_method_is_IrregularSamplingMedian_without_medianEvenReflection(self, nesting):
            parametersDict = DictUtils.make_nested({
                "continuumRemoval": {
                    "method": "IrregularSamplingMedian",
                    "medianKernelWidth": 1
                }
            }, nesting)
            accessor = ParametersAccessor(parametersDict)

            with pytest.raises(APIException, match="Missing parameter continuumRemoval medianEvenReflection"):
                ParametersChecker(accessor).custom_check()

    class TestMedianEvenWidth:
        def test_warning_present_and_method_is_not_IrregularSamplingMedian(self, zflag, nesting):
            parametersDict = DictUtils.make_nested({
                "continuumRemoval": {
                    "method": "sth",
                    "medianKernelWidth": 1,
                }
            }, nesting)
            accessor = ParametersAccessor(parametersDict)
            ParametersChecker(accessor).custom_check()
            assert WarningUtils.has_warning(zflag)

    class TestMedianEvenReflection:
        def test_warning_present_and_method_is_not_IrregularSamplingMedian(self, zflag, nesting):
            parametersDict = DictUtils.make_nested({
                "continuumRemoval": {
                    "method": "sth",
                    "medianEvenReflection": 1
                }
            }, nesting)
            accessor = ParametersAccessor(parametersDict)
            ParametersChecker(accessor).custom_check()
            assert WarningUtils.has_warning(zflag)
