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
from tests.python.utils import WarningUtils, check_from_parameter_dict


class TestLSF:
    class TestLSFTypeGaussianConstantWidth:
        def test_raises_an_error_if_GaussianConstantWidth_without_width_defined(self):
            parametersDict = {
                "LSF": {
                    "LSFType": "GaussianConstantWidth"
                }
            }
            accessor = ParametersAccessor(parametersDict)
            with pytest.raises(APIException, match=r"Missing parameter LSF width"):
                ParametersChecker(accessor).custom_check()

        def test_OK_if_GaussianConstantWidth_with_width_defined(self, zflag):
            parametersDict = {
                "LSF": {
                    "LSFType": "GaussianConstantWidth",
                    "width": "1"
                }
            }
            check_from_parameter_dict(parametersDict)
            assert not WarningUtils.has_warning(zflag)

    class TestLSFTypeGaussianConstantResolution:
        def test_raises_an_error_if_GaussianConstantResolution_without_width_defined(self):
            parametersDict = {
                "LSF": {
                    "LSFType": "GaussianConstantResolution"
                }
            }
            with pytest.raises(APIException, match=r"Missing parameter LSF resolution"):
                check_from_parameter_dict(parametersDict)

        def test_OK_if_GaussianConstantResolution_with_width_defined(self, zflag):
            parametersDict = {
                "LSF": {
                    "LSFType": "GaussianConstantResolution",
                    "resolution": "1"
                }
            }
            check_from_parameter_dict(parametersDict)
            assert not WarningUtils.has_warning(zflag)

    class TestLSFTypeGaussianNISPSIM2016:
        def test_raises_an_error_if_GaussianConstantResolution_without_width_defined(self):
            parametersDict = {
                "LSF": {
                    "LSFType": "GaussianNISPSIM2016"
                }
            }
            with pytest.raises(APIException, match=r"Missing parameter LSF sourcesize"):
                check_from_parameter_dict(parametersDict)

        def test_OK_if_GaussianConstantResolution_with_width_defined(self, zflag):
            parametersDict = {
                "LSF": {
                    "LSFType": "GaussianNISPSIM2016",
                    "sourcesize": "1"
                }
            }
            check_from_parameter_dict(parametersDict)
            assert not WarningUtils.has_warning(zflag)

    class TestLSFTypeGaussianNISPSIM201707:
        def test_raises_an_error_if_GaussianConstantResolution_without_width_defined(self):
            parametersDict = {
                "LSF": {
                    "LSFType": "GaussianNISPSIM201707"
                }
            }
            with pytest.raises(APIException, match=r"Missing parameter LSF sourcesize"):
                check_from_parameter_dict(parametersDict)

        def test_OK_if_GaussianConstantResolution_with_width_defined(self, zflag):
            parametersDict = {
                "LSF": {
                    "LSFType": "GaussianNISPSIM201707",
                    "sourcesize": "1"
                }
            }
            check_from_parameter_dict(parametersDict)
            assert not WarningUtils.has_warning(zflag)

    class TestLSFTypeGaussianVariablewidth:
        def test_raises_an_error_if_GaussianConstantResolution_without_width_defined(self):
            parametersDict = {
                "LSF": {
                    "LSFType": "GaussianVariablewidth"
                }
            }
            with pytest.raises(APIException, match=r"Missing parameter LSF GaussianVariablewidthFileName"):
                check_from_parameter_dict(parametersDict)

        def test_OK_if_GaussianConstantResolution_with_width_defined(self, zflag):
            parametersDict = {
                "LSF": {
                    "LSFType": "GaussianVariablewidth",
                    "GaussianVariablewidthFileName": "someFileName"
                }
            }
            check_from_parameter_dict(parametersDict)
            assert not WarningUtils.has_warning(zflag)

    class TestWidth:
        def test_warning_if_width_defined_with_other_LSF_type_than_GaussianConstantWidth(self, zflag):
            parametersDict = {
                "LSF": {
                    "LSFType": "sth",
                    "width": "1"
                }
            }
            check_from_parameter_dict(parametersDict)
            assert WarningUtils.has_warning(zflag)

        def test_OK_if_width_not_defined_with_other_LSF_type_than_GaussianConstantWidth(self, zflag):
            parametersDict = {
                "LSF": {
                    "LSFType": "sth",
                }
            }
            check_from_parameter_dict(parametersDict)
            assert not WarningUtils.has_warning(zflag)

    class TestResolution:
        def test_warning_if_resolution_defined_with_other_LSF_type_than_GaussianConstantResolution(
                self, zflag):
            parametersDict = {
                "LSF": {
                    "LSFType": "sth",
                    "resolution": "1"
                }
            }
            check_from_parameter_dict(parametersDict)
            assert WarningUtils.has_warning(zflag)

        def test_OK_if_resolution_not_defined_with_other_LSF_type_than_GaussianConstantResolution(
                self, zflag):
            parametersDict = {
                "LSF": {
                    "LSFType": "sth",
                }
            }
            check_from_parameter_dict(parametersDict)
            assert not WarningUtils.has_warning(zflag)

    class TestSourceSize:
        def test_warning_if_sourcesize_defined_with_other_LSF_type_than_NISP(self, zflag):
            parametersDict = {
                "LSF": {
                    "LSFType": "sth",
                    "sourcesize": "1"
                }
            }
            check_from_parameter_dict(parametersDict)
            assert WarningUtils.has_warning(zflag)

        def test_OK_if_sourcesize_not_defined_with_other_LSF_type_than_NISP(self, zflag):
            parametersDict = {
                "LSF": {
                    "LSFType": "sth",
                }
            }
            check_from_parameter_dict(parametersDict)
            assert not WarningUtils.has_warning(zflag)

    class TestFileName:
        def test_warning_if_filename_defined_with_other_LSF_type_than_NISP(self, zflag):
            parametersDict = {
                "LSF": {
                    "LSFType": "sth",
                    "GaussianVariablewidthFileName": "someFileName"
                }
            }
            check_from_parameter_dict(parametersDict)
            assert WarningUtils.has_warning(zflag)

        def test_OK_if_sourcesize_not_defined_with_other_LSF_type_than_NISP(self, zflag):
            parametersDict = {
                "LSF": {
                    "LSFType": "sth",
                }
            }
            check_from_parameter_dict(parametersDict)
            assert not WarningUtils.has_warning(zflag)
