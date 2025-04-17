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
from pylibamazed.CustomParametersChecker import CustomParametersChecker
from pylibamazed.redshift import WarningCode
from tests.python.utils import (
    WarningUtils,
    check_from_parameter_dict,
    make_parameter_dict_at_linemodelsolve_level,
)


class TestLSFUtils:
    def _make_parameter_dict(self, **kwargs):
        param_dict = make_parameter_dict_at_linemodelsolve_level()
        param_dict["lsf"] = kwargs
        return param_dict


class TestLSF:
    class TestLSFSectionPresence(TestLSFUtils):
        def test_raises_an_error_if_linemodelsolve_without_lsf(self):
            param_dict = make_parameter_dict_at_linemodelsolve_level()
            del param_dict["lsf"]
            with pytest.raises(APIException, match="Missing parameter lsf"):
                check_from_parameter_dict(param_dict)

        def test_raises_a_warning_if_lsf_without_linemodelsolve(self, zflag):
            param_dict = {"lsf": {}}
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_warning(WarningCode.UNUSED_PARAMETER)

        def test_ok_if_linemodelsolve_and_lsf(self, zflag):
            param_dict = make_parameter_dict_at_linemodelsolve_level()
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_warning(WarningCode.UNUSED_PARAMETER)

        def test_ok_if_only_one_object_with_linemodelsolve_and_lsf(self, zflag):
            param_dict = make_parameter_dict_at_linemodelsolve_level()
            param_dict["spectrumModels"].append("star")
            param_dict["star"] = {"redshiftSolve": {"method": "sth"}}
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_warning(WarningCode.UNUSED_PARAMETER)

    class TestLSFTypeGaussianConstantWidth(TestLSFUtils):
        def test_raises_an_error_if_GaussianConstantWidth_without_width_defined(self):
            param_dict = self._make_parameter_dict(**{"lsfType": "gaussianConstantWidth"})
            with pytest.raises(APIException, match=r"Missing parameter lsf width"):
                CustomParametersChecker(param_dict).check()

        def test_OK_if_GaussianConstantWidth_with_width_defined(self, zflag):
            param_dict = self._make_parameter_dict(**{"lsfType": "gaussianConstantWidth", "width": "1"})
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

    class TestLSFTypeGaussianConstantResolution(TestLSFUtils):
        def test_raises_an_error_if_GaussianConstantResolution_without_width_defined(self):
            param_dict = self._make_parameter_dict(**{"lsfType": "gaussianConstantResolution"})
            with pytest.raises(APIException, match=r"Missing parameter lsf resolution"):
                check_from_parameter_dict(param_dict)

        def test_OK_if_GaussianConstantResolution_with_width_defined(self, zflag):
            param_dict = self._make_parameter_dict(
                **{"lsfType": "gaussianConstantResolution", "resolution": "1"}
            )
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

    class TestLSFTypeGaussianNISPSIM201707(TestLSFUtils):
        def test_raises_an_error_if_GaussianConstantResolution_without_width_defined(self):
            param_dict = self._make_parameter_dict(**{"lsfType": "GaussianNISPSIM201707"})
            with pytest.raises(APIException, match=r"Missing parameter lsf sourceSize"):
                check_from_parameter_dict(param_dict)

        def test_OK_if_GaussianConstantResolution_with_width_defined(self, zflag):
            param_dict = self._make_parameter_dict(**{"lsfType": "GaussianNISPSIM201707", "sourceSize": "1"})
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

    class TestLSFTypeGaussianVariablewidth(TestLSFUtils):
        def test_raises_an_error_if_GaussianConstantResolution_without_width_defined(self):
            param_dict = self._make_parameter_dict(**{"lsfType": "gaussianVariableWidth"})
            with pytest.raises(APIException, match=r"Missing parameter lsf gaussianVariableWidthFileName"):
                check_from_parameter_dict(param_dict)

        def test_OK_if_GaussianConstantResolution_with_width_defined(self, zflag):
            param_dict = self._make_parameter_dict(
                **{"lsfType": "gaussianVariableWidth", "gaussianVariableWidthFileName": "someFileName"}
            )
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

    class TestWidth(TestLSFUtils):
        def test_warning_if_width_defined_with_other_LSF_type_than_GaussianConstantWidth(self, zflag):
            param_dict = self._make_parameter_dict(**{"lsfType": "sth", "width": "1"})
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_any_warning()

        def test_OK_if_width_not_defined_with_other_LSF_type_than_GaussianConstantWidth(self, zflag):
            param_dict = self._make_parameter_dict(**{"lsfType": "sth"})
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

    class TestResolution(TestLSFUtils):
        def test_warning_if_resolution_defined_with_other_LSF_type_than_GaussianConstantResolution(
            self, zflag
        ):
            param_dict = self._make_parameter_dict(**{"lsfType": "sth", "resolution": "1"})
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_any_warning()

        def test_OK_if_resolution_not_defined_with_other_LSF_type_than_GaussianConstantResolution(
            self, zflag
        ):
            param_dict = self._make_parameter_dict(
                **{
                    "lsfType": "sth",
                }
            )
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

    class TestSourceSize(TestLSFUtils):
        def test_warning_if_sourcesize_defined_with_other_LSF_type_than_NISP(self, zflag):
            param_dict = self._make_parameter_dict(**{"lsfType": "sth", "sourceSize": "1"})
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_any_warning()

        def test_OK_if_sourcesize_not_defined_with_other_LSF_type_than_NISP(self, zflag):
            param_dict = self._make_parameter_dict(**{"lsfType": "sth"})
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

    class TestFileName(TestLSFUtils):
        def test_warning_if_filename_defined_with_other_LSF_type_than_NISP(self, zflag):
            param_dict = self._make_parameter_dict(
                **{"lsfType": "sth", "gaussianVariableWidthFileName": "someFileName"}
            )
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_any_warning()

        def test_OK_if_sourcesize_not_defined_with_other_LSF_type_than_NISP(self, zflag):
            param_dict = self._make_parameter_dict(**{"lsfType": "sth"})
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()
