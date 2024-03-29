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
from tests.python.utils import (WarningUtils, check_from_parameter_dict,
                                make_parameter_dict_at_object_level)


class TestLineMeasSolve:

    class TestLineRatioType:

        def _make_parameter_dict(self, **kwargs):
            param_dict = make_parameter_dict_at_object_level(**{
                "LineMeasSolve": {"linemodel": kwargs}
            })
            return param_dict

        def test_error_if_lineRatioType_is_rules_but_rules_is_absent(self):
            param_dict = self._make_parameter_dict(**{"lineRatioType": "rules"})
            with pytest.raises(APIException, match=r"Missing parameter LineMeasSolve rules"):
                check_from_parameter_dict(param_dict)

        def test_ok_if_lineRatioType_is_rules_and_rules_is_present(self, zflag):
            param_dict = self._make_parameter_dict(**{"lineRatioType": "rules", "rules": "sth"})
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning(zflag)

        def test_warning_if_lineRatioType_is_not_rules_but_rules_is_present(self, zflag):
            param_dict = self._make_parameter_dict(**{"rules": "sth"})
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_any_warning(zflag)

    class TestFittingMethod:
        def _make_parameter_dict(self, **kwargs):
            param_dict = make_parameter_dict_at_object_level(**{
                "LineMeasSolve": {"linemodel": kwargs}
            })
            return param_dict

        def test_error_if_fittingmethod_is_lbfgsb_but_velocityfit_is_absent(self):
            param_dict = self._make_parameter_dict(**{"fittingmethod": "lbfgsb"})
            with pytest.raises(APIException, match=r"Missing parameter LineMeasSolve velocityfit"):
                check_from_parameter_dict(param_dict)

        def test_ok_if_fittingmethod_is_lbfgsb_and_velocityfit_is_present(self, zflag):
            param_dict = self._make_parameter_dict(**{"fittingmethod": "lbfgsb", "velocityfit": False})
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning(zflag)

        def test_warning_if_fitting_method_is_not_lbfgsb_but_velocityfit_is_present(self, zflag):
            param_dict = self._make_parameter_dict(**{"velocityfit": False})
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_any_warning(zflag)

    class TestVelocityFit:
        def _make_parameter_dict(self, **kwargs):
            param_dict = make_parameter_dict_at_object_level(**{
                "LineMeasSolve": {"linemodel": {"fittingmethod": "lbfgsb", **kwargs}}
            })
            return param_dict

        def test_error_if_velocityfit_is_true_but_a_velocity_param_is_absent(self):
            param_dict = self._make_parameter_dict(**{
                "velocityfit": True,
                "emvelocityfitmin": 1,
            })
            with pytest.raises(APIException, match=r"Missing parameter LineMeasSolve"):
                check_from_parameter_dict(param_dict)

        def test_ok_if_velocityfit_is_true_and_all_velocity_params_are_present(self, zflag):
            param_dict = self._make_parameter_dict(**{
                "velocityfit": True,
                "emvelocityfitmin": 1,
                "emvelocityfitmax": 1,
                "absvelocityfitmin": 1,
                "absvelocityfitmax": 1,
            })
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning(zflag)

        def test_warning_if_velocityfit_is_false_but_some_velocity_params_are_present(self, zflag):
            param_dict = self._make_parameter_dict(**{
                "velocityfit": False,
                "emvelocityfitmin": 1,
            })
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_any_warning(zflag)
