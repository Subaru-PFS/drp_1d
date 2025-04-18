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
from pylibamazed.redshift import WarningCode
from tests.python.utils import (
    WarningUtils,
    check_from_parameter_dict,
    make_parameter_dict_at_linemeas_solve_level,
    make_parameter_dict_linemeas_solve_piped_linemodel,
    make_parameter_dict_at_object_level,
)


class TestLineMeasSolve:
    class TestLineMeasSolverSection:
        def test_error_if_stage_but_not_section(self):
            param_dict = make_parameter_dict_at_object_level(**{"stages": ["lineMeasSolver"]})
            param_dict["lsf"] = {}
            with pytest.raises(APIException, match=r"Missing parameter galaxy lineMeasSolver"):
                check_from_parameter_dict(param_dict)

        def test_warning_if_not_stage_but_section(self, zflag):
            param_dict = make_parameter_dict_at_object_level(**{"lineMeasSolver": {}})
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_warning(WarningCode.UNUSED_PARAMETER)

        def test_ok_if_both(self, zflag):
            param_dict = make_parameter_dict_at_object_level(
                **{"stages": ["lineMeasSolver"], "lineMeasSolver": {}}
            )
            param_dict["lsf"] = {}
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_warning(WarningCode.UNUSED_PARAMETER)

    class TestLineRatioType:
        def _make_parameter_dict(self, **kwargs):
            param_dict = make_parameter_dict_at_linemeas_solve_level(**{"lineModel": kwargs})
            return param_dict

        def test_error_if_lineRatioType_is_rules_but_rules_is_absent(self):
            param_dict = self._make_parameter_dict(**{"lineRatioType": "rules"})
            with pytest.raises(APIException, match=r"Missing parameter lineMeasSolve rules"):
                check_from_parameter_dict(param_dict)

        def test_ok_if_lineRatioType_is_rules_and_rules_is_present(self, zflag):
            param_dict = self._make_parameter_dict(**{"lineRatioType": "rules", "rules": "sth"})
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

        def test_warning_if_lineRatioType_is_not_rules_but_rules_is_present(self, zflag):
            param_dict = self._make_parameter_dict(**{"rules": "sth"})
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_any_warning()

    class TestFittingMethod:
        def _make_parameter_dict(self, **kwargs):
            param_dict = make_parameter_dict_at_linemeas_solve_level(**{"lineModel": kwargs})
            return param_dict

        def test_error_if_fittingmethod_is_lbfgsb_but_velocityfit_is_absent(self):
            param_dict = self._make_parameter_dict(**{"fittingMethod": "lbfgsb"})
            with pytest.raises(APIException, match=r"Missing parameter lineMeasSolve velocityFit"):
                check_from_parameter_dict(param_dict)

        def test_ok_if_fittingmethod_is_lbfgsb_and_velocityfit_is_present(self, zflag):
            param_dict = self._make_parameter_dict(**{"fittingMethod": "lbfgsb", "velocityFit": False})
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

        def test_warning_if_fitting_method_is_not_lbfgsb_but_velocityfit_is_present(self, zflag):
            param_dict = self._make_parameter_dict(**{"velocityFit": False})
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_any_warning()

    class TestVelocityFit:
        def _make_parameter_dict(self, **kwargs):
            param_dict = make_parameter_dict_at_linemeas_solve_level(
                **{"lineModel": {"fittingMethod": "lbfgsb", **kwargs}}
            )
            return param_dict

        def _make_parameter_pipe_dict(self, linemodel_params, linemeas_params):
            param_dict = make_parameter_dict_linemeas_solve_piped_linemodel(
                linemodel_level_params=linemodel_params,
                linemeas_level_params={"lineModel": {"fittingMethod": "lbfgsb", **linemeas_params}},
            )
            return param_dict

        def test_error_if_velocityfit_is_true_but_a_velocity_param_is_absent(self):
            param_dict = self._make_parameter_dict(
                **{
                    "velocityFit": True,
                    "emVelocityFitMin": 1,
                }
            )
            with pytest.raises(APIException, match=r"Missing parameter lineMeasSolve"):
                check_from_parameter_dict(param_dict)

        def test_ok_if_velocityfit_is_true_and_all_velocity_params_are_present(self, zflag):
            param_dict = self._make_parameter_dict(
                **{
                    "velocityFit": True,
                    "velocityEmission": 100,
                    "velocityAbsorption": 100,
                    "emVelocityFitMin": 50,
                    "emVelocityFitMax": 300,
                    "absVelocityFitMin": 100,
                    "absVelocityFitMax": 500,
                }
            )
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

        def test_warning_if_velocityfit_is_false_but_some_velocity_params_are_present(self, zflag):
            param_dict = self._make_parameter_dict(
                **{
                    "velocityFit": False,
                    "emVelocityFitMin": 1,
                }
            )
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_any_warning()

        def test_error_if_velocity_not_in_range(self, zflag):
            param_dict = self._make_parameter_dict(
                **{
                    "velocityFit": True,
                    "velocityEmission": 100,
                    "velocityAbsorption": 100,
                    "emVelocityFitMin": 200,
                    "emVelocityFitMax": 300,
                    "absVelocityFitMin": 100,
                    "absVelocityFitMax": 500,
                }
            )
            with pytest.raises(APIException):
                check_from_parameter_dict(param_dict)

        def test_ok_piped_velocity(self, zflag):
            param_dict = self._make_parameter_pipe_dict(
                linemodel_params={
                    "velocityFit": True,
                    "velocityEmission": 100,
                    "velocityAbsorption": 100,
                    "emVelocityFitMin": 50,
                    "emVelocityFitMax": 300,
                    "emVelocityFitStep": 50,
                    "absVelocityFitMin": 100,
                    "absVelocityFitMax": 500,
                    "absVelocityFitStep": 50,
                },
                linemeas_params={
                    "velocityFit": True,
                    "velocityEmission": 100,
                    "velocityAbsorption": 100,
                    "emVelocityFitMin": 50,
                    "emVelocityFitMax": 300,
                    "absVelocityFitMin": 100,
                    "absVelocityFitMax": 500,
                },
            )
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

        def test_error1_piped_velocity(self, zflag):
            param_dict = self._make_parameter_pipe_dict(
                linemodel_params={
                    "velocityFit": False,
                    "velocityEmission": 50,
                    "velocityAbsorption": 100,
                },
                linemeas_params={
                    "velocityFit": True,
                    "velocityEmission": 100,
                    "velocityAbsorption": 100,
                    "emVelocityFitMin": 100,
                    "emVelocityFitMax": 300,
                    "absVelocityFitMin": 100,
                    "absVelocityFitMax": 500,
                },
            )
            with pytest.raises(APIException):
                check_from_parameter_dict(param_dict)

        def test_error2_piped_velocity(self, zflag):
            param_dict = self._make_parameter_pipe_dict(
                linemodel_params={
                    "velocityFit": True,
                    "velocityEmission": 100,
                    "velocityAbsorption": 100,
                    "emVelocityFitMin": 50,
                    "emVelocityFitMax": 300,
                    "emVelocityFitStep": 50,
                    "absVelocityFitMin": 100,
                    "absVelocityFitMax": 500,
                    "absVelocityFitStep": 50,
                },
                linemeas_params={
                    "velocityFit": True,
                    "velocityEmission": 100,
                    "velocityAbsorption": 100,
                    "emVelocityFitMin": 100,
                    "emVelocityFitMax": 300,
                    "absVelocityFitMin": 100,
                    "absVelocityFitMax": 500,
                },
            )
            with pytest.raises(APIException):
                check_from_parameter_dict(param_dict)

    class TestLyaFit:
        def _make_parameter_dict(self, **kwargs) -> dict:
            kwargs["method"] = "lineModelSolve"
            param_dict = make_parameter_dict_at_linemeas_solve_level(**kwargs)
            return param_dict

        def test_error_if_lya_profile_is_asym_but_asym_section_absent(self):
            param_dict = self._make_parameter_dict(**{"lineModel": {"lya": {"profile": "asym"}}})
            with pytest.raises(
                APIException, match=r"Missing parameter lineMeasSolve linemodel lya asymProfile section"
            ):
                check_from_parameter_dict(param_dict)

        def test_warning_if_lya_profile_is_not_asym_but_asym_section_present(self, zflag):
            param_dict = self._make_parameter_dict(
                **{"lineModel": {"lya": {"profile": "igm", "asymProfile": {}}}}
            )
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_any_warning()

        def test_ok_if_lya_profile_is_asym_and_asym_section_present(self, zflag):
            param_dict = self._make_parameter_dict(
                **{"lineModel": {"lya": {"profile": "asym", "asymProfile": {}}}}
            )
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

        def test_ok_if_lya_profile_is_not_asym_and_asym_section_absent(self, zflag):
            param_dict = self._make_parameter_dict(**{"lineModel": {"lya": {"profile": "igm"}}})
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

        class TestEmVelocityFitMin:
            def _make_parameter_dict(self, **kwargs) -> dict:
                kwargs["method"] = "lineModelSolve"
                param_dict = make_parameter_dict_at_linemeas_solve_level(**kwargs)
                return param_dict

            def test_error_if_em_velocity_fit_min_is_negative_with_velocity_driven(self):
                param_dict = self._make_parameter_dict(
                    **{"lineModel": {"lineWidthType": "velocityDriven", "emVelocityFitMin": -0.5}}
                )
                with pytest.raises(
                    APIException,
                    match=r"lineMeasSolve lineModel emVelocityFit min must be > 0 when lineWidthType is velocityDriven",
                ):
                    check_from_parameter_dict(param_dict)

            def test_OK_if_em_velocity_fit_min_is_positive_with_velocity_driven(self):
                param_dict = self._make_parameter_dict(
                    **{"lineModel": {"lineWidthType": "velocityDriven", "emVelocityFitMin": 0.5}}
                )
                check_from_parameter_dict(param_dict)
