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
from tests.python.utils import (WarningUtils, check_from_parameter_dict,
                                make_parameter_dict_at_linemeas_solve_level,
                                make_parameter_dict_at_object_level,
                                make_parameter_dict_at_redshift_solver_level)


class TestParametersCheckerObject:

    class TestObjectLinemeasDzhalf:

        def _make_param_dict(self, **kwargs):
            return make_parameter_dict_at_linemeas_solve_level(**kwargs)

        def test_error_if_linemeas_stage_and_dzhalf_not_defined(self):
            param_dict = self._make_param_dict(**{})
            del param_dict["galaxy"]["lineMeasDzHalf"]
            with pytest.raises(APIException, match=r"Missing parameter lineMeasDzHalf for object"):
                check_from_parameter_dict(param_dict)

        def test_ok_if_linemeas_stage_and_dzhalf_defined(self, zflag):
            param_dict = self._make_param_dict(**{})
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

        def test_warning_if_not_linemeas_stage_and_dzhalf_defined(self, zflag):
            param_dict = self._make_param_dict(**{})
            param_dict["galaxy"]["stages"] = []
            del param_dict["galaxy"]["lineMeasRedshiftStep"]
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_any_warning()

    class TestObjectLinemeasredshiftStep:

        def _make_param_dict(self, **kwargs):
            return make_parameter_dict_at_linemeas_solve_level(**kwargs)

        def test_error_if_linemeas_stage_and_redshift_step_not_defined(self):
            param_dict = self._make_param_dict(**{})
            del param_dict["galaxy"]["lineMeasRedshiftStep"]
            with pytest.raises(APIException, match=r"Missing parameter lineMeasRedshiftStep for object"):
                check_from_parameter_dict(param_dict)

        def test_ok_if_linemeas_stage_and_redshift_step_defined(self, zflag):
            param_dict = self._make_param_dict(**{})
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

        def test_warning_if_not_linemeas_stage_and_redshift_step_defined(self, zflag):
            param_dict = self._make_param_dict(**{})
            param_dict["galaxy"]["stages"] = []
            del param_dict["galaxy"]["lineMeasDzHalf"]
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_any_warning()

    class TestObjectReliability:

        def _make_parameter_dict(self, **kwargs) -> dict:
            return make_parameter_dict_at_object_level(**kwargs)

        def test_error_if_reliability_enabled_without_corresponding_section(self):
            param_dict = self._make_parameter_dict(**{
                "stages": ["reliabilitySolver"]
            })
            with pytest.raises(APIException, match=r"Missing parameter galaxy reliabilitySolver"):
                check_from_parameter_dict(param_dict)

        def test_OK_if_reliability_enabled_with_corresponding_section(self, zflag):
            param_dict = self._make_parameter_dict(**{
                "stages": ["reliabilitySolver"],
                "reliabilitySolver": {}
            })
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

        def test_warning_if_reliability_disabled_but_section_is_present(self, zflag):
            param_dict = self._make_parameter_dict(**{
                "stages": ["redshiftSolver"],
                "reliabilitySolver": {}
            })
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_warning(WarningCode.UNUSED_PARAMETER)

    class TestTemplateDir:

        def test_error_if_templateFittingSolve_without_template_dir(self):
            param_dict = make_parameter_dict_at_redshift_solver_level(**{
                "method": "templateFittingSolve", "templateFittingSolve": {}
            })
            with pytest.raises(APIException, match=r"Missing parameter galaxy templateDir"):
                check_from_parameter_dict(param_dict)

        def test_error_if_templateCombinationSolve_without_template_dir(self):
            param_dict = make_parameter_dict_at_redshift_solver_level(**{
                "method": "tplCombinationSolve", "tplCombinationSolve": {}
            })
            with pytest.raises(APIException, match=r"Missing parameter galaxy templateDir"):
                check_from_parameter_dict(param_dict)

        def test_error_if_lineModelSolve_and_continuum_tplfit_without_template_dir(self):
            param_dict = make_parameter_dict_at_redshift_solver_level(**{
                "method": "lineModelSolve",
                "lineModelSolve": {
                    "lineModel": {
                        "continuumComponent": "tplFit",  # This is the important line
                        "continuumFit": {},
                        "secondPass": {
                            "continuumFit": "sth"
                        }
                    }
                }
            })
            param_dict["continuumRemoval"] = {}
            with pytest.raises(APIException, match=r"Missing parameter galaxy templateDir"):
                check_from_parameter_dict(param_dict)

        def test_error_if_lineModelSolve_and_continuum_tplfitauto_without_template_dir(self):
            param_dict = make_parameter_dict_at_redshift_solver_level(**{
                "method": "lineModelSolve",
                "lineModelSolve": {
                    "lineModel": {
                        "continuumComponent": "tplFitAuto",  # This is the important line
                        "continuumFit": {},
                        "secondPass": {
                            "continuumFit": "sth"
                        },
                    }
                }
            })
            param_dict["continuumRemoval"] = {}
            with pytest.raises(APIException, match=r"Missing parameter galaxy templateDir"):
                check_from_parameter_dict(param_dict)

        def test_warning_if_template_dir_is_present_but_no_solve_method(self, zflag):
            param_dict = make_parameter_dict_at_object_level(**{
                "stages": [],
                "templateDir": "sth"
            })
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_warning(WarningCode.UNUSED_PARAMETER)

        def test_warning_if_template_dir_is_present_but_linemodelsolve_without_tplfit(self, zflag):
            param_dict = make_parameter_dict_at_redshift_solver_level(**{
                "method": "lineModelSolve",
                "lineModelSolve": {
                    "lineModel": {
                        "continuumComponent": "sth",  # This is the important line
                        "continuumFit": {},
                        "secondPass": {
                            "continuumFit": "sth"
                        },
                    }
                }
            })
            param_dict["continuumRemoval"] = {}
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_warning(WarningCode.UNUSED_PARAMETER)
