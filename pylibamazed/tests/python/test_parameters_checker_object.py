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


class TestParametersCheckerObject:

    class TestObjectLinemeasDzhalf:

        def _make_param_dict(self, **kwargs):
            if kwargs.get("linemeas_method") not in ["", None]:
                kwargs["lineMeasRedshiftStep"] = 1
            return make_parameter_dict_at_object_level(**kwargs)

        def test_error_if_method_not_null_and_dzhalf_not_defined(self):
            param_dict = self._make_param_dict(**{
                "linemeas_method": "sth"
            })
            with pytest.raises(APIException, match=r"Missing parameter lineMeasDzHalf for object"):
                check_from_parameter_dict(param_dict)

        def test_ok_if_method_null_and_dzhalf_not_defined(self, zflag):
            param_dict = self._make_param_dict(**{
                "linemeas_method": ""
            })
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning(zflag)

        def test_ok_if_method_not_null_and_dzhalf_defined(self, zflag):
            param_dict = self._make_param_dict(**{
                "linemeas_method": "sth",
                "lineMeasDzHalf": 1
            })
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning(zflag)

        def test_warning_if_method_null_and_dzhalf_defined(self, zflag):
            param_dict = self._make_param_dict(**{
                "linemeas_method": "",
                "lineMeasDzHalf": 1
            })
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_any_warning(zflag)

    class TestObjectLinemeasredshiftStep:

        def _make_param_dict(self, **kwargs):
            if kwargs.get("linemeas_method") not in ["", None]:
                kwargs["lineMeasDzHalf"] = 1
            return make_parameter_dict_at_object_level(**kwargs)

        def test_error_if_method_not_null_and_redshift_not_defined(self):
            param_dict = self._make_param_dict(**{
                "linemeas_method": "sth"
            })
            with pytest.raises(APIException, match=r"Missing parameter lineameas_redshiftstep for object"):
                check_from_parameter_dict(param_dict)

        def test_ok_if_method_null_and_redshift_not_defined(self, zflag):
            param_dict = self._make_param_dict(**{
                "linemeas_method": ""
            })
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning(zflag)

        def test_ok_if_method_not_null_and_redshift_defined(self, zflag):
            param_dict = self._make_param_dict(**{
                "linemeas_method": "sth",
                "lineMeasRedshiftStep": 1
            })
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning(zflag)

        def test_warning_if_method_null_and_redshift_defined(self, zflag):
            param_dict = self._make_param_dict(**{
                "linemeas_method": "",
                "lineMeasRedshiftStep": 1
            })
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_any_warning(zflag)

    class TestObjectReliability:

        def _make_parameter_dict(self, **kwargs) -> dict:
            kwargs["linemeas_method"] = ""
            return make_parameter_dict_at_object_level(**kwargs)

        def test_error_if_reliability_enabled_without_reliability_model(self):
            param_dict = self._make_parameter_dict(**{
                "enable_reliability": True
            })
            with pytest.raises(APIException, match=r"Missing parameter reliabilityModel"):
                check_from_parameter_dict(param_dict)

        def test_OK_if_reliability_enabled_with_reliability_model(self, zflag):

            param_dict = self._make_parameter_dict(**{
                "enable_reliability": True,
                "reliabilityModel": "sth"
            })
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning(zflag)

        def test_warning_if_reliability_model_without_enable_reliability(self, zflag):
            param_dict = self._make_parameter_dict(**{
                "enable_reliability": False,
                "reliabilityModel": "sth"
            })
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_any_warning(zflag)

        def test_OK_if_reliability_disabled_without_reliability_model(self, zflag):
            param_dict = self._make_parameter_dict(**{
                "enable_reliability": False,
            })
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning(zflag)

    class TestTemplateDir:

        def test_error_if_templateFittingSolve_without_template_dir(self):
            param_dict = make_parameter_dict_at_object_level(**{
                "method": "templateFittingSolver", "templateFittingSolver": {}
            })
            with pytest.raises(APIException, match=r"Missing parameter galaxy templateDir"):
                check_from_parameter_dict(param_dict)

        def test_error_if_templateCombinationSolve_without_template_dir(self):
            param_dict = make_parameter_dict_at_object_level(**{
                "method": "tplCombinationSolver", "tplCombinationSolver": {}
            })
            with pytest.raises(APIException, match=r"Missing parameter galaxy templateDir"):
                check_from_parameter_dict(param_dict)

        def test_error_if_lineModelSolve_and_continuum_tplfit_without_template_dir(self):
            param_dict = make_parameter_dict_at_object_level(**{
                "method": "lineModelSolver",
                "lineModelSolver": {
                    "lineModel": {
                        "continuumComponent": "tplFit",  # This is the important line
                        "continuumFit": {},
                        "secondPass": {
                            "continuumFit": "sth"
                        }
                    }
                }
            })
            with pytest.raises(APIException, match=r"Missing parameter galaxy templateDir"):
                check_from_parameter_dict(param_dict)

        def test_error_if_lineModelSolve_and_continuum_tplfitauto_without_template_dir(self):
            param_dict = make_parameter_dict_at_object_level(**{
                "method": "lineModelSolver",
                "lineModelSolver": {
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
                "method": "",
                "templateDir": "sth"
            })
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_warning(zflag, "UNUSED_PARAMETER")

        def test_warning_if_template_dir_is_present_but_linemodelsolve_without_tplfit(self, zflag):
            param_dict = make_parameter_dict_at_object_level(**{
                "method": "lineModelSolver",
                "lineModelSolver": {
                    "lineModel": {
                        "continuumComponent": "sth",  # This is the important line
                        "continuumFit": {},
                        "secondPass": {
                            "continuumFit": "sth"
                        },
                    }
                }
            })
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_warning(zflag, "UNUSED_PARAMETER")
