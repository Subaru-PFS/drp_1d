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
                kwargs["linemeas_redshiftstep"] = 1
            return make_parameter_dict_at_object_level(**kwargs)

        def test_error_if_method_not_null_and_dzhalf_not_defined(self):
            param_dict = self._make_param_dict(**{
                "linemeas_method": "sth"
            })
            with pytest.raises(APIException, match=r"Missing parameter lineameas_dzhalf for object"):
                check_from_parameter_dict(param_dict)

        def test_ok_if_method_null_and_dzhalf_not_defined(self, zflag):
            param_dict = self._make_param_dict(**{
                "linemeas_method": ""
            })
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_warning(zflag)

        def test_ok_if_method_not_null_and_dzhalf_defined(self, zflag):
            param_dict = self._make_param_dict(**{
                "linemeas_method": "sth",
                "linemeas_dzhalf": 1
            })
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_warning(zflag)

        def test_warning_if_method_null_and_dzhalf_defined(self, zflag):
            param_dict = self._make_param_dict(**{
                "linemeas_method": "",
                "linemeas_dzhalf": 1
            })
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_warning(zflag)

    class TestObjectLinemeasredshiftStep:

        def _make_param_dict(self, **kwargs):
            if kwargs.get("linemeas_method") not in ["", None]:
                print("in if")
                kwargs["linemeas_dzhalf"] = 1
            else:
                print("in else")
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
            assert not WarningUtils.has_warning(zflag)

        def test_ok_if_method_not_null_and_redshift_defined(self, zflag):
            param_dict = self._make_param_dict(**{
                "linemeas_method": "sth",
                "linemeas_redshiftstep": 1
            })
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_warning(zflag)

        def test_warning_if_method_null_and_redshift_defined(self, zflag):
            param_dict = self._make_param_dict(**{
                "linemeas_method": "",
                "linemeas_redshiftstep": 1
            })
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_warning(zflag)

    class TestObjectReliability:

        def _make_parameter_dict(self, **kwargs) -> dict:
            kwargs["linemeas_method"] = ""
            return make_parameter_dict_at_object_level(**kwargs)

        def test_error_if_reliability_enabled_without_reliability_model(self):
            param_dict = self._make_parameter_dict(**{
                "enable_reliability": True
            })
            with pytest.raises(APIException, match=r"Missing parameter reliability_model"):
                check_from_parameter_dict(param_dict)

        def test_OK_if_reliability_enabled_with_reliability_model(self, zflag):

            param_dict = self._make_parameter_dict(**{
                "enable_reliability": True,
                "reliability_model": "sth"
            })
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_warning(zflag)

        def test_warning_if_reliability_model_without_enable_reliability(self, zflag):
            param_dict = self._make_parameter_dict(**{
                "enable_reliability": False,
                "reliability_model": "sth"
            })
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_warning(zflag)

        def test_OK_if_reliability_disabled_without_reliability_model(self, zflag):
            param_dict = self._make_parameter_dict(**{
                "enable_reliability": False,
            })
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_warning(zflag)

    class TestTemplateFittingSolve:

        def _make_parameter_dict(self, **kwargs) -> dict:
            kwargs["linemeas_method"] = kwargs.get("linemeas_method", "")
            kwargs["method"] = kwargs.get("method", "TemplateFittingSolve")
            param_dict = make_parameter_dict_at_object_level(**kwargs)
            if kwargs.get("TemplateFittingSolve", {}).get("enablephotometry"):
                param_dict["photometryTransmissionDir"] = "sth"
                param_dict["photometryBand"] = []
            return param_dict

        def test_error_if_method_is_templateFittingSolve_and_section_is_absent(self):
            param_dict = self._make_parameter_dict(**{})
            with pytest.raises(APIException, match=r"Missing parameter TemplateFittingSolve"):
                check_from_parameter_dict(param_dict)

        def test_OK_if_method_is_templateFittingSolve_and_section_is_present(self, zflag):
            param_dict = self._make_parameter_dict(**{
                "TemplateFittingSolve": {}
            })
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_warning(zflag)

        def test_warning_if_method_is_not_templateFittingSolve_but_section_is_present(self, zflag):
            param_dict = self._make_parameter_dict(**{
                "method": "sth",
                "TemplateFittingSolve": {}
            })
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_warning(zflag)

        def test_ok_if_method_is_not_templateFittingSolve_and_section_is_absent(self, zflag):
            param_dict = self._make_parameter_dict(**{
                "method": "sth",
            })
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_warning(zflag)

        def test_error_if_ismfit_enabled_and_ebmv_section_is_not_present(self):
            param_dict = self._make_parameter_dict(**{
                "TemplateFittingSolve": {"ismfit": True}
            })
            with pytest.raises(APIException, match=r"Missing parameter ebmv"):
                check_from_parameter_dict(param_dict)

        def test_ok_if_ismfit_enabled_and_ebmv_is_present(self, zflag):
            param_dict = self._make_parameter_dict(**{
                "TemplateFittingSolve": {"ismfit": True}
            })
            param_dict["ebmv"] = {}
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_warning(zflag)

        def test_error_if_photometry_is_enabled_but_photometry_weight_is_absent(self):
            param_dict = self._make_parameter_dict(**{
                "TemplateFittingSolve": {"enablephotometry": True}
            })
            with pytest.raises(
                APIException,
                match=r"Missing parameter object galaxy TemplateFittingSolve photometry weight"
            ):
                check_from_parameter_dict(param_dict)

        def test_ok_if_photometry_is_enabled_and_photometry_weight_too(self, zflag):
            param_dict = self._make_parameter_dict(**{
                "TemplateFittingSolve": {
                    "enablephotometry": True,
                    "photometry": {"weight": 1}
                }
            })
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_warning(zflag)

        def test_warning_if_photometry_is_disabled_but_photometry_weight_is_present(self, zflag):
            param_dict = self._make_parameter_dict(**{
                "TemplateFittingSolve": {
                    "enablephotometry": False,
                    "photometry": {"weight": 1}
                }
            })
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_warning(zflag)

        def test_ok_if_photometry_is_disabled_and_photometry_weight_is_absent(self, zflag):
            param_dict = self._make_parameter_dict(**{
                "TemplateFittingSolve": {
                    "enablephotometry": False,
                    "photometry": {}
                }
            })
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_warning(zflag)

    class TestTemplateContinuationSolve:

        def _make_parameter_dict(self, **kwargs) -> dict:
            kwargs["linemeas_method"] = kwargs.get("linemeas_method", "")
            kwargs["method"] = kwargs.get("method", "TplcombinationSolve")
            return make_parameter_dict_at_object_level(**kwargs)

        def test_error_if_method_is_TplcombinationSolve_and_section_is_absent(self):
            param_dict = self._make_parameter_dict(**{})
            with pytest.raises(APIException, match=r"Missing parameter TplcombinationSolve"):
                check_from_parameter_dict(param_dict)

        def test_OK_if_method_is_templateCombinationSolve_and_section_is_present(self, zflag):
            param_dict = self._make_parameter_dict(**{
                "TplcombinationSolve": {}
            })
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_warning(zflag)

        def test_warning_if_method_is_not_templateCombinationSolve_but_section_is_present(self, zflag):
            param_dict = self._make_parameter_dict(**{
                "method": "sth",
                "TplcombinationSolve": {}
            })
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_warning(zflag)

        def test_ok_if_method_is_not_templateCombinationSolve_and_section_is_absent(self, zflag):
            param_dict = self._make_parameter_dict(**{
                "method": "sth",
            })
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_warning(zflag)

        def test_error_if_ismfit_enabled_and_ebmv_section_is_not_present(self):
            param_dict = self._make_parameter_dict(**{
                "TplcombinationSolve": {"ismfit": True}
            })
            with pytest.raises(APIException, match=r"Missing parameter ebmv"):
                check_from_parameter_dict(param_dict)

        def test_ok_if_ismfit_enabled_and_ebmv_is_present(self, zflag):
            param_dict = self._make_parameter_dict(**{
                "TplcombinationSolve": {"ismfit": True}
            })
            param_dict["ebmv"] = {}
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_warning(zflag)

    class TestContinuumFit:

        def _make_parameter_dict(self, **kwargs) -> dict:
            kwargs["method"] = kwargs.get("method", "LineModelSolve")
            param_dict = make_parameter_dict_at_object_level(**kwargs)
            return param_dict

        def test_error_if_ismfit_enabled_and_ebmv_section_is_not_present(self):
            param_dict = self._make_parameter_dict(**{
                "LineModelSolve": {"linemodel": {"continuumfit": {"ismfit": True}}}
            })
            print("param dict", param_dict)
            with pytest.raises(APIException, match=r"Missing parameter ebmv"):
                check_from_parameter_dict(param_dict)

        def test_ok_if_ismfit_enabled_and_ebmv_is_present(self, zflag):
            param_dict = self._make_parameter_dict(**{
                "LineModelSolve": {"ismfit": True}
            })
            param_dict["ebmv"] = {}
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_warning(zflag)

    class TestFirstPass:
        def _make_parameter_dict(self, **kwargs) -> dict:
            kwargs["method"] = "LineModelSolve"
            if kwargs["LineModelSolve"]["linemodel"]["lineRatioType"] in ["tplratio", "tplcorr"]:
                kwargs["LineModelSolve"]["linemodel"]["tplratio_catalog"] = "sth"
                kwargs["LineModelSolve"]["linemodel"]["tplratio_ismfit"] = False

            param_dict = make_parameter_dict_at_object_level(**kwargs)
            return param_dict

        def test_error_tplratio_ismfit_not_defined_but_lineratiotype_is_tpl(self):
            param_dict = self._make_parameter_dict(**{
                "LineModelSolve": {"linemodel": {
                    "lineRatioType": "tplratio",
                }}
            })
            with pytest.raises(
                APIException,
                match=r"Missing parameter object galaxy LineModelSolve firstpass tplratio_ismfit"
            ):
                check_from_parameter_dict(param_dict)

        def test_ok_tplratio_ismfit_defined_and_lineratiotype_is_tpl(self, zflag):
            param_dict = self._make_parameter_dict(**{
                "LineModelSolve": {"linemodel": {
                    "lineRatioType": "tplratio",
                    "firstpass": {"tplratio_ismfit": True}
                }}
            })
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_warning(zflag)

        def test_warning_tplratio_imsfit_defined_but_lineratiotype_is_not_tpl(self, zflag):
            param_dict = self._make_parameter_dict(**{
                "LineModelSolve": {"linemodel": {
                    "lineRatioType": "sth",
                    "firstpass": {"tplratio_ismfit": True}
                }}
            })
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_warning(zflag)

    class TestSkipSecondPass:
        def _make_parameter_dict(self, **kwargs) -> dict:
            kwargs["method"] = "LineModelSolve"
            param_dict = make_parameter_dict_at_object_level(**kwargs)
            return param_dict

        def test_error_skipsecondpass_false_but_secondpass_absent(self):
            param_dict = self._make_parameter_dict(**{
                "LineModelSolve": {"linemodel": {
                    "skipsecondpass": False,
                }}
            })
            with pytest.raises(
                APIException,
                match=r"Missing parameter object galaxy LineModelSolve secondpass"
            ):
                check_from_parameter_dict(param_dict)

        def test_ok_skipsecondpass_false_and_secondpass_present(self, zflag):
            param_dict = self._make_parameter_dict(**{
                "LineModelSolve": {"linemodel": {
                    "skipsecondpass": False,
                    "secondpass": {},
                }}
            })
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_warning(zflag)

        def test_warning_skipsecondpass_true_but_secondpass_present(self, zflag):
            param_dict = self._make_parameter_dict(**{
                "LineModelSolve": {"linemodel": {
                    "skipsecondpass": True,
                    "secondpass": {},
                }}
            })
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_warning(zflag)

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
                assert not WarningUtils.has_warning(zflag)

            def test_warning_if_lineRatioType_is_not_rules_but_rules_is_present(self, zflag):
                param_dict = self._make_parameter_dict(**{"rules": "sth"})
                check_from_parameter_dict(param_dict)
                assert WarningUtils.has_warning(zflag)

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
                assert not WarningUtils.has_warning(zflag)

            def test_warning_if_fitting_method_is_not_lbfgsb_but_velocityfit_is_present(self, zflag):
                param_dict = self._make_parameter_dict(**{"velocityfit": False})
                check_from_parameter_dict(param_dict)
                assert WarningUtils.has_warning(zflag)

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
                assert not WarningUtils.has_warning(zflag)

            def test_warning_if_velocityfit_is_false_but_some_velocity_params_are_present(self, zflag):
                param_dict = self._make_parameter_dict(**{
                    "velocityfit": False,
                    "emvelocityfitmin": 1,
                })
                check_from_parameter_dict(param_dict)
                assert WarningUtils.has_warning(zflag)
