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
from tests.python.utils import (
    WarningUtils,
    check_from_parameter_dict,
    make_parameter_dict_at_redshift_solver_level,
)


class TestLineModelSolve:
    class TestMethod:
        def _make_parameter_dict(self, **kwargs):
            kwargs["linemeas_method"] = ""
            kwargs["method"] = "lineModelSolve"
            return make_parameter_dict_at_redshift_solver_level(**kwargs)

        def test_error_if_method_is_lineModelSolve_and_section_is_absent(self):
            param_dict = self._make_parameter_dict(**{})
            with pytest.raises(APIException, match=r"Missing parameter lineModelSolve"):
                check_from_parameter_dict(param_dict)

        def test_OK_if_method_is_lineModelSolve_and_section_is_present(self, zflag):
            param_dict = self._make_parameter_dict(**{"lineModelSolve": {}})
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

        def test_OK_if_lineRatioType_is_rules_and_rules_section_is_present(self, zflag):
            param_dict = self._make_parameter_dict(
                **{"lineModelSolve": {"lineRatioType": "rules", "rules": {}}}
            )

            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

    class TestLineRatioType:
        def _make_parameter_dict(self, **kwargs):
            kwargs["linemeas_method"] = ""
            kwargs["method"] = "lineModelSolve"
            if kwargs.get("lineModelSolve", {}).get("lineModel", {}).get("lineRatioType") in [
                "tplRatio",
                "tplCorr",
            ]:
                kwargs["lineModelSolve"]["lineModel"]["firstPass"] = {"tplRatioIsmFit": False}
            param_dict = make_parameter_dict_at_redshift_solver_level(**kwargs)
            return param_dict

        def test_error_if_lineRatioType_is_rules_but_rules_section_is_absent(self):
            param_dict = self._make_parameter_dict(
                **{
                    "lineModelSolve": {
                        "lineModel": {
                            "lineRatioType": "rules",
                        }
                    }
                }
            )
            with pytest.raises(APIException, match=r"Missing parameter lineModelSolve rules"):
                check_from_parameter_dict(param_dict)

        def test_warning_if_lineRatioType_is_not_rules_and_rules_section_is_present(self, zflag):
            param_dict = self._make_parameter_dict(
                **{"lineModelSolve": {"lineModel": {"lineRatioType": "sth", "rules": {}}}}
            )
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_any_warning()

        def test_OK_if_lineRatioType_is_not_rules_and_rules_section_is_absent(self, zflag):
            param_dict = self._make_parameter_dict(**{"lineModelSolve": {}})
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

        @pytest.mark.parametrize("tpl_ratio", ["tplRatio", "tplCorr"])
        def test_OK_if_lineRatioType_is_tplratio_and_tplratio_params_are_present(self, zflag, tpl_ratio):
            param_dict = self._make_parameter_dict(
                **{
                    "lineModelSolve": {
                        "lineModel": {
                            "lineRatioType": tpl_ratio,
                            "tplRatioCatalog": "",
                            "tplRatioIsmFit": False,
                        }
                    }
                }
            )
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

        @pytest.mark.parametrize("tpl_ratio", ["tplRatio", "tplCorr"])
        def test_error_if_lineRatioType_is_tplratio_and_missing_tplratio_catalog(self, tpl_ratio):
            param_dict = self._make_parameter_dict(
                **{
                    "lineModelSolve": {
                        "lineModel": {
                            "lineRatioType": tpl_ratio,
                            "tplRatioIsmFit": False,
                        }
                    }
                }
            )
            with pytest.raises(APIException, match=r"Missing parameter lineModelSolve tplRatioCatalog"):
                check_from_parameter_dict(param_dict)

        @pytest.mark.parametrize("tpl_ratio", ["tplRatio", "tplCorr"])
        def test_error_if_lineRatioType_is_tplratio_and_missing_tplratio_ismfit(self, tpl_ratio):
            param_dict = self._make_parameter_dict(
                **{
                    "lineModelSolve": {
                        "lineModel": {
                            "lineRatioType": tpl_ratio,
                            "tplRatioCatalog": "",
                        }
                    }
                }
            )
            with pytest.raises(APIException, match=r"Missing parameter lineModelSolve tplRatioIsmFit"):
                check_from_parameter_dict(param_dict)

        @pytest.mark.parametrize("tplparam", ["tplRatioCatalog", "tplRatioIsmFit"])
        def test_warning_if_lineRatioType_is_rules_and_tplratio_catalog_is_present(self, zflag, tplparam):
            param_dict = self._make_parameter_dict(
                **{"lineModelSolve": {"lineModel": {"lineRatioType": "rules", "rules": {}, tplparam: "sth"}}}
            )
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_any_warning()

    class TestContinuumComponent:
        def _make_parameter_dict(self, object_level_params=None, **kwargs):
            kwargs["linemeas_method"] = ""
            kwargs["method"] = "lineModelSolve"
            if kwargs.get("lineModelSolve", {}).get("lineModel", {}).get("secondPass") is not None:
                kwargs["lineModelSolve"]["lineModel"]["skipSecondPass"] = False
            param_dict = make_parameter_dict_at_redshift_solver_level(object_level_params, **kwargs)
            return param_dict

        def test_error_if_continuumcomponent_is_fromspectrum_but_continuumreestimation_absent(self):
            param_dict = self._make_parameter_dict(
                **{
                    "lineModelSolve": {
                        "lineModel": {
                            "continuumComponent": "fromSpectrum",
                        }
                    }
                }
            )
            param_dict["continuumRemoval"] = {}

            with pytest.raises(
                APIException, match=r"Missing parameter object galaxy lineModelSolve continuumReestimation"
            ):
                check_from_parameter_dict(param_dict)

        def test_warning_if_continuumcomponent_is_not_fromspectrum_but_continuumreestimation_present(
            self, zflag
        ):
            param_dict = self._make_parameter_dict(
                **{
                    "lineModelSolve": {
                        "lineModel": {"continuumComponent": "sth", "continuumReestimation": "sth"}
                    }
                }
            )
            param_dict["continuumRemoval"] = {}
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_any_warning()

        def test_OK_if_continuumcomponent_is_fromspectrum_and_mandatory_fields_present(self, zflag):
            param_dict = self._make_parameter_dict(
                **{
                    "lineModelSolve": {
                        "lineModel": {"continuumComponent": "fromSpectrum", "continuumReestimation": "sth"}
                    }
                }
            )
            param_dict["continuumRemoval"] = {}
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

        @pytest.mark.parametrize("continuum_component", ["tplFit", "tplFitAuto"])
        def test_error_if_continuumcomponent_is_tplfit_but_continuumfit_is_absent(self, continuum_component):
            param_dict = self._make_parameter_dict(
                **{
                    "lineModelSolve": {
                        "lineModel": {
                            "continuumComponent": continuum_component,
                            "secondPass": {"continuumFit": ""},
                        }
                    }
                }
            )
            param_dict["continuumRemoval"] = {}
            with pytest.raises(
                APIException, match=r"Missing parameter object galaxy lineModelSolve continuumFit"
            ):
                check_from_parameter_dict(param_dict)

        @pytest.mark.parametrize("continuum_component", ["tplFit", "tplFitAuto"])
        def test_error_if_continuumcomponent_is_tplfit_but_secondpass_continuumfit_is_absent(
            self, continuum_component
        ):
            param_dict = self._make_parameter_dict(
                **{
                    "lineModelSolve": {
                        "lineModel": {
                            "continuumComponent": continuum_component,
                            "continuumFit": {},
                            "secondPass": {},
                        }
                    },
                }
            )
            param_dict["continuumRemoval"] = {}
            param_dict["galaxy"]["templateDir"] = "sth"

            with pytest.raises(
                APIException, match=r"Missing parameter object galaxy lineModelSolve secondpass continuumFit"
            ):
                check_from_parameter_dict(param_dict)

        @pytest.mark.parametrize("continuum_component", ["tplFit", "tplFitAuto"])
        def test_ok_if_continuumcomponent_is_tplfit_and_continuumfit_is_present(
            self, zflag, continuum_component
        ):
            param_dict = self._make_parameter_dict(
                {"templateDir": "sth"},
                **{
                    "lineModelSolve": {
                        "lineModel": {
                            "continuumComponent": continuum_component,
                            "continuumFit": {},
                            "secondPass": {"continuumFit": ""},
                        }
                    }
                },
            )
            if continuum_component == "tplFitAuto":
                param_dict["continuumRemoval"] = {}
            param_dict["galaxy"]["templateDir"] = "sth"
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()
