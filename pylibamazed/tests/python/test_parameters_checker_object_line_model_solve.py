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
    make_parameter_dict_at_redshift_solver_level,
)


class TestLineModelSolve:
    class TestImproveBalmerFit:
        def _make_parameter_dict(self, **kwargs) -> dict:
            kwargs = {"method": "lineModelSolve", "lineModelSolve": {"lineModel": kwargs}}
            param_dict = make_parameter_dict_at_redshift_solver_level(**kwargs)
            return param_dict

        def test_warning_if_improveBalmerFit_True_but_lineRatioType_is_not_rules(self, zflag):
            param_dict = self._make_parameter_dict(**{"improveBalmerFit": True, "lineRatioType": "sth"})
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_any_warning()

        def test_OK_if_improveBalmerFit_True_and_lineRatioType_is_rules(self, zflag):
            param_dict = self._make_parameter_dict(
                **{"improveBalmerFit": True, "lineRatioType": "rules", "rules": ""}
            )
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

    class TestContinuumFit:
        def _make_parameter_dict(self, **kwargs) -> dict:
            kwargs["method"] = kwargs.get("method", "lineModelSolve")
            param_dict = make_parameter_dict_at_redshift_solver_level(**kwargs)
            return param_dict

        def test_error_if_ismfit_enabled_and_ebmv_section_is_not_present(self):
            param_dict = self._make_parameter_dict(
                **{"lineModelSolve": {"lineModel": {"continuumFit": {"ismFit": True}}}}
            )
            with pytest.raises(APIException, match=r"Missing parameter ebmv"):
                check_from_parameter_dict(param_dict)

        def test_ok_if_ismfit_enabled_and_ebmv_is_present(self, zflag):
            param_dict = self._make_parameter_dict(**{"lineModelSolve": {"ismFit": True}})
            param_dict["ebmv"] = {}
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

        def test_error_if_fftprocessing_enabled_but_photometry_enabled(self):
            param_dict = self._make_parameter_dict(
                **{
                    "lineModelSolve": {
                        "lineModel": {"continuumFit": {"fftProcessing": True}, "enablePhotometry": True}
                    }
                }
            )
            param_dict["photometryTransmissionDir"] = "sth"
            param_dict["photometryBand"] = "sth"
            with pytest.raises(APIException, match=r"cannot activate both fft and photometry"):
                check_from_parameter_dict(param_dict)

    class TestFirstPass:
        def _make_parameter_dict(self, **kwargs) -> dict:
            kwargs["method"] = "lineModelSolve"
            if kwargs["lineModelSolve"]["lineModel"]["lineRatioType"] in ["tplRatio", "tplCorr"]:
                kwargs["lineModelSolve"]["lineModel"]["tplRatioCatalog"] = "sth"
                kwargs["lineModelSolve"]["lineModel"]["tplRatioIsmFit"] = False
            param_dict = make_parameter_dict_at_redshift_solver_level(**kwargs)
            return param_dict

        def test_error_tplratio_ismfit_not_defined_but_lineratiotype_is_tpl(self):
            param_dict = self._make_parameter_dict(
                **{
                    "lineModelSolve": {
                        "lineModel": {
                            "lineRatioType": "tplRatio",
                        }
                    }
                }
            )
            with pytest.raises(
                APIException, match=r"Missing parameter object galaxy lineModelSolve firstpass tplRatioIsmFit"
            ):
                check_from_parameter_dict(param_dict)

        def test_ok_tplratio_ismfit_defined_and_lineratiotype_is_tpl(self, zflag):
            param_dict = self._make_parameter_dict(
                **{
                    "lineModelSolve": {
                        "lineModel": {"lineRatioType": "tplRatio", "firstPass": {"tplRatioIsmFit": True}}
                    }
                }
            )
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

        def test_warning_tplratio_imsfit_defined_but_lineratiotype_is_not_tpl(self, zflag):
            param_dict = self._make_parameter_dict(
                **{
                    "lineModelSolve": {
                        "lineModel": {"lineRatioType": "sth", "firstPass": {"tplRatioIsmFit": True}}
                    }
                }
            )
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_any_warning()

        def test_ok_extremacount(self):
            param_dict = self._make_parameter_dict(
                **{
                    "lineModelSolve": {
                        "lineModel": {
                            "extremaCount": 5,
                            "lineRatioType": "sth",
                            "firstPass": {"extremaCount": 5},
                        }
                    }
                }
            )
            assert check_from_parameter_dict(param_dict) is None

        def test_error_extremacount(self):
            param_dict = self._make_parameter_dict(
                **{
                    "lineModelSolve": {
                        "lineModel": {
                            "extremaCount": 5,
                            "lineRatioType": "sth",
                            "firstPass": {"extremaCount": 3},
                        }
                    }
                }
            )
            with pytest.raises(
                APIException,
                match=r"INVALID_PARAMETER_FILE: linemodel.firstpass.extremaCount is lower than "
                "linemodel.extremaCount for object galaxy",
            ):
                check_from_parameter_dict(param_dict)

    class TestSkipSecondPass:
        def _make_parameter_dict(self, **kwargs) -> dict:
            kwargs["method"] = "lineModelSolve"
            param_dict = make_parameter_dict_at_redshift_solver_level(**kwargs)
            return param_dict

        def test_error_skipsecondpass_false_but_secondpass_absent(self):
            param_dict = self._make_parameter_dict(
                **{
                    "lineModelSolve": {
                        "lineModel": {
                            "skipSecondPass": False,
                        }
                    }
                }
            )
            with pytest.raises(
                APIException, match=r"Missing parameter object galaxy lineModelSolve secondPass"
            ):
                check_from_parameter_dict(param_dict)

        def test_ok_skipsecondpass_false_and_secondpass_present(self, zflag):
            param_dict = self._make_parameter_dict(
                **{
                    "lineModelSolve": {
                        "lineModel": {
                            "skipSecondPass": False,
                            "secondPass": {},
                        }
                    }
                }
            )
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

        def test_warning_skipsecondpass_true_but_secondpass_present(self, zflag):
            param_dict = self._make_parameter_dict(
                **{
                    "lineModelSolve": {
                        "lineModel": {
                            "skipSecondPass": True,
                            "secondPass": {},
                        }
                    }
                }
            )
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_any_warning()

    class TestLyaFit:
        def _make_parameter_dict(self, **kwargs) -> dict:
            kwargs["method"] = "lineModelSolve"
            param_dict = make_parameter_dict_at_redshift_solver_level(**kwargs)
            return param_dict

        def test_error_if_lya_profile_is_asym_but_asym_section_absent(self):
            param_dict = self._make_parameter_dict(
                **{"lineModelSolve": {"lineModel": {"lya": {"profile": "asym"}}}}
            )
            with pytest.raises(
                APIException, match=r"Missing parameter lineModelSolve linemodel lya asymProfile section"
            ):
                check_from_parameter_dict(param_dict)

        def test_warning_if_lya_profile_is_not_asym_but_asym_section_present(self, zflag):
            param_dict = self._make_parameter_dict(
                **{"lineModelSolve": {"lineModel": {"lya": {"profile": "igm", "asymProfile": {}}}}}
            )
            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_any_warning()

        def test_ok_if_lya_profile_is_asym_and_asym_section_present(self, zflag):
            param_dict = self._make_parameter_dict(
                **{"lineModelSolve": {"lineModel": {"lya": {"profile": "asym", "asymProfile": {}}}}}
            )
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

        def test_ok_if_lya_profile_is_not_asym_and_asym_section_absent(self, zflag):
            param_dict = self._make_parameter_dict(
                **{"lineModelSolve": {"lineModel": {"lya": {"profile": "igm"}}}}
            )
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_any_warning()

    class TestUseLogLambdaSampling:
        def make_parameter_dict(self, **kwargs) -> dict:
            kwargs["method"] = "lineModelSolve"
            param_dict = make_parameter_dict_at_redshift_solver_level(**kwargs)
            param_dict["galaxy"]["templateDir"] = "sth"
            return param_dict

        def test_must_be_present_if_fftprocessing_and_tplfit(self):
            param_dict = self.make_parameter_dict(
                **{
                    "lineModelSolve": {
                        "lineModel": {
                            "continuumComponent": "tplFit",
                            "continuumFit": {
                                "fftProcessing": True,
                            },
                            "secondPass": {"continuumFit": {}},
                        }
                    }
                }
            )
            param_dict["continuumRemoval"] = {}
            with pytest.raises(
                APIException, match=r"Missing parameter galaxy lineModelSolve lineModel useLogLambdaSampling"
            ):
                check_from_parameter_dict(param_dict)

        def test_must_be_present_if_fftprocessing_and_tplfitauto(self):
            param_dict = self.make_parameter_dict(
                **{
                    "lineModelSolve": {
                        "lineModel": {
                            "continuumComponent": "tplFitAuto",
                            "continuumFit": {
                                "fftProcessing": True,
                            },
                            "secondPass": {"continuumFit": {}},
                        }
                    }
                }
            )
            param_dict["continuumRemoval"] = {}

            with pytest.raises(
                APIException, match=r"Missing parameter galaxy lineModelSolve lineModel useLogLambdaSampling"
            ):
                check_from_parameter_dict(param_dict)

        def test_unused_warning_if_present_but_no_fftprocessing(self, zflag):
            param_dict = self.make_parameter_dict(
                **{
                    "lineModelSolve": {
                        "lineModel": {
                            "continuumComponent": "tplFit",
                            "continuumFit": {
                                "fftProcessing": False,
                            },
                            "secondPass": {"continuumFit": {}},
                            "useLogLambdaSampling": True,
                        }
                    }
                }
            )
            param_dict["continuumRemoval"] = {}

            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_warning(WarningCode.UNUSED_PARAMETER)

        def test_unused_warning_if_present_but_no_tplfit(self, zflag):
            param_dict = self.make_parameter_dict(
                **{
                    "lineModelSolve": {
                        "lineModel": {
                            "continuumComponent": "sth",
                            "continuumFit": {
                                "fftProcessing": True,
                            },
                            "useLogLambdaSampling": True,
                        }
                    }
                }
            )
            param_dict["continuumRemoval"] = {}

            check_from_parameter_dict(param_dict)
            assert WarningUtils.has_warning(WarningCode.UNUSED_PARAMETER)

        def test_ok_if_all_present(self, zflag):
            param_dict = self.make_parameter_dict(
                **{
                    "lineModelSolve": {
                        "lineModel": {
                            "continuumComponent": "tplFit",
                            "continuumFit": {
                                "fftProcessing": True,
                            },
                            "useLogLambdaSampling": True,
                        }
                    }
                }
            )
            check_from_parameter_dict(param_dict)
            assert not WarningUtils.has_warning(WarningCode.UNUSED_PARAMETER)
