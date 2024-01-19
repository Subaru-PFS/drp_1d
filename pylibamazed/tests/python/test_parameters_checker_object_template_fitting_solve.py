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
                                make_parameter_dict_at_redshift_solver_level)


class TestTemplateFittingSolve:

    def _make_parameter_dict(self, **kwargs) -> dict:
        kwargs["linemeas_method"] = kwargs.get("linemeas_method", "")
        kwargs["method"] = kwargs.get("method", "templateFittingSolve")
        object_level_params = None
        if kwargs["method"] == "templateFittingSolve":
            object_level_params = {"templateDir": "sth"}
        param_dict = make_parameter_dict_at_redshift_solver_level(object_level_params, **kwargs)
        if kwargs.get("templateFittingSolve", {}).get("enablePhotometry"):
            param_dict["photometryTransmissionDir"] = "sth"
            param_dict["photometryBand"] = ["someBand"]
        return param_dict

    def test_error_if_method_is_templateFittingSolve_and_section_is_absent(self):
        param_dict = self._make_parameter_dict(**{})
        with pytest.raises(APIException, match=r"Missing parameter templateFittingSolve"):
            check_from_parameter_dict(param_dict)

    def test_OK_if_method_is_templateFittingSolve_and_section_is_present(self, zflag):
        param_dict = self._make_parameter_dict(**{
            "templateFittingSolve": {}
        })
        check_from_parameter_dict(param_dict)
        assert not WarningUtils.has_any_warning(zflag)

    def test_warning_if_method_is_not_templateFittingSolve_but_section_is_present(self, zflag):
        param_dict = self._make_parameter_dict(**{
            "method": "sth",
            "templateFittingSolve": {}
        })
        check_from_parameter_dict(param_dict)
        assert WarningUtils.has_any_warning(zflag)

    def test_ok_if_method_is_not_templateFittingSolve_and_section_is_absent(self, zflag):
        param_dict = self._make_parameter_dict(**{
            "method": "sth",
        })
        check_from_parameter_dict(param_dict)
        assert not WarningUtils.has_any_warning(zflag)

    def test_error_if_ismfit_enabled_and_ebmv_section_is_not_present(self):
        param_dict = self._make_parameter_dict(**{
            "templateFittingSolve": {"ismFit": True}
        })
        with pytest.raises(APIException, match=r"Missing parameter ebmv"):
            check_from_parameter_dict(param_dict)

    def test_ok_if_ismfit_enabled_and_ebmv_is_present(self, zflag):
        param_dict = self._make_parameter_dict(**{
            "templateFittingSolve": {"ismFit": True}
        })
        param_dict["ebmv"] = {}
        check_from_parameter_dict(param_dict)
        assert not WarningUtils.has_any_warning(zflag)

    def test_error_if_photometry_is_enabled_but_photometry_weight_is_absent(self):
        param_dict = self._make_parameter_dict(**{
            "templateFittingSolve": {"enablePhotometry": True}
        })
        with pytest.raises(
            APIException,
            match=r"Missing parameter object galaxy TemplateFittingSolve photometry weight"
        ):
            check_from_parameter_dict(param_dict)

    def test_ok_if_photometry_is_enabled_and_photometry_weight_too(self, zflag):
        param_dict = self._make_parameter_dict(**{
            "templateFittingSolve": {
                "enablePhotometry": True,
                "photometry": {"weight": 1}
            }
        })
        check_from_parameter_dict(param_dict)
        assert not WarningUtils.has_any_warning(zflag)

    def test_warning_if_photometry_is_disabled_but_photometry_weight_is_present(self, zflag):
        param_dict = self._make_parameter_dict(**{
            "templateFittingSolve": {
                "enablePhotometry": False,
                "photometry": {"weight": 1}
            }
        })
        check_from_parameter_dict(param_dict)
        assert WarningUtils.has_any_warning(zflag)

    def test_ok_if_photometry_is_disabled_and_photometry_weight_is_absent(self, zflag):
        param_dict = self._make_parameter_dict(**{
            "templateFittingSolve": {
                "enablePhotometry": False,
                "photometry": {}
            }
        })
        check_from_parameter_dict(param_dict)
        assert not WarningUtils.has_any_warning(zflag)
