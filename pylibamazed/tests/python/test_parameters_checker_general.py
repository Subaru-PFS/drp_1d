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
from tests.python.utils import WarningUtils, make_parameter_dict_at_object_level


class TestParametersCheckGeneral:
    class TestFilters:
        class TestFiltersFormat:
            def test_ok_if_filter_format_is_correct(self, zflag):
                parametersDict = {"filters": [{"key": "errors", "instruction": "^", "value": "8"}]}
                CustomParametersChecker(parametersDict).check()
                assert not WarningUtils.has_any_warning()

            def test_error_if_filters_is_not_a_list(self):
                parametersDict = {"filters": {"key": "errors", "instruction": "^", "value": "8"}}
                with pytest.raises(APIException, match=r"Input filters json must be a list"):
                    CustomParametersChecker(parametersDict).check()

            def test_error_if_filters_is_missing_a_key(self):
                parametersDict = {"filters": [{"key": "errors", "instruction": "^"}]}
                with pytest.raises(APIException, match=r"Filters"):
                    CustomParametersChecker(parametersDict).check()

            def test_error_if_filters_has_an_additional_key(self):
                parametersDict = {
                    "filters": [{"key": "errors", "instruction": "^", "value": "8", "errorKey": "123"}]
                }
                with pytest.raises(APIException, match=r"Filters"):
                    CustomParametersChecker(parametersDict).check()

        class TestFiltersKeys:
            def test_no_error_if_no_filter(self, zflag):
                parametersDict = {}
                CustomParametersChecker(parametersDict).check()
                assert not WarningUtils.has_any_warning()

            def test_error_if_filter_uses_an_unknown_column(self):
                parametersDict = {"filters": [{"key": "zzz", "instruction": "^", "value": "8"}]}

                with pytest.raises(APIException, match=r"Unknown filter key zzz"):
                    CustomParametersChecker(parametersDict).check()

            def test_ok_if_filter_uses_a_default_or_additional_column(self, zflag):
                parametersDict = {
                    "filters": [
                        {"key": "errors", "instruction": "^", "value": "8"},
                        {"key": "zzz", "instruction": "^", "value": "8"},
                    ],
                    "additionalCols": ["zzz"],
                }

                CustomParametersChecker(parametersDict).check()
                assert not WarningUtils.has_any_warning()

    class TestPhotometryTransmissionDir:
        def _make_param_dict(self, **kwargs):
            new_kwargs = kwargs.copy()
            if kwargs.get("enablePhotometry"):
                new_kwargs["photometry"] = {"weight": 1}
            new_kwargs = {
                "stages": ["redshiftSolver"],
                "templateDir": "sth",
                "redshiftSolver": {
                    "method": "templateFittingSolve",
                    "templateFittingSolve": new_kwargs,
                },
            }
            param_dict = make_parameter_dict_at_object_level(**new_kwargs)
            if kwargs.get("enablePhotometry"):
                param_dict["photometryBand"] = "sth"

            return param_dict

        def test_photometry_enabled_without_transmission_dir_raises_an_error(self):
            param_dict = self._make_param_dict(**{"enablePhotometry": True})
            with pytest.raises(APIException, match=r"Missing parameter"):
                CustomParametersChecker(param_dict).check()

        def test_photometry_enabled_with_transmission_dir_ok(self, zflag):
            param_dict = self._make_param_dict(**{"enablePhotometry": True})
            param_dict["photometryTransmissionDir"] = "sth"
            CustomParametersChecker(param_dict).check()
            assert not WarningUtils.has_any_warning()

        def test_photometry_disabled_with_transmission_dir_raises_warning(self, zflag):
            param_dict = self._make_param_dict(**{"enablePhotometry": False})
            param_dict["photometryTransmissionDir"] = "sth"
            CustomParametersChecker(param_dict).check()
            assert WarningUtils.has_any_warning()

    class TestPhotometryBand:
        def _make_param_dict(self, **kwargs):
            new_kwargs = kwargs.copy()
            if kwargs.get("enablePhotometry"):
                new_kwargs["photometry"] = {"weight": 1}

            new_kwargs = {
                "stages": ["redshiftSolver"],
                "templateDir": "sth",
                "redshiftSolver": {
                    "method": "templateFittingSolve",
                    "templateFittingSolve": new_kwargs,
                },
            }

            param_dict = make_parameter_dict_at_object_level(**new_kwargs)
            if kwargs.get("enablePhotometry"):
                param_dict["photometryTransmissionDir"] = "sth"

            return param_dict

        def test_photometry_enabled_without_transmission_dir_raises_an_error(self):
            param_dict = self._make_param_dict(**{"enablePhotometry": True})
            with pytest.raises(APIException, match=r"Missing parameter"):
                CustomParametersChecker(param_dict).check()

        def test_photometry_enabled_with_transmission_dir_ok(self, zflag):
            param_dict = self._make_param_dict(**{"enablePhotometry": True})
            param_dict["photometryBand"] = "sth"
            CustomParametersChecker(param_dict).check()
            assert not WarningUtils.has_any_warning()

        def test_photometry_disabled_with_transmission_dir_raises_warning(self, zflag):
            param_dict = self._make_param_dict(**{"enablePhotometry": False})
            param_dict["photometryBand"] = "sth"
            CustomParametersChecker(param_dict).check()
            assert WarningUtils.has_any_warning()
