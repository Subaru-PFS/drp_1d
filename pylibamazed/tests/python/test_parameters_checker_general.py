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
from pylibamazed.ParametersAccessor import ParametersAccessor
from pylibamazed.ParametersChecker import ParametersChecker
from tests.python.utils import (WarningUtils,
                                make_parameter_dict_at_object_level)


class TestParametersCheckGeneral:
    class TestFilters:
        class TestFiltersFormat:
            def test_ok_if_filter_format_is_correct(self, zflag):
                parametersDict = {
                    "filters": [{"key": "Err", "instruction": "^", "value": "8"}]
                }

                accessor = ParametersAccessor(parametersDict)
                ParametersChecker(accessor)._check_filters()
                assert not WarningUtils.has_warning(zflag)

            def test_error_if_filters_is_not_a_list(self):
                parametersDict = {
                    "filters": {"key": "Err", "instruction": "^", "value": "8"}
                }

                accessor = ParametersAccessor(parametersDict)
                with pytest.raises(APIException, match=r"Input filters json must be a list"):
                    ParametersChecker(accessor)._check_filters()

            def test_error_if_filters_is_missing_a_key(self):
                parametersDict = {
                    "filters": [{"key": "Err", "instruction": "^"}]
                }
                accessor = ParametersAccessor(parametersDict)
                with pytest.raises(APIException, match=r"Filters"):
                    ParametersChecker(accessor)._check_filters()

            def test_error_if_filters_has_an_additional_key(self):
                parametersDict = {
                    "filters": [{"key": "Err", "instruction": "^", "value": "8", "errorKey": "123"}]
                }
                accessor = ParametersAccessor(parametersDict)
                with pytest.raises(APIException, match=r"Filters"):
                    ParametersChecker(accessor)._check_filters()

        class TestFiltersKeys:
            def test_no_error_if_no_filter(self, zflag):
                parametersDict = {}
                accessor = ParametersAccessor(parametersDict)

                ParametersChecker(accessor)._check_filters()
                assert not WarningUtils.has_warning(zflag)

            def test_error_if_filter_uses_an_unknown_column(self):
                parametersDict = {
                    "filters": [{"key": "zzz", "instruction": "^", "value": "8"}]
                }
                accessor = ParametersAccessor(parametersDict)

                with pytest.raises(APIException, match=r"Unknown filter key zzz"):
                    ParametersChecker(accessor)._check_filters()

            def test_ok_if_filter_uses_a_default_or_additional_column(self, zflag):
                parametersDict = {
                    "filters": [
                        {"key": "Err", "instruction": "^", "value": "8"},
                        {"key": "zzz", "instruction": "^", "value": "8"}
                    ],
                    "additional_cols": ["zzz"]
                }

                accessor = ParametersAccessor(parametersDict)
                ParametersChecker(accessor)._check_filters()
                assert not WarningUtils.has_warning(zflag)

    class TestPhotometryTransmissionDir:

        def _make_param_dict(self, **kwargs):
            new_kwargs = kwargs.copy()
            if kwargs.get("enablephotometry"):
                new_kwargs["photometry"] = {"weight": 1}
            new_kwargs = {"TemplateFittingSolve": new_kwargs}
            new_kwargs["method"] = "TemplateFittingSolve"

            param_dict = make_parameter_dict_at_object_level(**new_kwargs)
            if kwargs.get("enablephotometry"):
                param_dict["photometryBand"] = "sth"

            return param_dict

        def test_photometry_enabled_without_transmission_dir_raises_an_error(self):
            param_dict = self._make_param_dict(**{"enablephotometry": True})
            accessor = ParametersAccessor(param_dict)
            with pytest.raises(APIException, match=r"Missing parameter"):
                ParametersChecker(accessor).custom_check()

        def test_photometry_enabled_with_transmission_dir_ok(self, zflag):
            param_dict = self._make_param_dict(**{"enablephotometry": True})
            param_dict["photometryTransmissionDir"] = "sth"
            accessor = ParametersAccessor(param_dict)
            ParametersChecker(accessor).custom_check()
            assert not WarningUtils.has_warning(zflag)

        def test_photometry_disabled_with_transmission_dir_raises_warning(self, zflag):
            param_dict = self._make_param_dict(**{"enablephotometry": False})
            param_dict["photometryTransmissionDir"] = "sth"
            accessor = ParametersAccessor(param_dict)
            ParametersChecker(accessor).custom_check()
            assert WarningUtils.has_warning(zflag)

    class TestPhotometryBand:

        def _make_param_dict(self, **kwargs):
            new_kwargs = kwargs.copy()
            if kwargs.get("enablephotometry"):
                new_kwargs["photometry"] = {"weight": 1}
            new_kwargs = {"TemplateFittingSolve": new_kwargs}
            new_kwargs["method"] = "TemplateFittingSolve"

            param_dict = make_parameter_dict_at_object_level(**new_kwargs)
            if kwargs.get("enablephotometry"):
                param_dict["photometryTransmissionDir"] = "sth"

            return param_dict

        def test_photometry_enabled_without_transmission_dir_raises_an_error(self):
            param_dict = self._make_param_dict(**{"enablephotometry": True})
            accessor = ParametersAccessor(param_dict)
            with pytest.raises(APIException, match=r"Missing parameter"):
                ParametersChecker(accessor).custom_check()

        def test_photometry_enabled_with_transmission_dir_ok(self, zflag):
            param_dict = self._make_param_dict(**{"enablephotometry": True})
            print("param dict", param_dict)
            param_dict["photometryBand"] = "sth"

            accessor = ParametersAccessor(param_dict)
            ParametersChecker(accessor).custom_check()
            assert not WarningUtils.has_warning(zflag)

        def test_photometry_disabled_with_transmission_dir_raises_warning(self, zflag):
            param_dict = self._make_param_dict(**{"enablephotometry": False})
            param_dict["photometryBand"] = "sth"

            accessor = ParametersAccessor(param_dict)
            ParametersChecker(accessor).custom_check()
            assert WarningUtils.has_warning(zflag)
