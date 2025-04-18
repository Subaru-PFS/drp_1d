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
    make_parameter_dict_at_reliability_solver_level,
)


class TestTemplateFittingSolve:
    def _make_parameter_dict(self, **kwargs) -> dict:
        param_dict = make_parameter_dict_at_reliability_solver_level(**kwargs)
        return param_dict

    def test_error_if_method_is_deeplLearningSolver_and_section_is_absent(self):
        param_dict = self._make_parameter_dict(**{"method": ["deepLearningSolver"]})
        with pytest.raises(APIException, match=r"Missing parameter galaxy deepLearningSolver"):
            check_from_parameter_dict(param_dict)

    def test_OK_if_method_is_deeplLearningSolver_and_section_is_present(self, zflag):
        param_dict = self._make_parameter_dict(**{"method": ["deepLearningSolver"], "deepLearningSolver": {}})
        check_from_parameter_dict(param_dict)
        assert not WarningUtils.has_any_warning()

    def test_warning_if_method_is_not_deeplLearningSolver_but_section_is_present(self, zflag):
        param_dict = self._make_parameter_dict(**{"method": ["lalala"], "deepLearningSolver": {}})
        check_from_parameter_dict(param_dict)
        assert WarningUtils.has_any_warning()

    def test_error_if_method_is_skLearnClassifier_and_section_is_absent(self):
        param_dict = self._make_parameter_dict(**{"method": ["skLearnClassifier"]})
        with pytest.raises(APIException, match=r"Missing parameter galaxy skLearnClassifier"):
            check_from_parameter_dict(param_dict)

    def test_OK_if_method_is_skLearnClassifier_and_section_is_present(self, zflag):
        param_dict = self._make_parameter_dict(**{"method": ["skLearnClassifier"], "skLearnClassifier": {}})
        check_from_parameter_dict(param_dict)
        assert not WarningUtils.has_any_warning()

    def test_warning_if_method_is_not_skLearnClassifier_but_section_is_present(self, zflag):
        param_dict = self._make_parameter_dict(**{"method": ["lalala"], "skLearnClassifier": {}})
        check_from_parameter_dict(param_dict)
        assert WarningUtils.has_any_warning()

    def test_error_if_one_section_absent(self):
        param_dict = self._make_parameter_dict(
            **{"method": ["skLearnClassifier", "deepLearningSolver"], "deepLearningSolver": {}}
        )
        with pytest.raises(APIException, match=r"Missing parameter galaxy skLearnClassifier"):
            check_from_parameter_dict(param_dict)

    def test_OK_if_both_methods_and_sections_are_present(self, zflag):
        param_dict = self._make_parameter_dict(
            **{
                "method": ["skLearnClassifier", "deepLearningSolver"],
                "deepLearningSolver": {},
                "skLearnClassifier": {},
            }
        )
        check_from_parameter_dict(param_dict)
        assert not WarningUtils.has_any_warning()
