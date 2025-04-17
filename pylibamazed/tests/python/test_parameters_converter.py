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
import json
import os

import pytest
from pylibamazed.Exception import APIException
from pylibamazed.ParametersConverter import (
    ParametersConverter,
    ParametersConverterSelector,
    ParametersConverterV2,
)
from tests.python.config import test_dir

# For investigating errors only
# from deepdiff import DeepDiff


class TestParametersConverterSelector:
    selector = ParametersConverterSelector()

    def test_indicated_version_is_used_ex_2(self):
        Converter: ParametersConverter = self.selector.get_converter(2)
        assert isinstance(Converter(), ParametersConverterV2)

    def test_error_if_unknown_version(self):
        with pytest.raises(APIException, match=r"Unexpected parameters version"):
            self.selector.get_converter(100)


class TestParametersConverterV2:
    converter = ParametersConverterV2()

    def test_returns_same_than_input(self):
        raw_params = {"aaa": "b", "cc": 2}
        assert self.converter.convert(raw_params) == raw_params

    def test_spectrum_models_renamed(self):
        raw_params = {"spectrumModel_galaxy": {}, "spectrumModel_other": {"lala": "yopyop"}}
        expected = {"galaxy": {}, "other": {"lala": "yopyop"}}

        assert self.converter.convert(raw_params) == expected
