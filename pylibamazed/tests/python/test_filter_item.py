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

import pandas as pd
import numpy as np
import pytest
from pylibamazed.Exception import APIException
from pylibamazed.Filter import FilterItem, FilterMorphology


class TestFilterItem:
    def test_init(self):
        # Init ok if known instruction
        FilterItem("T", "<", 1)

        # Init raises error if unkown instruction
        with pytest.raises(APIException, match=r"INVALID_FILTER_INSTRUCTION"):
            FilterItem("T", "unkown", 1)

    def test_apply(self):
        df = pd.DataFrame({"col1": [1, 2, 3], "col2": [4, 5, 6]})

        filter = FilterItem("col1", "<", 2)
        assert filter.apply(df).equals(pd.Series([True, False, False]))

        filter.instruction = ">"
        assert filter.apply(df).equals(pd.Series([False, False, True]))

        filter.instruction = "<="
        assert filter.apply(df).equals(pd.Series([True, True, False]))

        filter.instruction = ">="
        assert filter.apply(df).equals(pd.Series([False, True, True]))

        filter.instruction = "="
        assert filter.apply(df).equals(pd.Series([False, True, False]))

        filter.instruction = "!="
        assert filter.apply(df).equals(pd.Series([True, False, True]))

        filter.instruction = "in"
        filter.value = [2, 3, 4]
        assert filter.apply(df).equals(pd.Series([False, True, True]))

        filter.instruction = "~in"
        filter.value = [2, 3, 4]
        assert filter.apply(df).equals(pd.Series([True, False, False]))

        filter.instruction = "&"
        filter.value = 1
        assert filter.apply(df).equals(pd.Series([True, False, True]))

        df2 = pd.DataFrame({"col1": [0, 1, 2, 3]})
        filter.instruction = "0&"
        assert filter.apply(df2).equals(pd.Series([True, True, False, True]))

        filter.instruction = "^"
        assert filter.apply(df2).equals(pd.Series([True, False, True, True]))

        filter.instruction = "~^"
        assert filter.apply(df2).equals(pd.Series([False, True, False, False]))

        filter.value = 1
        filter.instruction = "~&"
        assert filter.apply(df2).equals(pd.Series([True, False, True, False]))

        filter.key = "unexistant col"
        with pytest.raises(APIException, match=r"INVALID_FILTER_KEY"):
            filter.apply(df)

        with pytest.raises(APIException, match=r"INTERNAL_ERROR"):
            filter.apply()

    def test_repr(self):
        filter = FilterItem("col1", "<", 2)
        assert filter.__repr__() == "Filter {'key': 'col1', 'instruction': '<', 'value': 2}"

    def test_check_instructions(self):
        # Instruction is registered -> OK
        FilterItem.check_instruction(">")

        # Instruction is not registered -> error
        with pytest.raises(APIException, match=r"INVALID_FILTER_INSTRUCTION"):
            FilterItem.check_instruction("unkown")


class TestFilterMorphology:
    def test_init(self):
        FilterMorphology("opening", [1, 1])

    # Init raises error if unkown instruction
    with pytest.raises(APIException, match=r"INVALID_FILTER_INSTRUCTION"):
        FilterMorphology("unkown", [1, 1])

    def test_apply_opening(self):
        filter = FilterMorphology("opening", [1, 1])
        mask = pd.Series(np.array([True, True, False], dtype=bool))

        assert np.array_equal(filter.apply(mask=mask), mask)  # no change

        mask = pd.Series([0, 1, 0])
        assert np.array_equal(filter.apply(mask=mask), [0, 0, 0])  # isolated removed

        mask = pd.Series([1, 0, 0])
        assert np.array_equal(filter.apply(mask=mask), [0, 0, 0])  # isolated on left border removed

        mask = pd.Series([0, 0, 1])
        assert np.array_equal(filter.apply(mask=mask), [0, 0, 0])  # isolated on right border removed

        # apply without mask raises error
        with pytest.raises(APIException, match=r"INTERNAL_ERROR"):
            filter.apply()

    def test_apply_closing(self):
        filter = FilterMorphology("closing", [1, 1])
        mask = pd.Series(np.array([False, False, True], dtype=bool))

        assert np.array_equal(filter.apply(mask=mask), mask)  # no change

        mask = pd.Series([0, 1, 0, 1])
        assert np.array_equal(filter.apply(mask=mask), [0, 1, 1, 1])  # isolated removed

        mask = pd.Series([0, 1, 1])
        assert np.array_equal(filter.apply(mask=mask), [0, 1, 1])  # isolated on left border not removed

        mask = pd.Series([0, 1, 1, 0])
        assert np.array_equal(filter.apply(mask=mask), [0, 1, 1, 0])  # isolated on right border not removed

    def test_apply_erosion(self):
        filter = FilterMorphology("erosion", [1, 1, 1])
        mask = pd.Series(np.array([False, True, True, True, False], dtype=bool))

        assert np.array_equal(filter.apply(mask=mask), [0, 0, 1, 0, 0])

    def test_apply_dilation(self):
        filter = FilterMorphology("dilation", [1, 1, 1])
        mask = pd.Series(np.array([False, False, True, False, False], dtype=bool))

        assert np.array_equal(filter.apply(mask=mask), [0, 1, 1, 1, 0])
