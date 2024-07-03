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
from pylibamazed.Filter import FilterList, SpectrumFilterItem
from tests.python.utils import ComparisonUtils


class TestFilterList:
    def test_apply(self):
        df = pd.DataFrame({"col1": [1, 22, 111], "col2": [30, 1, 1]})

        # Result is as expected
        filter = FilterList()
        filter.add_filter(SpectrumFilterItem("col1", ">", 12))
        filter.add_filter(SpectrumFilterItem("col2", "<", 10))

        filtered = filter.apply(df)
        expected = pd.DataFrame({"col1": [22, 111], "col2": [1, 1]})
        ComparisonUtils.compare_dataframe_without_index(filtered, expected)

    def test_apply_without_items(self):
        # If filter list is empty, returns the initial data
        df = pd.DataFrame({"col1": [1, 22, 111], "col2": [30, 1, 1]})
        filter = FilterList()
        filtered = filter.apply(df)
        ComparisonUtils.compare_dataframe_without_index(filtered, df)

    def test_repr(self):
        filter = FilterList()
        assert filter.__repr__() == "FilterList []"

    class TestEquality:
        default_filter_item_1 = SpectrumFilterItem("col1", ">", 12)
        default_filter_item_2 = SpectrumFilterItem("col2", "<", 10)

        def test_eq_true_for_empty_filters(self):
            filter1 = FilterList()
            filter2 = FilterList()

            assert filter1 == filter2
            assert filter2 == filter1

        def test_eq_true_for_filters_containing_same_elements_in_same_order(self):
            filter1 = FilterList([self.default_filter_item_1, self.default_filter_item_2])
            filter2 = FilterList([self.default_filter_item_1, self.default_filter_item_2])

            assert filter1 == filter2
            assert filter2 == filter1

        def test_eq_false_for_filters_containing_same_elements_in_different_order(self):
            filter1 = FilterList([self.default_filter_item_1, self.default_filter_item_2])
            filter2 = FilterList([self.default_filter_item_2, self.default_filter_item_1])

            assert filter1 != filter2
            assert filter2 != filter1

        def test_eq_false_for_filters_with_different_length(self):
            filter1 = FilterList([self.default_filter_item_1])
            filter2 = FilterList([self.default_filter_item_1, self.default_filter_item_2])

            assert filter1 != filter2
            assert filter2 != filter1
