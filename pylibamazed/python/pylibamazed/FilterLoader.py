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
from abc import abstractmethod
from typing import List

from pylibamazed.Filter import FilterItem, FilterList, SpectrumFilterItem
from pylibamazed.Parameters import Parameters


class AbstractFilterLoader:

    def __init__(
        self,
        json_keys_to_filter_attrs,
        FilterItemClass,
    ):
        # dictionary name of the key in FilterItem object: name of the key in json
        self.json_keys_to_filter_attrs = json_keys_to_filter_attrs
        self.FilterItemClass = FilterItemClass

    @abstractmethod
    def get_filters(self):
        pass


class ParamJsonFilterLoader(AbstractFilterLoader):
    """Reads & loads filter from input param file.

    Loaded param must have the following format:
    [
        {"key": "col1", "instruction": "<" , "value": 2},
        {"key": "col2", "instruction": ">=", "value": 2},
        ...
    ]
    """

    default_json_keys_to_filter_attrs = {
        "key": "key",
        "instruction": "instruction",
        "value": "value"
    }

    def __init__(
        self,
        json_keys_to_filter_attrs=default_json_keys_to_filter_attrs,
        FitlerItemClass=SpectrumFilterItem
    ):
        super().__init__(json_keys_to_filter_attrs, FitlerItemClass)

    def get_filters(self, params: Parameters) -> FilterList:
        filters: List[FilterItem] = []
        json_filters = params.get_filters()
        if json_filters is None:
            return None
        self._check_json_format(json_filters)
        for filter in json_filters:
            filters.append(self.FilterItemClass(filter["key"], filter["instruction"], filter["value"]))
        return FilterList(filters)

    def _check_json_format(self, json: json) -> None:
        if type(json) != list:
            raise ValueError("Input filters json must be a list")

        filter: List
        for filter in json:
            json_keys = self.json_keys_to_filter_attrs.values()
            different_keys = set(filter.keys()) != set(json_keys)
            different_length = len(filter.keys()) != len(json_keys)
            if different_keys or different_length:
                raise ValueError(
                    f"Each dictionary in json list must have exactly the following keys {json_keys}"
                )
