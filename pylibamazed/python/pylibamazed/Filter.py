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
from typing import Callable, List, Any, Optional
from abc import ABCMeta, abstractmethod
from scipy.ndimage import binary_closing, binary_opening, binary_dilation, binary_erosion
import pandas as pd
import numpy as np

from pylibamazed.Exception import APIException
from pylibamazed.redshift import ErrorCode


class AbstractFilterItem(metaclass=ABCMeta):
    """Creates a filter object.

    This filter object is composed of 3 items:
    - key: name of the column we want to filter on
    - instruction: what type of comparison we want to make
    - value: with which value we want to make the comparison
    """

    allowed_instructions = []

    def __init__(self, key: str, instruction: str, value: Any):
        self.check_instruction(instruction)
        self.key = key
        self.instruction = instruction
        self.value = value

    def __repr__(self):
        return "Filter " + str(
            {
                "key": self.key,
                "instruction": self.instruction,
                "value": self.value,
            }
        )

    def __eq__(self, __value: object) -> bool:
        return (
            type(self) is type(__value)
            and self.key == __value.key
            and self.instruction == __value.instruction
            and self.value == __value.value
        )

    @abstractmethod
    def apply(self, df: Optional[pd.DataFrame] = None, mask: Optional[pd.Series] = None) -> pd.Series:
        raise NotImplementedError("Implement in derived class")

    @abstractmethod
    def _action_from_instruction(self) -> Callable:
        raise NotImplementedError("Implement in derived class")

    @classmethod
    def check_instruction(cls, instruction: str):
        if instruction not in cls.allowed_instructions:
            raise APIException(
                ErrorCode.IE_INVALID_FILTER_INSTRUCTION,
                f"Instruction {instruction} is not registered."
                f"Allowed instructions are: {cls.allowed_instructions}",
            )


def filterFactory(filterDict: dict) -> AbstractFilterItem:
    filterType = filterDict.get("type", "byValue")
    if filterType == "byValue":
        return FilterItem(filterDict["key"], filterDict["instruction"], filterDict["value"])
    elif filterType == "morphology":
        return FilterMorphology(filterDict["instruction"], filterDict["value"])
    else:
        raise APIException(ErrorCode.INTERNAL_ERROR, f"Wrong filter type: {filterType}")


class FilterItem(AbstractFilterItem):
    allowed_instructions = ["<", ">", "<=", ">=", "=", "in", "~in", "!=", "&", "~&", "0&", "^"]

    def apply(self, df: Optional[pd.DataFrame] = None, mask: Optional[pd.Series] = None) -> pd.Series:
        if df is None:
            raise APIException(ErrorCode.INTERNAL_ERROR, "df parameter is None, should be a pandas DataFrame")
        if self.key not in df:
            raise APIException(ErrorCode.IE_INVALID_FILTER_KEY, f"Column {self.key} does not exist")
        action = self._action_from_instruction()
        newMask = action(df[self.key])
        if mask is None:
            return newMask
        else:
            return newMask & mask

    def _action_from_instruction(self) -> Callable:
        str_to_action = {
            "<": self._inf,
            ">": self._sup,
            ">=": self._sup_equal,
            "<=": self._inf_equal,
            "=": self._equal,
            "!=": self._different,
            "in": self._is_in,
            "~in": self._is_not_in,
            "&": self._bitwise_and,
            "~&": self._bitwise_not_and,
            "0&": self._bitwise_and_or_0,
            "^": self._bitwise_xor,
            "~^": self._bitwise_not_xor,
        }
        return str_to_action[self.instruction]

    def _inf(self, a):
        return a < self.value

    def _sup(self, a):
        return a > self.value

    def _inf_equal(self, a):
        return a <= self.value

    def _sup_equal(self, a):
        return a >= self.value

    def _equal(self, a):
        return a == self.value

    def _different(self, a):
        return a != self.value

    def _is_in(self, a):
        return a.isin(self.value)

    def _is_not_in(self, a):
        return ~a.isin(self.value)

    def _bitwise_and(self, a):
        return (a & self.value).astype(bool)

    def _bitwise_not_and(self, a):
        return ~((a & self.value)).astype(bool)

    def _bitwise_and_or_0(self, a):
        return (a == 0) | (a & self.value).astype(bool)

    def _bitwise_xor(self, a):
        return (a ^ self.value).astype(bool)

    def _bitwise_not_xor(self, a):
        return ~(a ^ self.value).astype(bool)


class FilterMorphology(AbstractFilterItem):
    allowed_instructions = ["opening", "closing", "erosion", "dilation"]

    def __init__(self, instruction: str, value: List):
        super().__init__("unused", instruction, value)

    def apply(self, df: pd.DataFrame = None, mask: pd.Series = None) -> pd.Series:
        if mask is None:
            raise APIException(ErrorCode.INTERNAL_ERROR, "mask parameter is None, should be a pandas Series")
        action = self._action_from_instruction()
        return action(mask)

    def _action_from_instruction(self) -> Callable:
        str_to_action = {
            "opening": self._opening,
            "closing": self._closing,
            "erosion": self._erosion,
            "dilation": self._dilation,
        }
        return str_to_action[self.instruction]

    def _opening(self, a):
        return binary_opening(a, self.value)

    def _closing(self, a):
        return binary_closing(a, self.value)

    def _erosion(self, a):
        return binary_erosion(a, self.value)

    def _dilation(self, a):
        return binary_dilation(a, self.value)


class FilterList:
    def __init__(self, filters=[]):
        self.items: List = filters[:]

    def __repr__(self):
        return "FilterList " + str(self.items)

    def __eq__(self, __value__):
        if len(__value__.items) != len(self.items):
            return False
        for i, filter in enumerate(self.items):
            if filter != __value__.items[i]:
                return False
        return True

    def add_filter(self, filter) -> None:
        self.items.append(filter)

    def apply(self, df: pd.DataFrame):
        if not self.items:
            return None
        currentMask = pd.Series(np.ones(len(df), dtype=bool), index=df.index)
        for filt in self.items:
            currentMask = filt.apply(df, currentMask)
        return currentMask
