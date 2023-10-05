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
from typing import Callable, List

import pandas as pd
from pylibamazed.Exception import APIException
from pylibamazed.redshift import ErrorCode
from pylibamazed.Utils import LogicUtils


class FilterItem:
    """Creates a filter object.

    This filter object is composed of 3 items:
    - key: name of the column we want to filter on
    - instruction: what type of comparison we want to make
    - value: with which value we want to make the comparison
    """

    allowed_instructions = ["<", ">", "<=", ">=", "=", "in", "~in", "!=", "&", "~&", "0&", "^"]

    def __init__(self, key: str, instruction: str, value: any):
        self.check_instruction(instruction)
        self.key = key
        self.instruction = instruction
        self.value = value

    def __repr__(self):
        return "Filter " + str({
            "key": self.key,
            "instruction": self.instruction,
            "value": self.value,
        })

    def __eq__(self, __value: object) -> bool:
        return type(self) == type(__value) \
            and self.key == __value.key \
            and self.instruction == __value.instruction \
            and self.value == __value.value

    def compliant_lines(self, df: pd.DataFrame) -> pd.Series:
        if self.key not in df:
            raise APIException(
                ErrorCode.INVALID_FILTER_KEY,
                f"Column {self.key} does not exist"
            )
        comparator = self._comparator_from_instruction()
        return comparator(df[self.key])

    def _comparator_from_instruction(self) -> Callable:
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

    @classmethod
    def check_instruction(cls, instruction: str):
        if instruction not in cls.allowed_instructions:
            raise APIException(
                ErrorCode.INVALID_FILTER_INSTRUCTION,
                f"Instruction {instruction} is not registered."
                f"Allowed instructions are: {cls.allowed_instructions}"
            )


class SpectrumFilterItem(FilterItem):
    pass


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

    def apply(self, df: pd.DataFrame) -> pd.DataFrame:
        if not self.items:
            return df
        condition_list = [filt.compliant_lines(df) for filt in self.items]
        return df[LogicUtils.cumulate_conditions(condition_list)]
