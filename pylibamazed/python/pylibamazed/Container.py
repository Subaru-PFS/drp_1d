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
from typing import Generic, TypeVar

T = TypeVar("T")
class Container(Generic[T]):
    def __init__(self, **kwargs):
        self.data: dict[str, T] = dict()
        if(kwargs):
            for key in kwargs:
                self.append(kwargs[key], key)

    def __repr__(self):
        repr = ''
        for key in self.data:
            repr += f"\n{key}: {self.data[key]}"
        return repr
    
    def __eq__(self, __value__):
        if type(self) != type(__value__):
            return False
        if len(self.data) != len(__value__.data):
            return False
        for key in self.data:
            if self.data.get(key) != __value__.get(key):
                return False
        return True

    def append(self, dataToAppend: T, obs_id=""):
        self.data[obs_id] = dataToAppend
    
    def get(self, obs_id="") -> T:
        return self.data[obs_id]

    def keys(self):
        return self.data.keys()

    def size(self) -> int:
        return len(self.data)
