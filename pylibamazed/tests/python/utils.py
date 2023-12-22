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
from pylibamazed.ParametersAccessor import ParametersAccessor
from pylibamazed.ParametersChecker import ParametersChecker
from pylibamazed.redshift import CFlagWarning
from pylibamazed.Warning import extract_warning_flags

default_object_type = "galaxy"


class ComparisonUtils:
    @staticmethod
    def compare_dataframe_without_index(df1, df2):
        df1_without_index = df1.reset_index(drop=True)
        df2_without_index = df2.reset_index(drop=True)
        assert df1_without_index.equals(df2_without_index)


class WarningUtils:
    @staticmethod
    def has_any_warning(zflag: CFlagWarning):
        return zflag.getBitMask() != 0

    def has_warning(zflag: CFlagWarning, warning_name):
        return warning_name in extract_warning_flags(zflag.getBitMask())


class DictUtils:
    @staticmethod
    def make_nested(dict: dict, nested: str) -> dict:
        if nested is None:
            return dict
        else:
            return {nested: dict}


def make_parameter_dict(**kwargs) -> dict:
    param_dict = {}

    for key, value in kwargs.items():
        param_dict[key] = value
    return param_dict


def make_parameter_dict_at_object_level(**kwargs) -> dict:
    param_dict = {
        "spectrumModels": [default_object_type],
        default_object_type: kwargs
    }
    return param_dict


def make_parameter_dict_at_redshift_solver_level(object_level_params=None, **redshift_kwargs) -> dict:
    param_dict = {
        "spectrumModels": [default_object_type],
        default_object_type: {
            "stages": ["redshiftSolver"],
            "redshiftSolver": redshift_kwargs,
        }
    }
    if object_level_params is not None:
        for key, val in object_level_params.items():
            param_dict[default_object_type][key] = val
    return param_dict


def make_parameter_dict_at_linemeas_solve_level(object_level_params=None, **kwargs) -> dict:
    param_dict = {
        "spectrumModels": [default_object_type],
        default_object_type: {
            "lineMeasDzHalf": 0.1,
            "lineMeasRedshiftStep": 0.1,
            "stages": ["lineMeasSolver"],
            "lineMeasSolver": {
                "method": "lineMeasSolve",
                "lineMeasSolve": kwargs,
            }
        }
    }
    if object_level_params is not None:
        for key, val in object_level_params.items():
            param_dict[default_object_type][key] = val
    return param_dict


def make_parameter_dict_at_reliability_deep_learning_level(object_level_params=None, **kwargs) -> dict:
    param_dict = {
        "spectrumModels": [default_object_type],
        default_object_type: {
            "stages": ["reliabilitySolver"],
            "reliabilitySolver": {
                "method": "deepLearningSolver",
                "deepLearningSolver": kwargs,
            }
        }
    }
    if object_level_params is not None:
        for key, val in object_level_params.items():
            param_dict[default_object_type][key] = val
    return param_dict


def check_from_parameter_dict(param_dict: dict):
    accessor = ParametersAccessor(param_dict)
    ParametersChecker().custom_check(accessor)
