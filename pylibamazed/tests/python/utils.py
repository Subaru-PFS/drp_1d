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
from pylibamazed.CustomParametersChecker import CustomParametersChecker
from pylibamazed.redshift import CFlagWarning, WarningCode
from typing import Optional, Union, Dict

default_object_type = "galaxy"


class ComparisonUtils:
    @staticmethod
    def compare_dataframe_without_index(df1, df2):
        df1_without_index = df1.reset_index(drop=True)
        df2_without_index = df2.reset_index(drop=True)
        assert df1_without_index.equals(df2_without_index)


class WarningUtils:
    @staticmethod
    def has_any_warning():
        return bool(CFlagWarning.GetInstance().getBitMask())

    def has_warning(warning: WarningCode):
        return warning in CFlagWarning.GetInstance().getBitMask()


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
        "lineMeasRunMode": "pipe",
        "spectrumModels": [default_object_type],
        default_object_type: kwargs,
    }
    return param_dict


def make_parameter_dict_at_redshift_solver_level(
    object_level_params=None, object_type: Optional[str] = None, **redshift_kwargs
) -> dict:
    if object_type is None:
        object_type = default_object_type
    param_dict: Dict[str, Union[list, Dict]] = {
        "version": 2,
        "spectrumModels": [object_type],
        object_type: {
            "stages": ["redshiftSolver"],
            "redshiftSolver": redshift_kwargs,
        },
    }
    if object_level_params is not None:
        for key, val in object_level_params.items():
            param_dict[object_type][key] = val
    if redshift_kwargs.get("method") == "lineModelSolve":
        if param_dict.get("lsf") is None:
            param_dict["lsf"] = {}
    tfs = redshift_kwargs.get("templateFittingSolve")
    if (
        redshift_kwargs.get("method") == "templateFittingSolve"
        and tfs is not None
        and tfs.get("singlePass") is not False
    ):
        param_dict[object_type]["redshiftSolver"]["templateFittingSolve"]["singlePass"] = True  # type: ignore
    return param_dict


def make_parameter_dict_at_linemodelsolve_level(**kwargs):
    param_dict = make_parameter_dict_at_redshift_solver_level()
    param_dict[default_object_type]["redshiftSolver"] = {
        "method": "lineModelSolve",
        "lineModelSolve": {"lineModel": kwargs},
    }
    param_dict["lsf"] = {}
    return param_dict


def make_parameter_dict_at_linemeas_solve_level(object_level_params=None, **kwargs) -> dict:
    param_dict: Dict = {
        "lineMeasRunMode": "pipe",
        "spectrumModels": [default_object_type],
        default_object_type: {
            "lineMeasDzHalf": 0.1,
            "lineMeasRedshiftStep": 0.1,
            "stages": ["lineMeasSolver"],
            "lineMeasSolver": {
                "method": "lineMeasSolve",
                "lineMeasSolve": kwargs,
            },
        },
    }
    param_dict["lsf"] = {}
    if object_level_params is not None:
        for key, val in object_level_params.items():
            param_dict[default_object_type][key] = val
    return param_dict


def make_parameter_dict_linemeas_solve_piped_linemodel(
    linemodel_level_params: Dict, linemeas_level_params: Dict
) -> dict:
    param_dict = make_parameter_dict_at_linemodelsolve_level(**linemodel_level_params)
    param_dict["lineMeasRunMode"] = "pipe"
    param_dict[default_object_type]["stages"] += ["lineMeasSolver"]
    linemeas_dict = make_parameter_dict_at_linemeas_solve_level(**linemeas_level_params)
    del linemeas_dict[default_object_type]["stages"]
    param_dict[default_object_type] = param_dict[default_object_type] | linemeas_dict[default_object_type]
    return param_dict


def make_parameter_dict_at_reliability_solver_level(object_level_params=None, **kwargs) -> dict:
    param_dict: Dict = {
        "spectrumModels": [default_object_type],
        default_object_type: {
            "stages": ["reliabilitySolver"],
            "reliabilitySolver": {**kwargs},
        },
    }
    if object_level_params is not None:
        for key, val in object_level_params.items():
            param_dict[default_object_type][key] = val
    return param_dict


def make_parameter_dict_at_reliability_deep_learning_level(object_level_params=None, **kwargs) -> dict:
    param_dict: Dict = {
        "spectrumModels": [default_object_type],
        default_object_type: {
            "stages": ["reliabilitySolver"],
            "reliabilitySolver": {
                "method": "deepLearningSolver",
                "deepLearningSolver": kwargs,
            },
        },
    }
    if object_level_params is not None:
        for key, val in object_level_params.items():
            param_dict[default_object_type][key] = val
    return param_dict


def check_from_parameter_dict(param_dict: dict):
    CustomParametersChecker(param_dict).check()
