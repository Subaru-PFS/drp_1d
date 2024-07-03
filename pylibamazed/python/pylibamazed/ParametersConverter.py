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
from abc import ABC, abstractmethod

import pandas as pd
from pylibamazed.Exception import APIException
from pylibamazed.Paths import v1_to_treed_filename
from pylibamazed.redshift import ErrorCode


class ParametersConverter(ABC):
    """Converts raw input parameters into treed parameters, used for custom checks and pylibamazed
    calculations
    """

    @abstractmethod
    def convert(self, raw_params: dict) -> dict:
        pass


class ParametersConverterSelector:
    def get_converter(self, version: int) -> ParametersConverter:
        if version == 1:
            Converter = ParametersConverterV1
        elif version == 2:
            Converter = ParametersConverterV2
        else:
            raise APIException(ErrorCode.INVALID_PARAMETER_FILE, f"Unexpected parameters version {version}")
        return Converter


class ParametersConverterV2(ParametersConverter):
    def convert(self, raw_params: dict) -> dict:
        params = raw_params.copy()

        for key in raw_params:
            if "spectrumModel_" in key:
                new_key = key.replace("spectrumModel_", "")
                params[new_key] = params.pop(key)

        return params


class ParametersConverterV1(ParametersConverter):
    def convert(self, raw_params: dict) -> dict:
        params = raw_params.copy()
        params_str = json.dumps(params)
        renaming_df = pd.read_csv(v1_to_treed_filename)
        for _, row in renaming_df.iterrows():
            params_str = params_str.replace(f"\"{row.loc['old_name']}\"", f"\"{row.loc['new_name']}\"")
        renamed_params = json.loads(params_str)

        self.update_redshift_part(renamed_params)
        self.update_linemeas_part(renamed_params)
        self.update_reliability_part(renamed_params)
        return renamed_params

    def update_redshift_part(self, renamed_params):
        redshift_methods = ["lineModelSolve", "templateFittingSolve", "tplCombinationSolve"]
        for spectrum_model in renamed_params.get("spectrumModels", []):
            method = renamed_params[spectrum_model].pop("method", None)
            if method in redshift_methods:
                renamed_params[spectrum_model].setdefault("stages", []).append("redshiftSolver")
                renamed_params[spectrum_model]["redshiftSolver"] = {"method": method}
            for redshift_method in redshift_methods:
                if redshift_method in renamed_params[spectrum_model]:
                    renamed_params[spectrum_model].setdefault("redshiftSolver", {})
                    renamed_params[spectrum_model]["redshiftSolver"][redshift_method] = renamed_params[
                        spectrum_model
                    ].pop(redshift_method)

    def update_linemeas_part(self, renamed_params):
        for spectrum_model in renamed_params.get("spectrumModels", []):
            if renamed_params[spectrum_model].pop("linemeas_method", None) == "lineMeasSolve":
                renamed_params[spectrum_model].setdefault("stages", []).append("lineMeasSolver")
                renamed_params[spectrum_model]["lineMeasSolver"] = {
                    "method": "lineMeasSolve",
                    "lineMeasSolve": renamed_params[spectrum_model].pop("lineMeasSolve"),
                }

    def update_reliability_part(self, renamed_params):
        for spectrum_model in renamed_params.get("spectrumModels", []):
            if renamed_params[spectrum_model].pop("enable_reliability", False):
                renamed_params[spectrum_model].setdefault("stages", []).append("reliabilitySolver")
                renamed_params[spectrum_model]["reliabilitySolver"] = {
                    "method": "deepLearningSolver",
                    "deepLearningSolver": {
                        "reliabilityModel": renamed_params[spectrum_model].pop("reliabilityModel")
                    },
                }
            else:
                renamed_params[spectrum_model].pop("reliabilityModel", None)
