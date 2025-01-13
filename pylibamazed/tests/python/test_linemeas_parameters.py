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
import pytest
from pylibamazed.AbstractOutput import AbstractOutput
from pylibamazed.Exception import APIException
from pylibamazed.redshift import WarningCode
from pylibamazed.Parameters import Parameters
from pylibamazed.ParametersAccessor import EVelocityType, EVelocityFitParam
from pylibamazed.LinemeasParameters import LinemeasParameters
from tests.python.test_parameters_utils import TestParametersUtils
from tests.python import utils


class TestLinemeasParameters:
    generic_parameters: Parameters = Parameters(
        TestParametersUtils().make_parameters_dict(), make_checks=False
    )
    object_type: str = TestParametersUtils.default_object_type

    def test_load_from_catalog(self, mocker):
        source_id = "some_source_id"
        mocker.patch(
            "pandas.read_csv",
            return_value=pd.DataFrame(
                {
                    "ProcessingID": [source_id],
                    "redshift_column_name": [0.0],
                    "velocity_absorption_column_name": [1.0],
                    "velocity_emission_column_name": [2.0],
                }
            ),
        )
        config = TestParametersUtils().make_config()
        lp = LinemeasParameters()
        lp.load_from_catalogs(source_id, config["linemeascatalog"], config["linemeas_catalog_columns"])
        assert lp.redshift_ref[self.object_type] == 0.0
        assert lp.velocity[EVelocityType.Absorption][self.object_type] == 1.0
        assert lp.velocity[EVelocityType.Emission][self.object_type] == 2.0

        # Source id absent from loaded dataframe
        with pytest.raises(APIException):
            lp.load_from_catalogs(
                "some absent source id", config["linemeascatalog"], config["linemeas_catalog_columns"]
            )

    def test_load_from_result_store(self, mocker):
        def get_attribute_from_source_mocked(
            object_type, stage, method, dataset, attribute, rank=None, band_name=None, obs_id=None
        ):
            map = {"Redshift": 0.0, "VelocityAbsorption": 1.0, "VelocityEmission": 2.0}
            return map[attribute]

        FakeOutput = mocker.Mock(spec=AbstractOutput)
        FakeOutput.get_attribute_from_source = get_attribute_from_source_mocked
        lp = LinemeasParameters()
        lp.load_from_result_store(self.generic_parameters, FakeOutput, self.object_type)
        assert lp.redshift_ref[self.object_type] == 0.0
        assert lp.velocity[EVelocityType.Absorption][self.object_type] == 1.0
        assert lp.velocity[EVelocityType.Emission][self.object_type] == 2.0

    def _make_parameter_dict_for_update_parameters(self) -> Parameters:
        param_dict = utils.make_parameter_dict_at_linemeas_solve_level(
            **{
                "lineModel": {
                    "velocityFit": True,
                    "emVelocityFitMin": 500,
                    "emVelocityFitMax": 2000,
                    "absVelocityFitMin": 150,
                    "absVelocityFitMax": 500,
                }
            }
        )
        param_dict["version"] = 2
        parameters = Parameters(param_dict, make_checks=False)
        return parameters

    def test_update_parameters(self):
        parameters = self._make_parameter_dict_for_update_parameters()
        lp = LinemeasParameters()
        lp.redshift_ref[utils.default_object_type] = 1.0
        lp.velocity[EVelocityType.Absorption][utils.default_object_type] = 200
        lp.velocity[EVelocityType.Emission][utils.default_object_type] = 1000
        lp.update_parameters(parameters)
        assert parameters.get_spectrum_model_section(utils.default_object_type)["redshiftref"] == 1.0
        assert (
            parameters.get_velocity(utils.default_object_type, "lineMeasSolve", EVelocityType.Absorption)
            == 200
        )
        assert (
            parameters.get_velocity(utils.default_object_type, "lineMeasSolve", EVelocityType.Emission)
            == 1000
        )

    def test_update_parameters_below_fit_min(self):
        parameters = self._make_parameter_dict_for_update_parameters()
        lp = LinemeasParameters()
        lp.redshift_ref[utils.default_object_type] = 1.0
        lp.velocity[EVelocityType.Absorption][utils.default_object_type] = 200
        lp.velocity[EVelocityType.Emission][utils.default_object_type] = 400
        lp.update_parameters(parameters)
        assert utils.WarningUtils.has_warning(WarningCode.VELOCITY_FIT_RANGE)

        assert parameters.get_spectrum_model_section(utils.default_object_type)["redshiftref"] == 1.0
        assert (
            parameters.get_velocity(utils.default_object_type, "lineMeasSolve", EVelocityType.Absorption)
            == 200
        )
        assert (
            parameters.get_velocity(utils.default_object_type, "lineMeasSolve", EVelocityType.Emission) == 400
        )
        assert (
            parameters.get_velocity_fit_param(
                utils.default_object_type, "lineMeasSolve", EVelocityType.Emission, EVelocityFitParam.Min
            )
            == 400
        )

    def test_update_parameters_above_fit_min(self):
        parameters = self._make_parameter_dict_for_update_parameters()
        lp = LinemeasParameters()
        lp.redshift_ref[utils.default_object_type] = 1.0
        lp.velocity[EVelocityType.Absorption][utils.default_object_type] = 200
        lp.velocity[EVelocityType.Emission][utils.default_object_type] = 2500
        lp.update_parameters(parameters)
        assert utils.WarningUtils.has_warning(WarningCode.VELOCITY_FIT_RANGE)

        assert parameters.get_spectrum_model_section(utils.default_object_type)["redshiftref"] == 1.0
        assert (
            parameters.get_velocity(utils.default_object_type, "lineMeasSolve", EVelocityType.Absorption)
            == 200
        )
        assert (
            parameters.get_velocity(utils.default_object_type, "lineMeasSolve", EVelocityType.Emission)
            == 2500
        )
        assert (
            parameters.get_velocity_fit_param(
                utils.default_object_type, "lineMeasSolve", EVelocityType.Emission, EVelocityFitParam.Max
            )
            == 2500
        )
