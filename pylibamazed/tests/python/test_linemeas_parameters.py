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
from pylibamazed.Parameters import Parameters
from pylibamazed.LinemeasParameters import LinemeasParameters
from tests.python.test_parameters_utils import TestParametersUtils


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
        assert lp.velocity_abs[self.object_type] == 1.0
        assert lp.velocity_em[self.object_type] == 2.0

        # Source id absent from loaded dataframe
        with pytest.raises(APIException):
            lp.load_from_catalogs(
                "some absent source id", config["linemeascatalog"], config["linemeas_catalog_columns"]
            )

    def test_load_from_result_store(self, mocker):
        # mocker.patch(
        #     "pylibamazed.AbstractOutput.AbstractOutput.get_attribute_from_source", return_value="anything"
        # )
        def get_attribute_from_source(
            object_type, stage, method, dataset, attribute, rank=None, band_name=None, obs_id=None
        ):
            map = {"Redshift": 0.0, "VelocityAbsorption": 1.0, "VelocityEmission": 2.0}
            return map[attribute]

        mocker.patch(
            "pylibamazed.AbstractOutput.AbstractOutput.get_attribute_from_source",
            side_effect=get_attribute_from_source,
        )
        mocker.patch("pylibamazed.AbstractOutput.AbstractOutput.__init__", return_value=None)
        lp = LinemeasParameters()
        lp.load_from_result_store(self.generic_parameters, AbstractOutput(), self.object_type)
        assert lp.redshift_ref[self.object_type] == 0.0
        assert lp.velocity_abs[self.object_type] == 1.0
        assert lp.velocity_em[self.object_type] == 2.0

    def test_update_parameters(self):
        lp = LinemeasParameters()
        lp.redshift_ref[self.object_type] = 0.0
        lp.velocity_abs[self.object_type] = 1.0
        lp.velocity_em[self.object_type] = 2.0
        lp.update_parameters(self.generic_parameters)
        assert self.generic_parameters.get_spectrum_model_section(self.object_type)["redshiftref"] == 0.0
        assert self.generic_parameters.get_velocity_absorption(self.object_type, "lineMeasSolve") == 1.0
        assert self.generic_parameters.get_velocity_emission(self.object_type, "lineMeasSolve") == 2.0
