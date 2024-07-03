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
from tests.python.test_parameters_utils import TestParametersUtils


class TestParameters:
    generic_parameters: Parameters = Parameters(
        TestParametersUtils().make_parameters_dict(), make_checks=False
    )

    def test_load_linemeas_parameters_from_catalog(self, mocker):
        source_id = "some_source_id"
        mocker.patch(
            "pandas.read_csv",
            return_value=pd.DataFrame(
                {
                    "ProcessingID": [source_id],
                    "redshift_column_name": [0],
                    "velocity_absorption_column_name": [0],
                    "velocity_emission_column_name": [0],
                }
            ),
        )
        config = TestParametersUtils().make_config()
        self.generic_parameters.load_linemeas_parameters_from_catalog(source_id, config)

        # Source id absent from loaded dataframe
        with pytest.raises(APIException):
            self.generic_parameters.load_linemeas_parameters_from_catalog("some absent source id", config)

    def test_load_linemeas_parameters_from_result_store(self, mocker):
        mocker.patch(
            "pylibamazed.AbstractOutput.AbstractOutput.get_attribute_from_source", return_value="anything"
        )
        mocker.patch("pylibamazed.AbstractOutput.AbstractOutput.__init__", return_value=None)
        self.generic_parameters.load_linemeas_parameters_from_result_store(
            AbstractOutput(),
            TestParametersUtils.default_object_type,
        )

    def test_to_json(self):
        self.generic_parameters.to_json()

    def test_lineratio_catalog_enabled(self):
        self.generic_parameters.is_tplratio_catalog_needed(TestParametersUtils.default_object_type)

        # With other solve method
        parameters = Parameters(
            TestParametersUtils().make_parameters_dict(**{"method": "other"}), make_checks=False
        )
        parameters.is_tplratio_catalog_needed(TestParametersUtils.default_object_type)

    def test_stage_enabled(self):
        for stage in [
            "redshiftSolver",
            "lineMeasSolver",
            "linemeas_catalog_load",
            "reliabilitySolver",
            "subClassifSolver",
        ]:
            self.generic_parameters.stage_enabled(TestParametersUtils.default_object_type, stage)

        with pytest.raises(Exception):
            self.generic_parameters.stage_enabled(TestParametersUtils.default_object_type, "unkown stage")

    class TestIsARedshiftSolverUsed:
        def test_false_if_no_object_has_a_redshift_solver_defined(self):
            parameters: Parameters = Parameters(
                {"galaxy": {"method": None}, "spectrumModels": ["galaxy"]}, make_checks=False
            )
            assert parameters.is_a_redshift_solver_used() is False

        def test_true_if_an_object_has_a_redshift_solver(self):
            parameters: Parameters = Parameters(
                {
                    "version": 2,
                    "galaxy": {"method": None},
                    "star": {"stages": ["redshiftSolver"], "redshiftSolver": {"method": "lineModelsolve"}},
                    "spectrumModels": ["galaxy", "star"],
                },
                make_checks=False,
            )
            assert parameters.is_a_redshift_solver_used() is True
