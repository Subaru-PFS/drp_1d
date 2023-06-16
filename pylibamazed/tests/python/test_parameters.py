import pandas as pd
import pytest
from pylibamazed.AbstractOutput import AbstractOutput
from pylibamazed.Exception import APIException
from pylibamazed.Parameters import Parameters
from .test_parameters_utils import TestParametersUtils
        
class TestParameters:
    generic_parameters: Parameters = TestParametersUtils().make_parameters();

    def test_get_solve_methods(self):
        self.generic_parameters.get_solve_methods(TestParametersUtils.default_object_type)
    
    def test_get_redshift_sampling(self):
        self.generic_parameters.get_redshift_sampling(TestParametersUtils.default_object_type)
    
    def test_get_linemodel_methods(self):
        self.generic_parameters.get_linemodel_methods(TestParametersUtils.default_object_type)

    def test_check_lmskipsecondpass(self):
        # With Linemodel solve method
        self.generic_parameters.check_lmskipsecondpass(TestParametersUtils.default_object_type)

        # With other solve method
        parameters = TestParametersUtils().make_parameters(**{"method": "other"})
        parameters.check_lmskipsecondpass(TestParametersUtils.default_object_type)

        # With None solve method
        parameters = TestParametersUtils().make_parameters(**{"method": None})
        parameters.check_lmskipsecondpass(TestParametersUtils.default_object_type)

    def test_get_solve_method(self):
        self.generic_parameters.get_solve_method(TestParametersUtils.default_object_type)

    def test_get_linemeas_method(self):
        self.generic_parameters.get_linemeas_method(TestParametersUtils.default_object_type)

    def test_get_objects(self):
        self.generic_parameters.get_objects()

    def test_load_linemeas_parameters_from_catalog(self, mocker):
        source_id = 'some_source_id'
        mocker.patch('pandas.read_csv', return_value= pd.DataFrame({
            "ProcessingID": [source_id],
            "redshift_column_name": [0],
            "velocity_absorption_column_name": [0],
            "velocity_emission_column_name": [0],
        }))
        self.generic_parameters.load_linemeas_parameters_from_catalog(source_id)

        # Source id absent from loaded dataframe
        with pytest.raises(APIException):
            self.generic_parameters.load_linemeas_parameters_from_catalog("some absent source id")

    def test_load_linemeas_parameters_from_result_store(self, mocker):
        mocker.patch(
            'pylibamazed.AbstractOutput.AbstractOutput.get_attribute_from_source',
            return_value="anything"
        )
        mocker.patch(
            'pylibamazed.AbstractOutput.AbstractOutput.__init__',
            return_value=None
        )
        self.generic_parameters.load_linemeas_parameters_from_result_store(
            AbstractOutput(),
            TestParametersUtils.default_object_type,
        )
        
    def test_get_json(self):
        self.generic_parameters.get_json()

    def test_reliability_enabled(self):
        self.generic_parameters.reliability_enabled(TestParametersUtils.default_object_type)

    def test_lineratio_catalog_enabled(self):
        self.generic_parameters.lineratio_catalog_enabled(TestParametersUtils.default_object_type)

        # With other solve method
        parameters = TestParametersUtils().make_parameters(**{"method": "other"})
        parameters.lineratio_catalog_enabled(TestParametersUtils.default_object_type)

    def test_stage_enabled(self):
        for stage in [
            "redshift_solver",
            "linemeas_solver",
            "linemeas_catalog_load",
            "reliability_solver",
            "sub_classif_solver"
        ]:
            self.generic_parameters.stage_enabled(
                TestParametersUtils.default_object_type,
                stage
            )
        
        with pytest.raises(Exception):
            self.generic_parameters.stage_enabled(
                TestParametersUtils.default_object_type,
                "unkown stage"
            )
        