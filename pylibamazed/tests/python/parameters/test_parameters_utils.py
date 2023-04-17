from pylibamazed.Parameters import Parameters


class TestParametersUtils:

    default_object_type = "object_type"

    def make_input_parameters(self, **kwargs) -> dict:
        input_parameters = {
            self.default_object_type: {
                "method": kwargs.get("method", "LineModelSolve"),
                "linemeas_method": "some lineas method",
                "redshiftsampling": "some redshift sampling",
                "LineModelSolve": {
                    "linemodel": {
                        "skipsecondpass": True,
                        "lineRatioType": "tplratio",
                    }
                },
                "LineMeasSolve": {
                    "linemodel": {
                        "velocityabsorption": None,
                        "velocityemission": None,
                    }
                }
            },
            "objects": [self.default_object_type]
        }
        return input_parameters

    def make_config(self) -> dict:
        config = {
            "linemeascatalog": {
                self.default_object_type: "somefile.csv"
            },
            "linemeas_catalog_columns": {
                self.default_object_type: {
                    "Redshift": "redshift_column_name",
                    "VelocityAbsorption": "velocity_absorption_column_name",
                    "VelocityEmission": "velocity_emission_column_name",
                },
            },
        }
        return config

    def make_parameters(self, config=None, **kwargs) -> Parameters:
        return Parameters(self.make_input_parameters(**kwargs), self.make_config())
    