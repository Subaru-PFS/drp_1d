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

from pylibamazed.ParametersExtender import ParametersExtender
from tests.python.utils import (
    make_parameter_dict_at_linemodelsolve_level,
    make_parameter_dict_at_redshift_solver_level,
)


class TestParametersExtender:
    class TestContinuumReestimationExtension:
        def _make_parameter_dict(self, **kwargs):
            param_dict = make_parameter_dict_at_linemodelsolve_level(**kwargs)
            return param_dict

        def test_no_change_if_present_ok(self):
            initial_dict = self._make_parameter_dict(
                **{"fittingMethod": "hybrid", "continuumReestimation": "always"}
            )
            extended_dict = ParametersExtender().extend(initial_dict)

            assert (
                extended_dict["galaxy"]["redshiftSolver"]["lineModelSolve"]["lineModel"][
                    "continuumReestimation"
                ]
                == "always"
            )
            # Assert no changes on initial dict
            assert (
                initial_dict["galaxy"]["redshiftSolver"]["lineModelSolve"]["lineModel"][
                    "continuumReestimation"
                ]
                == "always"
            )

        def test_overwrites_if_present_but_should_not(self):
            initial_dict = self._make_parameter_dict(
                **{"fittingMethod": "sth", "continuumReestimation": "always"}
            )
            extended_dict = ParametersExtender().extend(initial_dict)

            assert (
                extended_dict["galaxy"]["redshiftSolver"]["lineModelSolve"]["lineModel"][
                    "continuumReestimation"
                ]
                == "no"
            )
            # Assert no change on initial dict
            assert (
                initial_dict["galaxy"]["redshiftSolver"]["lineModelSolve"]["lineModel"][
                    "continuumReestimation"
                ]
                == "always"
            )

        def test_adds_if_absent(self):
            initial_dict = self._make_parameter_dict(**{"fittingMethod": "sth"})
            extended_dict = ParametersExtender().extend(initial_dict)

            assert (
                extended_dict["galaxy"]["redshiftSolver"]["lineModelSolve"]["lineModel"][
                    "continuumReestimation"
                ]
                == "no"
            )
            # Assert no change on initial dict
            assert (
                initial_dict["galaxy"]["redshiftSolver"]["lineModelSolve"]["lineModel"].get(
                    "continuumReestimation"
                )
                is None
            )

    class TestUseLogLambdaSamplingExtension:
        def _make_parameter_dict(self, **kwargs):
            param_dict = make_parameter_dict_at_linemodelsolve_level(**kwargs)
            return param_dict

        def test_no_change_if_present_ok(self):
            initial_dict = self._make_parameter_dict(
                **{
                    "continuumComponent": "tplFit",
                    "continuumFit": {"fftProcessing": True},
                    "useLogLambdaSampling": True,
                }
            )
            extended_dict = ParametersExtender().extend(initial_dict)

            assert extended_dict["galaxy"]["redshiftSolver"]["lineModelSolve"]["lineModel"][
                "useLogLambdaSampling"
            ]
            assert initial_dict["galaxy"]["redshiftSolver"]["lineModelSolve"]["lineModel"][
                "useLogLambdaSampling"
            ]

        def test_overwrites_if_present_but_should_not(self):
            initial_dict = self._make_parameter_dict(
                **{
                    "continuumComponent": "tplFit",
                    "continuumFit": {"fftProcessing": False},
                    "useLogLambdaSampling": True,
                }
            )
            extended_dict = ParametersExtender().extend(initial_dict)

            assert not extended_dict["galaxy"]["redshiftSolver"]["lineModelSolve"]["lineModel"][
                "useLogLambdaSampling"
            ]

            # Assert no change on initial dict
            assert initial_dict["galaxy"]["redshiftSolver"]["lineModelSolve"]["lineModel"][
                "useLogLambdaSampling"
            ]

        def test_adds_if_absent(self):
            initial_dict = self._make_parameter_dict(
                **{
                    "continuumComponent": "tplFit",
                    "continuumFit": {"fftProcessing": False},
                }
            )
            extended_dict = ParametersExtender().extend(initial_dict)

            assert not extended_dict["galaxy"]["redshiftSolver"]["lineModelSolve"]["lineModel"][
                "useLogLambdaSampling"
            ]

            # Assert no change on initial dict
            assert (
                initial_dict["galaxy"]["redshiftSolver"]["lineModelSolve"]["lineModel"].get(
                    "useLogLambdaSampling"
                )
                is None
            )

    class TestUseMedianKernelExtension:
        def _make_parameter_dict(self, **kwargs):
            param_dict = make_parameter_dict_at_linemodelsolve_level(**{"continuumComponent": "fromSpectrum"})
            return param_dict | kwargs

        def test_no_change_if_present_ok(self):
            initial_dict = self._make_parameter_dict(
                **{
                    "continuumRemoval": {
                        "method": "irregularSamplingMedian",
                        "medianKernelWidth": 1,
                        "medianEvenReflection": 2,
                    }
                }
            )
            extended_dict = ParametersExtender().extend(initial_dict)

            assert extended_dict["continuumRemoval"]["medianKernelWidth"] == 1
            assert extended_dict["continuumRemoval"]["medianEvenReflection"] == 2

            assert initial_dict["continuumRemoval"]["medianKernelWidth"] == 1
            assert initial_dict["continuumRemoval"]["medianEvenReflection"] == 2

        def test_overwrites_if_present_but_should_not(self):
            initial_dict = self._make_parameter_dict(
                **{
                    "continuumRemoval": {
                        "method": "sth",
                        "medianKernelWidth": 1,
                        "medianEvenReflection": False,
                    }
                }
            )
            extended_dict = ParametersExtender().extend(initial_dict)

            assert extended_dict["continuumRemoval"].get("method") == "sth"
            assert extended_dict["continuumRemoval"]["medianKernelWidth"] == -1
            assert extended_dict["continuumRemoval"]["medianEvenReflection"] is True

            assert initial_dict["continuumRemoval"]["method"] == "sth"
            assert initial_dict["continuumRemoval"]["medianKernelWidth"] == 1
            assert initial_dict["continuumRemoval"]["medianEvenReflection"] is False

        def test_adds_if_absent(self):
            initial_dict = self._make_parameter_dict(**{"continuumRemoval": {"method": "sth"}})
            extended_dict = ParametersExtender().extend(initial_dict)

            assert extended_dict["continuumRemoval"].get("method") == "sth"
            assert extended_dict["continuumRemoval"]["medianKernelWidth"] == -1
            assert extended_dict["continuumRemoval"]["medianEvenReflection"] is True

            assert initial_dict["continuumRemoval"].get("method") == "sth"
            assert initial_dict["continuumRemoval"].get("medianKernelWidth") is None
            assert initial_dict["continuumRemoval"].get("medianEvenReflection") is None

    class TestTemplateCatalogContinuumRemovalExtension:
        def _make_parameter_dict(self, **kwargs):
            param_dict = make_parameter_dict_at_redshift_solver_level(
                **{"templateFittingSolve": {"spectrum": {"component": "sth"}}}
            )
            return param_dict | kwargs

        def test_no_change_if_present_ok(self):
            initial_dict = self._make_parameter_dict(
                **{"templateCatalog": {"continuumRemoval": {"some key": "some value"}}}
            )

            print("initial dict", initial_dict)
            extended_dict = ParametersExtender().extend(initial_dict)
            assert extended_dict["templateCatalog"]["continuumRemoval"].get("some key") == "some value"

        def test_overwrites_if_present_but_should_not(self):
            initial_dict = {"templateCatalog": {"continuumRemoval": {"some key": "some value"}}}
            extended_dict = ParametersExtender().extend(initial_dict)
            assert extended_dict["templateCatalog"]["continuumRemoval"].get("some key") is None

        def test_adds_if_absent(self):
            initial_dict = {"templateCatalog": {}}
            extended_dict = ParametersExtender().extend(initial_dict)
            assert extended_dict["templateCatalog"].get("continuumRemoval") is not None
