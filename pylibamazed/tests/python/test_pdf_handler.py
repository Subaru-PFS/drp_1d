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

import pytest
from pylibamazed.PdfHandler import BuilderPdfHandler, buildPdfParams, get_final_regular_z_grid
from pylibamazed.Parameters import Parameters
from tests.python.test_pdf_handler_utils import PdfHandlerTestUtils
from tests.python.utils import make_parameter_dict_at_linemodelsolve_level
from tests.python.utils import (
    WarningUtils,
    check_from_parameter_dict,
    make_parameter_dict_at_redshift_solver_level,
)


class TestBuildPdfParams:
    def test_first_pass(self):
        buildPdfParams(PdfHandlerTestUtils.pdf_params())

    def test_other_pass(self):
        buildPdfParams(PdfHandlerTestUtils.pdf_params(), True)


class TestBuildPdfhandler:
    def test_first_pass(self):
        abstract_output = PdfHandlerTestUtils.abstract_output()
        abstract_output.object_results = {
            "some_object_type": {
                "firstpass_pdf_params": PdfHandlerTestUtils.pdf_params(),
                "firstpass_pdf": {"FirstpassLogZPdfNative": ""},
            }
        }
        BuilderPdfHandler().add_params(abstract_output, "some_object_type", True, True).build()

    def test_other_pass(self):
        abstract_output = PdfHandlerTestUtils.abstract_output()
        abstract_output.object_results = {
            "some_object_type": {"pdf_params": PdfHandlerTestUtils.pdf_params(), "pdf": {"LogZPdfNative": ""}}
        }
        BuilderPdfHandler().add_params(abstract_output, "some_object_type", True).build()


class TestPdfHandlerClass:
    pdf_handler = PdfHandlerTestUtils.pdf_handler()

    def test_redshifts_property(self):
        self.pdf_handler.redshifts

    def test_valPropaLog_property(self):
        self.pdf_handler.valProbaLog

    def test_isRegular(self):
        self.pdf_handler.isRegular()

    def test_convertToRegular(self):
        self.pdf_handler.convertToRegular()
        self.pdf_handler.convertToRegular(True, 10)

    def test_isPdfValid(self):
        with pytest.raises(Exception):
            self.pdf_handler.isPdfValid()

    def test_getSumTrapez(self):
        self.pdf_handler.getSumTrapez()


class TestGetFinalRegularZgrid:
    def _make_parameter_dict(self, **kwargs) -> dict:
        kwargs["method"] = "lineModelSolve"
        param_dict = make_parameter_dict_at_redshift_solver_level(**kwargs)
        return param_dict

    def test(self):
        param_dict = self._make_parameter_dict(
            object_level_params={"redshiftRange": [0, 1], "redshiftStep": 0.1, "redshiftSampling": "log"},
            **{
                "lineModelSolve": {
                    "lineModel": {
                        "firstPass": {"largeGridStepRatio": 5},
                        "skipSecondPass": False,
                    }
                }
            },
        )
        get_final_regular_z_grid("galaxy", Parameters(param_dict, make_checks=False))
