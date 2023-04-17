import pytest
from pylibamazed.PdfHandler import PdfHandler
from pylibamazed.PdfHandler import buildPdfParams, buildPdfHandler
from .pdf_handler_utils import PdfHandlerTestUtils

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
                "firstpass_pdf": {
                    "FirstpassPDFProbaLog": ""
                }
            }
        }
        buildPdfHandler(abstract_output, "some_object_type", True, True)

    def test_other_pass(self):
        abstract_output = PdfHandlerTestUtils.abstract_output()
        abstract_output.object_results = {
            'some_object_type': {
                "pdf_params": PdfHandlerTestUtils.pdf_params(),
                "pdf": {
                    "PDFProbaLog": ""
                }
            }
        }
        buildPdfHandler(abstract_output, "some_object_type", True)

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


