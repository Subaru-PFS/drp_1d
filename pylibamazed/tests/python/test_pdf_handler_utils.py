from pylibamazed.AbstractOutput import AbstractOutput
from pylibamazed.Parameters import Parameters
from pylibamazed.PdfHandler import buildPdfHandler

class PdfHandlerTestUtils:
    @staticmethod
    def pdf_params():
        return {
            "FPZmin": [0],
            "FPZmax": [100],
            "FPZstep": [2],
            "zmin": [0],
            "zmax": [1],
            "zstep": [0.1],
            "zcenter": [0.5]
        }

    @staticmethod
    def parameters():
        return Parameters({"objects": []})
    
    @staticmethod
    def abstract_output():
        return AbstractOutput(PdfHandlerTestUtils.parameters())
    
    @staticmethod
    def pdf_handler():
        abstract_output = PdfHandlerTestUtils.abstract_output()
        abstract_output.object_results = {
            'some_object_type': {
                "pdf_params": PdfHandlerTestUtils.pdf_params(),
                "pdf": {
                    "PDFProbaLog": ""
                }
            }
        }
        pdf_handler = buildPdfHandler(abstract_output, "some_object_type", True)
        return pdf_handler
