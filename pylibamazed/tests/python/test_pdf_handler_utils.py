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
