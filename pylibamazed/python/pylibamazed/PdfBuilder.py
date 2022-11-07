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
import json
from pylibamazed.redshift import CLogZPdfResult_getLogZPdf_fine, CLogZPdfResult_buildLogMixedZPdfGrid
from pylibamazed.redshift import ZGridParameters, PC_Get_Float64Array
import numpy as np

def _extract_pdf_params(pdf_params, first_pass = False):
    v = [dict(zip(pdf_params, t)) for t in zip(*pdf_params.values())]
    if first_pass:
        return [ZGridParameters(p["FPZmin"], p["FPZmax"], p["FPZstep"], np.nan) for p in v]
    return [ZGridParameters(p["zmin"], p["zmax"], p["zstep"], p["zcenter"]) for p in v]


class PdfBuilder:
    def __init__(self, abstract_output):
        self.output = abstract_output

    def get_fine_zgrid_pdf(self, object_type, logsampling, first_pass=False):
        pdf_params = self.output.get_dataset(object_type, "pdf_params")
        pdf_proba = list(self.output.get_dataset(object_type, "pdf")["PDFProbaLog"])
        ret = CLogZPdfResult_getLogZPdf_fine(
            logsampling, _extract_pdf_params(pdf_params), pdf_proba
        )
        return {
            "zgrid": PC_Get_Float64Array(ret.zgrid),
            "probaLog": PC_Get_Float64Array(ret.probaLog),
        }

    def get_mixed_zgrid_pdf(self, object_type, logsampling, first_pass=False):
        if first_pass:
            pdf_params = self.output.get_dataset(object_type,"firstpass_pdf_params")
        else:
            pdf_params = self.output.get_dataset(object_type,"pdf_params")
        ret = CLogZPdfResult_buildLogMixedZPdfGrid(
            logsampling,_extract_pdf_params(pdf_params,first_pass)
        )
        return np.array(ret)

    def get_coarse_zgrid_pdf(self, object_type):
        pass
