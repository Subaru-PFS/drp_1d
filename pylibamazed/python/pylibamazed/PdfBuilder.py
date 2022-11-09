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
from pylibamazed.redshift import CLogZPdfResult_getLogZPdf_fine
from pylibamazed.redshift import TZGridParameters
from pylibamazed.redshift import CZGridListParams
from pylibamazed.redshift import TFloat64Range
import numpy as np

def _extract_pdf_params(pdf_params, first_pass = False):
    v = [dict(zip(pdf_params, t)) for t in zip(*pdf_params.values())]
    if first_pass:
        return CZGridListParams([TZGridParameters(TFloat64Range(p["FPZmin"], p["FPZmax"]), p["FPZstep"], np.nan) for p in v])
    return CZGridListParams([TZGridParameters(TFloat64Range(p["zmin"], p["zmax"]), p["zstep"], p["zcenter"]) for p in v])


class PdfBuilder:
    def __init__(self, abstract_output):
        self.output = abstract_output

    def interpolate_pdf_on_regular_grid(self, object_type, logsampling, c_zgrid_max=np.nan):
        pdf_params = self.output.get_dataset(object_type, "pdf_params")
        #tmp code to fix reliability, the aim is to create 
        extrapolate = False
        if not np.isnan(c_zgrid_max) and pdf_params['zmax'][0] != c_zgrid_max:
            extrapolate = True
            pdf_params['zmax'][0] = c_zgrid_max #overwrite read value
          
        pdf_proba = list(self.output.get_dataset(object_type, "pdf")["PDFProbaLog"])
        ret = CLogZPdfResult_getLogZPdf_fine(
            logsampling, _extract_pdf_params(pdf_params), pdf_proba
        )

        probalog = np.array(ret.probaLog)
        if extrapolate :
            probalog =  self.extrapolate_pdf_onborder(probalog)
        return {
            "zgrid": np.array(ret.zgrid),
            "probaLog": probalog
        }

    #temporary function, meant to disappear
    #search for nan at the end of probalog and replace with the last value
    #interpolation function adds up nan probalog at the end of the array
    def extrapolate_pdf_onborder(self, probalog):
        nan_indices = np.where(np.isnan(probalog))[0]
        ref_value = probalog[nan_indices[0]-1]#last valid value
        probalog = np.nan_to_num(probalog, nan=ref_value)
        return probalog

    def get_zgrid(self, object_type, logsampling, first_pass=False):
        if first_pass:
            pdf_params = self.output.get_dataset(object_type,"firstpass_pdf_params")
        else:
            pdf_params = self.output.get_dataset(object_type,"pdf_params")
        ret = _extract_pdf_params(pdf_params,first_pass).buildLogMixedZGrid(
            logsampling,
        )
        return np.array(ret)

    def get_coarse_zgrid_pdf(self, object_type):
        pass
