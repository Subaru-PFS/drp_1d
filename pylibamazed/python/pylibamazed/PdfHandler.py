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
import numpy as np
from pylibamazed.redshift import (CLogZPdfResult, CZGridListParams,
                                  CZGridParam, TFloat64Range)

from pylibamazed import AbstractOutput


def buildPdfParams(pdf_params, first_pass=False):
    v = [dict(zip(pdf_params, t)) for t in zip(*pdf_params.values())]
    if first_pass:
        return CZGridListParams(
            [CZGridParam(TFloat64Range(p["FPZmin"], p["FPZmax"]), p["FPZstep"], np.nan) for p in v]
        )
    return CZGridListParams(
        [CZGridParam(TFloat64Range(p["zmin"], p["zmax"]), p["zstep"], p["zcenter"]) for p in v]
    )


class BuilderPdfHandler:
    def add_params(
        self,
        abstract_output: AbstractOutput,
        object_type,
        logsampling,
        first_pass=False
    ):
        self.abstract_output = abstract_output
        self.object_type = object_type
        self.logsampling = logsampling
        self.first_pass = first_pass
        return self

    def build(self):
        dataset_prefix = ""
        name_prefix = ""

        if self.first_pass:
            dataset_prefix = "firstpass_"
            name_prefix = "Firstpass"

        pdf_params = self.abstract_output.get_dataset(
            self.object_type, dataset_prefix + "pdf_params"
        )
        c_pdf_params = buildPdfParams(pdf_params, self.first_pass)
        pdf_proba = self.abstract_output.get_dataset(
            self.object_type, dataset_prefix + "pdf"
        )[name_prefix + "PDFProbaLog"]

        return PdfHandler(c_pdf_params, self.logsampling, pdf_proba)


class PdfHandler:
    def __init__(self, pdf_params, logsampling: bool, pdf_proba):
        self._pdf = CLogZPdfResult(pdf_params, logsampling, pdf_proba)

    @property
    def redshifts(self):
        return self._pdf.redshifts.to_numpy()

    @property
    def valProbaLog(self):
        return self._pdf.valProbaLog.to_numpy()

    def isRegular(self):
        return self._pdf.isRegular()

    def convertToRegular(self, fine=True, zgrid_max=None):
        self._pdf.convertToRegular(fine)

        if zgrid_max is not None and self._pdf.zmax[0] < zgrid_max:
            self._pdf.extrapolate_on_right_border(zgrid_max)

    def isPdfValid(self):
        return self._pdf.isPdfValid()

    def getSumTrapez(self):
        return self._pdf.getSumTrapez()
