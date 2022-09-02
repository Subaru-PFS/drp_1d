// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
//
// https://www.lam.fr/
//
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
//
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use,
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info".
//
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability.
//
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or
// data to be ensured and,  more generally, to use and operate it in the
// same conditions as regards security.
//
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#include "RedshiftLibrary/spectrum/rebin/rebinNgp.h"
#include "RedshiftLibrary/common/indexing.h"
#include "RedshiftLibrary/log/log.h"

using namespace NSEpic;
using namespace std;

bool CRebinNgp::rebin(CSpectrumFluxAxis &rebinedFluxAxis,
                      const TFloat64Range &range,
                      const CSpectrumSpectralAxis &targetSpectralAxis,
                      CSpectrum &rebinedSpectrum, CMask &rebinedMask,
                      const std::string m_opt_error_interp,
                      const TAxisSampleList &Xsrc, const TAxisSampleList &Ysrc,
                      const TAxisSampleList &Xtgt, const TFloat64List &Error) {

  TAxisSampleList &Yrebin = rebinedFluxAxis.GetSamplesVector();
  TFloat64List &ErrorRebin = rebinedFluxAxis.GetError().GetSamplesVector();

  CSpectrumSpectralAxis spectralAxis = m_spectrum.GetSpectralAxis();

  // nearest sample, lookup
  Int32 k = 0;
  Int32 n = spectralAxis.GetSamplesCount();
  while (m_cursor < targetSpectralAxis.GetSamplesCount() &&
         Xtgt[m_cursor] <= range.GetEnd()) {
    // k = gsl_interp_bsearch
    // (Xsrc.data(), Xtgt[j], kprev,
    // n);
    k = CIndexing<Float64>::getCloserIndex(Xsrc, Xtgt[m_cursor]);
    Float64 xSrcStep = NAN;
    if (k == Xsrc.size() - 1)
      xSrcStep = Xsrc[k] - Xsrc[k - 1];
    else
      xSrcStep = Xsrc[k + 1] - Xsrc[k];

    // closest value
    Yrebin[m_cursor] = Ysrc[k];

    if (m_opt_error_interp != "no") {
      ErrorRebin[m_cursor] = Error[k];

      if (m_opt_error_interp == "rebinVariance") {
        Float64 xDestStep = NAN;
        Float64 xStepCompensation = 1.;
        if (m_cursor < targetSpectralAxis.GetSamplesCount() - 1) {
          xDestStep = Xtgt[m_cursor + 1] - Xtgt[m_cursor];
          xStepCompensation = xSrcStep / xDestStep;
        } else if (m_cursor > 0) {
          xDestStep = Xtgt[m_cursor] - Xtgt[m_cursor - 1];
          xStepCompensation = xSrcStep / xDestStep;
        }
        ErrorRebin[m_cursor] *= sqrt(xStepCompensation);
      }
    }

    rebinedMask[m_cursor] = 1;
    m_cursor++;
  }

  return true;
}