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
#include "RedshiftLibrary/spectrum/rebin/rebinLinearFull.h"

#include "RedshiftLibrary/log/log.h"

using namespace NSEpic;
using namespace std;

void CRebinLinearFull::rebin(CSpectrumFluxAxis &rebinedFluxAxis,
                             const TFloat64Range &range,
                             const CSpectrumSpectralAxis &targetSpectralAxis,
                             CMask &rebinedMask,
                             const std::string opt_error_interp,
                             const TAxisSampleList &Xtgt,
                             TFloat64List &error_tmp, Int32 &cursor) {
  bool const handle_error = handleError(opt_error_interp);

  const CFullSpectrum &origin = dynamic_cast<const CFullSpectrum &>(m_spectrum);
  Int32 n = origin.GetSampleCount();
  const TAxisSampleList &Xsrc = origin.GetSpectralAxis().GetSamplesVector();
  const TAxisSampleList &Ysrc = origin.GetFluxAxis().GetSamplesVector();
  const TFloat64List &Error =
      handle_error ? origin.GetErrorAxis().GetSamplesVector() : TFloat64List{};

  Int32 k = 0;
  // For each sample in the valid lambda range interval.
  while (k < n - 1 && Xsrc[k] <= range.GetEnd()) {
    // For each sample in the target spectrum that are in between two
    // continous source sample
    while (cursor < targetSpectralAxis.GetSamplesCount() &&
           Xtgt[cursor] <= Xsrc[k + 1]) {
      // perform linear interpolation of the flux
      Float64 xSrcStep = (Xsrc[k + 1] - Xsrc[k]);
      Float64 t = (Xtgt[cursor] - Xsrc[k]) / xSrcStep;

      // If both samples arround subsample are valid, compute flux, otherwise,
      // set mask & flux to 0
      if (origin.getMask()[k] && origin.getMask()[k + 1]) {
        rebinedFluxAxis[cursor] = Ysrc[k] + (Ysrc[k + 1] - Ysrc[k]) * t;
        rebinedMask[cursor] = 1;
      } else {
        rebinedMask[cursor] = 0;
        rebinedFluxAxis[cursor] = 0;
      }

      // Same for error
      if (handle_error) {
        if (opt_error_interp == "rebin" && origin.getMask()[k] &&
            origin.getMask()[k + 1])
          error_tmp[cursor] = Error[k] + (Error[k + 1] - Error[k]) * t;
        else if (opt_error_interp == "rebinVariance" && origin.getMask()[k] &&
                 origin.getMask()[k + 1]) {
          error_tmp[cursor] = sqrt(Error[k] * Error[k] * (1 - t) * (1 - t) +
                                   Error[k + 1] * Error[k + 1] * t * t);
          Float64 xStepCompensation = computeXStepCompensation(
              targetSpectralAxis, Xtgt, cursor, xSrcStep);
          error_tmp[cursor] = error_tmp[cursor] * sqrt(xStepCompensation);
        } else {
          Log.LogDetail(Formatter()
                        << "set error to dbl_min at " << cursor
                        << " mask before=" << (Int32)origin.getMask()[k]
                        << " after=" << (Int32)origin.getMask()[k + 1] << " at "
                        << targetSpectralAxis[cursor] << " between" << Xsrc[k]
                        << " and " << Xsrc[k + 1]);
          error_tmp[cursor] = DBL_MAX;
        }
      }
      cursor++;
    }

    k++;
  }
}
