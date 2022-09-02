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
#include "RedshiftLibrary/spectrum/rebin/rebinSpline.h"
#include "RedshiftLibrary/log/log.h"

#include <gsl/gsl_fit.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

using namespace NSEpic;
using namespace std;

bool CRebinSpline::rebin(
    CSpectrumFluxAxis &rebinedFluxAxis, const TFloat64Range &range,
    const CSpectrumSpectralAxis &targetSpectralAxis, CSpectrum &rebinedSpectrum,
    CMask &rebinedMask, const std::string m_opt_error_interp,
    const TAxisSampleList &Xsrc, const TAxisSampleList &Ysrc,
    const TAxisSampleList &Xtgt, const TFloat64List &Error) {

  TAxisSampleList &Yrebin = rebinedFluxAxis.GetSamplesVector();
  TFloat64List &ErrorRebin = rebinedFluxAxis.GetError().GetSamplesVector();

  CSpectrumSpectralAxis spectralAxis = m_spectrum.GetSpectralAxis();

  // GSL method spline
  // Initialize and allocate the gsl
  // objects
  Int32 n = spectralAxis.GetSamplesCount();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, n);
  gsl_spline_init(spline, Xsrc.data(), Ysrc.data(), n);
  gsl_interp_accel *accelerator = gsl_interp_accel_alloc();

  // For each sample in the valid
  // lambda range interval.
  while (m_cursor < targetSpectralAxis.GetSamplesCount() &&
         Xtgt[m_cursor] <= range.GetEnd()) {
    Yrebin[m_cursor] = gsl_spline_eval(spline, Xtgt[m_cursor], accelerator);
    rebinedMask[m_cursor] = 1;

    // note: error rebin not
    // implemented for spline interp
    if (m_opt_error_interp != "no") {
      gsl_spline_free(spline);
      gsl_interp_accel_free(accelerator);
      return false;
    }

    m_cursor++;
  }
  gsl_spline_free(spline);
  gsl_interp_accel_free(accelerator);

  return true;
}