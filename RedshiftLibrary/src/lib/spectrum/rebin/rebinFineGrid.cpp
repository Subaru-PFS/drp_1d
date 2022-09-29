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
#include "RedshiftLibrary/spectrum/rebin/rebinFineGrid.h"
#include "RedshiftLibrary/log/log.h"
#include <gsl/gsl_fit.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

using namespace NSEpic;
using namespace std;

void CRebinFineGrid::rebinFineGrid() {
  // Precalculate a fine grid template to be used for the 'closest value' rebin
  // method
  Int32 n = m_spectrum.GetSampleCount();
  if (!n)
    THROWG(INVALID_SPECTRUM, "Invalid spectrum : spectral axis is empty");

  CSpectrumSpectralAxis spectralAxis = m_spectrum.GetSpectralAxis();

  Float64 lmin = spectralAxis[0]; // template wavelength never starts at 0
  Float64 lmax = spectralAxis[n - 1];
  Int32 nTgt = (lmax - lmin) / m_dLambdaFineGrid + 2.0 / m_dLambdaFineGrid;

  m_pfgFlux.resize(nTgt);

  const TAxisSampleList &Ysrc = m_spectrum.GetFluxAxis().GetSamplesVector();
  const TAxisSampleList &Xsrc = spectralAxis.GetSamplesVector();

  // Initialize and allocate the gsl objects
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, n);
  gsl_spline_init(spline, Xsrc.data(), Ysrc.data(), n);
  gsl_interp_accel *accelerator = gsl_interp_accel_alloc();

  Int32 k = 0;
  Float64 x = 0.0;
  for (k = 0; k < nTgt; k++) {
    x = lmin + k * m_dLambdaFineGrid;
    if (x < spectralAxis[0] || x > spectralAxis[n - 1]) {
      m_pfgFlux[k] = 0.0;
    } else {
      m_pfgFlux[k] = gsl_spline_eval(spline, x, accelerator);
    }
  }
  gsl_spline_free(spline);
  gsl_interp_accel_free(accelerator);
  m_FineGridInterpolated = true;
}

void CRebinFineGrid::clearFineGrid() {
  m_FineGridInterpolated = false;
  m_pfgFlux.clear();
}

void CRebinFineGrid::rebin(
    CSpectrumFluxAxis &rebinedFluxAxis, const TFloat64Range &range,
    const CSpectrumSpectralAxis &targetSpectralAxis, CSpectrum &rebinedSpectrum,
    CMask &rebinedMask, const std::string m_opt_error_interp,
    const TAxisSampleList &Xsrc, const TAxisSampleList &Ysrc,
    const TAxisSampleList &Xtgt, const TFloat64List &Error, Int32 &cursor) {

  if (m_FineGridInterpolated == false)
    rebinFineGrid();

  if (m_pfgFlux.size() == 0 && m_FineGridInterpolated == true) {
    THROWG(INTERNAL_ERROR,
           "Problem finegrid interpolatted buffer couldnt be computed");
  }

  TAxisSampleList &Yrebin = rebinedFluxAxis.GetSamplesVector();
  CSpectrumSpectralAxis spectralAxis = m_spectrum.GetSpectralAxis();

  // Precomputed FINE GRID nearest
  // sample, 20150801
  Int32 k = 0;
  Float64 lmin = spectralAxis[0];
  // For each sample in the target
  // spectrum
  while (cursor < targetSpectralAxis.GetSamplesCount() &&
         Xtgt[cursor] <= range.GetEnd()) {
    k = int((Xtgt[cursor] - lmin) / m_dLambdaFineGrid + 0.5);
    Yrebin[cursor] = m_pfgFlux[k];
    rebinedMask[cursor] = 1;

    // note: error rebin not
    // implemented for
    // precomputedfinegrid
    if (m_opt_error_interp != "no")
      THROWG(INTERNAL_ERROR,
             "noise rebining not implemented for precomputedfinegrid");

    cursor++;
  }
}