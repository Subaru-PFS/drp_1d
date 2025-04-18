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
#include "RedshiftLibrary/spectrum/rebin/rebin.h"
#include "RedshiftLibrary/spectrum/rebin/rebinFineGrid.h"
#include "RedshiftLibrary/spectrum/rebin/rebinLinear.h"
#include "RedshiftLibrary/spectrum/rebin/rebinNgp.h"
#include "RedshiftLibrary/spectrum/rebin/rebinSpline.h"

using namespace NSEpic;
using namespace std;

void CRebin::compute(const TFloat64Range &range,
                     const CSpectrumSpectralAxis &targetSpectralAxis,
                     CSpectrum &rebinedSpectrum, CMask &rebinedMask,
                     const std::string opt_error_interp) {

  Int32 s = targetSpectralAxis.GetSamplesCount();
  const CSpectrumSpectralAxis &spectralAxis = m_spectrum.GetSpectralAxis();

  // find start/end indexs for both axes
  if (spectralAxis[0] > range.GetBegin() ||
      spectralAxis[spectralAxis.GetSamplesCount() - 1] < range.GetEnd()) {
    THROWG(ErrorCode::INTERNAL_ERROR, "input spectral range is not "
                                      "included in spectral axis");
  }

  CSpectrumFluxAxis rebinedFluxAxis(s);
  rebinedMask.SetSize(s);

  const TAxisSampleList &Xtgt = targetSpectralAxis.GetSamplesVector();
  TFloat64List error_tmp = rebinedFluxAxis.GetError().GetSamplesVector();

  // Move cursors up to lambda range start
  Int32 cursor = 0;
  while (cursor < targetSpectralAxis.GetSamplesCount() &&
         Xtgt[cursor] < range.GetBegin()) {
    rebinedMask[cursor] = 0;
    rebinedFluxAxis[cursor] = 0.0;
    if (opt_error_interp == "rebin" || opt_error_interp == "rebinVariance")
      error_tmp[cursor] = INFINITY;
    cursor++;
  }

  rebin(rebinedFluxAxis, range, targetSpectralAxis, rebinedMask,
        opt_error_interp, Xtgt, error_tmp, cursor);

  while (cursor < targetSpectralAxis.GetSamplesCount()) {
    rebinedMask[cursor] = 0;
    rebinedFluxAxis[cursor] = 0.0;
    if (opt_error_interp == "rebin" || opt_error_interp == "rebinVariance")
      error_tmp[cursor] = INFINITY;
    cursor++;
  }

  rebinedFluxAxis.setError(CSpectrumNoiseAxis(error_tmp));
  rebinedSpectrum.ResetContinuum();
  rebinedSpectrum.SetType(CSpectrum::EType::raw);
  rebinedSpectrum.SetSpectralAndFluxAxes(targetSpectralAxis,
                                         std::move(rebinedFluxAxis));
}

Float64 CRebin::computeXStepCompensation(
    const CSpectrumSpectralAxis &targetSpectralAxis,
    const TAxisSampleList &Xtgt, Int32 cursor, Float64 xSrcStep) {
  Float64 xDestStep = NAN;
  Float64 xStepCompensation = 1.;
  if (cursor < targetSpectralAxis.GetSamplesCount() - 1) {
    xDestStep = Xtgt[cursor + 1] - Xtgt[cursor];
    xStepCompensation = xSrcStep / xDestStep;
  } else {
    xDestStep = Xtgt[cursor] - Xtgt[cursor - 1];
    xStepCompensation = xSrcStep / xDestStep;
  }
  return xStepCompensation;
}

std::unique_ptr<CRebin> CRebin::convert(const std::string opt_interp) && {
  if (opt_interp == "lin")
    return std::unique_ptr<CRebin>(new CRebinLinear(std::move(*this)));
  if (opt_interp == "preComputedFineGrid")
    return std::unique_ptr<CRebin>(new CRebinFineGrid(std::move(*this)));
  if (opt_interp == "spline")
    return std::unique_ptr<CRebin>(new CRebinSpline(std::move(*this)));
  if (opt_interp == "ngp")
    return std::unique_ptr<CRebin>(new CRebinNgp(std::move(*this)));

  THROWG(ErrorCode::INVALID_PARAMETER,
         "Only {lin, precomputedfinegrid, ngp, spline} values are "
         "supported for TemplateFittingSolver.interpolation");
}

std::unique_ptr<CRebin> CRebin::create(const std::string &opt_interp,
                                       const CSpectrum &spectrum) {
  if (opt_interp == "lin")
    return std::unique_ptr<CRebin>(new CRebinLinear(spectrum));
  if (opt_interp == "preComputedFineGrid")
    return std::unique_ptr<CRebin>(new CRebinFineGrid(spectrum));
  if (opt_interp == "spline")
    return std::unique_ptr<CRebin>(new CRebinSpline(spectrum));
  if (opt_interp == "ngp")
    return std::unique_ptr<CRebin>(new CRebinNgp(spectrum));

  THROWG(ErrorCode::INVALID_PARAMETER,
         "Only {lin, precomputedfinegrid, ngp, spline} values are "
         "supported for TemplateFittingSolver.interpolation");
}
