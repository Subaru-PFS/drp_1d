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

using namespace NSEpic;
using namespace std;

// CRebin::CRebin(std::unique_ptr<CRebin> &&other) {
//   m_pfgFlux = std::move(other->m_pfgFlux);
//   m_FineGridInterpolated = std::move(other->m_FineGridInterpolated);
//   // m_spectrum = std::move(other->m_spectrum);
// }

// CRebin &CRebin::operator=(std::unique_ptr<CRebin> &&other) {
//   m_pfgFlux = std::move(other->m_pfgFlux);
//   m_FineGridInterpolated = std::move(other->m_FineGridInterpolated);
//   return *this;
// }

CRebin::CRebin(std::unique_ptr<CRebin> &&other)
    : m_pfgFlux(std::move(other->m_pfgFlux)),
      m_FineGridInterpolated(std::move(other->m_FineGridInterpolated)),
      m_spectrum(other->m_spectrum) {}

bool CRebin::compute(const TFloat64Range &range,
                     const CSpectrumSpectralAxis &targetSpectralAxis,
                     CSpectrum &rebinedSpectrum, CMask &rebinedMask,
                     const std::string m_opt_error_interp) {

  Int32 s = targetSpectralAxis.GetSamplesCount();
  CSpectrumSpectralAxis spectralAxis = m_spectrum.GetSpectralAxis();

  // find start/end indexs for both axes
  if (spectralAxis[0] > range.GetBegin() ||
      spectralAxis[spectralAxis.GetSamplesCount() - 1] < range.GetEnd()) {
    THROWG(INTERNAL_ERROR, "input spectral range is not "
                           "included in spectral axis");
  }

  CSpectrumFluxAxis rebinedFluxAxis =
      std::move(rebinedSpectrum.GetRawFluxAxis());
  rebinedFluxAxis.SetSize(s); // does not re-allocate if already allocated
  rebinedMask.SetSize(s);

  const TAxisSampleList &Xsrc = spectralAxis.GetSamplesVector();
  const TAxisSampleList &Ysrc = m_spectrum.GetFluxAxis().GetSamplesVector();
  const TAxisSampleList &Xtgt = targetSpectralAxis.GetSamplesVector();
  TAxisSampleList &Yrebin = rebinedFluxAxis.GetSamplesVector();
  const TFloat64List &Error =
      m_spectrum.GetFluxAxis().GetError().GetSamplesVector();
  TFloat64List &ErrorRebin = rebinedFluxAxis.GetError().GetSamplesVector();

  // Move cursors up to lambda range start
  m_cursor = 0;
  while (m_cursor < targetSpectralAxis.GetSamplesCount() &&
         Xtgt[m_cursor] < range.GetBegin()) {
    rebinedMask[m_cursor] = 0;
    Yrebin[m_cursor] = 0.0;
    if (m_opt_error_interp == "rebin" || m_opt_error_interp == "rebinVariance")
      ErrorRebin[m_cursor] = INFINITY;
    m_cursor++;
  }

  bool status =
      rebin(rebinedFluxAxis, range, targetSpectralAxis, rebinedSpectrum,
            rebinedMask, m_opt_error_interp, Xsrc, Ysrc, Xtgt, Error);

  while (m_cursor < targetSpectralAxis.GetSamplesCount()) {
    rebinedMask[m_cursor] = 0;
    Yrebin[m_cursor] = 0.0;
    if (m_opt_error_interp == "rebin" || m_opt_error_interp == "rebinVariance")
      ErrorRebin[m_cursor] = INFINITY;
    m_cursor++;
  }

  rebinedSpectrum.ResetContinuum();
  rebinedSpectrum.SetType(CSpectrum::EType::nType_raw);
  rebinedSpectrum.SetSpectralAndFluxAxes(targetSpectralAxis,
                                         std::move(rebinedFluxAxis));

  return status;
}
