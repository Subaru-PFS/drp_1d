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
#include "RedshiftLibrary/spectrum/fullspectrum.h"
#include "RedshiftLibrary/spectrum/rebin/rebinFineGrid.h"
#include "RedshiftLibrary/spectrum/rebin/rebinLinear.h"
#include "RedshiftLibrary/spectrum/rebin/rebinNgp.h"
#include "RedshiftLibrary/spectrum/rebin/rebinSpline.h"

using namespace NSEpic;
using namespace std;

CFullSpectrum::CFullSpectrum() : CSpectrum() {
  m_rebin = CRebin::create("linFull", *this);
}

CFullSpectrum::CFullSpectrum(const CFullSpectrum &other)
    : CSpectrum(other), m_mask(other.m_mask) {}

CFullSpectrum::CFullSpectrum(CFullSpectrum &&other)
    : CSpectrum(other), m_mask(std::move(other.m_mask)) {}

CFullSpectrum::CFullSpectrum(CSpectrumSpectralAxis spectralAxis,
                             CSpectrumFluxAxis fluxAxis,
                             TMaskList invalidPixels)
    : CSpectrum(std::move(spectralAxis), std::move(fluxAxis)),
      m_mask(std::move(invalidPixels)) {
  m_rebin = CRebin::create("linFull", *this);
}

CFullSpectrum::CFullSpectrum(const CFullSpectrum &other,
                             const TFloat64List &mask)
    : CSpectrum(other, mask) {

  auto new_mask = CSpectrumAxis::maskVector(
      mask, TFloat64List(other.m_mask.getMaskList().begin(),
                         other.m_mask.getMaskList().end()));
  m_mask = CMask(TMaskList(new_mask.begin(), new_mask.end()));
}

CFullSpectrum::CFullSpectrum(const std::string &name, const std::string &obsId)
    : CSpectrum(name, obsId) {}

CFullSpectrum::CFullSpectrum(CSpectrumSpectralAxis spectralAxis,
                             CSpectrumFluxAxis fluxAxis)
    : CSpectrum(spectralAxis, fluxAxis) {}

std::shared_ptr<CSpectrum> CFullSpectrum::getUnmaskedSpectrum() {
  return std::make_shared<CSpectrum>(
      *this,
      TFloat64List(m_mask.getMaskList().begin(), m_mask.getMaskList().end()));
}

/**
 * targetSpectralAxis should be expressed in same frame as source
 * SpetralAxis
 */
void CFullSpectrum::Rebin(const TFloat64Range &range,
                          const CSpectrumSpectralAxis &targetSpectralAxis,
                          CFullSpectrum &rebinedSpectrum, CMask &rebinedMask,
                          const std::string &opt_error_interp) const {
  ASSERT_CSpectrum_IS_VALID(*this);

  m_rebin->compute(range, targetSpectralAxis, rebinedSpectrum, rebinedMask,
                   opt_error_interp);

  rebinedSpectrum.setMask(rebinedMask);
}

bool CFullSpectrum::checkCorrectness(bool valid, Int32 index) const {

  if (!valid)
    Log.LogDebug(Formatter()
                 << "Invalid pixel with mask =" << (Int32)m_mask[index]);
  if (!valid && m_mask[index] == 0)
    return true;
  else
    return valid;
}
