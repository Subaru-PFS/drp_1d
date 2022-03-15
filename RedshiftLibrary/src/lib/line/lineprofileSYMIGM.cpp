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
#include "RedshiftLibrary/line/lineprofileSYMIGM.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/flag.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/line/lineprofileSYM.h"
#include "RedshiftLibrary/log/log.h"
using namespace NSEpic;
using namespace std;

CLineProfileSYMIGM::CLineProfileSYMIGM(
    const std::shared_ptr<CSpectrumFluxCorrectionMeiksin> &igmcorrectionMeiksin,
    const Float64 nsigmasupport)
    : CLineProfileSYM(SYMIGM, nsigmasupport),
      m_igmCorrectionMeiksin(igmcorrectionMeiksin) {}

void CLineProfileSYMIGM::MeiksinInitFailed() const {
  if (!m_igmCorrectionMeiksin) {
    THROWG(INTERNAL_ERROR, "m_igmCorrectionMeiksin is nullptr");
  }
}

Float64 CLineProfileSYMIGM::GetLineProfile(Float64 x, Float64 x0, Float64 sigma,
                                           Float64 redshift,
                                           Int32 igmIdx) const {
  MeiksinInitFailed();

  Float64 igm_fluxcorr_lbda = getIGMCorrection(x, redshift, igmIdx);

  Float64 val =
      CLineProfileSYM::GetLineProfile(x, x0, sigma) * igm_fluxcorr_lbda;
  return val;
}

/**
 * @brief Integrate over sigma range centered around mu
 * Similar to LineProfileSYM, we assume that amplitude is constant?
 *
 * @param A
 * @param sigma
 * @param redshift
 * @param mu
 * @param igmIdx
 * @return Float64
 */
Float64 CLineProfileSYMIGM::GetLineFlux(Float64 A, Float64 sigma,
                                        Float64 redshift, Float64 mu,
                                        Int32 igmIdx) const {
  MeiksinInitFailed();
  TFloat64Range range(mu - sigma / 2, mu + sigma / 2);
  Float64 step = IGM_INTERPOLATION_RATIO;
  TFloat64List list = range.SpreadOver(step);

  Float64 flux = 0.;
  for (Float64 x : list) {
    Float64 igm_fluxcorr_lbda = getIGMCorrection(x, redshift, igmIdx);
    flux += igm_fluxcorr_lbda;
  }
  return A * flux;
}

// we should derivate the IgmCorrection  at z_binwith respect to Z
Float64 CLineProfileSYMIGM::GetLineProfileDerivZ(Float64 x, Float64 lambda0,
                                                 Float64 redshift,
                                                 Float64 sigma,
                                                 Int32 igmIdx) const {
  THROWG(INTERNAL_ERROR,
         "GetLineProfileDerivZ is not yet defined for SYMIGM lineprofile");
}

Float64 CLineProfileSYMIGM::GetLineProfileDerivSigma(Float64 x, Float64 x0,
                                                     Float64 sigma,
                                                     Float64 redshift,
                                                     Int32 igmIdx) const {

  MeiksinInitFailed();

  Float64 igm_fluxcorr_lbda = getIGMCorrection(x, redshift, igmIdx);

  Float64 val = CLineProfileSYMIGM::GetLineProfileDerivSigma(x, x0, sigma,
                                                             redshift, igmIdx) *
                igm_fluxcorr_lbda;
  return val;
}

Int32 CLineProfileSYMIGM::getIGMIdxCount() const {
  return m_igmCorrectionMeiksin->getIdxCount();
}

Int32 CLineProfileSYMIGM::getIGMCorrection(Float64 x, Float64 redshift,
                                           Int32 igmIdx) const {
  return m_igmCorrectionMeiksin->getCorrection(redshift, igmIdx,
                                               x / (1 + redshift));
}