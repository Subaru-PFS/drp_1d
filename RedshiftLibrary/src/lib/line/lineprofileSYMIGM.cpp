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
    : CLineProfileSYM(nsigmasupport, SYMIGM),
      m_igmCorrectionMeiksin(igmcorrectionMeiksin) {}

void CLineProfileSYMIGM::CheckMeiksinInit() const {
  if (!m_igmCorrectionMeiksin) {
    THROWG(INTERNAL_ERROR, "m_igmCorrectionMeiksin is nullptr");
  }
}

Float64 CLineProfileSYMIGM::GetLineProfileVal(Float64 x, Float64 x0,
                                              Float64 sigma) const {
  CheckMeiksinInit();

  Float64 igm_fluxcorr_lbda = getIGMCorrection(x);

  Float64 val =
      CLineProfileSYM::GetLineProfileVal(x, x0, sigma) * igm_fluxcorr_lbda;
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
Float64 CLineProfileSYMIGM::GetLineFlux(Float64 x0, Float64 sigma,
                                        Float64 A) const {
  CheckMeiksinInit();

  const Float64 x0_rest = x0 / (1 + m_redshift);
  const Float64 winsize_rest = sigma / (1 + m_redshift) * N_SIGMA_SUPPORT;
  TFloat64Range range_rest(x0_rest - winsize_rest / 2,
                           x0_rest + winsize_rest / 2);

  // check igm, if not applicable return Gaussian flux
  if (m_igmidx < 0 || range_rest.GetBegin() > RESTLAMBDA_LYA)
    return A * sigma * sqrt(2 * M_PI);

  // Line support aligned to lambda rest grid of convolved igm
  TFloat64List xlist_rest = m_igmCorrectionMeiksin->getWaveVector(range_rest);
  if (xlist_rest.empty())
    THROWG(INTERNAL_ERROR,
           "not enough IGM extinction samples inside line profile to compute "
           "flux, increase IGM extinction curves resolution.");
  Float64 flux = 0.;
  for (Float64 x_rest : xlist_rest)
    flux += GetLineProfileVal(x_rest * (1 + m_redshift), x0, sigma);
  flux *= (xlist_rest.back() - xlist_rest.front()) * (1 + m_redshift);

  return A * flux;
}

// we should derive the IgmCorrection  at z_binwith respect to Z
Float64 CLineProfileSYMIGM::GetLineProfileDerivX0(Float64 x, Float64 x0,
                                                  Float64 sigma) const {
  CheckMeiksinInit();

  Float64 igm_fluxcorr = NAN;
  Float64 igm_fluxcorr_derivX0 = NAN;
  std::tie(igm_fluxcorr, igm_fluxcorr_derivX0) =
      getIGMCorrectionAndDerivX0(x, x0);

  Float64 val =
      CLineProfileSYM::GetLineProfileDerivX0(x, x0, sigma) * igm_fluxcorr;
  val +=
      CLineProfileSYM::GetLineProfileVal(x, x0, sigma) * igm_fluxcorr_derivX0;
  return val;
}

Float64 CLineProfileSYMIGM::GetLineProfileDerivSigma(Float64 x, Float64 x0,
                                                     Float64 sigma) const {
  CheckMeiksinInit();

  Float64 igm_fluxcorr_lbda = getIGMCorrection(x);

  Float64 val = CLineProfileSYM::GetLineProfileDerivSigma(x, x0, sigma) *
                igm_fluxcorr_lbda;
  return val;
}

TSymIgmParams CLineProfileSYMIGM::GetSymIgmParams() const {
  return TSymIgmParams(m_igmidx, m_redshift);
}

void CLineProfileSYMIGM::SetSymIgmParams(const TSymIgmParams &params) {
  m_redshift = params.m_redshift;
  m_igmidx = params.m_igmidx;
}

void CLineProfileSYMIGM::resetParams() {
  if (m_igmFit)
    SetSymIgmParams({});
}

Int32 CLineProfileSYMIGM::getIGMIdxCount() const {
  return m_igmCorrectionMeiksin->getIdxCount();
}

Float64 CLineProfileSYMIGM::getIGMCorrection(Float64 x) const {
  if (m_igmidx < 0)
    return 1.0;
  return m_igmCorrectionMeiksin->getCorrection(m_redshift, m_igmidx,
                                               x / (1 + m_redshift));
}

std::tuple<Float64, Float64>
CLineProfileSYMIGM::getIGMCorrectionAndDerivX0(Float64 x, Float64 x0) const {
  if (m_igmidx < 0)
    return std::make_tuple(1.0, 0.0);

  Float64 xRest = x / (1 + m_redshift);
  Float64 derivLbdaRest = NAN;
  Float64 corr = NAN;
  std::tie(corr, derivLbdaRest) =
      m_igmCorrectionMeiksin->getCorrectionAndDerivLbdaRest(m_redshift,
                                                            m_igmidx, xRest);
  Float64 derivX0 = -xRest / x0 * derivLbdaRest;
  return std::make_tuple(corr, derivX0);
}