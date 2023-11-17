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
#include "RedshiftLibrary/spectrum/fluxcorrectionmeiksin.h"
#include "RedshiftLibrary/common/indexing.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/spectrum/LSF.h"

#include <algorithm>

using namespace NSEpic;

CSpectrumFluxCorrectionMeiksin::CSpectrumFluxCorrectionMeiksin(
    std::vector<MeiksinCorrection> meiksinCorrectionCurves, TFloat64List zbins)
    : m_rawCorrections(std::move(meiksinCorrectionCurves)),
      m_LambdaMin(m_rawCorrections.at(0).lbda.at(0)),
      m_LambdaMax(RESTLAMBDA_LYA),
      m_LambdaSize(m_rawCorrections.at(0).lbda.size()),
      m_zbins(std::move(zbins)) {}

Float64 CSpectrumFluxCorrectionMeiksin::getCorrection(
    Float64 redshift, Int32 meiksinIdx, Float64 lambdaRest) const {

  if (lambdaRest > RESTLAMBDA_LYA) // && meiksinIdx == -1)
    return 1.;

  if (meiksinIdx == undefIdx)
    THROWG(INTERNAL_ERROR, "igmIdx undefined");
  Int32 zIdx = getRedshiftIndex(redshift);
  if (zIdx == undefIdx)
    return 1.;
  Int32 lbdaIdx = getWaveIndex(lambdaRest);

  return getCorrection(zIdx, meiksinIdx, lbdaIdx);
}

std::tuple<Float64, Float64>
CSpectrumFluxCorrectionMeiksin::getCorrectionAndDerivLbdaRest(
    Float64 redshift, Int32 meiksinIdx, Float64 lambdaRest) const {

  if (lambdaRest > RESTLAMBDA_LYA) // && meiksinIdx == -1)
    return std::make_tuple(1., 0.);

  if (meiksinIdx == -1)
    THROWG(INTERNAL_ERROR, "igmIdx cannot be negative");
  Int32 zIdx = getRedshiftIndex(redshift);
  if (zIdx == undefIdx)
    return std::make_tuple(1., 0.);
  Int32 lbdaIdx = getWaveIndex(lambdaRest);

  return std::make_tuple(getCorrection(zIdx, meiksinIdx, lbdaIdx),
                         getCorrectionDerivLbdaRest(zIdx, meiksinIdx, lbdaIdx));
}

Float64 CSpectrumFluxCorrectionMeiksin::getCorrectionDerivLbdaRest(
    Int32 zIdx, Int32 meiksinIdx, Int32 lbdaIdx) const {
  Int32 lbdaI0 = std::max(0, lbdaIdx - 1);
  Int32 lbdaI1 = std::min(m_fineLambdaSize - 1, lbdaIdx + 1);
  Float64 corr_diff = m_corrections[zIdx].fluxcorr[meiksinIdx].at(lbdaI1) -
                      m_corrections[zIdx].fluxcorr[meiksinIdx].at(lbdaI0);
  Float64 lbda_diff =
      m_corrections[zIdx].lbda.at(lbdaI1) - m_corrections[zIdx].lbda.at(lbdaI0);
  return corr_diff / lbda_diff;
}

/**
 * @brief CSpectrumFluxCorrectionMeiksin::getRedshiftIndex
 * @param z
 *
 * Returns the index corresponding to the table to be used for that redshift
 * value
 *
 * @return
 */
Int32 CSpectrumFluxCorrectionMeiksin::getRedshiftIndex(Float64 z) const {
  Int32 index = undefIdx;
  TFloat64Index::getClosestLowerIndex(m_zbins, z, index);

  // keep last curves above last bin
  if (index == m_zbins.size() - 1)
    --index;
  return index;
}

// NGP interp
Int32 CSpectrumFluxCorrectionMeiksin::getWaveIndex(Float64 w) const {
  Int32 LambdaSize = m_fineLambdaSize;
  Int32 wIdx = Int32(std::round((w - getLambdaMin()) / m_finegridstep));
  return std::min(LambdaSize - 1, std::max(0, wIdx));
}

// get wave vector inside range and aligned to fine lambda grid
TFloat64List CSpectrumFluxCorrectionMeiksin::getWaveVector(
    const TFloat64Range &wrange) const {
  return getWaveVector(wrange, false);
}

TFloat64List
CSpectrumFluxCorrectionMeiksin::getWaveVector(const TFloat64Range &wrange,
                                              bool raw) const {
  TFloat64List waves;
  Float64 step = raw ? IGM_RAW_STEP : m_finegridstep;
  TInt32Range indices = getWaveRangeIndices(wrange, raw);
  if (indices.GetLength() < 0)
    return waves;
  waves.resize(indices.GetLength() + 1);
  for (std::size_t i = 0; i < waves.size(); ++i)
    waves[i] = getLambdaMin() + step * (indices.GetBegin() + i);

  return waves;
}

TInt32Range
CSpectrumFluxCorrectionMeiksin::getWaveRangeIndices(const TFloat64Range &wrange,
                                                    bool raw) const {
  Float64 step = raw ? IGM_RAW_STEP : m_finegridstep;
  Int32 LambdaSize = raw ? m_LambdaSize : m_fineLambdaSize;
  Int32 imin = std::min(
      LambdaSize - 1,
      std::max(0,
               Int32(std::ceil((wrange.GetBegin() - getLambdaMin()) / step))));
  Int32 imax = std::min(
      LambdaSize - 1,
      std::max(0,
               Int32(std::floor((wrange.GetEnd() - getLambdaMin()) / step))));
  return TInt32Range(imin, imax);
}

// loop over lambda values
// for each lambda0, compute kernel_lambda0 and then multiply it by the igm
// curve
TFloat64List CSpectrumFluxCorrectionMeiksin::ConvolveByLSFOneCurve(
    const TFloat64List &arr, const TFloat64List &lambdas,
    const TFloat64List &fineLambdas, const TFloat64Range &zbin,
    const std::shared_ptr<const CLSF> &lsf) const {
  if (!arr.size()) {
    THROWG(INTERNAL_ERROR, "Cannot convolve: array is empty. ");
  }

  Int32 n = arr.size();
  TFloat64List convolvedArr(fineLambdas.size());

  // determine the restframe convolution range, i.e., convolRange/(1+z_center)
  TFloat64Range convRange_rest(m_convolRange.GetBegin() / (1 + zbin.GetEnd()),
                               m_convolRange.GetEnd() / (1 + zbin.GetBegin()));
  bool overlap = convRange_rest.IntersectWith(TFloat64Range(fineLambdas));
  if (!overlap)
    return convolvedArr;

  TInt32Range indices = getWaveRangeIndices(convRange_rest, false);
  Float64 z_center = (zbin.GetBegin() + zbin.GetEnd()) / 2.;
  Float64 sigmaSupport =
      lsf->GetProfile()->GetNSigmaSupport() / 2. / (1.0 + z_center);
  for (std::size_t i = indices.GetBegin(); i <= indices.GetEnd(); i++) {
    Float64 lambda0 = fineLambdas[i]; // lambda restframe

    // compute the LSF kernel centered at lambda0 (restframe)
    Float64 lambda0_obs = lambda0 * (1 + z_center);
    // clip lambda0 to lambdaRange for z != z_center
    Float64 half_lsf_range =
        lsf->GetWidth(lambda0_obs, true) * sigmaSupport; // restframe
    TFloat64Range lsf_range(lambda0 - half_lsf_range,
                            lambda0 + half_lsf_range); // restframe
    TInt32Range lsf_indices = getWaveRangeIndices(lsf_range, true);
    if (lsf_indices.GetBegin() < 0 || lsf_indices.GetEnd() >= m_LambdaSize)
      THROWG(INTERNAL_ERROR,
             "lsf kernel does not overlap meiksin curve samples");
    if (lsf_indices.GetLength() < 0)
      THROWG(INTERNAL_ERROR,
             "lsf kernel smaller than Meiksin wavelength steps");
    TFloat64List lambdas_obs(lsf_indices.GetLength() + 1);
    std::transform(lambdas.cbegin() + lsf_indices.GetBegin(),
                   lambdas.cbegin() + lsf_indices.GetEnd() + 1,
                   lambdas_obs.begin(),
                   [z_center](Float64 v) { return v * (1.0 + z_center); });
    TFloat64List kernel =
        lsf->getNormalizedProfileVector(lambdas_obs, lambda0_obs);

    for (std::size_t j = 0; j < kernel.size(); ++j)
      convolvedArr[i] += kernel[j] * arr[lsf_indices.GetBegin() + j];
  }
  return convolvedArr;
}

/**
 * convolve only on m_convolRange/(1+zbin_meiksin), while keeping vector size
 */
void CSpectrumFluxCorrectionMeiksin::convolveByLSF(
    const std::shared_ptr<const CLSF> &lsf, const TFloat64Range &convolRange) {

  m_convolRange = convolRange;

  TFloat64Range range(m_LambdaMin, m_rawCorrections[0].lbda.back());
  TFloat64List finelbdaGrid = range.SpreadOverEpsilon(m_finegridstep);
  m_fineLambdaSize = finelbdaGrid.size();

  // std::vector<MeiksinCorrection> corrections(m_rawCorrections.size());
  m_corrections.resize(m_rawCorrections.size());

  for (std::size_t i = 0; i < m_rawCorrections.size(); i++) {
    // z_center = (m_zbins[i + 1] + m_zbins[i]) / 2.;
    TFloat64Range zbin(m_zbins[i], m_zbins[i + 1]);
    m_corrections[i].lbda = finelbdaGrid;

    for (std::size_t j = 0; j < m_rawCorrections[i].fluxcorr.size(); j++) {
      TFloat64List interpolatedConvolvedArr = ConvolveByLSFOneCurve(
          m_rawCorrections[i].fluxcorr[j], m_rawCorrections[i].lbda,
          finelbdaGrid, zbin, lsf);
      m_corrections[i].fluxcorr.push_back(std::move(interpolatedConvolvedArr));
    }
  }

  m_convolved = true;
  return;
}
