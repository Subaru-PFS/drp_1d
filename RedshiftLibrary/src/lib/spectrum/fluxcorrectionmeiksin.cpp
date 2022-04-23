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
using namespace NSEpic;

CSpectrumFluxCorrectionMeiksin::CSpectrumFluxCorrectionMeiksin(
    std::vector<MeiksinCorrection> meiksinCorrectionCurves, TFloat64List zbins)
    : m_rawCorrections(std::move(meiksinCorrectionCurves)),
      m_LambdaMin(m_rawCorrections[0].lbda.front()), m_LambdaMax(1215.),
      m_zbins(std::move(zbins)) {}

Float64 CSpectrumFluxCorrectionMeiksin::getCorrection(
    Float64 redshift, Int32 meiksinIdx, Float64 lambdaRest) const {

  if (lambdaRest > RESTLAMBDA_LYA) // && meiksinIdx == -1)
    return 1.;

  if (meiksinIdx == -1)
    THROWG(INTERNAL_ERROR, "igmIdx cannot be negative");
  Int32 zIdx = getRedshiftIndex(redshift);
  if (zIdx == -1)
    return 1.;
  Int32 lbdaIdx = getWaveIndex(lambdaRest);

  return getCorrection(zIdx, meiksinIdx, lbdaIdx);
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
  Int32 index = -1;
  TFloat64Index::getClosestLowerIndex(m_zbins, z, index);

  // keep last curves above last bin
  if (index == m_zbins.size() - 1)
    --index;
  return index;
}

Int32 CSpectrumFluxCorrectionMeiksin::getWaveIndex(Float64 w) const {
  return Int32((w - getLambdaMin()) / m_finegridstep + 0.5);
}

// loop over lambda values
// for each lambda0, compute kernel_lambda0 and then multiply it by the igm
// curve
TFloat64List CSpectrumFluxCorrectionMeiksin::applyLSFKernel(
    const TFloat64List &arr, const TFloat64List &lambdas,
    const TFloat64Range &zbin, const std::shared_ptr<const CLSF> &lsf) {
  if (!arr.size()) {
    THROWG(INTERNAL_ERROR,
           "Cannot convolve: either kernel or array is empty. ");
  }

  Int32 n = arr.size();
  TFloat64List convolvedArr(lambdas.size());

  // determine the restframe convolution range, i.e., convolRange/(1+z_center)
  TFloat64Range convRange_rest(m_convolRange.GetBegin() / (1 + zbin.GetEnd()),
                               m_convolRange.GetEnd() / (1 + zbin.GetBegin()));
  Int32 i_min = -1, i_max = -1;
  bool ret =
      convRange_rest.getClosedIntervalIndices(lambdas, i_min, i_max, false);
  if (!ret) {
    return convolvedArr;
  }
  Float64 tmp;
  for (Int32 i = i_min; i <= i_max; i++) {
    Float64 lambda0 = lambdas[i]; // lambda restframe
    // compute the adpative kernel at lambda0
    Float64 z_center = (zbin.GetBegin() + zbin.GetEnd()) / 2.;
    TFloat64List kernel = lsf->getRestFrameProfileVector(lambda0, z_center);
    Int32 Nhalf = int(kernel.size() / 2);
    if (!kernel.size()) {
      THROWG(INTERNAL_ERROR,
             "Cannot convolve: either kernel or array is empty. ");
    }

    tmp = 0.0;                              // sum over the intersection area
    for (Int32 j = -Nhalf; j <= Nhalf; j++) // center kernel at arr[j]
    {
      if (i + j >= 0) {
        if (i + j >= n) {
          tmp += kernel[Nhalf + j];
        } else {
          tmp += arr[i + j] * kernel[Nhalf + j];
        }
      } else {
        tmp += arr[0] * kernel[Nhalf + j];
      }
    }
    convolvedArr[i] = tmp;
  }
  return convolvedArr;
}
/**
 * convolve only on m_convolRange/(1+zbin_meiksin), while keeping vector size
 */
void CSpectrumFluxCorrectionMeiksin::convolveByLSF(
    const std::shared_ptr<const CLSF> &lsf, const TFloat64Range &convolRange) {

  m_convolRange = convolRange;

  TFloat64Range range(m_LambdaMin, m_LambdaMax);
  TFloat64List finelbdaGrid = range.SpreadOver(m_finegridstep);

  std::vector<MeiksinCorrection> corrections(m_rawCorrections.size());
  m_corrections.resize(m_rawCorrections.size());

  for (Int32 i = 0; i < m_rawCorrections.size(); i++) {
    // z_center = (m_zbins[i + 1] + m_zbins[i]) / 2.;
    TFloat64Range zbin(m_zbins[i], m_zbins[i + 1]);
    m_corrections[i].lbda = finelbdaGrid;

    for (Int32 j = 0; j < m_rawCorrections[i].fluxcorr.size(); j++) {
      TFloat64List interpolatedConvolvedArr = applyLSFKernel(
          m_rawCorrections[i].fluxcorr[j], m_corrections[i].lbda, zbin, lsf);
      m_corrections[i].fluxcorr.push_back(std::move(interpolatedConvolvedArr));
    }
    m_corrections[i].lbda = finelbdaGrid;
  }

  m_convolved = true;
  return;
}
