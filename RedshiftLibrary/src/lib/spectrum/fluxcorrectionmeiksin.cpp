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

using namespace NSEpic;

CSpectrumFluxCorrectionMeiksin::CSpectrumFluxCorrectionMeiksin(
    std::vector<MeiksinCorrection> meiksinCorrectionCurves)
    : m_rawCorrections(std::move(meiksinCorrectionCurves)),
      m_LambdaMin(m_rawCorrections[0].lbda.front()),
      m_LambdaMax(m_rawCorrections[0].lbda.back()) {}

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

  Float64 zStart = 2.0;
  Float64 zStep = 0.5;
  Float64 zStop = 7.0;
  if (z < zStart) {
    index = 0;
  } else if (z >= zStart && z < zStop) {
    index = (Int32)((z - zStart) / zStep + 1.0);
  } else if (z >= zStop) {
    index = (Int32)((zStop - zStart) / zStep);
  }
  return index;
}

// loop over lambda values
// for each lambda0, compute kernel_lambda0 and then multiply it by the igm
// curve
TFloat64List CSpectrumFluxCorrectionMeiksin::applyAdaptativeKernel(
    const TFloat64List &arr, const Float64 z_center,
    const std::shared_ptr<const CLSF> &lsf, const TFloat64List &lambdas) {
  if (!arr.size()) {
    throw GlobalException(ErrorCode::INTERNAL_ERROR,
                          "Cannot convolve: either kernel or array is empty. ");
  }

  Int32 n = arr.size(), Nhalf = -1;
  TFloat64List convolvedArr(arr);

  // determine the restframe convolution range, i.e., convolRange/(1+z_center)
  TFloat64Range convRange_rest(m_convolRange.GetBegin() / (1 + z_center),
                               m_convolRange.GetEnd() /
                                   (1 + z_center)); // conv range in restframe
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
    TFloat64List kernel = lsf->getRestFrameProfileVector(lambda0, z_center);
    Nhalf = int(kernel.size() / 2);
    if (!kernel.size()) {
      throw GlobalException(
          ErrorCode::INTERNAL_ERROR,
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

  m_corrections.resize(m_rawCorrections.size());

  TFloat64List meiksin_Bins = getSegmentsStartRedshiftList();
  Float64 zstep = meiksin_Bins[2] - meiksin_Bins[1];

  // iterate over the redshift list
  Float64 z_center;
  for (Int32 i = 0; i < m_rawCorrections.size(); i++) {
    if (i == 0)
      z_center = meiksin_Bins[i + 1] - zstep / 2.;
    else
      z_center = meiksin_Bins[i] + zstep / 2.;

    m_corrections[i].lbda = m_rawCorrections[i].lbda;

    for (Int32 j = 0; j < m_rawCorrections[i].fluxcorr.size();
         j++) // iterating over the different curves
    {
      TFloat64List a =
          applyAdaptativeKernel(m_rawCorrections[i].fluxcorr[j], z_center, lsf,
                                m_corrections[i].lbda);
      m_corrections[i].fluxcorr.push_back(std::move(a));
    }
  }

  m_convolved = true;
  return;
}
