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
#include "RedshiftLibrary/operator/templatefittingBase.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/operator/modelspectrumresult.h"
#include "RedshiftLibrary/processflow/context.h"

using namespace NSEpic;
using namespace std;

COperatorTemplateFittingBase::COperatorTemplateFittingBase(
    const TFloat64List &redshifts)
    : m_spectra(Context.getSpectra()),
      m_lambdaRanges(Context.getClampedLambdaRanges()), m_redshifts(redshifts),
      m_templateRebined_bf(m_spectra.size()),
      m_spcSpectralAxis_restframe(m_spectra.size()),
      m_mskRebined_bf(m_spectra.size()),
      m_maskBuilder(std::make_shared<CMaskBuilder>()){};

/**
 * \brief this function estimates the likelihood_cstLog term withing the
 * wavelength range
 **/
Float64 COperatorTemplateFittingBase::EstimateLikelihoodCstLog() const {

  Float64 cstLog = 0.0;
  for (auto it = std::make_tuple(m_spectra.begin(), m_lambdaRanges.begin());
       std::get<0>(it) != m_spectra.end();
       ++std::get<0>(it), ++std::get<1>(it)) {
    const auto &spectrum = **std::get<0>(it);
    const auto &lambdaRange = **std::get<1>(it);
    const CSpectrumSpectralAxis &spcSpectralAxis = spectrum.GetSpectralAxis();
    const TFloat64List &error =
        spectrum.GetFluxAxis().GetError().GetSamplesVector();
    ;

    Int32 numDevs = 0;

    Float64 sumLogNoise = 0.0;

    Int32 imin;
    Int32 imax;
    lambdaRange.getClosedIntervalIndices(spcSpectralAxis.GetSamplesVector(),
                                         imin, imax);
    for (Int32 j = imin; j <= imax; j++) {
      numDevs++;
      sumLogNoise += log(error[j]);
    }
    cstLog += -numDevs * 0.5 * log(2 * M_PI) - sumLogNoise;
  }
  return cstLog;
}

void COperatorTemplateFittingBase::InitIsmIgmConfig(
    const TInt32List kStart, const TInt32List kEnd, Float64 redshift,
    const std::shared_ptr<const CSpectrumFluxCorrectionCalzetti>
        &ismCorrectionCalzetti,
    const std::shared_ptr<const CSpectrumFluxCorrectionMeiksin>
        &igmCorrectionMeiksin) {

  for (Int32 spcIndex = 0; spcIndex < m_spectra.size(); spcIndex++)
    m_templateRebined_bf[spcIndex].InitIsmIgmConfig(
        kStart[spcIndex], kEnd[spcIndex], redshift, ismCorrectionCalzetti,
        igmCorrectionMeiksin);
}

// return tuple with photmetric values
std::tuple<std::shared_ptr<CModelSpectrumResult>, const TPhotVal>
COperatorTemplateFittingBase::ComputeSpectrumModel(
    const std::shared_ptr<const CTemplate> &tpl, Float64 redshift,
    Float64 EbmvCoeff, Int32 meiksinIdx, Float64 amplitude,
    const Float64 overlapThreshold, Int32 spcIndex) {
  Log.LogDetail("  Operator-COperatorTemplateFitting: building spectrum model "
                "templateFitting for candidate Zcand=%f",
                redshift);

  Float64 overlapFraction = 0.0;
  TFloat64Range currentRange;
  RebinTemplate(tpl, redshift, currentRange, overlapFraction, overlapThreshold,
                spcIndex);

  const TAxisSampleList &Xspc =
      m_spcSpectralAxis_restframe[spcIndex].GetSamplesVector();

  if ((EbmvCoeff > 0.) || (meiksinIdx > -1)) {
    TInt32List kstart(1);
    TInt32List kend(1);
    bool ret = currentRange.getClosedIntervalIndices(
        m_templateRebined_bf[spcIndex].GetSpectralAxis().GetSamplesVector(),
        kstart[0], kend[0]);
    if (!ret)
      THROWG(INTERNAL_ERROR, "Impossible to get valid kstart or kend");
    InitIsmIgmConfig(kstart, kend, redshift, tpl->m_ismCorrectionCalzetti,
                     tpl->m_igmCorrectionMeiksin);
  }

  if (EbmvCoeff > 0.) {
    if (m_templateRebined_bf[spcIndex].CalzettiInitFailed()) {
      THROWG(INTERNAL_ERROR, "ISM is not initialized");
    }
    Int32 idxEbmv = -1;
    idxEbmv =
        m_templateRebined_bf[spcIndex].m_ismCorrectionCalzetti->GetEbmvIndex(
            EbmvCoeff);

    if (idxEbmv != -1)
      ApplyDustCoeff(idxEbmv, spcIndex);
  }

  if (meiksinIdx > -1) {
    if (m_templateRebined_bf[spcIndex].MeiksinInitFailed()) {
      THROWG(INTERNAL_ERROR, "IGM in not initialized");
    }
    ApplyMeiksinCoeff(meiksinIdx, spcIndex);
  }
  ApplyAmplitude(amplitude, spcIndex);

  // shift the spectralaxis to sync with the spectrum lambdaAxis
  const CSpectrumFluxAxis &modelflux =
      m_templateRebined_bf[spcIndex].GetFluxAxis();
  CSpectrumSpectralAxis modelwav =
      m_templateRebined_bf[spcIndex]
          .GetSpectralAxis(); // needs a copy to be shifted
  modelwav.ShiftByWaveLength((1.0 + redshift),
                             CSpectrumSpectralAxis::nShiftForward);

  return std::make_tuple(std::make_shared<CModelSpectrumResult>(
                             CSpectrum(std::move(modelwav), modelflux)),
                         getIntegratedFluxes());
}

void COperatorTemplateFittingBase::RebinTemplate(
    const std::shared_ptr<const CTemplate> &tpl, Float64 redshift,
    TFloat64Range &currentRange, Float64 &overlapFraction,
    const Float64 overlapThreshold, Int32 spcIndex) {
  Float64 onePlusRedshift = 1.0 + redshift;

  // shift lambdaRange backward to be in restframe
  TFloat64Range spcLambdaRange_restframe;
  TFloat64Range lambdaRange_restframe(
      m_lambdaRanges[spcIndex]->GetBegin() / onePlusRedshift,
      m_lambdaRanges[spcIndex]->GetEnd() / onePlusRedshift);

  // redshift in restframe the tgtSpectralAxis, i.e., division by (1+Z)
  // Shouldn't ShiftByWaveLength be static ?
  m_spcSpectralAxis_restframe[spcIndex].ShiftByWaveLength(
      Context.getSpectra()[spcIndex]->GetSpectralAxis(), onePlusRedshift,
      CSpectrumSpectralAxis::nShiftBackward);
  m_spcSpectralAxis_restframe[spcIndex].ClampLambdaRange(
      lambdaRange_restframe, spcLambdaRange_restframe);

  // the spectral and tpl axis should be in the same scale
  const CSpectrumSpectralAxis &tplSpectralAxis = tpl->GetSpectralAxis();

  // Compute clamped lambda range over template in restframe
  TFloat64Range tplLambdaRange;
  tplSpectralAxis.ClampLambdaRange(lambdaRange_restframe, tplLambdaRange);
  // Compute the intersected range
  TFloat64Range intersectedLambdaRange(0.0, 0.0);
  TFloat64Range::Intersect(tplLambdaRange, spcLambdaRange_restframe,
                           intersectedLambdaRange);

  tpl->Rebin(intersectedLambdaRange, m_spcSpectralAxis_restframe[spcIndex],
             m_templateRebined_bf[spcIndex], m_mskRebined_bf[spcIndex]);

  // overlapFraction
  overlapFraction = m_spcSpectralAxis_restframe[spcIndex]
                        .IntersectMaskAndComputeOverlapFraction(
                            lambdaRange_restframe, m_mskRebined_bf[spcIndex]);

  // Check for overlap rate
  if (overlapFraction < overlapThreshold || overlapFraction <= 0.0) {
    // status = nStatus_NoOverlap;
    THROWG(OVERLAPFRACTION_NOTACCEPTABLE,
           Formatter() << "tpl overlap rate is too small: " << overlapFraction);
  }

  // the spectral axis should be in the same scale
  currentRange = intersectedLambdaRange;
}

// get z at which igm starts given that LyA starts at lbda_rest=1216
Float64
COperatorTemplateFittingBase::GetIGMStartingRedshiftValue(Float64 spcLbda0) {
  return spcLbda0 / RESTLAMBDA_LYA - 1; // check the rounding thing
}

void COperatorTemplateFittingBase::applyPositiveAndNonNullConstraint(
    const Float64 amp_sigma, Float64 &ampl) const {
  if (amp_sigma < m_continuum_null_amp_threshold)
    ampl = 0.;
  return;
}
