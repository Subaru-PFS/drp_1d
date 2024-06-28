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

#include "RedshiftLibrary/operator/templatefittingwithphot.h"
#include "RedshiftLibrary/common/datatypes.h"

using namespace NSEpic;
using namespace std;

COperatorTemplateFittingPhot::COperatorTemplateFittingPhot(
    const std::shared_ptr<const CPhotBandCatalog> &photbandcat,
    const Float64 weight, const TFloat64List &redshifts)
    : COperatorTemplateFitting(redshifts), m_photBandCat(photbandcat),
      m_weight(weight) {

  if (m_spectra.size() > 1)
    THROWG(ErrorCode::MULTIOBS_WITH_PHOTOMETRY_NOTIMPLEMENTED,
           "Photometry not supported with multiobs");
  // check availability and coherence of photometric bands & data
  checkInputPhotometry();

  // initialize restframe photometric axis and rebined photmetric template
  for (auto const &band : *m_photBandCat) {
    const string &bandName = band.first;
    m_photSpectralAxis_restframe[bandName]; // default initialisation
    m_templateRebined_phot[bandName];       // default initialisation
  }

  // sort photometric band by increasing lambda (by lower bound of the bands)
  //  for IGM processing
  m_sortedBandNames = m_photBandCat->GetNameListSortedByLambda();
}

void COperatorTemplateFittingPhot::checkInputPhotometry() const {
  if (m_photBandCat == nullptr)
    THROWG(ErrorCode::MISSING_PHOTOMETRIC_TRANSMISSION,
           "Photometric band "
           "transmision not available");
  if (m_photBandCat->empty())
    THROWG(ErrorCode::MISSING_PHOTOMETRIC_TRANSMISSION,
           "Empty photometric band transmission");

  if (m_spectra[0]->GetPhotData() == nullptr)
    THROWG(ErrorCode::MISSING_PHOTOMETRIC_DATA,
           "photometric data not available in spectrum");

  const auto &dataNames = m_spectra[0]->GetPhotData()->GetNameList();
  for (const auto &bandName : m_photBandCat->GetNameList())
    if (std::find(dataNames.cbegin(), dataNames.cend(), bandName) ==
        dataNames.cend())
      THROWG(ErrorCode::MISSING_PHOTOMETRIC_DATA,
             Formatter() << " "
                            "photometry point for band name: "
                         << bandName << " is not available in the spectrum");
}

void COperatorTemplateFittingPhot::RebinTemplate(
    const std::shared_ptr<const CTemplate> &tpl, Float64 redshift,
    TFloat64Range &currentRange, Float64 &overlapFraction,
    const Float64 overlapThreshold, Int32 spcIndex) {

  COperatorTemplateFittingBase::RebinTemplate(
      tpl, redshift, currentRange, overlapFraction, overlapThreshold, spcIndex);

  if (spcIndex == 0)
    RebinTemplateOnPhotBand(tpl, redshift);
}

void COperatorTemplateFittingPhot::RebinTemplateOnPhotBand(
    const std::shared_ptr<const CTemplate> &tpl, Float64 redshift) {

  Float64 onePlusRedshift = 1.0 + redshift;

  for (auto const &band : *m_photBandCat) {
    const string &bandName = band.first;
    const TFloat64List &bandLambda = band.second.GetWavelength();

    CSpectrumSpectralAxis &photSpectralAxis_restframe =
        m_photSpectralAxis_restframe[bandName];
    CTemplate &templateRebined_phot = m_templateRebined_phot[bandName];

    CSpectrumSpectralAxis photSpectralaxis = bandLambda;

    photSpectralAxis_restframe.ShiftByWaveLength(
        photSpectralaxis, onePlusRedshift,
        CSpectrumSpectralAxis::nShiftBackward);

    CMask mskRebined;
    const TFloat64Range lambdaRange_restframe =
        photSpectralAxis_restframe.GetLambdaRange();
    tpl->Rebin(lambdaRange_restframe, photSpectralAxis_restframe,
               templateRebined_phot, mskRebined);

    const Float64 overlapFraction =
        photSpectralAxis_restframe.IntersectMaskAndComputeOverlapFraction(
            lambdaRange_restframe, mskRebined);

    if (overlapFraction < 1.0) {
      THROWG(ErrorCode::OVERLAPFRACTION_NOTACCEPTABLE,
             Formatter() << "tpl overlap too small: " << overlapFraction);
    }
  }
}

void COperatorTemplateFittingPhot::InitIsmIgmConfig(
    Float64 redshift, Int32 kstart, Int32 kend,
    const std::shared_ptr<const CSpectrumFluxCorrectionCalzetti>
        &ismCorrectionCalzetti,
    const std::shared_ptr<const CSpectrumFluxCorrectionMeiksin>
        &igmCorrectionMeiksin,
    Int32 spcIndex) {

  COperatorTemplateFitting::InitIsmIgmConfig(redshift, kstart, kend,
                                             ismCorrectionCalzetti,
                                             igmCorrectionMeiksin, spcIndex);

  if (spcIndex > 0)
    return;

  // init ism & igm on all rebined photometric templates
  for (auto &band : m_templateRebined_phot)
    band.second.InitIsmIgmConfig(redshift, ismCorrectionCalzetti,
                                 igmCorrectionMeiksin);
}

void COperatorTemplateFittingPhot::init_fast_igm_processing(
    Int32 EbmvListSize) {

  COperatorTemplateFitting::init_fast_igm_processing(EbmvListSize);

  // determine last band not impacted by IGM
  m_BandIgmEnd = std::find_if(
      m_sortedBandNames.cbegin(), m_sortedBandNames.cend(),
      [this](std::string const &s) {
        return m_templateRebined_phot.at(s).GetIgmEndIndex() == undefIdx;
      });

  m_sumCross_outsideIGM_phot.assign(EbmvListSize, 0.0);
  m_sumS_outsideIGM_phot.assign(EbmvListSize, 0.0);
  m_sumT_outsideIGM_phot.assign(EbmvListSize, 0.0);
}

bool COperatorTemplateFittingPhot::igmIsInRange(
    const TFloat64RangeList &ranges) const {

  if (COperatorTemplateFitting::igmIsInRange(ranges))
    return true;

  for (const auto &band : m_templateRebined_phot)
    if (COperatorTemplateFitting::igmIsInRange({band.second.GetLambdaRange()}))
      return true;

  return false;
}

bool COperatorTemplateFittingPhot::ApplyMeiksinCoeff(Int32 meiksinIdx,
                                                     Int32 spcIndex) {
  bool ret = COperatorTemplateFitting::ApplyMeiksinCoeff(meiksinIdx, spcIndex);

  if (spcIndex > 0)
    return ret;

  for (auto &band : m_templateRebined_phot)
    if (band.second.ApplyMeiksinCoeff(meiksinIdx))
      ret = true;

  return ret;
}

bool COperatorTemplateFittingPhot::ApplyDustCoeff(Int32 kEbmv, Int32 spcIndex) {
  bool ret = COperatorTemplateFitting::ApplyDustCoeff(kEbmv, spcIndex);

  if (spcIndex > 0)
    return ret;

  for (auto &band : m_templateRebined_phot)
    if (band.second.ApplyDustCoeff(kEbmv))
      ret = true;

  return ret;
}

TCrossProductResult COperatorTemplateFittingPhot::ComputeCrossProducts(
    Int32 kM, Int32 kEbmv_, Float64 redshift, Int32 spcIndex) {
  TCrossProductResult crossResult =
      COperatorTemplateFitting::ComputeCrossProducts(kM, kEbmv_, redshift,
                                                     spcIndex);
  if (spcIndex == 0)
    ComputePhotCrossProducts(kM, kEbmv_, crossResult);

  return crossResult;
}

void COperatorTemplateFittingPhot::ComputeAmplitudeAndChi2(
    TFittingResult &fitResult,
    const CPriorHelper::SPriorTZE &logpriorTZ) const {
  COperatorTemplateFitting::ComputeAmplitudeAndChi2(fitResult, logpriorTZ);
  // save photometric chi2 part
  const Float64 &ampl = fitResult.ampl;
  fitResult.chiSquare_phot = fitResult.cross_result.sumS_phot +
                             fitResult.cross_result.sumT_phot * ampl * ampl -
                             2. * ampl * fitResult.cross_result.sumCross_phot;
}

void COperatorTemplateFittingPhot::ComputePhotCrossProducts(
    Int32 kM, Int32 kEbmv_, TCrossProductResult &fitResult) {

  Float64 &sumCross_phot = fitResult.sumCross_phot;
  Float64 &sumT_phot = fitResult.sumT_phot;
  Float64 &sumS_phot = fitResult.sumS_phot;

  Float64 sumCross_IGM = 0.0;
  Float64 sumT_IGM = 0.0;
  Float64 sumS_IGM = 0.0;
  bool sumsIgmSaved = false;

  const auto &photData = m_spectra[0]->GetPhotData();

  // determine bands impacted by IGM
  const auto &Start = m_sortedBandNames.cbegin();
  const auto &End = m_sortedBandNames.cend();

  auto EndLoop = m_option_igmFastProcessing && kM > 0 ? m_BandIgmEnd : End;

  // compute photometric Leastquare
  for (auto it = Start; it != EndLoop; it++) {

    if (m_option_igmFastProcessing && !sumsIgmSaved && it == m_BandIgmEnd) {
      // store intermediate sums for IGM range
      sumCross_IGM = sumCross_phot;
      sumT_IGM = sumT_phot;
      sumS_IGM = sumS_phot;
      sumsIgmSaved = true;
    }

    const auto &bandName = *it;
    const auto &photBand = m_photBandCat->at(bandName);
    const auto &flux =
        m_templateRebined_phot.at(bandName).GetFluxAxis().GetSamplesVector();

    // integrate flux
    const auto integ_flux = photBand.IntegrateFlux(flux);

    // add to leastsquare sum
    const auto d = photData->GetFlux(bandName);
    const auto oneOverErr2 =
        photData->GetOneOverErr2(bandName) * m_weight * m_weight;
    const auto fluxOverErr2 =
        photData->GetFluxOverErr2(bandName) * m_weight * m_weight;

    if (std::isinf(oneOverErr2) || std::isnan(oneOverErr2))
      THROWG(ErrorCode::INTERNAL_ERROR,
             Formatter() << "found invalid inverse variance : err2="
                         << oneOverErr2 << ", for band=" << bandName);

    sumCross_phot += d * integ_flux * oneOverErr2;
    sumT_phot += integ_flux * integ_flux * oneOverErr2;
    sumS_phot += fluxOverErr2;
  }

  if (m_option_igmFastProcessing) {
    if (kM == 0) {
      m_sumCross_outsideIGM_phot[kEbmv_] = sumCross_phot - sumCross_IGM;
      m_sumT_outsideIGM_phot[kEbmv_] = sumT_phot - sumT_IGM;
      m_sumS_outsideIGM_phot[kEbmv_] = sumS_phot - sumS_IGM;
    } else {
      sumCross_phot += m_sumCross_outsideIGM_phot[kEbmv_];
      sumT_phot += m_sumT_outsideIGM_phot[kEbmv_];
      sumS_phot += m_sumS_outsideIGM_phot[kEbmv_];
    }
  }

  // add photometric terms to spectroscopic ones
  fitResult.sumCross += sumCross_phot;
  fitResult.sumT += sumT_phot;
  fitResult.sumS += sumS_phot;
}

Float64 COperatorTemplateFittingPhot::EstimateLikelihoodCstLog() const {

  Float64 cstlog = COperatorTemplateFitting::EstimateLikelihoodCstLog();
  auto spectrum = m_spectra[0];
  Float64 sumLogNoise = 0.0;
  const auto &photData = spectrum->GetPhotData();
  for (const auto &b : *m_photBandCat) {
    const std::string &bandName = b.first;
    sumLogNoise += log(photData->GetFluxErr(bandName) * m_weight);
  }
  cstlog -= m_photBandCat->size() * 0.5 * log(2 * M_PI) - sumLogNoise;

  return cstlog;
}

// used to output the model photometric values
TPhotVal COperatorTemplateFittingPhot::getIntegratedFluxes(Float64 amp) const {

  TPhotVal modelPhotValue;

  // compute photometric Leastquare
  for (auto bandName : m_sortedBandNames) {
    const auto &photBand = m_photBandCat->at(bandName);
    const auto &flux =
        m_templateRebined_phot.at(bandName).GetFluxAxis().GetSamplesVector();

    // integrate flux
    modelPhotValue[bandName] = photBand.IntegrateFlux(flux) * amp;
  }
  return modelPhotValue;
}
