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
    THROWG(INTERNAL_ERROR, "Photometry not supported with multiobs");
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
    THROWG(INTERNAL_ERROR, "Photometric band "
                           "transmision not available");
  if (m_photBandCat->empty())
    THROWG(INTERNAL_ERROR, "Empty photometric band transmission");

  if (m_spectra[0]->GetPhotData() == nullptr)
    THROWG(INTERNAL_ERROR, "photometric data not available in spectrum");

  const auto &dataNames = m_spectra[0]->GetPhotData()->GetNameList();
  for (const auto &bandName : m_photBandCat->GetNameList())
    if (std::find(dataNames.cbegin(), dataNames.cend(), bandName) ==
        dataNames.cend())
      THROWG(INTERNAL_ERROR,
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
      // status = nStatus_NoOverlap;
      THROWG(OVERLAPFRACTION_NOTACCEPTABLE,
             Formatter() << "tpl overlap too small: " << overlapFraction);
    }
  }
}

void COperatorTemplateFittingPhot::InitIsmIgmConfig(
    Float64 redshift,
    const std::shared_ptr<const CSpectrumFluxCorrectionCalzetti>
        &ismCorrectionCalzetti,
    const std::shared_ptr<const CSpectrumFluxCorrectionMeiksin>
        &igmCorrectionMeiksin,
    Int32 EbmvListSize) {

  COperatorTemplateFitting::InitIsmIgmConfig(
      redshift, ismCorrectionCalzetti, igmCorrectionMeiksin, EbmvListSize);

  // init ism on all rebined photometric templates
  for (auto &band : m_templateRebined_phot)
    band.second.InitIsmIgmConfig(redshift, ismCorrectionCalzetti,
                                 igmCorrectionMeiksin);

  m_sumCross_outsideIGM_phot.resize(EbmvListSize);
  m_sumS_outsideIGM_phot.resize(EbmvListSize);
  m_sumT_outsideIGM_phot.resize(EbmvListSize);
}

bool COperatorTemplateFittingPhot::CheckLyaIsInCurrentRange(
    const TFloat64Range &currentRange) const {

  bool ret = COperatorTemplateFitting::CheckLyaIsInCurrentRange(currentRange);
  for (const auto &band : m_templateRebined_phot)
    ret = ret || COperatorTemplateFitting::CheckLyaIsInCurrentRange(
                     band.second.GetLambdaRange());
  return ret;
}

bool COperatorTemplateFittingPhot::ApplyMeiksinCoeff(Int32 meiksinIdx,
                                                     Int32 spcIndex) {
  bool ret = COperatorTemplateFitting::ApplyMeiksinCoeff(meiksinIdx, spcIndex);
  for (auto &band : m_templateRebined_phot)
    ret = ret || band.second.ApplyMeiksinCoeff(
                     meiksinIdx); // will return true if at list igm is applied
                                  // on one template
  return ret;
}

bool COperatorTemplateFittingPhot::ApplyDustCoeff(Int32 kEbmv, Int32 spcIndex) {
  bool ret = COperatorTemplateFitting::ApplyDustCoeff(kEbmv, spcIndex);
  for (auto &band : m_templateRebined_phot)
    ret = ret || band.second.ApplyDustCoeff(kEbmv);
  return ret;
}

void COperatorTemplateFittingPhot::ApplyAmplitude(Float64 amplitude,
                                                  Int32 spcIndex) {
  COperatorTemplateFitting::ApplyAmplitude(amplitude, spcIndex);
  for (auto &band : m_templateRebined_phot)
    band.second.ApplyAmplitude(amplitude);
  return;
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

  Float64 sumCross_phot = 0.0;
  Float64 sumT_phot = 0.0;
  Float64 sumS_phot = 0.0;

  Float64 sumCross_IGM = 0.0;
  Float64 sumT_IGM = 0.0;
  Float64 sumS_IGM = 0.0;
  bool sumsIgmSaved = false;

  const auto &photData = m_spectra[0]->GetPhotData();

  // determine bands impacted by IGM
  const auto &Start = m_sortedBandNames.cbegin();
  const auto &End = m_sortedBandNames.cend();
  const auto IgmEnd =
      std::lower_bound(m_sortedBandNames.cbegin(), m_sortedBandNames.cend(),
                       RESTLAMBDA_LYA, [this](const std::string &s, Float64 v) {
                         return m_photBandCat->at(s).GetMinLambda() < v;
                       });

  auto EndLoop = m_option_igmFastProcessing && kM > 0 ? IgmEnd : End;

  // compute photometric Leastquare
  for (auto it = Start; it < EndLoop; it++) {

    if (m_option_igmFastProcessing && !sumsIgmSaved && it > IgmEnd) {
      // store intermediate sums for IGM range
      sumCross_IGM = fitResult.sumCross_phot;
      sumT_IGM = fitResult.sumT_phot;
      sumS_IGM = fitResult.sumS_phot;
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
      THROWG(INTERNAL_ERROR, Formatter()
                                 << "found invalid inverse variance : err2="
                                 << oneOverErr2 << ", for band=" << bandName);

    fitResult.sumCross_phot += d * integ_flux * oneOverErr2;
    fitResult.sumT_phot += integ_flux * integ_flux * oneOverErr2;
    fitResult.sumS_phot += fluxOverErr2;
  }

  if (m_option_igmFastProcessing) {
    if (kM == 0) {
      m_sumCross_outsideIGM_phot[kEbmv_] = sumCross_phot - sumCross_IGM;
      m_sumT_outsideIGM_phot[kEbmv_] = sumT_phot - sumT_IGM;
      m_sumS_outsideIGM_phot[kEbmv_] = sumS_phot - sumS_IGM;
    } else {
      fitResult.sumCross_phot += m_sumCross_outsideIGM_phot[kEbmv_];
      fitResult.sumT_phot += m_sumT_outsideIGM_phot[kEbmv_];
      fitResult.sumS_phot += m_sumS_outsideIGM_phot[kEbmv_];
    }
  }
  fitResult.sumCross += fitResult.sumCross_phot;
  fitResult.sumT += fitResult.sumT_phot;
  fitResult.sumS += fitResult.sumS_phot;
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
const TPhotVal COperatorTemplateFittingPhot::getIntegratedFluxes() {

  TPhotVal modelPhotValue;

  // compute photometric Leastquare
  for (auto bandName : m_sortedBandNames) {
    const auto &photBand = m_photBandCat->at(bandName);
    const auto &flux =
        m_templateRebined_phot.at(bandName).GetFluxAxis().GetSamplesVector();

    // integrate flux
    modelPhotValue[bandName] = photBand.IntegrateFlux(flux);
  }
  return modelPhotValue;
}