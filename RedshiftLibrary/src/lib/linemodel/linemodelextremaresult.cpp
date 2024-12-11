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
#include "RedshiftLibrary/linemodel/linemodelextremaresult.h"
#include "RedshiftLibrary/common/flag.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/linemodel/continuummanager.h"
#include "RedshiftLibrary/linemodel/linemodelfitting.h"
#include "RedshiftLibrary/linemodel/tplratiomanager.h"
#include "RedshiftLibrary/statistics/pdfcandidatesz.h"
using namespace NSEpic;

void TLineModelResult::updateFromContinuumModelSolution(
    std::shared_ptr<const CContinuumModelSolution> cms) {
  fittedContinuum = *cms;
  fittedContinuum.name =
      fittedContinuum.tplAmplitude ? fittedContinuum.name : "noContinuum";
}

void TLineModelResult::updateFromContinuumModelSolution(
    const CContinuumModelSolution &cms) {
  fittedContinuum = cms;
  fittedContinuum.name =
      fittedContinuum.tplAmplitude ? fittedContinuum.name : "noContinuum";
}

void TLineModelResult::updateFromLineModelSolution(
    const CLineModelSolution &cms) {
  Elv = cms.EmissionVelocity;
  Alv = cms.AbsorptionVelocity;
}

void TLineModelResult::updateTplRatioFromModel(
    const std::shared_ptr<const CTplratioManager> &ratioMgr) {
  FittedTplratioName = ratioMgr->getTplratio_bestTplName();
  FittedTplratioIsmCoeff = ratioMgr->getTplratio_bestTplIsmCoeff();
  FittedTplratioAmplitudeEm = ratioMgr->getTplratio_bestAmplitudeEm();
  FittedTplratioAmplitudeAbs = ratioMgr->getTplratio_bestAmplitudeAbs();
  FittedTplratioAmplitudeUncertaintyEm =
      ratioMgr->getTplratio_bestAmplitudeUncertaintyEm();
  FittedTplratioAmplitudeUncertaintyAbs =
      ratioMgr->getTplratio_bestAmplitudeUncertaintyAbs();
  FittedTplratioSNREm = std::abs(FittedTplratioAmplitudeEm) /
                        FittedTplratioAmplitudeUncertaintyEm;
  FittedTplratioSNRAbs = std::abs(FittedTplratioAmplitudeAbs) /
                         FittedTplratioAmplitudeUncertaintyAbs;
}

void TLineModelResult::updateFromModel(
    const std::shared_ptr<const CLineModelFitting> &lmel,
    const std::shared_ptr<const CLineModelResult> &lmresult,
    bool estimateLeastSquareFast, int idx) {
  Merit = lmresult->ChiSquare[idx];

  // LineModelSolutions
  Elv = lmel->GetVelocityEmission();
  Alv = lmel->GetVelocityAbsorption();

  if (!estimateLeastSquareFast) {
    MeritContinuum = lmel->getLeastSquareContinuumMerit();
  } else {
    MeritContinuum = lmel->getLeastSquareContinuumMeritFast();
  }

  // store model Ha SNR & Flux
  lfHa = lmresult->LineModelSolutions[idx].lfHa;
  if (lmel->getLineRatioType() == "rules")
    snrHa = lmresult->LineModelSolutions[idx].snrHa;
  lfHa_DI = lmresult->LineModelSolutions[idx].lfHa_DI;
  snrHa_DI = lmresult->LineModelSolutions[idx].snrHa_DI;

  // store model OII SNR & Flux
  lfOII = lmresult->LineModelSolutions[idx].lfOII;
  if (lmel->getLineRatioType() == "rules")
    snrOII = lmresult->LineModelSolutions[idx].snrOII;
  lfOII_DI = lmresult->LineModelSolutions[idx].lfOII_DI;
  snrOII_DI = lmresult->LineModelSolutions[idx].snrOII_DI;

  // store Lya fitting parameters
  LyaWidthCoeff = lmresult->LineModelSolutions[idx].LyaWidthCoeff;
  LyaAlpha = lmresult->LineModelSolutions[idx].LyaAlpha;
  LyaDelta = lmresult->LineModelSolutions[idx].LyaDelta;
  LyaIgm = lmresult->LineModelSolutions[idx].LyaIgm;

  // scale marginalization correction
  Float64 corrScaleMarg = lmel->getScaleMargCorrection(); //
  CorrScaleMarg = corrScaleMarg;

  Float64 static const cutThres = SNR_THRESHOLD_FOR_NLINESOVER;
  if (lmel->getLineRatioType() == "rules")
    NLinesOverThreshold =
        lmresult->getNLinesOverCutThreshold(idx, cutThres, cutThres);

  std::tie(ELSNR, StrongELSNR) = lmel->getCumulSNRStrongEL();

  StrongELSNRAboveCut = lmel->getLinesAboveSNR(3.5);

  NDof = lmel->getNonZeroElementsNDdl();

  Int32 const nddl = lmresult->LineModelSolutions[idx].nDDL;
  bic = lmresult->ChiSquare[idx] + nddl * log(lmresult->nSpcSamples); // BIC
  // Float64 aic = m + 2*nddl; //AIC

  // lmresult->bic = aic + (2*nddl*(nddl+1) )/(nsamples-nddl-1);
  // //AICc, better when nsamples small

  // compute continuum indexes
  // TODO VB is this useful/necessary now ? if there is a computation it should
  // be done before
  // NB AA commented to avoid adding spectrum to getFromModel arguments
  /*
  CContinuumIndexes continuumIndexes;
  ContinuumIndexes =
    continuumIndexes.getIndexes(spectrum, z);
  */

  // save the outsideLinesMask
  OutsideLinesMask = lmel->getOutsideLinesMask();

  std::tie(OutsideLinesResidualRMS, OutsideLinesInputStDevRMS) =
      lmel->getOutsideLinesRMS(OutsideLinesMask);

  if (std::isnan(OutsideLinesResidualRMS) ||
      std::isnan(OutsideLinesInputStDevRMS)) {
    Flag.warning(WarningCode::STD_ESTIMATION_FAILED,
                 "StDev estimation outside lines failed");
    return;
  }

  Float64 ratioSTD = OutsideLinesResidualRMS / OutsideLinesInputStDevRMS;
  Float64 ratio_thres = 1.5;
  if (abs(ratioSTD) > ratio_thres || abs(ratioSTD) < 1. / ratio_thres) {
    Flag.warning(WarningCode::ESTIMATED_STD_FAR_FROM_INPUT,
                 Formatter()
                     << "StDev estimation outside lines do not match: ratio = "
                     << ratioSTD
                     << ", residual RMS = " << OutsideLinesResidualRMS
                     << ", Input Noise RMS = " << OutsideLinesInputStDevRMS);
  } else {
    Log.LogInfo(Formatter()
                << "StDev estimations outside lines match: ratio = " << ratioSTD
                << ", residual RMS = " << OutsideLinesResidualRMS
                << ", Input Noise RMS = " << OutsideLinesInputStDevRMS);
  }
}

std::shared_ptr<const COperatorResult> LineModelExtremaResult::getCandidate(
    const int &rank, const std::string &dataset, bool firstpassResult) const {

  if (firstpassResult) {
    return getCandidateParent(rank, dataset);
  }
  if (dataset == "model_parameters" || dataset == "line_mask" ||
      dataset == "continuum_polynom")
    return this->m_ranked_candidates[rank].second;
  else if (dataset == "fitted_lines" || dataset == "fp_fitted_lines")
    return m_savedModelFittingResults[rank];
  else if (dataset == "model")
    return this->m_savedModelSpectrumResults[rank];
  else if (dataset == "continuum")
    return this->m_savedModelContinuumSpectrumResults[rank];
  else if (dataset == "PhotometricModel")
    return this->m_modelPhotValues[rank];
  else
    THROWG(ErrorCode::UNKNOWN_ATTRIBUTE, "Unknown dataset");
}

const std::string &LineModelExtremaResult::getCandidateDatasetType(
    const std::string &dataset) const {
  if (dataset == "model_parameters" || dataset == "line_mask" ||
      dataset == "continuum_polynom")
    return this->m_ranked_candidates[0].second->getType();
  else if (dataset == "fitted_lines" || dataset == "fp_fitted_lines")
    return this->m_savedModelFittingResults[0]->getType();
  else if (dataset == "model")
    return this->m_savedModelSpectrumResults[0]->getType();
  else if (dataset == "continuum")
    return this->m_savedModelContinuumSpectrumResults[0]->getType();
  else if (dataset == "PhotometricModel")
    return this->m_modelPhotValues[0]->getType();
  else
    THROWG(ErrorCode::UNKNOWN_ATTRIBUTE, "Unknown dataset");
}

bool LineModelExtremaResult::HasCandidateDataset(
    const std::string &dataset) const {
  return (dataset == "model_parameters" || dataset == "model" ||
          dataset == "continuum" || dataset == "fitted_lines" ||
          dataset == "fp_fitted_lines" || dataset == "line_mask" ||
          dataset == "continuum_polynom" || dataset == "PhotometricModel");
}

std::shared_ptr<const COperatorResult>
LineModelExtremaResult::getCandidateParent(const int &rank,
                                           const std::string &dataset) const {
  if (dataset == "model_parameters") {
    return m_ranked_candidates[rank].second->ParentObject;
  }

  else
    THROWG(ErrorCode::UNKNOWN_ATTRIBUTE, "Unknown dataset for parentObject");
}

TCandidateZbyRank LineModelExtremaResult::getCandidatesZByRank() {
  TCandidateZbyRank ret;
  for (auto &cand : m_ranked_candidates) {
    ret.push_back(std::make_pair(
        cand.first, std::dynamic_pointer_cast<TCandidateZ>(cand.second)));
  }
  return ret;
}
