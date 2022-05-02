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
#include "RedshiftLibrary/linemodel/linemodelfitting.h"
#include "RedshiftLibrary/operator/spectraFluxResult.h"
using namespace NSEpic;

// TODO this should be a TExtremaResult constructor, using member initialization
// list
TLineModelResult::TLineModelResult(const CContinuumModelSolution &cms) {
  FittedTplName = cms.tplName;
  FittedTplAmplitude = cms.tplAmplitude;
  FittedTplAmplitudeError = cms.tplAmplitudeError;
  FittedTplMerit = cms.tplMerit;
  FittedTplMeritPhot = cms.tplMeritPhot;
  FittedTplEbmvCoeff = cms.tplEbmvCoeff;
  FittedTplMeiksinIdx = cms.tplMeiksinIdx;
  FittedTplRedshift = cms.tplRedshift;
  FittedTplDtm = cms.tplDtm;
  FittedTplMtm = cms.tplMtm;
  FittedTplLogPrior = cms.tplLogPrior;
  FittedTplpCoeffs = cms.pCoeffs;
}

void TLineModelResult::updateFromContinuumModelSolution(
    const CContinuumModelSolution &cms, bool all) {
  if (all) {
    FittedTplName = cms.tplName;
    FittedTplAmplitude = cms.tplAmplitude;
    FittedTplAmplitudeError = cms.tplAmplitudeError;
    FittedTplMerit = cms.tplMerit;
    FittedTplMeritPhot = cms.tplMeritPhot;
    FittedTplEbmvCoeff = cms.tplEbmvCoeff;
    FittedTplMeiksinIdx = cms.tplMeiksinIdx;
  }
  FittedTplRedshift = cms.tplRedshift;
  FittedTplDtm = cms.tplDtm;
  FittedTplMtm = cms.tplMtm;
  FittedTplLogPrior = cms.tplLogPrior;
  FittedTplpCoeffs = cms.pCoeffs;
}

void TLineModelResult::updateFromLineModelSolution(
    const CLineModelSolution &cms) {
  Elv = cms.EmissionVelocity;
  Alv = cms.AbsorptionVelocity;
}

void TLineModelResult::updateContinuumFromModel(
    const std::shared_ptr<const CLineModelFitting> &lmel) {
  FittedTplName = lmel->getFitContinuum_tplName();
  FittedTplAmplitude = lmel->getFitContinuum_tplAmplitude();
  FittedTplAmplitudeError = lmel->getFitContinuum_tplAmplitudeError();
  FittedTplMerit = lmel->getFitContinuum_tplMerit();
  FittedTplMeritPhot = lmel->getFitContinuum_tplMeritPhot();
  FittedTplEbmvCoeff = lmel->getFitContinuum_tplIsmEbmvCoeff();
  FittedTplMeiksinIdx = lmel->getFitContinuum_tplIgmMeiksinIdx();
}

void TLineModelResult::updateTplRatioFromModel(
    const std::shared_ptr<const CLineModelFitting> &lmel) {
  FittedTplratioName = lmel->getTplshape_bestTplName();
  FittedTplratioIsmCoeff = lmel->getTplshape_bestTplIsmCoeff();
  FittedTplratioAmplitude = lmel->getTplshape_bestAmplitude();
  FittedTplratioDtm = lmel->getTplshape_bestDtm();
  FittedTplratioMtm = lmel->getTplshape_bestMtm();
}

void TLineModelResult::updateFromModel(
    const std::shared_ptr<const CLineModelFitting> &lmel,
    const std::shared_ptr<const CLineModelResult> &lmresult,
    bool estimateLeastSquareFast, int idx, int i_2pass) {
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
  snrHa = lmresult->LineModelSolutions[idx].snrHa;
  lfHa = lmresult->LineModelSolutions[idx].lfHa;

  // store model OII SNR & Flux
  snrOII = lmresult->LineModelSolutions[idx].snrOII;
  lfOII = lmresult->LineModelSolutions[idx].lfOII;

  // store the model norm
  mTransposeM = lmel->EstimateMTransposeM();

  // scale marginalization correction
  Float64 corrScaleMarg = lmel->getScaleMargCorrection(); //
  CorrScaleMarg = corrScaleMarg;

  static Float64 cutThres = 3.0;
  Int32 nValidLines =
      lmresult->getNLinesOverCutThreshold(i_2pass, cutThres, cutThres);
  NLinesOverThreshold = nValidLines; // m/Float64(1+nValidLines);
  Float64 cumulStrongELSNR =
      lmel->getCumulSNRStrongEL(); // getStrongerMultipleELAmpCoeff(); //
  StrongELSNR = cumulStrongELSNR;

  TStringList strongELSNRAboveCut = lmel->getLinesAboveSNR(3.5);
  StrongELSNRAboveCut = strongELSNRAboveCut;

  Int32 nddl = lmel->GetNElements(); // get the total number of
  // elements in the model
  nddl = lmresult->LineModelSolutions[idx]
             .nDDL; // override nddl by the actual number of elements in
  // the fitted model
  NDof = lmel->m_Elements.GetModelNonZeroElementsNDdl();

  Float64 _bic =
      lmresult->ChiSquare[idx] + nddl * log(lmresult->nSpcSamples); // BIC
  // Float64 aic = m + 2*nddl; //AIC
  bic = _bic;
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

  OutsideLinesSTDFlux = lmel->getOutsideLinesSTD(1);
  OutsideLinesSTDError = lmel->getOutsideLinesSTD(2);

  if (OutsideLinesSTDError > 0.0) {
    Float64 ratioSTD = OutsideLinesSTDFlux / OutsideLinesSTDError;
    Float64 ratio_thres = 1.5;
    if (abs(ratioSTD) > ratio_thres || abs(ratioSTD) < 1. / ratio_thres) {
      Flag.warning(Flag.STDESTIMATION_NO_MATCHING,
                   Formatter()
                       << "  TLineModelResult::" << __func__
                       << ": STD estimations outside lines do not match: ratio="
                       << ratioSTD << ", flux-STD=" << OutsideLinesSTDFlux
                       << ", error-std=" << OutsideLinesSTDError);
    } else {
      Log.LogInfo("  Operator-Linemodel: STD estimations outside lines found "
                  "matching: ratio=%e, flux-STD=%e, error-std=%e",
                  ratioSTD, OutsideLinesSTDFlux, OutsideLinesSTDError);
    }
  } else {
    Flag.warning(Flag.STDESTIMATION_FAILED,
                 Formatter() << "  TLineModelResult::" << __func__
                             << ": unable to get STD estimations...");
  }
}

std::shared_ptr<const COperatorResult> LineModelExtremaResult::getCandidate(
    const int &rank, const std::string &dataset, bool firstpassResult) const {

  if (firstpassResult) {
    return getCandidateParent(rank, dataset);
  }
  if (dataset == "model_parameters") {
    return this->m_ranked_candidates[rank].second;
  } else if (dataset == "fitted_lines" || dataset == "fp_fitted_lines") {
    std::shared_ptr<const COperatorResult> cop =
        this->m_savedModelFittingResults[rank];
    return cop;
  } else if (dataset == "model")
    return this->m_savedModelSpectrumResults[rank];
  else if (dataset == "continuum")
    return this->m_savedModelContinuumSpectrumResults[rank];

  else
    throw GlobalException(UNKNOWN_ATTRIBUTE, "Unknown dataset");
}

const std::string &LineModelExtremaResult::getCandidateDatasetType(
    const std::string &dataset) const {
  if (dataset == "model_parameters")
    return this->m_ranked_candidates[0].second->getType();
  else if (dataset == "fitted_lines" || dataset == "fp_fitted_lines")
    return this->m_savedModelFittingResults[0]->getType();
  else if (dataset == "model")
    return this->m_savedModelSpectrumResults[0]->getType();
  else if (dataset == "continuum")
    return this->m_savedModelContinuumSpectrumResults[0]->getType();
  else
    throw GlobalException(UNKNOWN_ATTRIBUTE, "Unknown dataset");
}

bool LineModelExtremaResult::HasCandidateDataset(
    const std::string &dataset) const {
  return (dataset == "model_parameters" || dataset == "model" ||
          dataset == "fp_model_parameters" || dataset == "continuum" ||
          dataset == "fitted_lines" || dataset == "fp_fitted_lines");
}

std::shared_ptr<const COperatorResult>
LineModelExtremaResult::getCandidateParent(const int &rank,
                                           const std::string &dataset) const {
  if (dataset == "model_parameters") {
    return m_ranked_candidates[rank].second->ParentObject;
  }

  else
    throw GlobalException(UNKNOWN_ATTRIBUTE,
                          "Unknown dataset for parentObject");
}