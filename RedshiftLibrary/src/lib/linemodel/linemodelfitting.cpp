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
#include <climits>
#include <cmath>
#include <memory>
#include <numeric>

#include <boost/chrono/thread_clock.hpp>
#include <boost/format.hpp>
#include <boost/numeric/conversion/bounds.hpp>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_vector.h>

#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/continuum/irregularsamplingmedian.h"
#include "RedshiftLibrary/extremum/extremum.h"
#include "RedshiftLibrary/line/catalogsTplRatio.h"
#include "RedshiftLibrary/line/line.h"
#include "RedshiftLibrary/line/regulament.h"
#include "RedshiftLibrary/linemodel/element.h"
#include "RedshiftLibrary/linemodel/linemodelfitting.h"
#include "RedshiftLibrary/linemodel/rulesmanager.h"
#include "RedshiftLibrary/linemodel/tplcorrmanager.h"
#include "RedshiftLibrary/linemodel/tplratiomanager.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/processflow/autoscope.h"
#include "RedshiftLibrary/processflow/context.h"
#include "RedshiftLibrary/spectrum/template/template.h"

using namespace NSEpic;
using namespace std;

/**
 * \brief Prepares the state for Linemodel operation.
 * Loads the catalog.
 * Sets many state variables.
 * Sets the continuum either as a nocontinuum or a fromSpectrum.
 **/
CLineModelFitting::CLineModelFitting(
    const std::shared_ptr<COperatorTemplateFittingBase> &TFOperator,
    ElementComposition element_composition)
    : m_RestLineList(Context.getCLineMap()) {
  initParameters();

  m_inputSpcs = std::make_shared<std::vector<std::shared_ptr<const CSpectrum>>>(
      Context.getSpectra(m_useloglambdasampling));
  m_lambdaRanges = // std::make_shared<std::vector<std::shared_ptr<const
                   // TLambdaRange>>>(
      Context.getClampedLambdaRanges(m_useloglambdasampling);

  initMembers(TFOperator, element_composition);
  setLineRatioType(m_lineRatioType);
  if (m_lineRatioType == "rules")
    dynamic_cast<CRulesManager *>(m_lineRatioManager.get())->setRulesOption();
}

CLineModelFitting::CLineModelFitting(
    const std::shared_ptr<const CSpectrum> &template_,
    const TLambdaRange &lambdaRange,
    const std::shared_ptr<COperatorTemplateFittingBase> &TFOperator)
    : m_RestLineList(Context.getCLineMap()) {
  m_inputSpcs =
      std::make_shared<std::vector<std::shared_ptr<const CSpectrum>>>();

  m_inputSpcs->push_back(template_);
  m_lambdaRanges.push_back(std::make_shared<const TLambdaRange>(lambdaRange));
  initParameters();
  // override ortho specific parameters
  m_fittingmethod = "hybrid";
  // temporary options override to be removed when full tpl ortho is implemented
  m_lineRatioType = "rules";

  initMembers(TFOperator, ElementComposition::Default);
  setLineRatioType(m_lineRatioType);

  dynamic_cast<CRulesManager *>(m_lineRatioManager.get())->setRulesOption("no");
  setContinuumComponent("fromSpectrum");
}

void CLineModelFitting::initParameters() {
  CAutoScope autoscope(Context.m_ScopeStack, "lineModel");

  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();
  m_fittingmethod = ps->GetScoped<std::string>("fittingMethod");
  m_enableAmplitudeOffsets = ps->GetScoped<bool>("ampOffsetFit");
  m_enableLbdaOffsets = ps->GetScoped<bool>("lbdaOffsetFit");

  m_lineRatioType = ps->GetScoped<std::string>("lineRatioType");

  if (Context.GetCurrentMethod() == "lineModelSolve") {
    m_opt_firstpass_fittingmethod =
        ps->GetScoped<std::string>("firstPass.fittingMethod");
    m_opt_secondpass_fittingmethod = m_fittingmethod;
    m_opt_firstpass_forcedisableMultipleContinuumfit =
        ps->GetScoped<bool>("firstPass.multipleContinuumFitDisable");
  }

  std::string continuumComponent =
      ps->GetScoped<std::string>("continuumComponent");
  if (continuumComponent == "tplFit" || continuumComponent == "tplFitAuto") {
    m_opt_firstpass_forcedisableMultipleContinuumfit =
        ps->GetScoped<bool>("firstPass.multipleContinuumFitDisable");
    m_useloglambdasampling = ps->GetScoped<bool>("useLogLambdaSampling");
    m_opt_fitcontinuum_maxN = ps->GetScoped<Int32>("continuumFit.count");
  }
}

void CLineModelFitting::initMembers(
    const std::shared_ptr<COperatorTemplateFittingBase> &TFOperator,
    ElementComposition element_composition) {
  m_nbObs = m_inputSpcs->size();
  m_curObs = std::make_shared<Int32>(0);
  m_nominalWidthDefault = 13.4; // euclid 1 px
  m_continuumFitValues = std::make_shared<CTplModelSolution>();
  m_models = std::make_shared<std::vector<CSpectrumModel>>();
  if (element_composition == ElementComposition::Default &&
      (m_lineRatioType == "tplRatio" || m_lineRatioType == "tplCorr"))
    element_composition = ElementComposition::EmissionAbsorption;
  m_ElementsVector = std::make_shared<CLMEltListVector>(
      m_lambdaRanges, m_curObs, m_RestLineList, element_composition);
  for (*m_curObs = 0; *m_curObs < m_nbObs; (*m_curObs)++) {
    Log.LogDetail(Formatter() << "    model: Continuum winsize found is "
                              << std::fixed << std::setprecision(2)
                              << getSpectrum().GetMedianWinsize() << " A");
    m_models->push_back(CSpectrumModel(
        std::make_shared<CLineModelElementList>(getElementList()),
        getSpectrumPtr(), m_RestLineList, m_continuumFitValues, TFOperator,
        *m_curObs));
  }

  m_continuumManager = std::make_shared<CContinuumManager>(
      m_models, m_continuumFitValues, m_curObs);

  SetFittingMethod(m_fittingmethod, m_enableAmplitudeOffsets,
                   m_enableLbdaOffsets);
  SetLSF();
  LogCatalogInfos();

  // TODO restore check the continuum flux axis for NaN
}
// hook

void CLineModelFitting::logParameters() {
  Log.LogDetail(Formatter() << "m_pass" << m_pass);
  Log.LogDetail(Formatter()
                << " m_enableAmplitudeOffsets" << m_enableAmplitudeOffsets);
  Log.LogDetail(Formatter() << " m_LambdaOffsetMin" << m_LambdaOffsetMin);
  Log.LogDetail(Formatter() << " m_LambdaOffsetMax" << m_LambdaOffsetMax);
  Log.LogDetail(Formatter() << " m_LambdaOffsetStep" << m_LambdaOffsetStep);
  Log.LogDetail(Formatter()
                << " m_opt_firstpass_forcedisableMultipleContinuumfit"
                << m_opt_firstpass_forcedisableMultipleContinuumfit);
  Log.LogDetail(Formatter() << "m_opt_firstpass_fittingmethod "
                            << m_opt_firstpass_fittingmethod);
  Log.LogDetail(Formatter() << "m_opt_secondpass_fittingMethod"
                            << m_opt_secondpass_fittingmethod);

  Log.LogDetail(Formatter() << "nominalWidthDefault=" << m_nominalWidthDefault);

  Log.LogDetail(Formatter() << "fittingMethod=" << m_fittingmethod);

  Log.LogDetail(Formatter() << "lineRatioType=" << m_lineRatioType);

  // Log.LogDetail(Formatter()<<"tplCatalog="<<m_tplCatalog);
  // Log.LogDetail(Formatter()<<"tplCategoryList="<<m_tplCategoryList);

  //  Log.LogDetail(Formatter()<<"fitContinuum_tplFitPolyCoeffs="<<m_fitContinuum_tplFitPolyCoeffs);
  //  // only used with
  // m_fitContinuum_option==2 for now
  Log.LogDetail(Formatter() << "forcedisableMultipleContinuumfit="
                            << m_forcedisableMultipleContinuumfit);
}

/**
 * @brief setPassMode
 * @param iPass
 * set the fitting parameters according the the iPass argument.
 * @return
 */
Int32 CLineModelFitting::setPassMode(Int32 iPass) {
  m_pass = iPass;
  if (iPass == 1) {

    m_forcedisableMultipleContinuumfit =
        m_opt_firstpass_forcedisableMultipleContinuumfit;
    SetFittingMethod(m_opt_firstpass_fittingmethod);
  }
  if (iPass == 2) {
    m_forcedisableMultipleContinuumfit = false;
    SetFittingMethod(m_opt_secondpass_fittingmethod, m_enableAmplitudeOffsets,
                     m_enableLbdaOffsets);
  }
  if (iPass == 3) {
    m_forcedisableMultipleContinuumfit = false;
  }

  m_lineRatioManager->setPassMode(iPass);

  return true;
}
Int32 CLineModelFitting::GetPassNumber() const { return m_pass; }

/**
 * \brief LogDetail the number of lines for each element, and their nominal
 *amplitudes.
 **/
void CLineModelFitting::LogCatalogInfos() {
  Log.LogDetail("\n");
  Log.LogDetail("LineModel Infos: %d elements", getElementParam().size());
  int iElts = 0;
  for (auto const &elt_param : getElementParam()) {

    Int32 nLines = elt_param->size();
    if (nLines < 1) {
      Log.LogDetail(Formatter()
                    << "LineModel ctlg: elt " << iElts << " ("
                    << CLine::ETypeString.at(elt_param->GetElementType())
                    << "): no lines");
    }
    for (Int32 index = 0; index != elt_param->size(); ++index) {
      std::string nominalAmpStr = "";
      nominalAmpStr = boost::str(boost::format("(nominal amp = %.4e)") %
                                 elt_param->GetNominalAmplitude(index));
      Log.LogDetail(Formatter()
                    << "LineModel ctlg: elt " << iElts << " ("
                    << CLine::ETypeString.at(elt_param->GetElementType())
                    << "): line" << index << "= "
                    << elt_param->GetLineName(index) << nominalAmpStr);
    }
    iElts++;
  }
  Log.LogDetail("\n");
}

/*
Change the actual value of redshift.
the continuum can be reinterpolate.
*/
void CLineModelFitting::setRedshift(Float64 redshift,
                                    bool reinterpolatedContinuum) {
  for (*m_curObs = 0; *m_curObs < m_nbObs; (*m_curObs)++) {
    getSpectrumModel().m_Redshift = redshift;
  }
  if (reinterpolatedContinuum) {
    m_continuumManager->reinterpolateContinuum(redshift);
  }
}

Int32 CLineModelFitting::getTplratio_count() const {
  return m_lineRatioManager->getTplratio_count();
}

TFloat64List CLineModelFitting::getTplratio_priors() const {
  return m_lineRatioManager->getTplratio_priors();
}

bool CLineModelFitting::initDtd() {
  //  m_dTransposeDLambdaRange = TLambdaRange(*(m_lambdaRange));
  *m_curObs = 0; // we choose arbitrarily first obs to check if dtd is already
                 // initialized
                 // TODO check statement above
  m_dTransposeDLambdaRange = getLambdaRange();
  if (isContinuumComponentTplfitxx())
    m_dTransposeD = EstimateDTransposeD("raw");
  else
    m_dTransposeD = EstimateDTransposeD("noContinuum");

  m_likelihood_cstLog = EstimateLikelihoodCstLog();
  return true;
}

void CLineModelFitting::prepareAndLoadContinuum(Int32 k, Float64 redshift) {
  if (getContinuumComponent() == "noContinuum")
    return;

  if (!isContinuumComponentTplfitxx()) {
    for (*m_curObs = 0; *m_curObs < m_nbObs; (*m_curObs)++) {
      getSpectrumModel().setContinuumToInputSpc();
    }
    return;
  }

  // the support has to be already computed
  // when LoadFitContinuum() is called
  for (*m_curObs = 0; *m_curObs < m_nbObs; (*m_curObs)++) {

    getSpectrumModel().initObserveGridContinuumFlux(
        getSpectrum().GetSampleCount());
    m_continuumManager->LoadFitContinuum(k, redshift);
  }
}

void CLineModelFitting::computeSpectrumFluxWithoutContinuum() {
  for (*m_curObs = 0; *m_curObs < m_nbObs; (*m_curObs)++) {
    getSpectrumModel().initModelWithContinuum();
  }
}

/**
 * \brief Prepares the context and fits the Linemodel to the spectrum,
 *returning the bestMerit of the fit. Prepare the continuum. Initialize the
 *model spectrum. Prepare the elements. Fit the amplitudes of each element
 *independently. Fit the amplitude of all elements together with iterative
 *solver: Nelder Mead Simplex. Fit the amplitude of all elements together with
 *linear solver: gsl_multifit_wlinear. Fit the amplitudes of each element
 *independently, unless there is overlap. Apply a continuum iterative
 *re-estimation with lines removed from the initial spectrum. Apply rules.
 *Create spectrum model. Return bestMerit.
 **/
Float64 CLineModelFitting::fit(Float64 redshift,
                               CLineModelSolution &modelSolution,
                               CTplModelSolution &continuumModelSolution,
                               Int32 contreest_iterations, bool enableLogging) {
  // initialize the model spectrum
  m_fitter->m_cont_reestim_iterations = contreest_iterations;

  setRedshift(redshift);

  *m_curObs = 0; // we choose arbitrarily first obs to check if dtd is already
                 // initialized
  if (m_dTransposeDLambdaRange != getLambdaRange())
    initDtd();

  Int32 ntplratio = m_lineRatioManager->prepareFit(
      redshift); // multiple fitting steps for lineRatioType=tplratio/tplratio
  Int32 nContinuum = 1;
  Int32 savedIdxContinuumFitted = -1; // for continuum tplfit
  if (isContinuumComponentTplfitxx() && !m_forcedisableMultipleContinuumfit)
    nContinuum = m_opt_fitcontinuum_maxN;
  // 'on the fly' initialization
  Float64 bestMerit = INFINITY;
  Float64 bestMeritPrior = 0.0;

  for (Int32 k = 0; k < nContinuum; k++) {

    Float64 _merit = INFINITY;
    Float64 _meritprior = 0.; // only relevant for "tplRatio"
    prepareAndLoadContinuum(k, redshift);
    if (getContinuumComponent() != "noContinuum")
      computeSpectrumFluxWithoutContinuum();

    for (Int32 itratio = 0; itratio < ntplratio; itratio++) {

      if (m_lineRatioManager->init(redshift, itratio))
        continue;

      m_fitter->fit(redshift);

      std::string bestTplratioName = undefStr;

      _merit = m_lineRatioManager->computeMerit(itratio);

      if (bestMerit + bestMeritPrior > _merit + _meritprior) {
        bestMerit = _merit;
        bestMeritPrior = _meritprior;
        savedIdxContinuumFitted = k;
        Int32 modelSolutionLevel =
            m_lineRatioType == "rules" ? Int32(enableLogging) : 0;
        modelSolution = GetModelSolution(modelSolutionLevel);
        continuumModelSolution =
            m_continuumManager->GetContinuumModelSolutionCopy();

        m_lineRatioManager->saveResults(itratio);
      }
      if (getContinuumComponent() == "noContinuum") {
        for (*m_curObs = 0; *m_curObs < m_models->size(); (*m_curObs)++)
          getSpectrumModel().reinitModel();
      }
    }
  }

  if (!enableLogging)
    return bestMerit;

  if (isContinuumComponentTplfitxx()) {
    if (m_fittingmethod != "svdlc" && nContinuum > 1) {
      // TODO savedIdxContinuumFitted=-1 if lineRatioType!=tplratio
      for (*m_curObs = 0; *m_curObs < m_models->size(); (*m_curObs)++) {
        m_continuumManager->LoadFitContinuum(savedIdxContinuumFitted, redshift);
      }
    }
  }
  if (m_lineRatioType == "tplRatio") {
    m_lineRatioManager->resetToBestRatio(redshift);
    Int32 modelSolutionLevel = Int32(enableLogging);
    modelSolution = GetModelSolution(modelSolutionLevel);
    continuumModelSolution =
        m_continuumManager->GetContinuumModelSolutionCopy();
  }
  return bestMerit;
}

void CLineModelFitting::SetFittingMethod(const std::string &fitMethod,
                                         bool enableAmplitudeOffsets,
                                         bool enableLambdaOffsetsFit) {
  m_fittingmethod = fitMethod;
  *m_curObs = 0; // dummy implementation for svdlc and svdlcp2
  m_fitter = CAbstractFitter::makeFitter(
      fitMethod, m_ElementsVector, m_inputSpcs, m_lambdaRanges, m_models,
      m_RestLineList, m_continuumManager, m_curObs, enableAmplitudeOffsets,
      enableLambdaOffsetsFit);
  for (*m_curObs = 0; *m_curObs < m_nbObs; (*m_curObs)++) {
    getSpectrumModel().m_enableAmplitudeOffsets = enableAmplitudeOffsets;
  }
}

void CLineModelFitting::setLineRatioType(const std::string &lineRatioType) {
  m_lineRatioManager = CLineRatioManager::makeLineRatioManager(
      lineRatioType, m_ElementsVector, m_models, m_inputSpcs, m_lambdaRanges,
      m_continuumManager, m_RestLineList, m_fitter, m_curObs);
}

void CLineModelFitting::SetAbsLinesLimit(Float64 limit) {
  for (auto &elt_param : getElementParam()) {
    elt_param->SetAbsLinesLimit(limit);
  }
}

/**
 * \brief Creates and returns a Mask with 0 in the lines support, 1 under the
 *lines
 **/
CMask CLineModelFitting::getOutsideLinesMask() const {
  // TODO temp basic impl -> #8796
  *m_curObs = 0;
  // initialize the model spectrum
  const CSpectrumSpectralAxis &spectralAxis = getSpectrum().GetSpectralAxis();
  CMask _mask(spectralAxis.GetSamplesCount(), 1);

  TInt32List validEltsIdx = getElementList().GetModelValidElementsIndexes();
  TInt32List supportIdxes = getElementList().getSupportIndexes(validEltsIdx);

  // setting masks
  for (auto i : supportIdxes)
    _mask[i] = 0;

  return _mask;
}

/**
 * \brief Estimates the STD outside the lines for the observed-model spectrum
 * NB: supposes the spectrum whithout continuum has a null mean value
 * input: which = 1: uses the spectrum flux continuum subtracted to compute
 *STD input: which = 2: uses the spectrum error to compute STD
 **/
Float64 CLineModelFitting::getOutsideLinesSTD(Int32 which) const {
  // TODO temp basic impl
  *m_curObs = 0;

  if (which != 1 && which != 2)
    THROWG(INTERNAL_ERROR, Formatter()
                               << "wrong argument, which (1 or 2): " << which);

  CMask _mask = getOutsideLinesMask();

  const CSpectrumSpectralAxis &spectralAxis = getSpectrum().GetSpectralAxis();
  Float64 sum2 = 0.0;
  Int32 nsum = 0;
  Int32 imin = spectralAxis.GetIndexAtWaveLength(getLambdaRange().GetBegin());
  Int32 imax = spectralAxis.GetIndexAtWaveLength(getLambdaRange().GetEnd());
  const auto &spcFluxAxisNoContinuum =
      getSpectrumModel().getSpcFluxAxisNoContinuum();
  const auto &ErrorNoContinuum = getSpectrum().GetErrorAxis();
  for (Int32 i = imin; i < imax; i++) {
    if (!_mask[i])
      continue;

    if (which == 1)
      sum2 += spcFluxAxisNoContinuum[i] * spcFluxAxisNoContinuum[i];
    else if (which == 2)
      sum2 += ErrorNoContinuum[i] * ErrorNoContinuum[i];
    nsum++;
  }

  if (!nsum)
    return NAN;
  return sqrt(sum2 / nsum);
}

Float64 CLineModelFitting::getLeastSquareContinuumMerit() const {

  Float64 fit = 0.0;
  for (*m_curObs = 0; *m_curObs < m_nbObs; (*m_curObs)++) {

    const CSpectrumSpectralAxis &spcSpectralAxis =
        getSpectrum().GetSpectralAxis();
    const CSpectrumFluxAxis &Yspc = getSpectrumModel().getSpcFluxAxis();
    const auto &ErrorNoContinuum = getSpectrum().GetErrorAxis();

    const CSpectrumFluxAxis &YCont = getSpectrumModel().getContinuumFluxAxis();
    Float64 diff = 0.0;

    Float64 imin =
        spcSpectralAxis.GetIndexAtWaveLength(getLambdaRange().GetBegin());
    Float64 imax =
        spcSpectralAxis.GetIndexAtWaveLength(getLambdaRange().GetEnd());
    for (Int32 j = imin; j < imax; j++) {
      diff = (Yspc[j] - YCont[j]);
      fit += (diff * diff) / (ErrorNoContinuum[j] * ErrorNoContinuum[j]);
    }
  }
  if (isContinuumComponentTplfitxx()) {
    fit += m_continuumManager->getFittedLogPrior();
  }
  return fit;
}

Float64 CLineModelFitting::getLeastSquareContinuumMeritFast() const {
  Float64 fit;

  fit = m_dTransposeD;

  if (!isContinuumComponentTplfitxx())
    return fit;

  Float64 term1 = m_continuumManager->getTerm1();
  Float64 term2 = m_continuumManager->getTerm2();

  fit += term1 + term2;
  fit += m_continuumManager->getFittedLogPrior();

  return fit;
}

/**
 * \brief Returns the number of spectral samples between lambdaRange.
 **/
Int32 CLineModelFitting::computeSpcNSamples() const {

  Int32 imin;
  Int32 imax;
  Int32 lambdaMinSpcIndex = -1;
  Int32 lambdaMaxSpcIndex = -1;

  Int32 nSamples = 0;

  for (*m_curObs = 0; *m_curObs < m_nbObs; (*m_curObs)++) {
    const CSpectrumSpectralAxis &spcSpectralAxis =
        getSpectrum().GetSpectralAxis();
    imin = spcSpectralAxis.GetIndexAtWaveLength(getLambdaRange().GetBegin());
    imax = spcSpectralAxis.GetIndexAtWaveLength(getLambdaRange().GetEnd());
    nSamples += abs(imax - imin);
  }

  return nSamples;
}

/**
 * \brief Returns the cumulative SNR under the Strong Emission Lines
 * 1. retrieve the lines support
 * 2. process each
 **/
std::pair<Float64, Float64> CLineModelFitting::getCumulSNRStrongEL() const {

  *m_curObs = 0; // TODO #8797
  // Retrieve all the strone emission lines supports in a list of range
  TInt32RangeList supportList;
  TBoolList isStrongList;
  TInt32List validEltsIdx = getElementList().GetModelValidElementsIndexes();
  for (Int32 iElts : validEltsIdx) {
    auto const &elt = getElementList()[iElts];
    auto const &elt_param = elt->getElementParam();
    if (!elt_param->IsEmission())
      continue;
    for (Int32 index = 0; index != elt->GetSize(); ++index) {
      if (elt->IsOutsideLambdaRangeLine(index))
        continue;
      auto const &line = elt_param->GetLines()[index];
      isStrongList.push_back(line.IsStrong());
      supportList.push_back(elt->getTheoreticalSupportSubElt(index));
    }
  }

  // merge overlapping ranges
  TInt32RangeList nonOverlappingSupportList;
  TBoolList nonOverlappingIsStrongList;
  TInt32List processedSupport;
  for (Int32 k = 0; k < supportList.size(); k++) {
    // skip if already fitted
    if (std::find(processedSupport.cbegin(), processedSupport.cend(), k) !=
        processedSupport.cend())
      continue;

    processedSupport.push_back(k);
    TInt32Range &support = supportList[k];

    for (Int32 l = k + 1; l < supportList.size(); l++) {
      // skip if already fitted
      if (std::find(processedSupport.cbegin(), processedSupport.cend(), l) !=
          processedSupport.cend())
        continue;

      // try if current range is bluer than l and overlaps ?
      Float64 const xinf = support.GetBegin();
      Float64 const xsup = support.GetEnd();
      Float64 const yinf = supportList[l].GetBegin();
      Float64 const ysup = supportList[l].GetEnd();
      Float64 const max = std::max(xinf, yinf);
      Float64 const min = std::min(xsup, ysup);
      if (max - min < 0) {
        processedSupport.push_back(l);
        support.SetBegin(std::min(xinf, yinf));
        support.SetEnd(std::max(xsup, ysup));
        isStrongList[k] = isStrongList[k] || isStrongList[l];
      }
    }

    nonOverlappingSupportList.push_back(support);
    nonOverlappingIsStrongList.push_back(isStrongList[k]);
  }

  // process SNR on the non overlapping ranges
  Float64 sumFlux = 0.0;
  Float64 sumSquaredErr = 0.0;
  Float64 sumStrongFlux = 0.0;
  Float64 sumStrongSquaredErr = 0.0;
  for (size_t k = 0; k != nonOverlappingSupportList.size(); ++k) {
    auto const &[sumFlux_onRange, sumSquaredErr_onRange] =
        getSNROnRange(nonOverlappingSupportList[k]);
    sumFlux += sumFlux_onRange;
    sumSquaredErr += sumSquaredErr_onRange;
    if (nonOverlappingIsStrongList[k]) {
      sumStrongFlux += sumFlux_onRange;
      sumStrongSquaredErr += sumSquaredErr_onRange;
    }
  }

  Float64 const SNR = sumFlux / std::sqrt(sumSquaredErr);
  Float64 const StrongSNR = sumStrongFlux / std::sqrt(sumStrongSquaredErr);

  return std::make_pair(SNR, StrongSNR);
}

/**
 * \brief Returns the SNR on the idxRange
 **/
std::pair<Float64, Float64>
CLineModelFitting::getSNROnRange(TInt32Range idxRange) const {
  Int32 n = idxRange.GetEnd() - idxRange.GetBegin() + 1;
  if (idxRange.GetLength() < 1)
    return std::make_pair(0.0, 0.0);

  const CSpectrumFluxAxis &Ymodel =
      getSpectrumModel().GetModelSpectrum().GetFluxAxis();
  const auto &ErrorNoContinuum = getSpectrum().GetErrorAxis();
  const auto &ContinuumFluxAxis = getSpectrumModel().getContinuumFluxAxis();

  Float64 sumF = 0.0;
  Float64 sumM = 0.0;
  for (Int32 idx = idxRange.GetBegin(); idx <= idxRange.GetEnd(); idx++) {
    Float64 const flux =
        Ymodel[idx] - ContinuumFluxAxis[idx]; // using only the no-continuum
                                              // component to estimate SNR
    sumF += flux;
    sumM += ErrorNoContinuum[idx] * ErrorNoContinuum[idx];
  }

  return std::make_pair(std::abs(sumF), sumM);
}

/*
Reset all the model value to the previous solution found.
only used by linemeas throug getFittedModelWithoutContinuum to get
linemeas_model
*/
void CLineModelFitting::LoadModelSolution(
    const CLineModelSolution &modelSolution) {

  setRedshift(modelSolution.Redshift, false);

  // reset before loading
  for (auto param_ptr : m_ElementsVector->getElementParam()) {
    param_ptr->resetFittingParams();
    param_ptr->m_VelocityAbsorption = NAN;
    param_ptr->m_VelocityEmission = NAN;
  }
  m_ElementsVector->resetLambdaOffsets();

  // should also reset nominal amplitudes...
  // but not touched without using template-ratio

  TBoolList element_done(getElementParam().size(), false);
  for (Int32 iRestLine = 0; iRestLine < m_RestLineList.size(); iRestLine++) {
    Int32 eIdx = modelSolution.ElementId[iRestLine];
    if (eIdx == undefIdx)
      continue; // TODO should throw exception here
    if (modelSolution.OutsideLambdaRange[iRestLine])
      continue;
    Int32 line_id = modelSolution.lineId[iRestLine];

    auto const &elt_param = getElementParam()[eIdx];
    Int32 elt_line_index = elt_param->getLineIndex(line_id);
    if (elt_line_index == undefIdx)
      continue; // or throw an exception ?

    elt_param->setFittedAmplitude(
        elt_line_index, modelSolution.Amplitudes[iRestLine],
        modelSolution.AmplitudesUncertainties[iRestLine]);
    m_ElementsVector->getElementParam()[eIdx]->setLambdaOffset(
        elt_line_index, modelSolution.Offset[iRestLine]);

    if (element_done[eIdx])
      continue;

    elt_param->setVelocity(modelSolution.Velocity[iRestLine]);
    elt_param->SetFittingGroupInfo(modelSolution.fittingGroupInfo[iRestLine]);
    if (m_enableAmplitudeOffsets) {
      TPolynomCoeffs contPolynomCoeffs = {
          modelSolution.continuum_pCoeff0[iRestLine],
          modelSolution.continuum_pCoeff1[iRestLine],
          modelSolution.continuum_pCoeff2[iRestLine]};
      elt_param->SetPolynomCoeffs(std::move(contPolynomCoeffs));
    }
    element_done[eIdx] = true;
  }

  if (!std::isnan(modelSolution.LyaWidthCoeff) or
      !std::isnan(modelSolution.LyaAlpha) or
      !std::isnan(modelSolution.LyaDelta)) {

    std::string lyaTag = linetags::lya_em;
    auto const [idxLyaE, _] = m_ElementsVector->findElementIndex(lyaTag);
    if (idxLyaE != undefIdx)
      getElementParam()[idxLyaE]->SetAsymfitParams({modelSolution.LyaWidthCoeff,
                                                    modelSolution.LyaAlpha,
                                                    modelSolution.LyaDelta});
  }

  if (modelSolution.LyaIgm != undefIdx) {
    auto const indices_Igm = m_ElementsVector->getIgmLinesIndices();
    if (!indices_Igm.empty())
      for (auto const &[elt_idx, _] : indices_Igm)
        getElementParam()[elt_idx]->SetSymIgmParams(
            {modelSolution.LyaIgm, modelSolution.Redshift});
  }

  for (*m_curObs = 0; *m_curObs < m_inputSpcs->size(); (*m_curObs)++) {
    const CSpectrumSpectralAxis &spectralAxis = getSpectrum().GetSpectralAxis();
    for (Int32 iElts = 0; iElts < getElementList().size(); iElts++) {
      getElementList()[iElts]->computeOutsideLambdaRange();

      if (!getElementList()[iElts]->IsOutsideLambdaRange())

        getElementList()[iElts]->prepareSupport(
            spectralAxis, modelSolution.Redshift, getLambdaRange());
    }
  }

  return;
}

/**
 * \brief Returns a CLineModelSolution object populated with the current
 *solutions.
 **/
// this is not really a const method as spectrum model(s) have to be modified
// (cf CSpectrumModel::getContinuumUncertainty)
CLineModelSolution CLineModelFitting::GetModelSolution(Int32 opt_level) {
  Int32 s = m_RestLineList.size();
  CLineModelSolution modelSolution(m_RestLineList);

  modelSolution.EmissionVelocity =
      m_ElementsVector->getElementParam()[0]->m_VelocityEmission;
  modelSolution.AbsorptionVelocity =
      m_ElementsVector->getElementParam()[0]->m_VelocityAbsorption;

  // For some quantities it is more simple to get them from the first
  // observation objects There could be refactor but it can be complicated for
  // no gain of clarity or robustness
  *m_curObs = 0;
  modelSolution.Redshift = getSpectrumModel().m_Redshift;
  const CLineModelElementList &firstEltList = getElementList();

  TInt32List eIdx_oii;
  TInt32List subeIdx_oii;
  TInt32List eIdx_ha;
  TInt32List subeIdx_ha;
  Float64 flux_oii = 0.0;
  Float64 fluxVar_oii = 0.0;
  Float64 flux_ha = 0.0;
  Float64 fluxVar_ha = 0.0;

  modelSolution.nDDL = m_ElementsVector->GetModelNonZeroElementsNDdl();

  for (Int32 iRestLine = 0; iRestLine < s; iRestLine++) {
    Int32 line_id = modelSolution.lineId[iRestLine];
    auto const &line = m_RestLineList.at(line_id);
    auto [eIdx, line_index] = m_ElementsVector->findElementIndex(line_id);
    modelSolution.ElementId[iRestLine] = eIdx;
    if (eIdx == undefIdx || line_index == undefIdx ||
        m_ElementsVector->isOutsideLambdaRangeLine(eIdx, line_index)) {
      continue; // data already set to its default values
    }

    Float64 amp = m_ElementsVector->getElementParam()[eIdx]
                      ->m_FittedAmplitudes[line_index];
    modelSolution.Amplitudes[iRestLine] = amp;
    Float64 ampError = m_ElementsVector->getElementParam()[eIdx]
                           ->m_FittedAmplitudesStd[line_index];
    modelSolution.AmplitudesUncertainties[iRestLine] = ampError;
    if (getLineRatioType() == "rules")
      modelSolution.SNR[iRestLine] = std::abs(amp) / ampError;

    modelSolution.LambdaObs[iRestLine] =
        firstEltList[eIdx]->GetObservedPosition(line_index,
                                                modelSolution.Redshift);
    modelSolution.Velocity[iRestLine] =
        m_ElementsVector->getElementParam()[eIdx]->getVelocity();
    modelSolution.Offset[iRestLine] =
        m_ElementsVector->getElementParam()[eIdx]->m_Offsets[line_index];

    if (opt_level) // brief, to save processing time, do not estimate fluxes
                   // and high level line properties
    {
      modelSolution.ResidualRMS[iRestLine] =
          m_fitter->getModelResidualRmsUnderElements({eIdx}, true);
      if (m_enableAmplitudeOffsets) {
        const auto &polynom_coeffs =
            m_ElementsVector->getElementParam()[eIdx]->m_ampOffsetsCoeffs;
        modelSolution.continuum_pCoeff0[iRestLine] = polynom_coeffs.a0;
        modelSolution.continuum_pCoeff1[iRestLine] = polynom_coeffs.a1;
        modelSolution.continuum_pCoeff2[iRestLine] = polynom_coeffs.a2;
      }

      auto const [cont, cont_std] =
          GetMeanContinuumUnderLine(eIdx, line_index, modelSolution.Redshift);
      modelSolution.CenterContinuumFlux[iRestLine] = cont;
      modelSolution.CenterContinuumFluxUncertainty[iRestLine] = cont_std;

      Float64 flux = NAN;
      Float64 fluxError = NAN;
      TInt32List eIdx_line(1, eIdx);
      TInt32List subeIdx_line(1, line_index);
      Int32 opt_cont_substract_abslinesmodel = 0;
      bool isEmission = false;
      if (line.GetType() == CLine::EType::nType_Emission) {
        opt_cont_substract_abslinesmodel = 1;
        isEmission = true;
      }
      if (!std::isnan(amp) && amp >= 0.0) {
        if (!isEmission) {
          Float64 const positive_cont = std::max(0.0, cont);
          ampError *= positive_cont;
          ampError += amp * cont_std; // add continuum error (will lead to a
                                      // non-null error for a null continuum)
          amp *= -positive_cont; // minus sign to get a negative flux for abs
        }

        for (*m_curObs = 0; *m_curObs < m_inputSpcs->size(); (*m_curObs)++) {
          const auto &eltList = getElementList();
          if (!eltList[eIdx]->IsOutsideLambdaRangeLine(line_index)) {
            auto const &[mu, sigma] =
                eltList[eIdx]->getObservedPositionAndLineWidth(
                    modelSolution.Redshift, line_index,
                    false); // do not apply Lya asym offset
            modelSolution.Sigmas[iRestLine] = sigma;
            const auto &profile =
                eltList[eIdx]->getElementParam()->getLineProfile(line_index);

            Float64 const lineFlux = profile->GetLineFlux(mu, sigma);

            flux = amp * lineFlux;
            fluxError = ampError * lineFlux;
            break;
          }
        }
        if (*m_curObs >= m_inputSpcs->size())
          THROWG(INTERNAL_ERROR,
                 Formatter() << "Failed finding a spectrum containing"
                             << line_index << " at " << modelSolution.Redshift);
      }
      auto [fluxDI, snrDI] = getFluxDirectIntegration(
          eIdx_line, subeIdx_line, opt_cont_substract_abslinesmodel);
      modelSolution.Flux[iRestLine] = flux;
      modelSolution.FluxUncertainty[iRestLine] = fluxError;
      modelSolution.FluxDirectIntegration[iRestLine] = fluxDI;
      modelSolution.FluxDirectIntegrationUncertainty[iRestLine] =
          std::abs(fluxDI) / snrDI;

      // sum Ha complex fluxes
      if (isEmission && (line.GetName() == linetags::halpha_em ||
                         line.GetName() == linetags::niia_em ||
                         line.GetName() == linetags::niib_em)) {
        eIdx_ha.push_back(eIdx);
        subeIdx_ha.push_back(line_index);
        if (flux > 0.0)
          flux_ha += flux;
        if (fluxError > 0.0)
          fluxVar_ha += fluxError * fluxError;

        if (eIdx_ha.size() == 3) {

          Int32 opt_cont_substract_abslinesmodel = 0;
          auto [fluxDI, snrDI] = getFluxDirectIntegration(
              eIdx_ha, subeIdx_ha, opt_cont_substract_abslinesmodel);
          modelSolution.snrHa_DI = snrDI;
          modelSolution.lfHa_DI = fluxDI > 0.0 ? log10(fluxDI) : -INFINITY;
          modelSolution.lfHa = flux_ha > 0.0 ? log10(flux_ha) : -INFINITY;
          modelSolution.snrHa = flux_ha / std::sqrt(fluxVar_ha);
        }
      }

      // sum OII doublet fluxes
      if (isEmission && (line.GetName() == linetags::oII3726_em ||
                         line.GetName() == linetags::oII3729_em)) {
        eIdx_oii.push_back(eIdx);
        subeIdx_oii.push_back(line_index);
        if (flux > 0.0)
          flux_oii += flux;
        if (fluxError > 0.0)
          fluxVar_oii += fluxError * fluxError;

        if (eIdx_oii.size() == 2) {
          Int32 opt_cont_substract_abslinesmodel = 0;
          auto [fluxDI, snrDI] = getFluxDirectIntegration(
              eIdx_oii, subeIdx_oii, opt_cont_substract_abslinesmodel);

          modelSolution.snrOII_DI = snrDI;
          modelSolution.lfOII_DI = fluxDI > 0 ? log10(fluxDI) : -INFINITY;
          modelSolution.lfOII = flux_oii > 0.0 ? log10(flux_oii) : -INFINITY;
          modelSolution.snrOII = flux_oii / std::sqrt(fluxVar_oii);
        }
      }
    }

    modelSolution.fittingGroupInfo[iRestLine] =
        m_ElementsVector->getElementParam()[eIdx]->m_fittingGroupInfo;
    modelSolution.OutsideLambdaRange[iRestLine] = false;
  }

  // retrieve Lya params if fitted
  std::string lyaTag = linetags::lya_em;
  auto const [idxLyaE, _] = m_ElementsVector->findElementIndex(lyaTag);
  if (idxLyaE != undefIdx) {
    TAsymParams params =
        m_ElementsVector->getElementParam()[idxLyaE]->GetAsymfitParams(0);
    modelSolution.LyaWidthCoeff = params.sigma;
    modelSolution.LyaAlpha = params.alpha;
    modelSolution.LyaDelta = params.delta;
    TSymIgmParams params_igm =
        m_ElementsVector->getElementParam()[idxLyaE]->GetSymIgmParams(0);
    modelSolution.LyaIgm = params_igm.m_igmidx;
  }

  std::unordered_set<std::string>
      strongELSNRAboveCut; // = getLinesAboveSNR(3.5);
  modelSolution.NLinesAboveSnrCut = strongELSNRAboveCut.size();

  return modelSolution;
}

void CLineModelFitting::SetLSF() {
  for (*m_curObs = 0; *m_curObs < m_nbObs; (*m_curObs)++) {

    const std::shared_ptr<const CLSF> &lsf = getSpectrum().GetLSF();

    if (lsf == nullptr) {
      THROWG(INTERNAL_ERROR,
             "Cannot enable LSF, LSF spectrum member is not initialized");
    } else if (!lsf->IsValid()) {
      THROWG(INTERNAL_ERROR,
             " Cannot enable LSF, LSF spectrum member is not valid");
    }
    for (Int32 j = 0; j < getElementList().size(); j++) {
      getElementList()[j]->SetLSF(
          lsf); // lsf has now a type to be used for width computations
    }
  }
  //*m_curObs = 0;
}

void CLineModelFitting::SetVelocityEmission(Float64 vel) {
  for (auto &lmep : m_ElementsVector->getElementParam()) {
    lmep->m_VelocityEmission = vel;
  }
}

void CLineModelFitting::setVelocityEmissionByGroup(Float64 vel,
                                                   const TInt32List &inds) {
  for (auto idxElt : inds)
    m_ElementsVector->getElementParam()[idxElt]->m_VelocityEmission = vel;
}

void CLineModelFitting::SetVelocityAbsorption(Float64 vel) {

  for (auto &lmep : m_ElementsVector->getElementParam()) {
    lmep->m_VelocityAbsorption = vel;
  }
}

void CLineModelFitting::setVelocityAbsorptionByGroup(Float64 vel,
                                                     const TInt32List &inds) {

  for (auto idxElt : inds)
    m_ElementsVector->getElementParam()[idxElt]->m_VelocityAbsorption = vel;
}

Float64 CLineModelFitting::GetVelocityEmission() const {
  return m_ElementsVector->getElementParam()[0]->m_VelocityEmission;
}

Float64 CLineModelFitting::GetVelocityAbsorption() const {
  return m_ElementsVector->getElementParam()[0]->m_VelocityAbsorption;
}

Float64 CLineModelFitting::GetRedshift() const {
  return getSpectrumModel().m_Redshift;
}

/**
 * \brief this function returns the dtd value withing the wavelength range for
 *a given spcComponent
 *
 **/
// TODO rename this ! not a simple getter
Float64 CLineModelFitting::getDTransposeD() {

  *m_curObs = 0; // we choose arbitrarily first obs to check if dtd is already
                 // initialized
  if (m_dTransposeDLambdaRange != getLambdaRange()) {
    initDtd();
  }

  return m_dTransposeD;
}

/**
 * \brief this function returns the dtd value withing the wavelength range for
 *a given spcComponent
 *
 **/
// TODO rename this ! not a simple getter
Float64 CLineModelFitting::getLikelihood_cstLog() {

  *m_curObs = 0; // we choose arbitrarily first obs to check if dtd is already
                 // initialized
  if (m_dTransposeDLambdaRange != getLambdaRange()) {
    initDtd();
  }

  return m_likelihood_cstLog;
}

// below code could be moved to CSpectrum
/**
 * \brief this function estimates the dtd value withing the wavelength range
 **/
Float64
CLineModelFitting::EstimateDTransposeD(const std::string &spcComponent) const {

  Float64 dtd = 0.0;
  Float64 flux = 0.0;
  for (*m_curObs = 0; *m_curObs < m_nbObs; (*m_curObs)++) {

    const CSpectrumSpectralAxis &spcSpectralAxis =
        getSpectrum().GetSpectralAxis();
    const CSpectrumFluxAxis &Yspc = getSpectrumModel().getSpcFluxAxis();
    const CSpectrumFluxAxis &YspcNoContinuum =
        getSpectrumModel().getSpcFluxAxisNoContinuum();
    const auto &ErrorNoContinuum = getSpectrum().GetErrorAxis();

    Float64 imin =
        spcSpectralAxis.GetIndexAtWaveLength(getLambdaRange().GetBegin());
    Float64 imax =
        spcSpectralAxis.GetIndexAtWaveLength(getLambdaRange().GetEnd());
    for (Int32 j = imin; j < imax; j++) {
      if (spcComponent == "noContinuum")
        flux = YspcNoContinuum[j];
      else
        flux = Yspc[j];

      dtd += (flux * flux) / (ErrorNoContinuum[j] * ErrorNoContinuum[j]);
    }
    Log.LogDebug(Formatter()
                 << "CLineModelFitting::EstimateDTransposeD val = " << dtd);
  }
  return dtd;
}

/**
 * \brief this function estimates the mtm value withing the wavelength range
 **/
Float64 CLineModelFitting::EstimateMTransposeM()
    const // duplicate with getMTranposeMCumulative, except for return values
{
  Float64 mtm = 0.0;
  for (*m_curObs = 0; *m_curObs < m_nbObs; (*m_curObs)++) {

    const CSpectrumSpectralAxis &spcSpectralAxis =
        getSpectrum().GetSpectralAxis();
    const CSpectrumFluxAxis &spcFluxAxis =
        getSpectrumModel().GetModelSpectrum().GetFluxAxis();
    const auto &ErrorNoContinuum = getSpectrum().GetErrorAxis();

    Float64 diff = 0.0;

    Float64 imin =
        spcSpectralAxis.GetIndexAtWaveLength(getLambdaRange().GetBegin());
    Float64 imax =
        spcSpectralAxis.GetIndexAtWaveLength(getLambdaRange().GetEnd());
    for (Int32 j = imin; j < imax; j++) {
      diff = spcFluxAxis[j];
      mtm += (diff * diff) / (ErrorNoContinuum[j] * ErrorNoContinuum[j]);
    }
  }
  return mtm;
}

void CLineModelFitting::setContinuumComponent(std::string component) {
  m_continuumManager->setContinuumComponent(std::move(component));
}
/**
 * \brief this function estimates the likelihood_cstLog term withing the
 *wavelength range
 **/
Float64 CLineModelFitting::EstimateLikelihoodCstLog() const {

  Float64 cstLog = 0.0;
  for (*m_curObs = 0; *m_curObs < m_nbObs; (*m_curObs)++) {

    const CSpectrumSpectralAxis &spcSpectralAxis =
        getSpectrum().GetSpectralAxis();
    const auto &ErrorNoContinuum = getSpectrum().GetErrorAxis();

    Float64 sumLogNoise = 0.0;

    Int32 imin;
    Int32 imax;
    getLambdaRange().getClosedIntervalIndices(
        spcSpectralAxis.GetSamplesVector(), imin, imax);
    Int32 numDevs = std::abs(imax - imin + 1);
    for (Int32 j = imin; j <= imax; j++)
      sumLogNoise += log(ErrorNoContinuum[j]);

    cstLog += -numDevs * 0.5 * log(2 * M_PI) - sumLogNoise;
  }
  return cstLog;
}

std::pair<Float64, Float64>
CLineModelFitting::GetMeanContinuumUnderLine(Int32 eltIdx, Int32 line_index,
                                             Float64 redshift) {

  Float64 error = NAN;

  Float64 sumContinuumAll = 0.0;
  Float64 sumWeightAll = 0.0;
  Float64 sumSquaredWeightAll = 0.0;
  Float64 sumResidualAll = 0.0;
  Int32 nsumAll = 0;
  auto const &polynomCoeffs =
      m_ElementsVector->getElementParam()[eltIdx]->m_ampOffsetsCoeffs;
  for (*m_curObs = 0; *m_curObs < m_nbObs; (*m_curObs)++) {
    // TODO check this : it has been added because it caused Line position does
    // not belong to LSF range with lsf variable width
    if (m_ElementsVector->getElementList()[eltIdx]->IsOutsideLambdaRangeLine(
            line_index))
      continue;
    auto &model = getSpectrumModel();
    auto const &[indexRange, weights] =
        model.GetLineRangeAndProfile(eltIdx, line_index, redshift);
    auto [continuum_weighted_sum, sum_weight, sum_squared_weight] =
        model.GetContinuumWeightedSumInRange(indexRange, weights,
                                             polynomCoeffs);

    sumContinuumAll += continuum_weighted_sum;
    sumWeightAll += sum_weight;
    sumSquaredWeightAll += sum_squared_weight;

    auto const &[residual_sum, nsum] =
        model.getContinuumSquaredResidualInRange(indexRange);
    sumResidualAll += residual_sum;
    nsumAll += nsum;
  }
  if (sumWeightAll == 0.0)
    return std::make_pair(NAN, NAN);

  Float64 const continuum = sumContinuumAll / sumWeightAll;

  Float64 const std = nsumAll >= MIN_SAMPLE_NUMBER_CONTINUUMM_UNCERTAINTY
                          ? sqrt(sumResidualAll / nsumAll) *
                                sqrt(sumSquaredWeightAll) / sumWeightAll
                          : NAN;

  return std::make_pair(continuum, std);
}

std::pair<Float64, Float64> CLineModelFitting::getFluxDirectIntegration(
    const TInt32List &eIdx_list, const TInt32List &subeIdx_list,
    bool substract_abslinesmodel) const {

  Float64 sumFlux = 0;
  Float64 sumErr = 0;
  for (*m_curObs = 0; *m_curObs < m_nbObs; (*m_curObs)++) {
    auto const [flux, err] = getSpectrumModel().getFluxDirectIntegration(
        eIdx_list, subeIdx_list, substract_abslinesmodel, getLambdaRange());
    sumFlux += flux;
    sumErr += err;
  }
  if (sumErr <= 0.)
    return std::make_pair(NAN, NAN);

  Float64 snrdi = std::abs(sumFlux) / sqrt(sumErr);
  return std::make_pair(sumFlux, snrdi);
}

void CLineModelFitting::refreshAllModels() {

  for (*m_curObs = 0; *m_curObs < m_nbObs; (*m_curObs)++) {
    getSpectrumModel().refreshModel();
  }
}
