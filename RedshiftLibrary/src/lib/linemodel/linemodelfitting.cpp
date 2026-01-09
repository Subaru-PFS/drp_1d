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

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/common/size.h"
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
#include "RedshiftLibrary/spectrum/LSFFactory.h"
#include "RedshiftLibrary/spectrum/template/template.h"

using namespace NSEpic;
using namespace std;

CLineModelFitting::CLineModelFitting(Int32 spectraIndex)
    : m_RestLineList(Context.getCLineMap()), m_spectraIndex(spectraIndex) {}

/**
 * \brief Prepares the state for Linemodel operation.
 * Loads the catalog.
 * Sets many state variables.
 * Sets the continuum either as a nocontinuum or a fromSpectrum.
 **/
CLineModelFitting::CLineModelFitting(
    const std::shared_ptr<COperatorContinuumFitting> &continuumFittingOperator,
    ElementComposition element_composition)
    : CLineModelFitting(Context.getSpectra().size()) {
  initParameters();

  m_inputSpcs = std::make_shared<std::vector<std::shared_ptr<const CSpectrum>>>(
      Context.getSpectra(m_useloglambdasampling));
  m_lambdaRanges = Context.getClampedLambdaRanges(m_useloglambdasampling);
  auto lineRatioType = CLineRatioManager::stringToType.at(
      Context.GetParameterStore()->GetScoped<std::string>("lineRatioType"));
  initMembers(continuumFittingOperator, lineRatioType, element_composition);
  setLineRatioManager(lineRatioType);
  if (isLineRatioRules())
    dynamic_cast<CRulesManager *>(m_lineRatioManager.get())->setRulesOption();
}

void CLineModelFitting::initParameters() {
  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();
  m_fittingmethod = ps->GetScoped<std::string>("fittingMethod");

  if (Context.GetCurrentMethod() == "lineModelSolve") {
    m_opt_firstpass_fittingmethod =
        ps->GetScoped<std::string>("firstPass.fittingMethod");
    m_opt_secondpass_fittingmethod = m_fittingmethod;
    m_opt_firstpass_forcedisableMultipleContinuumfit =
        ps->GetScoped<bool>("firstPass.multipleContinuumFitDisable");
  }

  std::set<std::string> const lbdaOffsetFitters{"svd", "hybrid", "lbfgsb"};
  if (lbdaOffsetFitters.find(m_fittingmethod) != lbdaOffsetFitters.end() ||
      lbdaOffsetFitters.find(m_opt_firstpass_fittingmethod) !=
          lbdaOffsetFitters.end()) {
    m_enableAmplitudeOffsets = ps->GetScoped<bool>("ampOffsetFit");
    m_enableLbdaOffsets = ps->GetScoped<bool>("lbdaOffsetFit");
  }

  TContinuumComponent continuumComponent(
      ps->GetScoped<std::string>("continuumComponent"));
  if (continuumComponent.isTplFitXXX()) {
    m_opt_firstpass_forcedisableMultipleContinuumfit =
        ps->GetScoped<bool>("firstPass.multipleContinuumFitDisable");
    m_useloglambdasampling = ps->GetScoped<bool>("useLogLambdaSampling");
    m_opt_fitcontinuum_maxN = ps->GetScoped<Int32>("continuumFit.count");
  }
}

void CLineModelFitting::initMembers(
    const std::shared_ptr<COperatorContinuumFitting> &continuumFittingOperator,
    CLineRatioManager::EType const &lineRatioType,
    ElementComposition element_composition) {

  m_nominalWidthDefault = 13.4; // euclid 1 px
  m_continuumFitValues = std::make_shared<CContinuumModelSolution>();
  m_models = std::make_shared<CSpcModelVector>(m_spectraIndex);
  if (element_composition == ElementComposition::Default &&
      (lineRatioType == CLineRatioManager::EType::tplRatio ||
       lineRatioType == CLineRatioManager::EType::tplCorr ||
       lineRatioType == CLineRatioManager::EType::ratioToFree))
    element_composition = ElementComposition::EmissionAbsorption;
  setElementsVector(lineRatioType, element_composition);
  for ([[maybe_unused]] auto &spcIndex : m_spectraIndex) {
    Log.LogDetail(Formatter() << "    model: Continuum winsize found is "
                              << std::fixed << std::setprecision(2)
                              << getSpectrum().GetMedianWinsize() << " A");
    m_models->push_back(CSpectrumModel(
        getElementList(), getSpectrumPtr(), m_RestLineList,
        m_continuumFitValues, continuumFittingOperator, m_spectraIndex.get()));
  }

  m_continuumManager = std::make_shared<CContinuumManager>(
      m_models, m_continuumFitValues, m_spectraIndex);

  SetFittingMethod(m_fittingmethod, m_enableAmplitudeOffsets,
                   m_enableLbdaOffsets);
  SetLSF();
  LogCatalogInfos();
}

void CLineModelFitting::reloadFor2ndPass() {

  auto lineRatioType = m_lineRatioManager->getStrictType();

  setElementsVector(lineRatioType, ElementComposition::Default);

  for ([[maybe_unused]] auto &spcIndex : m_spectraIndex) {
    m_models->setModelsElements(m_ElementsVector->getElementList());
  }

  SetLSF();
  LogCatalogInfos();
  setLineRatioManager(lineRatioType);
  if (isLineRatioRules())
    dynamic_cast<CRulesManager *>(m_lineRatioManager.get())->setRulesOption();
}

void CLineModelFitting::setElementsVector(
    CLineRatioManager::EType const &lineRatioType,
    ElementComposition const &element_composition) {
  // Here must pass lineRatioType as arg because is used before
  // m_lineRatioManager initialization
  m_ElementsVector = std::make_shared<CLMEltListVector>(
      m_spectraIndex, m_RestLineList, element_composition,
      m_enableAmplitudeOffsets);
}

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
  Log.LogDetail("LineModel Infos: %d elements", getElementsParams().size());
  int iElts = 0;
  for (auto const &elt_param : getElementsParams()) {

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
void CLineModelFitting::setRedshift(Float64 redshift) {
  for ([[maybe_unused]] auto &spcIndex : m_spectraIndex) {
    getSpectrumModel().m_Redshift = redshift;
  }
}

Int32 CLineModelFitting::getTplratio_count() const {
  return m_lineRatioManager->getTplratio_count();
}

TFloat64List CLineModelFitting::getTplratio_priors() const {
  return m_lineRatioManager->getTplratio_priors();
}

void CLineModelFitting::initDtd() {
  m_spectraIndex.setAtBegining();
  m_dTransposeDLambdaRange = getLambdaRange();
  auto const &component = isContinuumComponentFitter() ? "raw" : "noContinuum";
  m_dTransposeD = EstimateDTransposeD(component);
  m_likelihood_cstLog = EstimateLikelihoodCstLog(component);
}

void CLineModelFitting::prepareAndLoadContinuum(Int32 k, Float64 redshift) {
  if (isContinuumComponentNoContinuum()) {
    m_ElementsVector->setAllAbsLinesNullContinuum();
    return;
  }

  if (isContinuumComponentFromSpectrum()) {
    for ([[maybe_unused]] auto &spcIndex : m_spectraIndex) {
      getSpectrumModel().setContinuumToInputSpc();
    }
  } else {
    // the support has to be already computed
    // when LoadFitContinuum() is called
    for ([[maybe_unused]] auto &spcIndex : m_spectraIndex)
      getSpectrumModel().initObserveGridContinuumFlux(
          getSpectrum().GetSampleCount());
    m_continuumManager->LoadFitContinuum(k, redshift);
  }

  computeSpectrumFluxWithoutContinuum();

  if (isContinuumFittedToNull())
    m_ElementsVector->setAllAbsLinesNullContinuum();
  else
    m_ElementsVector->unsetAllAbsLinesNullContinuum();
}

void CLineModelFitting::computeSpectrumFluxWithoutContinuum() {
  for ([[maybe_unused]] auto &spcIndex : m_spectraIndex) {
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
                               CContinuumModelSolution &continuumModelSolution,
                               Int32 contreest_iterations, bool fullSolution) {
  // initialize the model spectrum
  m_fitter->m_cont_reestim_iterations = contreest_iterations;

  setRedshift(redshift);

  m_spectraIndex.setAtBegining(); // we choose arbitrarily first obs to check if
                                  // dtd is already initialized

  if (m_dTransposeDLambdaRange != getLambdaRange())
    initDtd();

  Int32 ntplratio = m_lineRatioManager->prepareFit(
      redshift); // multiple fitting steps for lineRatioType=tplratio
  Int32 nContinuum = 1;
  Int32 savedIdxContinuumFitted = -1; // for continuum tplfit
  if (isContinuumComponentTplFitXXX() && !m_forcedisableMultipleContinuumfit)
    nContinuum = m_opt_fitcontinuum_maxN;
  // 'on the fly' initialization
  Float64 bestMerit = INFINITY;
  Float64 bestMeritPrior = 0.0;

  for (Int32 k = 0; k < nContinuum; k++) {

    Float64 _merit = INFINITY;
    Float64 _meritprior = 0.; // only relevant for "tplRatio" and "ratioToFree"

    prepareAndLoadContinuum(k, redshift);

    for (Int32 itratio = 0; itratio < ntplratio; itratio++) {

      if (m_lineRatioManager->init(redshift, itratio))
        continue;

      m_fitter->fit(redshift);

      std::string bestTplratioName = undefStr;

      std::tie(_merit, _meritprior) = m_lineRatioManager->computeMerit(itratio);

      if (bestMerit + bestMeritPrior > _merit + _meritprior) {
        bestMerit = _merit;
        bestMeritPrior = _meritprior;
        savedIdxContinuumFitted = k;
        bool modelSolutionLevel = isLineRatioRules() ? fullSolution : false;
        m_lineRatioManager->saveResults(itratio);
        modelSolution = GetModelSolution(modelSolutionLevel);
        continuumModelSolution =
            m_continuumManager->GetContinuumModelSolutionCopy();
      }
      if (isContinuumComponentNoContinuum()) {
        m_models->reinitAllModels();
      }
    }
  }

  if (!fullSolution)
    return bestMerit;

  if (isContinuumComponentFitter()) {
    if (m_fittingmethod != "svdlc" && nContinuum > 1) {
      m_continuumManager->LoadFitContinuum(savedIdxContinuumFitted, redshift);
    }
  }
  if (isLineRatioTplRatio()) {
    m_lineRatioManager->resetToBestRatio(redshift);
    modelSolution = GetModelSolution(fullSolution);
    continuumModelSolution =
        m_continuumManager->GetContinuumModelSolutionCopy();
  }
  return bestMerit;
}

void CLineModelFitting::SetFittingMethod(const std::string &fitMethod,
                                         bool enableAmplitudeOffsets,
                                         bool enableLambdaOffsetsFit) {
  // NB dummy multiobs implementation (functional for one obs only) for svdlc
  // and svdlcp2
  m_fittingmethod = fitMethod;
  m_spectraIndex.setAtBegining(); // temporary multiobs implementation
  m_fitter = CAbstractFitter::makeFitter(
      fitMethod, m_ElementsVector, m_inputSpcs, m_lambdaRanges, m_models,
      m_RestLineList, m_continuumManager, m_spectraIndex,
      enableAmplitudeOffsets, enableLambdaOffsetsFit);
  m_models->setEnableAmplitudeOffsets(enableAmplitudeOffsets);
}

void CLineModelFitting::setLineRatioManager(
    CLineRatioManager::EType lineRatioType) {
  m_lineRatioManager = CLineRatioManager::makeLineRatioManager(
      lineRatioType, m_ElementsVector, m_models, m_inputSpcs, m_lambdaRanges,
      m_continuumManager, m_RestLineList, m_fitter, m_spectraIndex);
}

void CLineModelFitting::SetAbsLinesLimit(Float64 limit) {
  for (auto &elt_param : getElementsParams()) {
    elt_param->SetAbsLinesLimit(limit);
  }
}

/**
 * \brief Creates and returns a Mask with 0 in the lines support, 1 under the
 *lines
 **/
CMask CLineModelFitting::getOutsideLinesMask() const {
  // NB temporary basic implementation
  m_spectraIndex.setAtBegining();
  // initialize the model spectrum
  const CSpectrumSpectralAxis &spectralAxis = getSpectrum().GetSpectralAxis();
  CMask _mask(spectralAxis.GetSamplesCount(), 1);

  TInt32List validEltsIdx =
      getElementList().GetElementsIndicesInsideLambdaRange();
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
std::pair<Float64, Float64>
CLineModelFitting::getOutsideLinesRMS(CMask const &_mask) const {
  // NB dummy multiobs implementation (functional for one obs only)
  m_spectraIndex.setAtBegining(); // temporary multiobs implementation

  const CSpectrumSpectralAxis &spectralAxis = getSpectrum().GetSpectralAxis();
  Float64 sum2_flux = 0.0;
  Float64 sum2_error = 0.0;
  Int32 nsum = 0;
  auto const &[imin, imax] =
      getLambdaRange().getClosestInnerIndices(spectralAxis.GetSamplesVector());

  const auto &spcFluxAxisNoContinuum =
      getSpectrumModel().getSpcFluxAxisNoContinuum();
  const auto &ErrorNoContinuum = getSpectrum().GetErrorAxis();
  for (Int32 i = imin; i <= imax; i++) {
    if (!_mask[i])
      continue;

    sum2_flux += spcFluxAxisNoContinuum[i] * spcFluxAxisNoContinuum[i];
    sum2_error += ErrorNoContinuum[i] * ErrorNoContinuum[i];
    nsum++;
  }

  if (nsum < RMS_MIN_SAMPLE_NUMBER)
    return std::make_pair(NAN, NAN);
  return std::make_pair(sqrt(sum2_flux / nsum), sqrt(sum2_error / nsum));
}

Float64 CLineModelFitting::getLeastSquareContinuumMerit() const {

  Float64 fit = 0.0;
  for ([[maybe_unused]] auto &spcIndex : m_spectraIndex) {

    const CSpectrumSpectralAxis &spcSpectralAxis =
        getSpectrum().GetSpectralAxis();
    const CSpectrumFluxAxis &Yspc = getSpectrumModel().getSpcFluxAxis();
    const auto &Error = getSpectrum().GetErrorAxis();

    const CSpectrumFluxAxis &YCont = getSpectrumModel().getContinuumFluxAxis();
    Float64 diff = 0.0;

    auto const &[imin, imax] = getLambdaRange().getClosestInnerIndices(
        spcSpectralAxis.GetSamplesVector());

    for (Int32 j = imin; j <= imax; j++) {
      diff = (Yspc[j] - YCont[j]);
      fit += (diff * diff) / (Error[j] * Error[j]);
    }
  }
  if (isContinuumComponentFitter()) {
    fit += m_continuumManager->getFittedLogPrior();
  }
  return fit;
}

Float64 CLineModelFitting::getLeastSquareContinuumMeritFast() const {
  Float64 fit;

  fit = m_dTransposeD;

  if (!isContinuumComponentFitter())
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

  Int32 nSamples = 0;

  for ([[maybe_unused]] auto &spcIndex : m_spectraIndex) {
    const CSpectrumSpectralAxis &spcSpectralAxis =
        getSpectrum().GetSpectralAxis();
    auto const &[imin, imax] = getLambdaRange().getClosestInnerIndices(
        spcSpectralAxis.GetSamplesVector());

    nSamples += abs(imax - imin + 1);
  }

  return nSamples;
}

/**
 * \brief Returns the cumulative SNR under the Strong Emission Lines
 * 1. retrieve the lines support
 * 2. process each
 **/
std::pair<Float64, Float64> CLineModelFitting::getCumulSNRStrongEL() const {
  // NB dummy multiobs implementation (functional for one obs only)
  m_spectraIndex.setAtBegining(); // temporary multiobs implementation
  // Retrieve all the strone emission lines supports in a list of range
  TInt32RangeList supportList;
  TBoolList isStrongList;
  TInt32List validEltsIdx =
      getElementList().GetElementsIndicesInsideLambdaRange();
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
      supportList.push_back(elt->getSupportSubElt(index));
    }
  }

  // merge overlapping ranges
  TInt32RangeList nonOverlappingSupportList;
  TBoolList nonOverlappingIsStrongList;
  TInt32List processedSupport;
  for (Int32 k = 0; k < ssize(supportList); k++) {
    // skip if already fitted
    if (std::find(processedSupport.cbegin(), processedSupport.cend(), k) !=
        processedSupport.cend())
      continue;

    processedSupport.push_back(k);
    TInt32Range &support = supportList[k];

    for (Int32 l = k + 1; l < ssize(supportList); l++) {
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

  setRedshift(modelSolution.Redshift);

  // reset before loading
  for (auto param_ptr : m_ElementsVector->getElementsParams()) {
    param_ptr->resetFittingParams();
    param_ptr->setVelocity(NAN);
  }
  m_ElementsVector->resetLambdaOffsets();

  // should also reset nominal amplitudes...
  // but not touched without using template-ratio

  TBoolList element_done(getElementsParams().size(), false);
  for (Int32 iRestLine = 0; iRestLine < ssize(m_RestLineList); iRestLine++) {
    Int32 eIdx = modelSolution.ElementId[iRestLine];
    if (eIdx == undefIdx)
      THROWG(ErrorCode::INTERNAL_ERROR,
             Formatter() << "Undefined element index, for rest line index "
                         << iRestLine << " in model solution");
    auto const &elt_param = getElementsParams()[eIdx];
    if (modelSolution.NotFitted[iRestLine]) {
      // set outsidelambdrarangeList
      for ([[maybe_unused]] auto const spcIndex : m_spectraIndex) {
        auto const &elt_ptr = getElementList()[eIdx];
        for (Int32 line_idx = 0; line_idx < elt_ptr->GetSize(); ++line_idx)
          elt_ptr->SetOutsideLambdaRangeList(line_idx);
      }
      continue;
    }
    Int32 line_id = modelSolution.lineId[iRestLine];

    Int32 elt_line_index = elt_param->getLineIndex(line_id);
    if (elt_line_index == undefIdx)
      continue; // or throw an exception ?

    elt_param->setFittedAmplitude(
        elt_line_index, modelSolution.Amplitudes[iRestLine],
        modelSolution.AmplitudesUncertainties[iRestLine]);
    elt_param->setLambdaOffset(elt_line_index, modelSolution.Offset[iRestLine]);

    if (element_done[eIdx])
      continue;

    elt_param->setVelocity(modelSolution.Velocity[iRestLine]);
    elt_param->SetFittingGroupInfo(modelSolution.fittingGroupInfo[iRestLine]);
    if (m_enableAmplitudeOffsets) {
      CPolynomCoeffs contPolynomCoeffs = {
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
      getElementsParams()[idxLyaE]->SetAsymfitParams(
          {modelSolution.LyaWidthCoeff, modelSolution.LyaAlpha,
           modelSolution.LyaDelta});
  }

  if (modelSolution.LyaIgm != undefIdx) {
    auto const indices_Igm = m_ElementsVector->getIgmLinesIndices();
    if (!indices_Igm.empty())
      for (auto const &[elt_idx, _] : indices_Igm)
        getElementsParams()[elt_idx]->SetSymIgmParams(
            {modelSolution.LyaIgm, modelSolution.Redshift});
  }

  for ([[maybe_unused]] auto &spcIndex : m_spectraIndex) {
    const CSpectrumSpectralAxis &spectralAxis = getSpectrum().GetSpectralAxis();
    for (auto const &elt_ptr : getElementList()) {
      elt_ptr->computeOutsideLambdaRange();
      if (!elt_ptr->IsOutsideLambdaRange())
        elt_ptr->prepareSupport(spectralAxis, modelSolution.Redshift,
                                getLambdaRange());
    }
  }
  m_ElementsVector->computeGlobalOutsideLambdaRange();

  return;
}

// should be called only on final candidates
void CLineModelFitting::ComputeAndAddOptionalLineProperties(
    CLineModelSolution &modelSolution) const {
  TInt32List eIdx_oii;
  TInt32List subeIdx_oii;
  TInt32List eIdx_ha;
  TInt32List subeIdx_ha;
  Float64 flux_oii = 0.0;
  Float64 fluxVar_oii = 0.0;
  Float64 flux_ha = 0.0;
  Float64 fluxVar_ha = 0.0;

  for (Int32 iRestLine = 0; iRestLine < ssize(m_RestLineList); ++iRestLine) {
    processSingleLine(iRestLine, modelSolution, eIdx_oii, subeIdx_oii, flux_oii,
                      fluxVar_oii, eIdx_ha, subeIdx_ha, flux_ha, fluxVar_ha);
  }

  addLyaParams(modelSolution);
}

void CLineModelFitting::processSingleLine(
    Int32 iRestLine, CLineModelSolution &modelSolution, TInt32List &eIdx_oii,
    TInt32List &subeIdx_oii, Float64 &flux_oii, Float64 &fluxVar_oii,
    TInt32List &eIdx_ha, TInt32List &subeIdx_ha, Float64 &flux_ha,
    Float64 &fluxVar_ha) const {
  Int32 line_id = modelSolution.lineId[iRestLine];
  Int32 eIdx = modelSolution.ElementId[iRestLine];
  auto const elt_param = getElementsParams()[eIdx];
  Int32 line_index = elt_param->getLineIndex(line_id);

  if (eIdx == undefIdx || line_index == undefIdx ||
      elt_param->isNotFittable() ||
      elt_param->isOutsideLambdaRangeLine(line_index))
    return;

  modelSolution.NSamples[iRestLine] =
      computeNSamplesUnderLine(eIdx, line_index);

  updateResidualsAndContinuum(iRestLine, modelSolution, eIdx, line_index);

  auto [flux, fluxError, isEmission] =
      computeLineFlux(iRestLine, modelSolution, eIdx, line_index);

  modelSolution.Flux[iRestLine] = flux;
  if (isLineRatioRules())
    modelSolution.FluxUncertainty[iRestLine] = fluxError;

  auto [fluxDI, snrDI] =
      getFluxDirectIntegration({eIdx}, {line_index}, isEmission ? 1 : 0);
  modelSolution.FluxDirectIntegration[iRestLine] = fluxDI;
  modelSolution.FluxDirectIntegrationUncertainty[iRestLine] =
      std::abs(fluxDI) / snrDI;

  if (isEmission)
    accumulateLineFluxes(flux, fluxError, eIdx, line_index, line_id,
                         modelSolution, eIdx_ha, subeIdx_ha, flux_ha,
                         fluxVar_ha, eIdx_oii, subeIdx_oii, flux_oii,
                         fluxVar_oii);

  modelSolution.fittingGroupInfo[iRestLine] = elt_param->m_fittingGroupInfo;
}

void CLineModelFitting::updateResidualsAndContinuum(
    Int32 iRestLine, CLineModelSolution &modelSolution, Int32 eIdx,
    Int32 line_index) const {
  modelSolution.ResidualRMS[iRestLine] =
      m_fitter->getModelResidualRmsUnderElements({eIdx});

  if (m_enableAmplitudeOffsets) {
    const auto &polynom_coeffs = getElementsParams()[eIdx]->m_ampOffsetsCoeffs;
    modelSolution.continuum_pCoeff0[iRestLine] = polynom_coeffs.m_a0;
    modelSolution.continuum_pCoeff1[iRestLine] = polynom_coeffs.m_a1;
    modelSolution.continuum_pCoeff2[iRestLine] = polynom_coeffs.m_a2;
  }

  Float64 cont, cont_std;
  if (m_fittingmethod == "svd" || m_fittingmethod == "hybrid" ||
      m_fittingmethod == "lbfgsb")
    std::tie(cont, cont_std) =
        GetContinuumAtCenterProfile(eIdx, line_index, modelSolution.Redshift);
  else
    std::tie(cont, cont_std) =
        GetMeanContinuumUnderLine(eIdx, line_index, modelSolution.Redshift);

  modelSolution.CenterContinuumFlux[iRestLine] = cont;
  modelSolution.CenterContinuumFluxUncertainty[iRestLine] = cont_std;
}

Int32 CLineModelFitting::computeNSamplesUnderLine(Int32 eIdx,
                                                  Int32 line_index) const {
  Int32 NSamples = 0;
  for ([[maybe_unused]] auto &spcIndex : m_spectraIndex) {
    const auto &elt = getElementList()[eIdx];
    if (elt->IsOutsideLambdaRangeLine(line_index))
      continue;
    NSamples += elt->getSupportSubElt(line_index).GetLength() + 1;
  }
  return NSamples;
}

std::tuple<Float64, Float64, bool>
CLineModelFitting::computeLineFlux(Int32 iRestLine,
                                   CLineModelSolution &modelSolution,
                                   Int32 eIdx, Int32 line_index) const {
  const auto &line = m_RestLineList.at(modelSolution.lineId[iRestLine]);
  const bool isEmission = line.GetType() == CLine::EType::nType_Emission;

  Float64 amp = modelSolution.Amplitudes[iRestLine];
  Float64 ampError = modelSolution.AmplitudesUncertainties[iRestLine];
  Float64 flux = NAN;
  Float64 fluxError = NAN;

  if (!std::isnan(amp) && amp >= 0.0) {
    if (!isEmission) {
      Float64 const positive_cont =
          std::max(0.0, modelSolution.CenterContinuumFlux[iRestLine]);
      ampError *= positive_cont;
      ampError += amp * modelSolution.CenterContinuumFluxUncertainty[iRestLine];
      amp *= -positive_cont;
    }

    for ([[maybe_unused]] auto &spcIndex : m_spectraIndex) {
      const auto &eltList = getElementList();
      if (!eltList[eIdx]->IsOutsideLambdaRangeLine(line_index)) {
        auto const &[mu, sigma] =
            eltList[eIdx]->getObservedPositionAndLineWidth(
                modelSolution.Redshift, line_index, false);
        modelSolution.Sigmas[iRestLine] = sigma;
        const auto &profile =
            eltList[eIdx]->getElementParam()->getLineProfile(line_index);
        auto const &rawFlux = profile->GetLineFlux(mu, sigma);
        flux = amp * rawFlux;
        fluxError = ampError * rawFlux;
        break;
      }
    }
  }

  return {flux, fluxError, isEmission};
}

void CLineModelFitting::accumulateLineFluxes(
    Float64 flux, Float64 fluxError, Int32 eIdx, Int32 line_index,
    Int32 line_id, CLineModelSolution &modelSolution, TInt32List &eIdx_ha,
    TInt32List &subeIdx_ha, Float64 &flux_ha, Float64 &fluxVar_ha,
    TInt32List &eIdx_oii, TInt32List &subeIdx_oii, Float64 &flux_oii,
    Float64 &fluxVar_oii) const {
  auto const &line = m_RestLineList.at(line_id);
  if (line.GetName() == linetags::halpha_em ||
      line.GetName() == linetags::niia_em ||
      line.GetName() == linetags::niib_em) {
    eIdx_ha.push_back(eIdx);
    subeIdx_ha.push_back(line_index);
    if (flux > 0.0)
      flux_ha += flux;
    if (fluxError > 0.0)
      fluxVar_ha += fluxError * fluxError;
    if (eIdx_ha.size() == 3) {
      auto [fluxDI, snrDI] = getFluxDirectIntegration(eIdx_ha, subeIdx_ha, 0);
      modelSolution.snrHa_DI = snrDI;
      modelSolution.lfHa_DI = fluxDI > 0.0 ? log10(fluxDI) : -INFINITY;
      modelSolution.lfHa = flux_ha > 0.0 ? log10(flux_ha) : -INFINITY;
      if (isLineRatioRules())
        modelSolution.snrHa = flux_ha / std::sqrt(fluxVar_ha);
    }
  } else if (line.GetName() == linetags::oII3726_em ||
             line.GetName() == linetags::oII3729_em) {
    eIdx_oii.push_back(eIdx);
    subeIdx_oii.push_back(line_index);
    if (flux > 0.0)
      flux_oii += flux;
    if (fluxError > 0.0)
      fluxVar_oii += fluxError * fluxError;
    if (eIdx_oii.size() == 2) {
      auto [fluxDI, snrDI] = getFluxDirectIntegration(eIdx_oii, subeIdx_oii, 0);
      modelSolution.snrOII_DI = snrDI;
      modelSolution.lfOII_DI = fluxDI > 0 ? log10(fluxDI) : -INFINITY;
      modelSolution.lfOII = flux_oii > 0 ? log10(flux_oii) : -INFINITY;
      if (isLineRatioRules())
        modelSolution.snrOII = flux_oii / std::sqrt(fluxVar_oii);
    }
  }
}

void CLineModelFitting::addLyaParams(CLineModelSolution &modelSolution) const {
  std::string lyaTag = linetags::lya_em;
  auto const [idxLyaE, _] = m_ElementsVector->findElementIndex(lyaTag);
  if (idxLyaE != undefIdx) {
    const auto &params =
        m_ElementsVector->getElementsParams()[idxLyaE]->GetAsymfitParams(0);
    modelSolution.LyaWidthCoeff = params.sigma;
    modelSolution.LyaAlpha = params.alpha;
    modelSolution.LyaDelta = params.delta;
    const auto &params_igm =
        m_ElementsVector->getElementsParams()[idxLyaE]->GetSymIgmParams(0);
    modelSolution.LyaIgm = params_igm.m_igmidx;
  }
}

/**
 * \brief Returns a CLineModelSolution object populated with the current
 *solutions.
 **/
// this is not really a const method as spectrum model(s) have to be modified
// (cf CSpectrumModel::getContinuumUncertainty)
CLineModelSolution
CLineModelFitting::GetModelSolution(bool fullSolution) const {
  Int32 s = m_RestLineList.size();
  CLineModelSolution modelSolution(m_RestLineList);

  auto const &elt_param_vect =
      static_pointer_cast<const CLMEltListVector>(m_ElementsVector)
          ->getElementsParams();

  modelSolution.EmissionVelocity = GetVelocityEmission();
  modelSolution.AbsorptionVelocity = GetVelocityAbsorption();

  // For some quantities it is more simple to get them from the first
  // observation objects There could be refactor but it can be complicated for
  // no gain of clarity or robustness
  m_spectraIndex.setAtBegining();
  modelSolution.Redshift = getSpectrumModel().m_Redshift;
  const CLineModelElementList &firstEltList = getElementList();

  modelSolution.nDDL = m_ElementsVector->getNonZeroElementsNDdl();

  for (Int32 iRestLine = 0; iRestLine < s; iRestLine++) {
    Int32 line_id = modelSolution.lineId[iRestLine];
    auto [eIdx, line_index] = m_ElementsVector->findElementIndex(line_id);
    modelSolution.ElementId[iRestLine] = eIdx;
    if (eIdx == undefIdx || line_index == undefIdx ||
        elt_param_vect[eIdx]->isNotFittable() ||
        elt_param_vect[eIdx]->isOutsideLambdaRangeLine(line_index))
      continue; // data already set to its default values

    Float64 amp = elt_param_vect[eIdx]->m_FittedAmplitudes[line_index];
    modelSolution.Amplitudes[iRestLine] = amp;
    Float64 ampError = elt_param_vect[eIdx]->m_FittedAmplitudesStd[line_index];
    modelSolution.AmplitudesUncertainties[iRestLine] = ampError;
    if (isLineRatioRules())
      modelSolution.SNR[iRestLine] = std::abs(amp) / ampError;

    modelSolution.LambdaObs[iRestLine] =
        firstEltList[eIdx]->GetObservedPosition(line_index,
                                                modelSolution.Redshift);
    modelSolution.Velocity[iRestLine] = elt_param_vect[eIdx]->getVelocity();
    modelSolution.VelocityUncertainty[iRestLine] =
        elt_param_vect[eIdx]->getVelocityStd();
    modelSolution.Offset[iRestLine] =
        elt_param_vect[eIdx]->getLambdaOffset(line_index);
    modelSolution.OffsetUncertainty[iRestLine] =
        elt_param_vect[eIdx]->getLambdaOffsetStd(line_index);
    modelSolution.NotFitted[iRestLine] = false;
  }

  std::unordered_set<std::string>
      strongELSNRAboveCut; // = getLinesAboveSNR(3.5);
  modelSolution.NLinesAboveSnrCut = strongELSNRAboveCut.size();

  // brief, to save processing time, do not estimate fluxes
  // and high level line properties
  if (fullSolution)
    ComputeAndAddOptionalLineProperties(modelSolution);

  return modelSolution;
}

void CLineModelFitting::SetLSF(std::shared_ptr<const CLSF> const &lsf_) {
  for ([[maybe_unused]] auto &spcIndex : m_spectraIndex) {

    const std::shared_ptr<const CLSF> &lsf =
        lsf_ ? lsf_ : getSpectrum().GetLSF();

    if (lsf == nullptr) {
      THROWG(ErrorCode::INTERNAL_ERROR,
             "Cannot enable LSF, LSF spectrum member is not initialized");
    }
    const auto &[valid, message] = lsf->IsValid();

    if (!valid) {
      THROWG(ErrorCode::INTERNAL_ERROR,
             " Cannot enable LSF, LSF spectrum member is not valid " + message);
    }

    getElementList().setLSF(lsf);
  }
}

void CLineModelFitting::SetVelocityEmission(Float64 vel) {
  for (auto &lmep : m_ElementsVector->getElementsParams())
    if (lmep->IsEmission())
      lmep->setVelocity(vel);
}

void CLineModelFitting::SetVelocityAbsorption(Float64 vel) {
  for (auto &lmep : m_ElementsVector->getElementsParams())
    if (lmep->IsAbsorption())
      lmep->setVelocity(vel);
}

void CLineModelFitting::setVelocityByGroup(Float64 vel,
                                           const TInt32List &inds) {
  for (auto idxElt : inds)
    m_ElementsVector->getElementsParams()[idxElt]->setVelocity(vel);
}

Float64 CLineModelFitting::GetVelocityEmission() const {

  // no global emission or absorption velocities
  if (isLineRatioRules() && m_fittingmethod == "lbfgsb")
    return NAN;

  // find 1st emission element
  auto elt_param_vect = m_ElementsVector->getElementsParams();
  auto const it = std::find_if(elt_param_vect.begin(), elt_param_vect.end(),
                               [](TLineModelElementParam_ptr const &p) {
                                 return p->IsEmission() && p->isFittable();
                               });
  if (it == elt_param_vect.end())
    return NAN;

  return (*it)->getVelocity();
}

Float64 CLineModelFitting::GetVelocityAbsorption() const {
  // no global emission or absorption velocities
  if (isLineRatioRules() && m_fittingmethod == "lbfgsb")
    return NAN;

  // find 1st emission element
  auto elt_param_vect = m_ElementsVector->getElementsParams();
  auto const it = std::find_if(elt_param_vect.begin(), elt_param_vect.end(),
                               [](TLineModelElementParam_ptr const &p) {
                                 return p->IsAbsorption() && p->isFittable();
                               });
  if (it == elt_param_vect.end())
    return NAN;

  return (*it)->getVelocity();
}

/**
 * \brief this function returns the dtd value withing the wavelength range for
 *a given spcComponent
 *
 **/
Float64 CLineModelFitting::getOrInitDtD() {

  m_spectraIndex.setAtBegining(); // we choose arbitrarily first obs to check if
                                  // dtd is already initialized
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
Float64 CLineModelFitting::getOrInitLikelihoodCstLog() {

  m_spectraIndex.setAtBegining(); // we choose arbitrarily first obs to check if
                                  // dtd is already initialized
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
  for ([[maybe_unused]] auto &spcIndex : m_spectraIndex) {

    const CSpectrumSpectralAxis &spcSpectralAxis =
        getSpectrum().GetSpectralAxis();
    const CSpectrumFluxAxis &Yspc =
        (spcComponent == "noContinuum")
            ? getSpectrumModel().getSpcFluxAxisNoContinuum()
            : getSpectrumModel().getSpcFluxAxis();

    auto const &iRange = TInt32Range(getLambdaRange().getClosestInnerIndices(
        spcSpectralAxis.GetSamplesVector()));
    dtd += std::transform_reduce(
        iRange.begin(), iRange.end(), 0., std::plus(),
        [&Yspc](Int32 j) { return Yspc[j] * Yspc[j] * Yspc.GetWeight(j); });
    Log.LogDebug(Formatter()
                 << "CLineModelFitting::EstimateDTransposeD val = " << dtd);
  }
  return dtd;
}

/**
 * \brief this function estimates the mtm value withing the wavelength range
 **/
Float64 CLineModelFitting::EstimateMTransposeM() const {
  Float64 mtm = 0.0;
  for ([[maybe_unused]] auto &spcIndex : m_spectraIndex) {

    const CSpectrumSpectralAxis &spcSpectralAxis =
        getSpectrum().GetSpectralAxis();
    const CSpectrumFluxAxis &spcFluxAxis =
        getSpectrumModel().GetModelSpectrum().GetFluxAxis();
    const auto &ErrorNoContinuum = getSpectrum().GetErrorAxis();

    Float64 diff = 0.0;

    auto const &[imin, imax] = getLambdaRange().getClosestInnerIndices(
        spcSpectralAxis.GetSamplesVector());

    for (Int32 j = imin; j <= imax; j++) {
      diff = spcFluxAxis[j];
      mtm += (diff * diff) / (ErrorNoContinuum[j] * ErrorNoContinuum[j]);
    }
  }
  return mtm;
}

void CLineModelFitting::setContinuumComponent(TContinuumComponent component) {
  m_continuumManager->setContinuumComponent(std::move(component));
}
/**
 * \brief this function estimates the likelihood_cstLog term withing the
 *wavelength range
 **/
Float64 CLineModelFitting::EstimateLikelihoodCstLog(
    const std::string &spcComponent) const {

  Float64 cstLog = 0.0;
  for ([[maybe_unused]] auto &spcIndex : m_spectraIndex) {

    const CSpectrumSpectralAxis &spcSpectralAxis =
        getSpectrum().GetSpectralAxis();

    auto const &flux_for_weight =
        (spcComponent == "noContinuum")
            ? getSpectrumModel().getSpcFluxAxisNoContinuum()
            : getSpectrumModel().getSpcFluxAxis();

    auto const &iRange = TInt32Range(getLambdaRange().getClosestInnerIndices(
        spcSpectralAxis.GetSamplesVector()));

    Int32 numDevs = iRange.GetLength() + 1;
    Float64 sumLogNoise =
        std::transform_reduce(iRange.begin(), iRange.end(), 0., std::plus(),
                              [&flux_for_weight](Int32 j) {
                                return log(flux_for_weight.GetWeight(j));
                              });

    cstLog += -numDevs * 0.5 * log(2 * M_PI) + 0.5 * sumLogNoise;
  }
  return cstLog;
}

// assumes model is refreshed and continuum uptodate
std::pair<Float64, Float64>
CLineModelFitting::GetMeanContinuumUnderLine(Int32 eltIdx, Int32 line_index,
                                             Float64 redshift) const {

  Float64 sumContinuumAll = 0.0;
  Float64 sumWeightAll = 0.0;
  Float64 sumSquaredWeightAll = 0.0;
  Float64 sumResidualAll = 0.0;
  Int32 npixResidualAll = 0;
  Float64 sumWeightResidualAll = 0;
  Float64 sumSquaredWeightResidualAll = 0;
  auto const &polynomCoeffs =
      m_ElementsVector->getElementsParams()[eltIdx]->m_ampOffsetsCoeffs;
  for ([[maybe_unused]] auto const &spcIndex : m_spectraIndex) {
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

    auto const &[residual_sum, nsum, nsum2] =
        model.getContinuumSquaredResidualInRange(indexRange);
    sumResidualAll += residual_sum;
    sumWeightResidualAll += nsum;
    sumSquaredWeightResidualAll += nsum2;
    npixResidualAll += indexRange.GetLength() + 1;
  }
  if (sumWeightAll == 0.0)
    return std::make_pair(NAN, NAN);

  Float64 const continuum = sumContinuumAll / sumWeightAll;

  Float64 const std = npixResidualAll >= RMS_MIN_SAMPLE_NUMBER
                          ? sqrt(sumResidualAll / sumSquaredWeightResidualAll) *
                                sqrt(sumSquaredWeightAll) / sumWeightAll
                          : NAN;

  return std::make_pair(continuum, std);
}

std::pair<Float64, Float64>
CLineModelFitting::GetContinuumAtCenterProfile(Int32 eltIdx, Int32 line_index,
                                               Float64 redshift) const {
  for ([[maybe_unused]] auto &spcIndex : m_spectraIndex) {
    auto const &elt = *m_ElementsVector->getElementList()[eltIdx];
    if (elt.IsOutsideLambdaRangeLine(line_index))
      continue;

    auto &model = getSpectrumModel();
    auto const &spectralAxis = model.GetModelSpectrum().GetSpectralAxis();
    auto const continuumFluxAxis = model.GetModelContinuum();
    return elt.GetContinuumAtCenterProfile(line_index, spectralAxis, redshift,
                                           continuumFluxAxis,
                                           m_enableAmplitudeOffsets);
  }
  return std::make_pair(NAN, NAN);
}

std::pair<Float64, Float64> CLineModelFitting::getFluxDirectIntegration(
    const TInt32List &eIdx_list, const TInt32List &subeIdx_list,
    bool substract_abslinesmodel) const {

  Float64 sumFlux = 0;
  Float64 sumErr = 0;
  for ([[maybe_unused]] auto &spcIndex : m_spectraIndex) {
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

void CLineModelFitting::refreshAllModels() { m_models->refreshAllModels(); }

void CLineModelFitting::setChiSquareRatioResult(
    const Int32 index_z, const std::shared_ptr<CLineModelResult> &lmResult) {
  return m_lineRatioManager->setChiSquareRatioResult(index_z, lmResult);
}

std::shared_ptr<const CLSF>
CLineModelFitting::buildEquivConstantResolLSF() const {
  // NB dummy multiobs implementation (functional for one obs only)
  getSpectraIndex().setAtBegining(); // temporary multiobs implementation

  Float64 lambda = getLambdaRange().GetMidRange();

  Float64 resolution = CLSFGaussianConstantResolution::computeResolution(
      lambda, getSpectrum().GetLSF()->GetWidth(lambda));
  std::shared_ptr<TLSFArguments> args =
      std::make_shared<TLSFGaussianConstantResolutionArgs>(resolution);

  return LSFFactory.Create("gaussianConstantResolution", args);
}
