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
#include "RedshiftLibrary/linemodel/linemodelfitting.h"
#include "RedshiftLibrary/line/catalogsTplRatio.h"
#include "RedshiftLibrary/line/line.h"
#include "RedshiftLibrary/line/regulament.h"
#include "RedshiftLibrary/linemodel/element.h"

#include "RedshiftLibrary/linemodel/rulesmanager.h"
#include "RedshiftLibrary/linemodel/tplcorrmanager.h"
#include "RedshiftLibrary/linemodel/tplratiomanager.h"

#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/continuum/irregularsamplingmedian.h"
#include "RedshiftLibrary/extremum/extremum.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/processflow/autoscope.h"
#include "RedshiftLibrary/processflow/context.h"
#include "RedshiftLibrary/spectrum/template/template.h"

#include <boost/chrono/thread_clock.hpp>
#include <boost/format.hpp>
#include <boost/numeric/conversion/bounds.hpp>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_vector.h>

#include <climits>
#include <cmath>
#include <memory>
#include <numeric>

using namespace NSEpic;
using namespace std;

/**
 * \brief Prepares the state for Linemodel operation.
 * Loads the catalog.
 * Sets many state variables.
 * Sets the continuum either as a nocontinuum or a fromspectrum.
 **/
CLineModelFitting::CLineModelFitting(
    const std::shared_ptr<COperatorTemplateFittingBase> &TFOperator)
    : m_RestLineList(Context.getCLineMap()) {
  initParameters();

  m_inputSpcs = std::make_shared<std::vector<std::shared_ptr<const CSpectrum>>>(
      Context.getSpectra(m_useloglambdasampling));
  m_lambdaRanges = // std::make_shared<std::vector<std::shared_ptr<const
                   // TLambdaRange>>>(
      Context.getClampedLambdaRanges(m_useloglambdasampling);

  initMembers(TFOperator);
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

  initMembers(TFOperator);
  setLineRatioType(m_lineRatioType);

  dynamic_cast<CRulesManager *>(m_lineRatioManager.get())->setRulesOption("no");
  setContinuumComponent("fromspectrum");
}

void CLineModelFitting::initParameters() {
  CAutoScope autoscope(Context.m_ScopeStack, "linemodel");

  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();
  m_fittingmethod = ps->GetScoped<std::string>("fittingmethod");
  m_enableAmplitudeOffsets = ps->GetScoped<bool>("ampoffsetfit");
  m_enableLbdaOffsets = ps->GetScoped<bool>("lbdaoffsetfit");
  m_LineWidthType = ps->GetScoped<std::string>("linewidthtype");

  m_lineRatioType = ps->GetScoped<std::string>("lineRatioType");

  if (Context.GetCurrentMethod() == "LineModelSolve") {
    m_opt_firstpass_fittingmethod =
        ps->GetScoped<std::string>("firstpass.fittingmethod");
    m_opt_secondpass_fittingmethod = m_fittingmethod;
    m_opt_firstpass_forcedisableMultipleContinuumfit =
        ps->GetScoped<bool>("firstpass.multiplecontinuumfit_disable");
  }

  std::string continuumComponent =
      ps->GetScoped<std::string>("continuumcomponent");
  if (continuumComponent == "tplfit" || continuumComponent == "tplfitauto") {
    m_opt_firstpass_forcedisableMultipleContinuumfit =
        ps->GetScoped<bool>("firstpass.multiplecontinuumfit_disable");
    m_useloglambdasampling = ps->GetScoped<bool>("useloglambdasampling");
    m_opt_fitcontinuum_maxN = ps->GetScoped<Int32>("continuumfit.count");
  }
}

void CLineModelFitting::initMembers(
    const std::shared_ptr<COperatorTemplateFittingBase> &TFOperator) {
  m_nbObs = m_inputSpcs->size();
  m_curObs = std::make_shared<Int32>(0);
  m_ElementsVector = std::make_shared<std::vector<CLineModelElementList>>();
  m_nominalWidthDefault = 13.4; // euclid 1 px

  for (*m_curObs = 0; *m_curObs < m_nbObs; (*m_curObs)++) {
    Log.LogDetail("    model: Continuum winsize found is %.2f A",
                  getSpectrum().GetMedianWinsize());
    m_ElementsVector->push_back(CLineModelElementList());
    if (m_lineRatioType == "rules") {
      // load the regular catalog
      LoadCatalog();
    } else { //"tplratio" and "tplcorr"
      // load the tplratio catalog with only 1 element for all lines
      // LoadCatalogOneMultiline(restLineList);
      // load the tplratio catalog with 2 elements: 1 for the Em lines + 1 for
      // the Abs lines
      LoadCatalogTwoMultilinesAE();
    }
  }
  *m_curObs = 0;
  m_continuumFitValues = std::make_shared<CTplModelSolution>();
  m_models = std::make_shared<std::vector<CSpectrumModel>>();
  m_models->push_back(CSpectrumModel(getElementList(), getSpectrumPtr(),
                                     m_RestLineList, m_continuumFitValues,
                                     TFOperator));
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
  Log.LogDetail(Formatter() << "m_opt_secondpass_fittingmethod"
                            << m_opt_secondpass_fittingmethod);

  Log.LogDetail(Formatter() << "LineWidthType=" << m_LineWidthType);

  Log.LogDetail(Formatter() << "nominalWidthDefault=" << m_nominalWidthDefault);

  Log.LogDetail(Formatter() << "fittingmethod=" << m_fittingmethod);

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

void CLineModelFitting::AddElement(CLineVector lines, Float64 velocityEmission,
                                   Float64 velocityAbsorption, Int32 ig) {
  if (*m_curObs == 0) {
    m_ElementParam.push_back(std::make_shared<TLineModelElementParam>(
        std::move(lines), velocityEmission, velocityAbsorption));
  }
  getElementList().push_back(
      std::make_shared<CLineModelElement>(m_ElementParam[ig], m_LineWidthType));
}

/**
 * \brief For each line in each group of the argument, finds the associated
 *line in the catalog and saves this information to getElementList(). Converts
 *the argument restLineList to a group list. For each entry in this list: For
 *each line in this entry: Finds the index in the catalog from the line name and
 *type. Saves the line, the catalog index and the nominal amplitude for the
 *line thusly associated to this line. If at least one line was found, save
 *this result in getElementList().
 **/
void CLineModelFitting::LoadCatalog() {
  CAutoScope autoscope(Context.m_ScopeStack, "linemodel");

  auto const ps = Context.GetParameterStore();
  Float64 const velocityEmission = ps->GetScoped<Float64>("velocityemission");
  Float64 const velocityAbsorption =
      ps->GetScoped<Float64>("velocityabsorption");
  Int32 lastEltIndex = 0;

  auto const groupList = CLineCatalog::ConvertToGroupList(m_RestLineList);

  for (auto [_, lines] : groupList) {
    AddElement(std::move(lines), velocityEmission, velocityAbsorption,
               lastEltIndex);
    lastEltIndex++;
  }
}

void CLineModelFitting::LoadCatalogOneMultiline() {
  CAutoScope const autoscope(Context.m_ScopeStack, "linemodel");

  auto const ps = Context.GetParameterStore();
  Float64 const velocityEmission = ps->GetScoped<Float64>("velocityemission");
  Float64 const velocityAbsorption =
      ps->GetScoped<Float64>("velocityabsorption");

  CLineVector RestLineVector;
  RestLineVector.reserve(m_RestLineList.size());
  for (auto const &[_, line] : m_RestLineList)
    RestLineVector.push_back(line);

  AddElement(std::move(RestLineVector), velocityEmission, velocityAbsorption,
             0);
}

void CLineModelFitting::LoadCatalogTwoMultilinesAE() {
  CAutoScope autoscope(Context.m_ScopeStack, "linemodel");

  auto const ps = Context.GetParameterStore();
  Float64 const velocityEmission = ps->GetScoped<Float64>("velocityemission");
  Float64 const velocityAbsorption =
      ps->GetScoped<Float64>("velocityabsorption");

  std::vector<CLine::EType> const types = {CLine::EType::nType_Absorption,
                                           CLine::EType::nType_Emission};

  Int32 lastEltIndex = 0;
  for (auto type : types) {
    CLineVector lines;
    for (auto const &[id, line] : m_RestLineList) {
      if (line.GetType() == type)
        lines.push_back(line);
    }

    if (lines.size() > 0) {
      AddElement(std::move(lines), velocityEmission, velocityAbsorption,
                 lastEltIndex);
      lastEltIndex++;
    }
  }
}

/**
 * \brief LogDetail the number of lines for each element, and their nominal
 *amplitudes.
 **/
void CLineModelFitting::LogCatalogInfos() {
  Log.LogDetail("\n");
  Log.LogDetail("LineModel Infos: %d elements", getElementList().size());
  for (Int32 iElts = 0; iElts < getElementList().size(); iElts++) {
    auto const &elt = getElementList()[iElts];
    Int32 nLines = elt->GetSize();
    if (nLines < 1) {
      Log.LogDetail("LineModel ctlg: elt %d (%s): no lines", iElts,
                    CLine::ETypeString.at(elt->GetElementType()).c_str());
    }
    for (Int32 index = 0; index != elt->GetSize(); ++index) {
      std::string nominalAmpStr = "";
      nominalAmpStr = boost::str(boost::format("(nominal amp = %.4e)") %
                                 elt->GetNominalAmplitude(index));
      Log.LogDetail("LineModel ctlg: elt %d (%s): line %d = %s %s", iElts,
                    CLine::ETypeString.at(elt->GetElementType()).c_str(), index,
                    elt->GetLineName(index).c_str(), nominalAmpStr.c_str());
    }
  }
  Log.LogDetail("\n");
}

/*
Change the actual value of redshift.
the continuum can be reinterpolate.
*/
void CLineModelFitting::setRedshift(Float64 redshift,
                                    bool reinterpolatedContinuum) {
  getSpectrumModel().m_Redshift = redshift;

  if (reinterpolatedContinuum) {
    m_continuumManager->reinterpolateContinuum(redshift);
  }
}

Int32 CLineModelFitting::getTplratio_count() const {
  return m_lineRatioManager->getTplratio_count();
}

TFloat64List CLineModelFitting::getTplratio_priors() {
  return m_lineRatioManager->getTplratio_priors();
}

bool CLineModelFitting::initDtd() {
  //  m_dTransposeDLambdaRange = TLambdaRange(*(m_lambdaRange));
  m_dTransposeDLambdaRange = getLambdaRange();
  if (isContinuumComponentTplfitxx())
    m_dTransposeD = EstimateDTransposeD("raw");
  else
    m_dTransposeD = EstimateDTransposeD("nocontinuum");

  m_likelihood_cstLog = EstimateLikelihoodCstLog();
  return true;
}

void CLineModelFitting::prepareAndLoadContinuum(Int32 k, Float64 redshift) {
  if (getContinuumComponent() == "nocontinuum")
    return;

  if (!isContinuumComponentTplfitxx()) {
    getSpectrumModel().setContinuumToInputSpc();
    return;
  }

  // the support has to be already computed
  // when LoadFitContinuum() is called
  getSpectrumModel().initObserveGridContinuumFlux(
      getSpectrum().GetSampleCount());
  Int32 autoselect = getContinuumComponent() == "tplfitauto";
  m_continuumManager->LoadFitContinuum(k, autoselect, redshift);
}

void CLineModelFitting::computeSpectrumFluxWithoutContinuum() {
  getSpectrumModel().initModelWithContinuum();
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
    Float64 _meritprior = 0.; // only relevant for "tplratio"
    prepareAndLoadContinuum(k, redshift);
    if (getContinuumComponent() != "nocontinuum")
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
      if (getContinuumComponent() == "nocontinuum")
        getSpectrumModel().reinitModel();
    }
  }

  if (!enableLogging)
    return bestMerit;

  if (isContinuumComponentTplfitxx()) {
    if (m_fittingmethod != "svdlc" && nContinuum > 1) {
      Int32 autoselect =
          m_continuumManager->getContinuumComponent() == "tplfitauto";
      // TODO savedIdxContinuumFitted=-1 if lineRatioType!=tplratio
      m_continuumManager->LoadFitContinuum(savedIdxContinuumFitted, autoselect,
                                           redshift);
    }
    /*
    Log.LogDetail("    model - Linemodel: fitcontinuum = %d (%s, with "
                  "ebmv=%.3f), and A=%e",
                  savedIdxContinuumFitted, m_fitContinuum_tplName.c_str(),
                  m_fitContinuum_tplFitEbmvCoeff,
                  m_fitContinuum_tplFitAmplitude);
    */
  }
  if (m_lineRatioType == "tplratio") {
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
  m_fitter = CAbstractFitter::makeFitter(
      fitMethod, m_ElementsVector, m_inputSpcs, m_lambdaRanges, m_models,
      m_RestLineList, m_continuumManager, m_ElementParam, m_curObs,
      enableAmplitudeOffsets, enableLambdaOffsetsFit);

  getSpectrumModel().m_enableAmplitudeOffsets = enableAmplitudeOffsets;
}

void CLineModelFitting::setLineRatioType(const std::string &lineRatioType) {
  m_lineRatioManager = CLineRatioManager::makeLineRatioManager(
      lineRatioType, m_ElementsVector, m_models, m_inputSpcs, m_lambdaRanges,
      m_continuumManager, m_RestLineList, m_fitter);
}

void CLineModelFitting::SetAbsLinesLimit(Float64 limit) {
  for (Int32 iElts = 0; iElts < getElementList().size(); iElts++) {
    getElementList()[iElts]->SetAbsLinesLimit(limit);
  }
}

/**
 * \brief Creates and returns a Mask with 0 in the lines support, 1 under the
 *lines
 **/
CMask CLineModelFitting::getOutsideLinesMask() const {

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
  const CSpectrumSpectralAxis &spcSpectralAxis =
      getSpectrum().GetSpectralAxis();
  const CSpectrumFluxAxis &Yspc = getSpectrumModel().getSpcFluxAxis();
  const auto &ErrorNoContinuum = getSpectrum().GetErrorAxis();

  Float64 fit = 0.0;
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
// TODO rename this ! not a simple getter
Int32 CLineModelFitting::getSpcNSamples() const {
  const CSpectrumSpectralAxis &spcSpectralAxis =
      getSpectrum().GetSpectralAxis();

  Float64 imin =
      spcSpectralAxis.GetIndexAtWaveLength(getLambdaRange().GetBegin());
  Float64 imax =
      spcSpectralAxis.GetIndexAtWaveLength(getLambdaRange().GetEnd());

  return abs(imax - imin);
}

/**
 * \brief Accumulates the squared differences between model and spectrum and
 *returns the sum.
 **/

Float64 CLineModelFitting::getLeastSquareMeritUnderElements() const {
  const CSpectrumFluxAxis &Yspc = getSpectrumModel().getSpcFluxAxis();
  const CSpectrumFluxAxis &Ymodel =
      getSpectrumModel().GetModelSpectrum().GetFluxAxis();
  const CSpectrumNoiseAxis &ErrorNoContinuum = getSpectrum().GetErrorAxis();

  Float64 fit = 0;
  Float64 diff = 0.0;

  TInt32RangeList support;
  for (Int32 iElts = 0; iElts < getElementList().size(); iElts++) {
    if (getElementList()[iElts]->IsOutsideLambdaRange())
      continue;

    TInt32RangeList s = getElementList()[iElts]->getSupport();
    for (Int32 iS = 0; iS < s.size(); iS++)
      support.push_back(s[iS]);
  }

  for (Int32 iS = 0; iS < support.size(); iS++) {
    for (Int32 j = support[iS].GetBegin(); j <= support[iS].GetEnd(); j++) {
      diff = (Yspc[j] - Ymodel[j]);
      fit += (diff * diff) / (ErrorNoContinuum[j] * ErrorNoContinuum[j]);
    }
  }
  return fit;
}

/**
 * \brief Returns the Stronger Multiple Emission Lines Amplitude Coefficient
 *(SMELAC)
 * 1. retrieve the lines amplitudes list for the Strong, and the Weak lines
 * 2. TODO: estimate the coefficient to be used as prior to penalise solutions
 *with less Strong lines
 **/
Float64 CLineModelFitting::getStrongerMultipleELAmpCoeff() const {
  TFloat64List AmpsStrong;
  TFloat64List AmpsWeak;
  Float64 sumAmps = 0.0;

  // Retrieve all the lines amplitudes in two lists (1 Strong, 1 weak)
  TInt32List validEltsIdx = getElementList().GetModelValidElementsIndexes();
  for (Int32 iElts : validEltsIdx) {
    auto const &elt = getElementList()[iElts];
    for (Int32 index = 0; index != elt->GetSize(); ++index) {
      auto const &line = elt->GetLines()[index];
      if (!line.IsEmission())
        continue;

      Float64 const amp = elt->GetFittedAmplitude(index);
      sumAmps += amp;
      if (line.IsStrong()) {
        AmpsStrong.push_back(amp);
      } else {
        AmpsWeak.push_back(amp);
      }
    }
  }

  Float64 sumAmpsStrong =
      std::reduce(AmpsStrong.cbegin(), AmpsStrong.cend(), 0.0);

  return sumAmpsStrong;
}

/**
 * \brief Returns the cumulative SNR under the Strong Emission Lines
 * 1. retrieve the lines support
 * 2. process each
 **/
std::pair<Float64, Float64> CLineModelFitting::getCumulSNRStrongEL() const {

  // Retrieve all the strone emission lines supports in a list of range
  TInt32RangeList supportList;
  TBoolList isStrongList;
  TInt32List validEltsIdx = getElementList().GetModelValidElementsIndexes();
  for (Int32 iElts : validEltsIdx) {
    auto const &elt = getElementList()[iElts];
    if (!elt->IsEmission())
      continue;
    for (Int32 index = 0; index != elt->GetSize(); ++index) {
      if (elt->IsOutsideLambdaRange(index))
        continue;
      auto const &line = elt->GetLines()[index];
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
return 0 if every was ok; else -1
*/
void CLineModelFitting::LoadModelSolution(
    const CLineModelSolution &modelSolution) {

  setRedshift(modelSolution.Redshift, false);

  auto &eltList = getElementList();

  if (m_enableAmplitudeOffsets)
    eltList.resetAmplitudeOffset();

  TBoolList element_done(eltList.size(), false);
  for (Int32 iRestLine = 0; iRestLine < m_RestLineList.size(); iRestLine++) {
    Int32 eIdx = modelSolution.ElementId[iRestLine];
    if (eIdx == undefIdx)
      continue;
    Int32 line_id = modelSolution.lineId[iRestLine];
    auto const &elt = eltList[eIdx];
    Int32 elt_line_index = elt->getLineIndex(line_id);
    if (elt_line_index == undefIdx)
      continue; // or throw an exception ?

    if (modelSolution.OutsideLambdaRange[iRestLine]) {
      elt->SetOutsideLambdaRangeList(elt_line_index);
      continue;
    }

    elt->SetFittedAmplitude(elt_line_index, modelSolution.Amplitudes[iRestLine],
                            modelSolution.AmplitudesUncertainties[iRestLine]);
    elt->SetOffset(elt_line_index, modelSolution.Offset[iRestLine]);

    if (element_done[eIdx])
      continue;

    elt->setVelocity(modelSolution.Velocity[iRestLine]);
    elt->SetFittingGroupInfo(modelSolution.fittingGroupInfo[iRestLine]);
    if (m_enableAmplitudeOffsets) {
      TPolynomCoeffs contPolynomCoeffs = {
          modelSolution.continuum_pCoeff0[iRestLine],
          modelSolution.continuum_pCoeff1[iRestLine],
          modelSolution.continuum_pCoeff2[iRestLine]};
      elt->SetPolynomCoeffs(std::move(contPolynomCoeffs));
    }

    element_done[eIdx] = true;
  }

  if (!std::isnan(modelSolution.LyaWidthCoeff) or
      !std::isnan(modelSolution.LyaAlpha) or
      !std::isnan(modelSolution.LyaDelta)) {

    std::string lyaTag = linetags::lya_em;
    auto const [idxLyaE, _] = eltList.findElementIndex(lyaTag);
    if (idxLyaE != undefIdx)
      eltList[idxLyaE]->SetAsymfitParams({modelSolution.LyaWidthCoeff,
                                          modelSolution.LyaAlpha,
                                          modelSolution.LyaDelta});
  }

  if (modelSolution.LyaIgm != undefIdx) {
    auto const indices_Igm = eltList.getIgmLinesIndices();
    if (!indices_Igm.empty())
      for (auto const &[elt_idx, _] : indices_Igm)
        eltList[elt_idx]->SetSymIgmParams(
            {modelSolution.LyaIgm, modelSolution.Redshift});
  }

  const CSpectrumSpectralAxis &spectralAxis = getSpectrum().GetSpectralAxis();
  for (Int32 iElts = 0; iElts < eltList.size(); iElts++) {
    eltList[iElts]->SetOutsideLambdaRange();

    if (!eltList[iElts]->IsOutsideLambdaRange())
      eltList[iElts]->prepareSupport(spectralAxis, modelSolution.Redshift,
                                     getLambdaRange());
  }

  return;
}

/**
 * \brief Returns a CLineModelSolution object populated with the current
 *solutions.
 **/
// this is not really a const method as spectrum model(s) have to be modified
// (cf CSpectrumModel::getContinuumError)
CLineModelSolution CLineModelFitting::GetModelSolution(Int32 opt_level) {
  Int32 s = m_RestLineList.size();
  CLineModelSolution modelSolution(m_RestLineList);

  auto const &eltList = getElementList();

  modelSolution.nDDL = eltList.GetModelNonZeroElementsNDdl();

  modelSolution.EmissionVelocity = m_ElementParam[0]->m_VelocityEmission;
  modelSolution.AbsorptionVelocity = m_ElementParam[0]->m_VelocityAbsorption;
  modelSolution.Redshift = getSpectrumModel().m_Redshift;

  TInt32List eIdx_oii;
  TInt32List subeIdx_oii;
  TInt32List eIdx_ha;
  TInt32List subeIdx_ha;
  Float64 flux_oii = 0.0;
  Float64 fluxVar_oii = 0.0;
  Float64 flux_ha = 0.0;
  Float64 fluxVar_ha = 0.0;

  for (Int32 iRestLine = 0; iRestLine < s; iRestLine++) {
    Int32 line_id = modelSolution.lineId[iRestLine];
    auto const &line = m_RestLineList.at(line_id);
    auto [eIdx, line_index] = eltList.findElementIndex(line_id);
    modelSolution.ElementId[iRestLine] = eIdx;
    if (eIdx == undefIdx || line_index == undefIdx ||
        eltList[eIdx]->IsOutsideLambdaRange(line_index)) {
      continue; // data already set to its default values
    }

    Float64 amp = eltList[eIdx]->GetFittedAmplitude(line_index);
    modelSolution.Amplitudes[iRestLine] = amp;
    Float64 ampError = eltList[eIdx]->GetFittedAmplitudeErrorSigma(line_index);
    modelSolution.AmplitudesUncertainties[iRestLine] = ampError;

    modelSolution.LambdaObs[iRestLine] =
        eltList[eIdx]->GetObservedPosition(line_index, modelSolution.Redshift);
    modelSolution.Velocity[iRestLine] = eltList[eIdx]->getVelocity();
    modelSolution.Offset[iRestLine] = eltList[eIdx]->GetOffset(line_index);

    if (opt_level) // brief, to save processing time, do not estimate fluxes
                   // and high level line properties
    {
      modelSolution.FittingError[iRestLine] =
          getSpectrumModel().getModelErrorUnderElement(
              eIdx, getSpectrumModel().getSpcFluxAxis());
      if (m_enableAmplitudeOffsets) {
        const auto &polynom_coeffs = eltList.getPolynomCoeffs(eIdx);
        modelSolution.continuum_pCoeff0[iRestLine] = polynom_coeffs.a0;
        modelSolution.continuum_pCoeff1[iRestLine] = polynom_coeffs.a1;
        modelSolution.continuum_pCoeff2[iRestLine] = polynom_coeffs.a2;
      }

      Float64 cont = eltList[eIdx]->GetContinuumAtCenterProfile(
          line_index, getSpectrum().GetSpectralAxis(), modelSolution.Redshift,
          getSpectrumModel().getContinuumFluxAxis(), m_enableAmplitudeOffsets);
      modelSolution.CenterContinuumFlux[iRestLine] = cont;
      modelSolution.ContinuumError[iRestLine] =
          getSpectrumModel().GetContinuumError(eIdx, line_index);
      Float64 mu = NAN;
      Float64 sigma = NAN;
      eltList[eIdx]->getObservedPositionAndLineWidth(
          line_index, modelSolution.Redshift, mu, sigma,
          false); // do not apply Lya asym offset

      Float64 flux = NAN;
      Float64 fluxError = NAN;
      Float64 fluxDI = NAN;
      Float64 snrDI = NAN;
      TInt32List eIdx_line(1, eIdx);
      TInt32List subeIdx_line(1, line_index);
      Int32 opt_cont_substract_abslinesmodel = 0;
      bool isEmission = false;
      if (line.GetType() == CLine::EType::nType_Emission) {
        opt_cont_substract_abslinesmodel = 1;
        isEmission = true;
      }
      getSpectrumModel().getFluxDirectIntegration(
          eIdx_line, subeIdx_line, opt_cont_substract_abslinesmodel, fluxDI,
          snrDI, getLambdaRange());
      if (!std::isnan(amp) && amp >= 0.0) {
        if (!isEmission) {
          amp *= cont;
          ampError *= cont;
        }
        const auto &profile = eltList[eIdx]->getLineProfile(line_index);

        Float64 lineFlux = profile->GetLineFlux(mu, sigma);
        flux = amp * lineFlux;
        fluxError = ampError * lineFlux;

        if (!isEmission)
          flux = -flux;
      }
      modelSolution.Sigmas[iRestLine] = sigma;
      modelSolution.Fluxs[iRestLine] = flux;
      modelSolution.FluxErrors[iRestLine] = fluxError;
      modelSolution.FluxDirectIntegration[iRestLine] = fluxDI;
      modelSolution.FluxDirectIntegrationError[iRestLine] =
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
          fluxDI = NAN;
          snrDI = NAN;
          Int32 opt_cont_substract_abslinesmodel = 0;
          getSpectrumModel().getFluxDirectIntegration(
              eIdx_ha, subeIdx_ha, opt_cont_substract_abslinesmodel, fluxDI,
              snrDI, getLambdaRange());
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
          fluxDI = NAN;
          snrDI = NAN;
          Int32 opt_cont_substract_abslinesmodel = 0;
          getSpectrumModel().getFluxDirectIntegration(
              eIdx_oii, subeIdx_oii, opt_cont_substract_abslinesmodel, fluxDI,
              snrDI, getLambdaRange());

          modelSolution.snrOII_DI = snrDI;
          modelSolution.lfOII_DI = fluxDI > 0 ? log10(fluxDI) : -INFINITY;
          modelSolution.lfOII = flux_oii > 0.0 ? log10(flux_oii) : -INFINITY;
          modelSolution.snrOII = flux_oii / std::sqrt(fluxVar_oii);
        }
      }
    }

    modelSolution.fittingGroupInfo[iRestLine] =
        eltList[eIdx]->GetFittingGroupInfo();
    modelSolution.OutsideLambdaRange[iRestLine] =
        eltList[eIdx]->IsOutsideLambdaRange(line_index);
  }

  // retrieve Lya params if fitted
  std::string lyaTag = linetags::lya_em;
  auto const [idxLyaE, _] = eltList.findElementIndex(lyaTag);
  if (idxLyaE != undefIdx) {
    TAsymParams params = eltList[idxLyaE]->GetAsymfitParams(0);
    modelSolution.LyaWidthCoeff = params.sigma;
    modelSolution.LyaAlpha = params.alpha;
    modelSolution.LyaDelta = params.delta;
    TSymIgmParams params_igm = eltList[idxLyaE]->GetSymIgmParams(0);
    modelSolution.LyaIgm = params_igm.m_igmidx;
  }

  std::unordered_set<std::string>
      strongELSNRAboveCut; // = getLinesAboveSNR(3.5);
  modelSolution.NLinesAboveSnrCut = strongELSNRAboveCut.size();

  return modelSolution;
}

/**
 * \brief Returns the size of getElementList().
 **/
Int32 CLineModelFitting::GetNElements() const {
  Int32 nddl = getElementList().size();
  return nddl;
}

void CLineModelFitting::SetLSF() {
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

void CLineModelFitting::SetVelocityEmission(Float64 vel) {

  for (Int32 j = 0; j < getElementList().size(); j++) {
    m_ElementParam[j]->m_VelocityEmission = vel;
  }
}

void CLineModelFitting::setVelocityEmissionByGroup(Float64 vel,
                                                   const TInt32List &inds) {

  for (auto idxElt : inds)
    m_ElementParam[idxElt]->m_VelocityEmission = vel;
}

void CLineModelFitting::SetVelocityAbsorption(Float64 vel) {

  for (Int32 j = 0; j < getElementList().size(); j++) {
    m_ElementParam[j]->m_VelocityAbsorption = vel;
  }
}

void CLineModelFitting::setVelocityAbsorptionByGroup(Float64 vel,
                                                     const TInt32List &inds) {

  for (auto idxElt : inds)
    m_ElementParam[idxElt]->m_VelocityAbsorption = vel;
}

Float64 CLineModelFitting::GetVelocityEmission() const {
  return m_ElementParam[0]->m_VelocityEmission;
}

Float64 CLineModelFitting::GetVelocityAbsorption() const {
  return m_ElementParam[0]->m_VelocityAbsorption;
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
  const CSpectrumSpectralAxis &spcSpectralAxis =
      getSpectrum().GetSpectralAxis();
  const CSpectrumFluxAxis &Yspc = getSpectrumModel().getSpcFluxAxis();
  const CSpectrumFluxAxis &YspcNoContinuum =
      getSpectrumModel().getSpcFluxAxisNoContinuum();
  const auto &ErrorNoContinuum = getSpectrum().GetErrorAxis();

  Float64 dtd = 0.0;
  Float64 flux = 0.0;

  Float64 imin =
      spcSpectralAxis.GetIndexAtWaveLength(getLambdaRange().GetBegin());
  Float64 imax =
      spcSpectralAxis.GetIndexAtWaveLength(getLambdaRange().GetEnd());
  for (Int32 j = imin; j < imax; j++) {
    if (spcComponent == "nocontinuum")
      flux = YspcNoContinuum[j];
    else
      flux = Yspc[j];

    dtd += (flux * flux) / (ErrorNoContinuum[j] * ErrorNoContinuum[j]);
  }
  Log.LogDebug("CLineModelFitting::EstimateDTransposeD val = %f", dtd);

  return dtd;
}

/**
 * \brief this function estimates the mtm value withing the wavelength range
 **/
Float64 CLineModelFitting::EstimateMTransposeM()
    const // duplicate with getMTranposeMCumulative, except for return values
{
  const CSpectrumSpectralAxis &spcSpectralAxis =
      getSpectrum().GetSpectralAxis();
  const CSpectrumFluxAxis &spcFluxAxis =
      getSpectrumModel().GetModelSpectrum().GetFluxAxis();
  const auto &ErrorNoContinuum = getSpectrum().GetErrorAxis();

  Float64 mtm = 0.0;
  Float64 diff = 0.0;

  Float64 imin =
      spcSpectralAxis.GetIndexAtWaveLength(getLambdaRange().GetBegin());
  Float64 imax =
      spcSpectralAxis.GetIndexAtWaveLength(getLambdaRange().GetEnd());
  for (Int32 j = imin; j < imax; j++) {
    diff = spcFluxAxis[j];
    mtm += (diff * diff) / (ErrorNoContinuum[j] * ErrorNoContinuum[j]);
  }
  // Log.LogDebug( "CLineModelFitting::EstimateMTransposeM val = %f", mtm );

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
  const CSpectrumSpectralAxis &spcSpectralAxis =
      getSpectrum().GetSpectralAxis();
  const auto &ErrorNoContinuum = getSpectrum().GetErrorAxis();

  Float64 cstLog = 0.0;
  Float64 sumLogNoise = 0.0;

  Int32 imin;
  Int32 imax;
  getLambdaRange().getClosedIntervalIndices(spcSpectralAxis.GetSamplesVector(),
                                            imin, imax);

  Int32 numDevs = std::abs(imax - imin + 1);
  for (Int32 j = imin; j <= imax; j++)
    sumLogNoise += log(ErrorNoContinuum[j]);

  // Log.LogDebug( "CLineModelFitting::EstimateMTransposeM val = %f", mtm );

  cstLog = -numDevs * 0.5 * log(2 * M_PI) - sumLogNoise;

  return cstLog;
}
