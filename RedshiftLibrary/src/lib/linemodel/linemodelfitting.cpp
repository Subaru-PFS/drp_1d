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

using namespace NSEpic;
using namespace std;

/**
 * \brief Prepares the state for Linemodel operation.
 * Loads the catalog.
 * Sets many state variables.
 * Sets the continuum either as a nocontinuum or a fromspectrum.
 **/
CLineModelFitting::CLineModelFitting()
    : m_RestLineList(Context.getLineVector()), m_enableAmplitudeOffsets(false) {
  initParameters();
  if (m_useloglambdasampling) {
    m_inputSpc = Context.GetRebinnedSpectrum();
    m_lambdaRange = Context.GetRebinnedClampedLambdaRange();
  } else {
    m_inputSpc = Context.GetSpectrum();
    m_lambdaRange = Context.GetClampedLambdaRange();
  }

  initMembers();
  setLineRatioType(m_lineRatioType);
  if (m_lineRatioType == "rules")
    dynamic_cast<CRulesManager *>(m_lineRatioManager.get())->setRulesOption();
}

CLineModelFitting::CLineModelFitting(
    const std::shared_ptr<const CSpectrum> &template_,
    const TLambdaRange &lambdaRange)
    : m_enableAmplitudeOffsets(false), m_RestLineList(Context.getLineVector()) {
  m_inputSpc = template_;
  m_lambdaRange = std::make_shared<const TLambdaRange>(lambdaRange);
  initParameters();
  // override ortho specific parameters
  m_fittingmethod = "hybrid";
  // temporary options override to be removed when full tpl ortho is implemented
  m_lineRatioType = "rules";

  initMembers();
  setLineRatioType(m_lineRatioType);

  dynamic_cast<CRulesManager *>(m_lineRatioManager.get())->setRulesOption("no");
  setContinuumComponent("fromspectrum");
}

void CLineModelFitting::initParameters() {
  CAutoScope autoscope(Context.m_ScopeStack, "linemodel");

  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();
  m_fittingmethod = ps->GetScoped<std::string>("fittingmethod");
  m_LineWidthType = ps->GetScoped<std::string>("linewidthtype");

  m_velocityEmission = ps->GetScoped<Float64>("velocityemission");
  m_velocityAbsorption = ps->GetScoped<Float64>("velocityabsorption");
  m_velocityEmissionInit = m_velocityEmission;
  m_velocityAbsorptionInit = m_velocityAbsorption;

  m_lineRatioType = ps->GetScoped<std::string>("lineRatioType");

  if (Context.GetCurrentMethod() == "LineModelSolve") {
    m_opt_firstpass_fittingmethod =
        ps->GetScoped<std::string>("firstpass.fittingmethod");
    m_opt_secondpass_fittingmethod = m_fittingmethod;
    m_opt_firstpass_forcedisableMultipleContinuumfit =
        ps->GetScoped<bool>("firstpass.multiplecontinuumfit_disable");
    // TODO dedicated function ?

    SetSecondpassContinuumFitPrms();
  }

  std::string continuumComponent =
      ps->GetScoped<std::string>("continuumcomponent");
  if (continuumComponent == "tplfit" || continuumComponent == "tplfitauto") {
    m_opt_firstpass_forcedisableMultipleContinuumfit =
        ps->GetScoped<bool>("firstpass.multiplecontinuumfit_disable");
    // useloglambdasampling param is relevant only if linemodel.continuumfit is
    // set to use fftprocessing below we explicit this check on this condition
    m_useloglambdasampling = ps->GetScoped<bool>("useloglambdasampling");
    m_useloglambdasampling &= ps->GetScoped<bool>("continuumfit.fftprocessing");
  }
}

void CLineModelFitting::initMembers() {
  // m_nominalWidthDefaultEmission = 1.15;// suited to new pfs simulations
  m_nominalWidthDefaultEmission = 13.4; // euclid 1 px
  m_nominalWidthDefaultAbsorption = m_nominalWidthDefaultEmission;

  m_enableLambdaOffsetsFit =
      true; // enable lambdaOffsetFit. Once enabled, the offset fixed value or
  // the fitting on/off switch is done through the offset calibration
  // file.

  Log.LogDetail("    model: Continuum winsize found is %.2f A",
                m_inputSpc->GetMedianWinsize());
  // Load the line catalog
  Log.LogDebug("About to load line catalog.");
  if (m_lineRatioType == "rules") {
    // load the regular catalog
    LoadCatalog(m_RestLineList);
  } else { //"tplratio" and "tplcorr"
    // load the tplratio catalog with only 1 element for all lines
    // LoadCatalogOneMultiline(restLineList);
    // load the tplratio catalog with 2 elements: 1 for the Em lines + 1 for the
    // Abs lines
    LoadCatalogTwoMultilinesAE(m_RestLineList);
  }

  m_model = std::make_shared<CSpectrumModel>(
      CSpectrumModel(m_Elements, m_inputSpc, m_RestLineList));
  m_continuumManager = std::make_shared<CContinuumManager>(
      CContinuumManager(m_inputSpc, m_model, *(m_lambdaRange)));

  SetFittingMethod(m_fittingmethod);
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
                << " m_enableLambdaOffsetsFit" << m_enableLambdaOffsetsFit);

  Log.LogDetail(Formatter()
                << " m_opt_firstpass_forcedisableMultipleContinuumfit"
                << m_opt_firstpass_forcedisableMultipleContinuumfit);
  Log.LogDetail(Formatter() << "m_opt_firstpass_fittingmethod "
                            << m_opt_firstpass_fittingmethod);
  Log.LogDetail(Formatter() << "m_opt_secondpass_fittingmethod"
                            << m_opt_secondpass_fittingmethod);

  Log.LogDetail(Formatter() << "LineWidthType=" << m_LineWidthType);

  Log.LogDetail(Formatter() << "velocityEmission=" << m_velocityEmission);
  Log.LogDetail(Formatter() << "velocityAbsorption=" << m_velocityAbsorption);
  Log.LogDetail(Formatter()
                << "velocityEmissionInit=" << m_velocityEmissionInit);
  Log.LogDetail(Formatter()
                << "velocityAbsorptionInit=" << m_velocityAbsorptionInit);

  Log.LogDetail(Formatter() << "nominalWidthDefaultEmission="
                            << m_nominalWidthDefaultEmission);
  Log.LogDetail(Formatter() << "nominalWidthDefaultAbsorption="
                            << m_nominalWidthDefaultAbsorption);

  Log.LogDetail(Formatter() << "fittingmethod=" << m_fittingmethod);

  Log.LogDetail(Formatter() << "lineRatioType=" << m_lineRatioType);

  // Log.LogDetail(Formatter()<<"tplCatalog="<<m_tplCatalog);
  // Log.LogDetail(Formatter()<<"tplCategoryList="<<m_tplCategoryList);

  Log.LogDetail(Formatter() << "secondpass_fitContinuum_dustfit="
                            << m_secondpass_fitContinuum_dustfit);
  Log.LogDetail(Formatter() << "secondpass_fitContinuum_igm="
                            << m_secondpass_fitContinuum_igm);

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
    SetFittingMethod(m_opt_secondpass_fittingmethod);
  }
  if (iPass == 3) {

    m_forcedisableMultipleContinuumfit = false;
    m_model->m_enableAmplitudeOffsets = true;
    m_enableAmplitudeOffsets = true;
    m_fitter->enableAmplitudeOffsets();
  }

  m_lineRatioManager->setPassMode(iPass);

  return true;
}
Int32 CLineModelFitting::GetPassNumber() const { return m_pass; }

/**
 * \brief For each line in each group of the argument, finds the associated
 *line in the catalog and saves this information to m_Elements. Converts the
 *argument restLineList to a group list. For each entry in this list: For each
 *line in this entry: Finds the index in the catalog from the line name and
 *type. Saves the line, the catalog index and the nominal amplitude for the
 *line thusly associated to this line. If at least one line was found, save
 *this result in m_Elements.
 **/
void CLineModelFitting::LoadCatalog(
    const CLineCatalog::TLineVector &restLineList) {

  std::vector<CLineCatalog::TLineVector> groupList =
      CLineCatalog::ConvertToGroupList(restLineList);
  for (Int32 ig = 0; ig < groupList.size(); ig++) {
    std::vector<CLine> lines;
    TFloat64List amps;
    TInt32List inds;
    for (Int32 i = 0; i < groupList[ig].size(); i++) {
      TInt32List idx = findLineIdxInCatalog(
          restLineList, groupList[ig][i].GetName(), groupList[ig][i].GetType());
      inds.push_back(idx[0]);
      amps.push_back(groupList[ig][i].GetNominalAmplitude());
      lines.push_back(groupList[ig][i]);
    }
    if (lines.size() > 0) {
      m_Elements.push_back(std::shared_ptr<CLineModelElement>(
          new CLineModelElement(lines, m_LineWidthType, m_velocityEmission,
                                m_velocityAbsorption, amps,
                                m_nominalWidthDefaultAbsorption, inds)));
    }
  }
}

void CLineModelFitting::LoadCatalogOneMultiline(
    const CLineCatalog::TLineVector &restLineList) {

  std::vector<CLine> lines;
  TFloat64List amps;
  TInt32List inds;
  for (Int32 ir = 0; ir < restLineList.size(); ir++) {
    inds.push_back(ir);
    amps.push_back(restLineList[ir].GetNominalAmplitude());
    lines.push_back(restLineList[ir]);
  }

  if (lines.size() > 0) {
    m_Elements.push_back(std::shared_ptr<CLineModelElement>(
        new CLineModelElement(lines, m_LineWidthType, m_velocityEmission,
                              m_velocityAbsorption, amps,
                              m_nominalWidthDefaultAbsorption, inds)));
  }
}

void CLineModelFitting::LoadCatalogTwoMultilinesAE(
    const CLineCatalog::TLineVector &restLineList) {
  std::vector<CLine::EType> types = {CLine::nType_Absorption,
                                     CLine::nType_Emission};

  for (Int32 iType = 0; iType < 2; iType++) {
    std::vector<CLine> lines;
    TFloat64List amps;
    TInt32List inds;
    for (Int32 ir = 0; ir < restLineList.size(); ir++) {
      if (restLineList[ir].GetType() == types[iType]) {
        inds.push_back(ir);
        amps.push_back(restLineList[ir].GetNominalAmplitude());
        lines.push_back(restLineList[ir]);
      }
    }

    if (lines.size() > 0) {
      m_Elements.push_back(std::shared_ptr<CLineModelElement>(
          new CLineModelElement(lines, m_LineWidthType, m_velocityEmission,
                                m_velocityAbsorption, amps,
                                m_nominalWidthDefaultAbsorption, inds)));
    }
  }
}

/**
 * \brief LogDetail the number of lines for each element, and their nominal
 *amplitudes.
 **/
void CLineModelFitting::LogCatalogInfos() {
  Log.LogDetail("\n");
  Log.LogDetail("LineModel Infos: %d elements", m_Elements.size());
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    Int32 nLines = m_Elements[iElts]->GetSize();
    if (nLines < 1) {
      Log.LogDetail("LineModel ctlg: elt %d (%s): no lines", iElts,
                    m_Elements[iElts]->GetElementTypeTag().c_str());
    }
    for (Int32 j = 0; j < nLines; j++) {
      std::string nominalAmpStr = "";
      if (nLines > 0) {
        nominalAmpStr = boost::str(boost::format("(nominal amp = %.4e)") %
                                   m_Elements[iElts]->GetNominalAmplitude(j));
      }
      Log.LogDetail("LineModel ctlg: elt %d (%s): line %d = %s %s", iElts,
                    m_Elements[iElts]->GetElementTypeTag().c_str(), j,
                    m_Elements[iElts]->GetLineName(j).c_str(),
                    nominalAmpStr.c_str());
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
  m_model->m_Redshift = redshift;

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

bool CLineModelFitting::initModelAtZ(
    Float64 redshift, const CSpectrumSpectralAxis &spectralAxis) {

  m_model->m_Redshift = redshift;

  // prepare the elements support
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    m_Elements[iElts]->resetAsymfitParams();
    m_Elements[iElts]->prepareSupport(spectralAxis, redshift, *(m_lambdaRange));
  }

  return true;
}

bool CLineModelFitting::initDtd() {
  //  m_dTransposeDLambdaRange = TLambdaRange(*(m_lambdaRange));
  m_dTransposeDLambdaRange = *(m_lambdaRange);
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
    m_model->setContinuumToInputSpc();
    return;
  }

  // the support has to be already computed
  // when LoadFitContinuum() is called
  m_continuumManager->initObserveGridContinuumFlux(
      m_inputSpc->GetSampleCount());
  Int32 autoselect = getContinuumComponent() == "tplfitauto";
  m_continuumManager->LoadFitContinuum(k, autoselect, redshift);
}

void CLineModelFitting::computeSpectrumFluxWithoutContinuum() {
  m_model->initModelWithContinuum();
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
                               Int32 contreest_iterations, bool enableLogging) {
  // initialize the model spectrum
  const CSpectrumSpectralAxis &spectralAxis = m_inputSpc->GetSpectralAxis();
  m_fitter->m_cont_reestim_iterations = contreest_iterations;
  initModelAtZ(redshift, spectralAxis);

  if (m_dTransposeDLambdaRange != *(m_lambdaRange))
    initDtd();

  Int32 ntplratio = m_lineRatioManager->prepareFit(
      redshift); // multiple fitting steps for lineRatioType=tplratio/tplratio

  Int32 nContinuum = 1;
  Int32 savedIdxContinuumFitted = -1; // for continuum tplfit
  if (isContinuumComponentTplfitxx() && !m_forcedisableMultipleContinuumfit)
    nContinuum = m_continuumManager->m_opt_fitcontinuum_maxCount;
  // 'on the fly' initialization
  Float64 bestMerit = INFINITY;
  Float64 bestMeritPrior = 0.0;

  for (Int32 k = 0; k < nContinuum; k++) {

    Float64 _merit = INFINITY;
    Float64 _meritprior = 0.; // only relevant for "tplratio"
    prepareAndLoadContinuum(k, redshift);
    if (getContinuumComponent() != "nocontinuum")
      computeSpectrumFluxWithoutContinuum();

    if (m_enableAmplitudeOffsets)
      m_Elements.prepareAmplitudeOffset();

    for (Int32 itratio = 0; itratio < ntplratio; itratio++) {

      if (m_lineRatioManager->init(redshift, itratio))
        continue;

      m_fitter->fit(redshift);

      std::string bestTplratioName = "undefined";

      _merit = m_lineRatioManager->computeMerit(itratio);

      if (bestMerit + bestMeritPrior > _merit + _meritprior) {
        bestMerit = _merit;
        bestMeritPrior = _meritprior;
        savedIdxContinuumFitted = k;
        Int32 modelSolutionLevel =
            m_lineRatioType == "rules" ? Int32(enableLogging) : 0;
        modelSolution = GetModelSolution(modelSolutionLevel);
        continuumModelSolution =
            m_continuumManager->GetContinuumModelSolution();

        m_lineRatioManager->saveResults(itratio);
      }
      if (getContinuumComponent() == "nocontinuum")
        m_model->reinitModel();
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
    /*    Log.LogDetail("    model - Linemodel: fitcontinuum = %d (%s, with "
                  "ebmv=%.3f), and A=%e",
                  savedIdxContinuumFitted, m_fitContinuum_tplName.c_str(),
                  m_fitContinuum_tplFitEbmvCoeff,
                  m_fitContinuum_tplFitAmplitude);
    */
  }
  if (m_lineRatioType == "tplratio") {
    m_lineRatioManager->finish(redshift);
    Int32 modelSolutionLevel = Int32(enableLogging);
    modelSolution = GetModelSolution(modelSolutionLevel);
    continuumModelSolution = m_continuumManager->GetContinuumModelSolution();
  }
  return bestMerit;
}

void CLineModelFitting::SetSecondpassContinuumFitPrms() {

  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();

  Int32 dustfit = -1;
  if (ps->GetScoped<bool>("continuumfit.ismfit"))
    dustfit = -10;

  Int32 meiksinfit = ps->GetScoped<bool>("continuumfit.igmfit");
  m_ignoreLinesSupport = ps->GetScoped<bool>("continuumfit.ignorelinesupport");
  m_secondpass_fitContinuum_dustfit = dustfit;
  m_secondpass_fitContinuum_igm = meiksinfit;

  Log.LogDetail("Elementlist: SetSecondpassContinuumFitPrms "
                "fitContinuum_dustfit = %d",
                m_secondpass_fitContinuum_dustfit);
  Log.LogDetail(
      "Elementlist: SetSecondpassContinuumFitPrms fitContinuum_igm = %d",
      m_secondpass_fitContinuum_igm);
}

void CLineModelFitting::SetFittingMethod(const std::string &fitMethod) {
  m_fittingmethod = fitMethod;
  m_fitter = CAbstractFitter::makeFitter(fitMethod, m_Elements, m_inputSpc,
                                         m_lambdaRange, m_model, m_RestLineList,
                                         m_continuumManager);
}

void CLineModelFitting::setLineRatioType(const std::string &lineRatioType) {
  m_lineRatioManager = CLineRatioManager::makeLineRatioManager(
      lineRatioType, m_Elements, m_model, m_inputSpc, m_lambdaRange,
      m_continuumManager, m_RestLineList);
}

void CLineModelFitting::SetAbsLinesLimit(Float64 limit) {
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    m_Elements[iElts]->SetAbsLinesLimit(limit);
  }
}

/**
 * \brief Creates and returns a Mask with 0 in the lines support, 1 under the
 *lines
 **/
CMask CLineModelFitting::getOutsideLinesMask() const {
  CMask _mask;
  // initialize the model spectrum
  const CSpectrumSpectralAxis &spectralAxis = m_inputSpc->GetSpectralAxis();
  _mask.SetSize(spectralAxis.GetSamplesCount());

  TInt32List validEltsIdx = m_Elements.GetModelValidElementsIndexes();
  TInt32List supportIdxes = m_Elements.getSupportIndexes(validEltsIdx);
  for (Int32 i = 0; i < spectralAxis.GetSamplesCount(); i++) {
    _mask[i] = 1;
  }
  // setting masks
  for (Int32 i = 0; i < supportIdxes.size(); i++) {
    _mask[supportIdxes[i]] = 0;
  }
  return _mask;
}

/**
 * \brief Estimates the STD outside the lines for the observed-model spectrum
 * NB: supposes the spectrum whithout continuum has a null mean value
 * input: which = 1: uses the spectrum flux continuum subtracted to compute
 *STD input: which = 2: uses the spectrum error to compute STD
 **/
Float64 CLineModelFitting::getOutsideLinesSTD(Int32 which) const {
  if (which != 1 && which != 2) {
    Log.LogError("    model: getOutsideLinesSTD - Failed to parse input "
                 "argument, which");
    return -1;
  }

  CMask _mask = getOutsideLinesMask();

  const CSpectrumSpectralAxis &spectralAxis = m_inputSpc->GetSpectralAxis();
  Float64 sum2 = 0.0;
  Int32 nsum = 0;
  Int32 imin = spectralAxis.GetIndexAtWaveLength(m_lambdaRange->GetBegin());
  Int32 imax = spectralAxis.GetIndexAtWaveLength(m_lambdaRange->GetEnd());
  const auto &spcFluxAxisNoContinuum = m_model->getSpcFluxAxisNoContinuum();
  const auto &ErrorNoContinuum = m_inputSpc->GetFluxAxis().GetError();
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
  const CSpectrumSpectralAxis &spcSpectralAxis = m_inputSpc->GetSpectralAxis();
  const CSpectrumFluxAxis &spcFluxAxis = m_model->getSpcFluxAxis();
  const auto &ErrorNoContinuum = m_inputSpc->GetFluxAxis().GetError();

  Int32 numDevs = 0;
  Float64 fit = 0.0;
  const Float64 *YCont = m_model->getContinuumFluxAxis().GetSamples();
  const Float64 *Yspc = spcFluxAxis.GetSamples();
  Float64 diff = 0.0;

  Float64 imin =
      spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange->GetBegin());
  Float64 imax = spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange->GetEnd());
  for (Int32 j = imin; j < imax; j++) {
    numDevs++;
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
  const CSpectrumSpectralAxis &spcSpectralAxis = m_inputSpc->GetSpectralAxis();

  Int32 numDevs = 0;

  Float64 imin =
      spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange->GetBegin());
  Float64 imax = spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange->GetEnd());

  for (Int32 j = imin; j < imax; j++) {
    numDevs++;
  }

  return numDevs;
}

/**
 * \brief Accumulates the squared differences between model and spectrum and
 *returns the sum.
 **/

Float64 CLineModelFitting::getLeastSquareMeritUnderElements() const {
  const CSpectrumFluxAxis &spcFluxAxis = m_model->getSpcFluxAxis();
  const CSpectrumFluxAxis &modelFluxAxis =
      m_model->GetModelSpectrum().GetFluxAxis();
  const CSpectrumNoiseAxis &ErrorNoContinuum =
      m_inputSpc->GetFluxAxis().GetError();

  Int32 numDevs = 0;
  Float64 fit = 0;
  const Float64 *Ymodel = modelFluxAxis.GetSamples();
  const Float64 *Yspc = spcFluxAxis.GetSamples();
  Float64 diff = 0.0;

  TInt32RangeList support;
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    if (m_Elements[iElts]->IsOutsideLambdaRange()) {
      continue;
    }
    TInt32RangeList s = m_Elements[iElts]->getSupport();
    for (Int32 iS = 0; iS < s.size(); iS++) {
      support.push_back(s[iS]);
    }
  }

  for (Int32 iS = 0; iS < support.size(); iS++) {
    for (Int32 j = support[iS].GetBegin(); j < support[iS].GetEnd(); j++) {
      numDevs++;
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
  TInt32List validEltsIdx = m_Elements.GetModelValidElementsIndexes();
  for (Int32 iValidElts = 0; iValidElts < validEltsIdx.size(); iValidElts++) {
    Int32 iElts = validEltsIdx[iValidElts];
    Int32 nlines = m_Elements[iElts]->GetLines().size();
    for (Int32 lineIdx = 0; lineIdx < nlines; lineIdx++) {
      if (!m_Elements[iElts]->m_Lines[lineIdx].GetIsEmission()) {
        continue;
      }

      Float64 amp = m_Elements[iElts]->GetFittedAmplitude(lineIdx);
      sumAmps += amp;
      if (m_Elements[iElts]->m_Lines[lineIdx].GetIsStrong()) {
        AmpsStrong.push_back(amp);
      } else {
        AmpsWeak.push_back(amp);
      }
    }
  }

  Float64 sumAmpsStrong = 0.0;
  for (Int32 k = 0; k < AmpsStrong.size(); k++) {
    sumAmpsStrong += AmpsStrong[k];
  }

  return sumAmpsStrong;
}

/**
 * \brief Returns the cumulative SNR under the Strong Emission Lines
 * 1. retrieve the lines support
 * 2. process each
 **/
Float64 CLineModelFitting::getCumulSNRStrongEL() const {

  // Retrieve all the liens supports in a list of range
  TInt32RangeList supportList;
  TInt32List validEltsIdx = m_Elements.GetModelValidElementsIndexes();
  for (Int32 iValidElts = 0; iValidElts < validEltsIdx.size(); iValidElts++) {
    Int32 iElts = validEltsIdx[iValidElts];
    Int32 nlines = m_Elements[iElts]->GetLines().size();
    for (Int32 lineIdx = 0; lineIdx < nlines; lineIdx++) {
      if (!m_Elements[iElts]->m_Lines[lineIdx].GetIsStrong()) {
        continue;
      }
      if (!m_Elements[iElts]->m_Lines[lineIdx].GetIsEmission()) {
        continue;
      }

      TInt32Range support =
          m_Elements[iElts]->getTheoreticalSupportSubElt(lineIdx);
      supportList.push_back(support);
    }
  }

  // merge overlapping ranges
  TInt32RangeList nonOverlappingSupportList;
  TInt32List processedSupport;
  for (Int32 k = 0; k < supportList.size(); k++) {
    // skip if already fitted
    bool alreadyProcessed = false;
    for (Int32 i = 0; i < processedSupport.size(); i++) {
      if (k == processedSupport[i]) {
        alreadyProcessed = true;
        break;
      }
    }
    if (alreadyProcessed) {
      continue;
    }

    processedSupport.push_back(k);
    TInt32Range support = supportList[k];

    for (Int32 l = 0; l < supportList.size(); l++) {
      // skip if already fitted
      bool alreadyProcessed = false;
      for (Int32 i = 0; i < processedSupport.size(); i++) {
        if (l == processedSupport[i]) {
          alreadyProcessed = true;
          break;
        }
      }
      if (alreadyProcessed) {
        continue;
      }

      // try if current range is bluer than l and overlaps ?
      Float64 xinf = support.GetBegin();
      Float64 xsup = support.GetEnd();
      Float64 yinf = supportList[l].GetBegin();
      Float64 ysup = supportList[l].GetEnd();
      Float64 max = std::max(xinf, yinf);
      Float64 min = std::min(xsup, ysup);
      if (max - min < 0) {
        processedSupport.push_back(l);
        support.SetBegin(std::min(xinf, yinf));
        support.SetEnd(std::max(xsup, ysup));
      }
    }
    nonOverlappingSupportList.push_back(support);
  }

  TFloat64List snrList;
  Float64 sumSNR = 0.0;
  // process SNR on the non overlapping ranges
  for (Int32 k = 0; k < nonOverlappingSupportList.size(); k++) {
    snrList.push_back(getCumulSNROnRange(nonOverlappingSupportList[k]));
    sumSNR += snrList[k];
  }
  std::sort(snrList.rbegin(), snrList.rend());

  // compute the snr metric
  TInt32List snrIsrelevantList;
  for (Int32 k = 0; k < snrList.size(); k++) {
    snrIsrelevantList.push_back(0);
  }
  Float64 thresRatio = 0.8;
  Float64 curRatio = 0.0;
  for (Int32 k = 0; k < snrList.size(); k++) {
    // relevant if snr>8.0
    if (snrList[k] > 8.0) {
      snrIsrelevantList[k] = 1;
    }

    // relevant if contributes to the 'thresRatio'*100 percent (ex. 80%) of
    // the SumSNR value
    curRatio += snrList[k] / sumSNR;
    if (curRatio <= thresRatio) {
      snrIsrelevantList[k] = 1;
    }
  }

  Float64 sumIsRelevant = 0.0;
  for (Int32 k = 0; k < snrIsrelevantList.size(); k++) {
    sumIsRelevant += snrIsrelevantList[k];
  }

  Float64 snrMetric = sumSNR * sumIsRelevant;
  return snrMetric;
}

/**
 * \brief Returns the cumulative SNR on the idxRange
 **/
Float64 CLineModelFitting::getCumulSNROnRange(TInt32Range idxRange) const {
  Int32 n = idxRange.GetEnd() - idxRange.GetBegin() + 1;
  if (n < 2) {
    return -1;
  }

  const CSpectrumFluxAxis &modelFluxAxis =
      m_model->GetModelSpectrum().GetFluxAxis();
  const Float64 *Ymodel = modelFluxAxis.GetSamples();
  const auto &ErrorNoContinuum = m_inputSpc->GetFluxAxis().GetError();
  const auto &ContinuumFluxAxis = m_model->getContinuumFluxAxis();

  Int32 idx = 0;
  Float64 sumF = 0.0;
  Float64 sumM = 0.0;
  for (Int32 i = 0; i < n; i++) {
    idx = i + idxRange.GetBegin();
    Float64 flux =
        Ymodel[idx] - ContinuumFluxAxis[idx]; // using only the no-continuum
                                              // component to estimate SNR
    sumF += flux;
    sumM += ErrorNoContinuum[idx] * ErrorNoContinuum[idx];
  }
  Float64 Err = std::sqrt(sumM);
  Float64 rangeSNR = sumF / Err;

  return rangeSNR;
}

/**
 * \brief Search the line catalog for lines whose name match the argument
 *strTag.
 **/
TInt32List CLineModelFitting::findLineIdxInCatalog(
    const CLineCatalog::TLineVector &restLineList, const std::string &strTag,
    Int32 type) const {
  TInt32List indexes;
  for (Int32 iRestLine = 0; iRestLine < restLineList.size(); iRestLine++) {
    if (restLineList[iRestLine].GetType() != type) {
      continue;
    }
    std::string name = restLineList[iRestLine].GetName();
    std::size_t foundstra = name.find(strTag.c_str());
    if (foundstra != std::string::npos) {
      indexes.push_back(iRestLine);
    }
  }
  return indexes;
}

/**
 * \brief Adds an entry to m_Elements as a CLineModelElement constructed from
 *the arguments.
 **/
void CLineModelFitting::addDoubleLine(const CLine &r1, const CLine &r2,
                                      Int32 index1, Int32 index2,
                                      Float64 nominalWidth, Float64 a1,
                                      Float64 a2) {
  std::vector<CLine> lines;
  lines.push_back(r1);
  lines.push_back(r2);
  TFloat64List amps;
  amps.push_back(a1);
  amps.push_back(a2);
  TInt32List a;
  a.push_back(index1);
  a.push_back(index2);
  m_Elements.push_back(std::shared_ptr<CLineModelElement>(
      new CLineModelElement(lines, m_LineWidthType, m_velocityEmission,
                            m_velocityAbsorption, amps, nominalWidth, a)));
}

/*
Reset all the model value to the previous solution found.
return 0 if every was ok; else -1
*/
void CLineModelFitting::LoadModelSolution(
    const CLineModelSolution &modelSolution) {

  SetVelocityEmission(modelSolution.EmissionVelocity);
  SetVelocityAbsorption(modelSolution.AbsorptionVelocity);

  const CSpectrumSpectralAxis &spectralAxis = m_inputSpc->GetSpectralAxis();
  for (Int32 iRestLine = 0; iRestLine < m_RestLineList.size(); iRestLine++) {
    Int32 eIdx = modelSolution.ElementId[iRestLine];
    if (eIdx == undefIdx)
      continue;

    m_Elements[eIdx]->prepareSupport(spectralAxis, modelSolution.Redshift,
                                     *(m_lambdaRange));
  }

  if (m_enableAmplitudeOffsets)
    m_Elements.prepareAmplitudeOffset();

  for (Int32 iRestLine = 0; iRestLine < m_RestLineList.size(); iRestLine++) {
    Int32 eIdx = modelSolution.ElementId[iRestLine];
    if (eIdx == undefIdx)
      continue;
    if (m_enableAmplitudeOffsets) {
      TPolynomCoeffs contPolynomCoeffs = {
          modelSolution.continuum_pCoeff0[iRestLine],
          modelSolution.continuum_pCoeff1[iRestLine],
          modelSolution.continuum_pCoeff2[iRestLine]};
      applyPolynomCoeffs(eIdx, contPolynomCoeffs);
    }

    Int32 subeIdx = m_Elements[eIdx]->findElementIndex(iRestLine);

    m_Elements[eIdx]->SetFittedAmplitude(
        subeIdx, modelSolution.Amplitudes[iRestLine],
        modelSolution.AmplitudesUncertainties[iRestLine]);
  }

  if (!std::isnan(modelSolution.LyaWidthCoeff) or
      !std::isnan(modelSolution.LyaAlpha) or
      !std::isnan(modelSolution.LyaDelta)) {

    std::string lyaTag = linetags::lya_em;
    Int32 idxLyaE = m_Elements.findElementIndex(lyaTag);
    if (idxLyaE != undefIdx)
      m_Elements[idxLyaE]->SetAsymfitParams({modelSolution.LyaWidthCoeff,
                                             modelSolution.LyaAlpha,
                                             modelSolution.LyaDelta});
  }

  if (modelSolution.LyaIgm != undefIdx) {
    auto const idxEltIGM = m_Elements.getIgmLinesIndices().front();
    if (!idxEltIGM.empty())
      for (auto const iElt : idxEltIGM)
        m_Elements[iElt]->SetSymIgmParams(
            {modelSolution.LyaIgm, modelSolution.Redshift});
  }
  return;
}

/**
 * \brief Returns a CLineModelSolution object populated with the current
 *solutions.
 **/
CLineModelSolution CLineModelFitting::GetModelSolution(Int32 opt_level) {
  Int32 s = m_RestLineList.size();
  CLineModelSolution modelSolution(m_RestLineList);
  modelSolution.nDDL = m_Elements.GetModelNonZeroElementsNDdl();

  modelSolution.EmissionVelocity = m_velocityEmission;
  modelSolution.AbsorptionVelocity = m_velocityAbsorption;
  modelSolution.Redshift = m_model->m_Redshift;

  TInt32List eIdx_oii;
  TInt32List subeIdx_oii;

  for (Int32 iRestLine = 0; iRestLine < s; iRestLine++) {
    Int32 subeIdx = undefIdx;
    Int32 eIdx = m_Elements.findElementIndex(iRestLine, subeIdx);
    modelSolution.ElementId[iRestLine] = eIdx;

    if (eIdx == undefIdx || subeIdx == undefIdx ||
        m_Elements[eIdx]->IsOutsideLambdaRange(subeIdx)) {
      continue; // data already set to its default values
    }
    Float64 amp = m_Elements[eIdx]->GetFittedAmplitude(subeIdx);
    modelSolution.Amplitudes[iRestLine] = amp;
    Float64 ampError = m_Elements[eIdx]->GetFittedAmplitudeErrorSigma(subeIdx);
    modelSolution.AmplitudesUncertainties[iRestLine] = ampError;

    modelSolution.LambdaObs[iRestLine] =
        m_Elements[eIdx]->GetObservedPosition(subeIdx, modelSolution.Redshift);
    modelSolution.Velocity[iRestLine] = m_Elements[eIdx]->GetVelocity();
    modelSolution.Offset[iRestLine] =
        m_Elements[eIdx]->m_Lines[subeIdx].GetOffset();

    if (opt_level) // brief, to save processing time, do not estimate fluxes
                   // and high level line properties
    {
      modelSolution.FittingError[iRestLine] =
          m_model->getModelErrorUnderElement(eIdx, m_model->getSpcFluxAxis());
      TPolynomCoeffs polynom_coeffs = m_Elements.getPolynomCoeffs(eIdx);
      // save polynom info to output them in hdf5, mainly to recontruct
      // linemeas model
      modelSolution.continuum_pCoeff0[iRestLine] = polynom_coeffs.x0;
      modelSolution.continuum_pCoeff1[iRestLine] = polynom_coeffs.x1;
      modelSolution.continuum_pCoeff2[iRestLine] = polynom_coeffs.x2;

      Float64 cont = m_Elements[eIdx]->GetContinuumAtCenterProfile(
          subeIdx, m_inputSpc->GetSpectralAxis(), modelSolution.Redshift,
          m_model->getContinuumFluxAxis(), polynom_coeffs);
      modelSolution.CenterContinuumFlux[iRestLine] = cont;
      modelSolution.ContinuumError[iRestLine] =
          m_model->GetContinuumError(eIdx, subeIdx);
      Float64 mu = NAN;
      Float64 sigma = NAN;
      m_Elements[eIdx]->getObservedPositionAndLineWidth(
          subeIdx, modelSolution.Redshift, mu, sigma,
          false); // do not apply Lya asym offset

      Float64 flux = NAN;
      Float64 fluxError = NAN;
      Float64 fluxDI = NAN;
      Float64 snrDI = NAN;
      TInt32List eIdx_line(1, eIdx);
      TInt32List subeIdx_line(1, subeIdx);
      Int32 opt_cont_substract_abslinesmodel = 0;
      bool isEmission = false;
      if (m_RestLineList[iRestLine].GetType() == CLine::nType_Emission) {
        opt_cont_substract_abslinesmodel = 1;
        isEmission = true;
      }
      m_model->getFluxDirectIntegration(eIdx_line, subeIdx_line,
                                        opt_cont_substract_abslinesmodel,
                                        fluxDI, snrDI, *(m_lambdaRange));
      if (!std::isnan(amp) && amp >= 0.0) {
        if (!isEmission) {
          amp *= cont;
          ampError *= cont;
        }
        const CLineProfile &profile = m_Elements[eIdx]->getLineProfile(subeIdx);

        Float64 lineFlux = profile.GetLineFlux(mu, sigma);
        flux = amp * lineFlux;
        fluxError = ampError * lineFlux;

        if (!isEmission)
          flux = -flux;
      }
      modelSolution.Sigmas[iRestLine] = sigma;
      modelSolution.Fluxs[iRestLine] = flux;
      modelSolution.FluxErrors[iRestLine] = fluxError;
      modelSolution.FluxDirectIntegration[iRestLine] = fluxDI;

      // rough estimation of SNR_Ha, using the given model and fitting method
      //(warning: Ha flux and error could have been obtained by global fitting
      // of the model, which leads to different results than fitted
      // individually...)
      bool directIntegration = true;
      if (isEmission &&
          m_RestLineList[iRestLine].GetName() == linetags::halpha_em) {
        if (directIntegration) {
          modelSolution.snrHa = snrDI;
          if (fluxDI > 0.0)
            modelSolution.lfHa = log10(fluxDI);
          else if (fluxDI == 0.)
            modelSolution.lfHa = -INFINITY;
        } else {
          if (fluxError > 0.0)
            modelSolution.snrHa = flux / fluxError;
          if (flux > 0.0)
            modelSolution.lfHa = log10(flux);
        }
      }
      if (isEmission &&
          (m_RestLineList[iRestLine].GetName() == linetags::oII3726_em ||
           m_RestLineList[iRestLine].GetName() == linetags::oII3729_em)) {
        // here we only cover the fluxDI case.
        eIdx_oii.push_back(eIdx);
        subeIdx_oii.push_back(subeIdx);
        if (directIntegration) {
          fluxDI = NAN;
          snrDI = NAN;
          Int32 opt_cont_substract_abslinesmodel = 0;
          m_model->getFluxDirectIntegration(eIdx_oii, subeIdx_oii,
                                            opt_cont_substract_abslinesmodel,
                                            fluxDI, snrDI, *(m_lambdaRange));

          modelSolution.snrOII = snrDI;
          if (fluxDI > 0.0)
            modelSolution.lfOII = log10(fluxDI);
          else if (fluxDI == 0.)
            modelSolution.lfOII = -INFINITY;
        }
      }
    }

    modelSolution.fittingGroupInfo[iRestLine] =
        m_Elements[eIdx]->m_fittingGroupInfo;
    modelSolution.OutsideLambdaRange[iRestLine] =
        m_Elements[eIdx]->IsOutsideLambdaRange(subeIdx);
  }

  // retrieve Lya params if fitted
  std::string lyaTag = linetags::lya_em;
  Int32 idxLyaE = m_Elements.findElementIndex(lyaTag);
  if (idxLyaE != undefIdx) {
    TAsymParams params = m_Elements[idxLyaE]->GetAsymfitParams(0);
    modelSolution.LyaWidthCoeff = params.sigma;
    modelSolution.LyaAlpha = params.alpha;
    modelSolution.LyaDelta = params.delta;
    TSymIgmParams params_igm = m_Elements[idxLyaE]->GetSymIgmParams(0);
    modelSolution.LyaIgm = params_igm.m_igmidx;
  }

  std::unordered_set<std::string>
      strongELSNRAboveCut; // = getLinesAboveSNR(3.5);
  modelSolution.NLinesAboveSnrCut = strongELSNRAboveCut.size();

  return modelSolution;
}

/**
 * @brief Add polynom coeffs into the class variables, for all correctly in
 * m_ampOffsetCoeffs
 *
 * @param eIdx
 * @param polynom_coeffs
 * @return
 */
void CLineModelFitting::applyPolynomCoeffs(
    Int32 eIdx, const TPolynomCoeffs &polynom_coeffs) {
  TInt32List xInds = m_Elements.getSupportIndexes({Int32(eIdx)});
  if (xInds.size() && m_Elements.m_ampOffsetsCoeffs.size()) {
    Int32 idxAmpOffset = m_Elements.getIndexAmpOffset(xInds[0]);
    m_Elements.m_ampOffsetsCoeffs[idxAmpOffset] = polynom_coeffs;
  }
  return;
}

/**
 * \brief Returns the size of m_Elements.
 **/
Int32 CLineModelFitting::GetNElements() const {
  Int32 nddl = m_Elements.size();
  return nddl;
}

void CLineModelFitting::SetLSF() {
  const std::shared_ptr<const CLSF> &lsf = m_inputSpc->GetLSF();

  if (lsf == nullptr) {
    THROWG(INTERNAL_ERROR,
           "Cannot enable LSF, LSF spectrum member is not initialized");
  } else if (!lsf->IsValid()) {
    THROWG(INTERNAL_ERROR,
           " Cannot enable LSF, LSF spectrum member is not valid");
  }

  for (Int32 j = 0; j < m_Elements.size(); j++) {
    m_Elements[j]->SetLSF(
        lsf); // lsf has now a type to be used for width computations
  }
}

void CLineModelFitting::SetVelocityEmission(Float64 vel) {
  m_velocityEmission = vel;
  for (Int32 j = 0; j < m_Elements.size(); j++) {
    m_Elements[j]->SetVelocityEmission(vel);
  }
}

void CLineModelFitting::setVelocity(Float64 vel, Int32 lineType) {
  if (lineType == CLine::nType_Absorption)
    m_velocityEmission = vel;
  else if (lineType == CLine::nType_Emission)
    m_velocityAbsorption = vel;
  else {
    m_velocityAbsorption = vel;
    m_velocityEmission = vel;
  }
  for (Int32 j = 0; j < m_Elements.size(); j++) {
    m_Elements[j]->setVelocity(vel);
  }
}

// TODO lineType may be useless, because element are homogeneous in lineType,
// but maybe m_velocityEmission and m_velocityAbsorption are useful ?
void CLineModelFitting::setVelocity(Float64 vel, Int32 idxElt, Int32 lineType) {
  if (lineType == CLine::nType_Absorption)
    m_velocityEmission = vel;
  else if (lineType == CLine::nType_Emission)
    m_velocityAbsorption = vel;
  else {
    m_velocityAbsorption = vel;
    m_velocityEmission = vel;
  }
  if (idxElt < m_Elements.size()) {
    m_Elements[idxElt]->setVelocity(vel);
  } else
    THROWG(INTERNAL_ERROR,
           Formatter() << "Wrong index for line model element " << idxElt);
}

void CLineModelFitting::SetVelocityEmissionOneElement(Float64 vel,
                                                      Int32 idxElt) {
  m_velocityEmission = vel;
  if (idxElt < m_Elements.size()) {
    m_Elements[idxElt]->SetVelocityEmission(vel);
  }
}

void CLineModelFitting::SetVelocityAbsorption(Float64 vel) {
  m_velocityAbsorption = vel;
  for (Int32 j = 0; j < m_Elements.size(); j++) {
    m_Elements[j]->SetVelocityAbsorption(vel);
  }
}

void CLineModelFitting::SetVelocityAbsorptionOneElement(Float64 vel,
                                                        Int32 idxElt) {
  m_velocityAbsorption = vel;
  if (idxElt < m_Elements.size()) {
    m_Elements[idxElt]->SetVelocityAbsorption(vel);
  }
}

Float64 CLineModelFitting::GetVelocityEmission() const {
  return m_velocityEmission;
}

Float64 CLineModelFitting::GetVelocityAbsorption() const {
  return m_velocityAbsorption;
}

Float64 CLineModelFitting::GetRedshift() const { return m_model->m_Redshift; }

Int32 CLineModelFitting::ApplyVelocityBound(Float64 inf, Float64 sup) {

  Int32 corrected = false;
  static Float64 velInfFromInstrument = inf;
  static Float64 velSupEmission = sup;
  static Float64 velSupAbsorption = sup;

  Float64 vel;
  vel = GetVelocityEmission();
  if (vel > velSupEmission || vel < velInfFromInstrument) {
    SetVelocityEmission(m_velocityEmissionInit);
    corrected = true;
    Log.LogInfo("\nLineModel Infos: Reset Velocity Emission, to v = %.1f",
                m_velocityEmissionInit);
  }
  vel = GetVelocityAbsorption();
  if (vel > velSupAbsorption || vel < velInfFromInstrument) {
    SetVelocityAbsorption(m_velocityAbsorptionInit);
    corrected = true;
    Log.LogInfo("\nLineModel Infos: Reset Velocity Absorption, to v = %.1f",
                m_velocityAbsorptionInit);
  }
  return corrected;
}

/**
 * \brief this function returns the dtd value withing the wavelength range for
 *a given spcComponent
 *
 **/
// TODO rename this ! not a simple getter
Float64 CLineModelFitting::getDTransposeD() {
  if (m_dTransposeDLambdaRange != *(m_lambdaRange)) {
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
  if (m_dTransposeDLambdaRange != *(m_lambdaRange)) {
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
  const CSpectrumSpectralAxis &spcSpectralAxis = m_inputSpc->GetSpectralAxis();
  const CSpectrumFluxAxis &spcFluxAxis = m_model->getSpcFluxAxis();
  const CSpectrumFluxAxis &spcFluxAxisNoContinuum =
      m_model->getSpcFluxAxisNoContinuum();
  const auto &ErrorNoContinuum = m_inputSpc->GetFluxAxis().GetError();

  Int32 numDevs = 0;
  Float64 dtd = 0.0;
  const Float64 *Yspc = spcFluxAxis.GetSamples();
  const Float64 *YspcNoContinuum = spcFluxAxisNoContinuum.GetSamples();
  Float64 flux = 0.0;

  Float64 imin =
      spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange->GetBegin());
  Float64 imax = spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange->GetEnd());
  for (Int32 j = imin; j < imax; j++) {
    numDevs++;
    if (spcComponent == "nocontinuum") {
      flux = YspcNoContinuum[j];
    } else {
      flux = Yspc[j];
    }
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
  const CSpectrumSpectralAxis &spcSpectralAxis = m_inputSpc->GetSpectralAxis();
  const CSpectrumFluxAxis &spcFluxAxis =
      m_model->GetModelSpectrum().GetFluxAxis();
  const auto &ErrorNoContinuum = m_inputSpc->GetFluxAxis().GetError();

  Int32 numDevs = 0;
  Float64 mtm = 0.0;
  const Float64 *Yspc = spcFluxAxis.GetSamples();
  Float64 diff = 0.0;

  Float64 imin =
      spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange->GetBegin());
  Float64 imax = spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange->GetEnd());
  for (Int32 j = imin; j < imax; j++) {
    numDevs++;
    { diff = Yspc[j]; }
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
  const CSpectrumSpectralAxis &spcSpectralAxis = m_inputSpc->GetSpectralAxis();
  const auto &ErrorNoContinuum = m_inputSpc->GetFluxAxis().GetError();

  Int32 numDevs = 0;
  Float64 cstLog = 0.0;
  Float64 sumLogNoise = 0.0;

  Float64 imin =
      spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange->GetBegin());
  Float64 imax = spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange->GetEnd());
  for (Int32 j = imin; j < imax; j++) {
    numDevs++;
    sumLogNoise += log(ErrorNoContinuum[j]);
  }
  // Log.LogDebug( "CLineModelFitting::EstimateMTransposeM val = %f", mtm );

  cstLog = -numDevs * 0.5 * log(2 * M_PI) - sumLogNoise;

  return cstLog;
}
