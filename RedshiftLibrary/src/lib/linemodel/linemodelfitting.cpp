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
#include "RedshiftLibrary/linemodel/hybridfitter.h"
#include "RedshiftLibrary/linemodel/individualfitter.h"
#include "RedshiftLibrary/linemodel/onesfitter.h"
#include "RedshiftLibrary/linemodel/randomfitter.h"
#include "RedshiftLibrary/linemodel/svdfitter.h"
#include "RedshiftLibrary/linemodel/svdlcfitter.h"

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

#include <float.h>
#include <gsl/gsl_blas.h>

#include <gsl/gsl_interp.h>

#include <gsl/gsl_spline.h>
#include <math.h>

#include <boost/numeric/conversion/bounds.hpp>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_vector.h>
using namespace NSEpic;
using namespace std;

/**
 * \brief Prepares the state for Linemodel operation.
 * Loads the catalog.
 * Sets many state variables.
 * Sets the continuum either as a nocontinuum or a fromspectrum.
 **/
CLineModelFitting::CLineModelFitting()
    : m_RestLineList(Context.getLineVector()), m_Regulament(),
      m_enableAmplitudeOffsets(false) {
  initParameters();
  if (m_useloglambdasampling) {
    m_inputSpc = Context.GetRebinnedSpectrum();
    m_lambdaRange = Context.GetRebinnedClampedLambdaRange();
  } else {
    m_inputSpc = Context.GetSpectrum();
    m_lambdaRange = Context.GetClampedLambdaRange();
  }

  initMembers();
}

CLineModelFitting::CLineModelFitting(
    const std::shared_ptr<const CSpectrum> &template_,
    const TLambdaRange &lambdaRange)
    : m_Regulament(), m_enableAmplitudeOffsets(false),
      m_RestLineList(Context.getLineVector()) {
  m_inputSpc = template_;
  m_lambdaRange = std::make_shared<const TLambdaRange>(lambdaRange);
  initParameters();
  // override ortho specific parameters
  m_ContinuumComponent = "fromspectrum";
  m_fittingmethod = "hybrid";
  // temporary options override to be removed when full tpl ortho is implemented
  m_rigidity = "rules";
  m_rulesoption = "no";
  initMembers();
}

void CLineModelFitting::initParameters() {
  CAutoScope autoscope(Context.m_ScopeStack, "linemodel");

  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();
  m_fittingmethod = ps->GetScoped<std::string>("fittingmethod");
  m_ContinuumComponent = ps->GetScoped<std::string>("continuumcomponent");

  m_LineWidthType = ps->GetScoped<std::string>("linewidthtype");
  m_NSigmaSupport = ps->GetScoped<Float64>("nsigmasupport");
  m_velocityEmission = ps->GetScoped<Float64>("velocityemission");
  m_velocityAbsorption = ps->GetScoped<Float64>("velocityabsorption");
  m_velocityEmissionInit = m_velocityEmission;
  m_velocityAbsorptionInit = m_velocityAbsorption;
  m_rulesoption = ps->GetScoped<std::string>("rules");
  m_rigidity = ps->GetScoped<std::string>("rigidity");

  TStringList rigidityValues = {"tplratio", "rules", "tplcorr"};
  if (std::find(rigidityValues.begin(), rigidityValues.end(), m_rigidity) ==
      rigidityValues.end())
    THROWG(INVALID_PARAMETER, "Only {tplratio, rules, tpcorr} values are "
                              "supported for linemodel.rigidity");

  m_opt_lya_forcedisablefit = ps->GetScoped<bool>("lyaforcedisablefit");
  if (m_rigidity == "tplratio")
    m_opt_haprior = ps->GetScoped<Float64>("haprior");
  // TODO this restriction should be checked elsewhere
  /*
  if (m_rigidity == "rules" && m_fittingmethod == "hybrid")
    m_opt_enable_improveBalmerFit = ps->GetScoped<bool>("improveBalmerFit");
  */
  if (Context.GetCurrentMethod() == "LineModelSolve") {
    m_opt_firstpass_fittingmethod =
        ps->GetScoped<std::string>("firstpass.fittingmethod");
    m_opt_secondpass_fittingmethod = m_fittingmethod;
    m_opt_firstpass_forcedisableMultipleContinuumfit =
        ps->GetScoped<bool>("firstpass.multiplecontinuumfit_disable");
    m_opt_firstpass_forcedisableTplratioISMfit =
        !ps->GetScoped<bool>("firstpass.tplratio_ismfit");
    // TODO dedicated function ?
    m_opt_lya_forcefit = ps->GetScoped<bool>("lyaforcefit");

    m_opt_lya_fit_asym_min = ps->GetScoped<Float64>("lyafit.asymfitmin");
    m_opt_lya_fit_asym_max = ps->GetScoped<Float64>("lyafit.asymfitmax");
    m_opt_lya_fit_asym_step = ps->GetScoped<Float64>("lyafit.asymfitstep");
    m_opt_lya_fit_width_min = ps->GetScoped<Float64>("lyafit.widthfitmin");
    m_opt_lya_fit_width_max = ps->GetScoped<Float64>("lyafit.widthfitmax");
    m_opt_lya_fit_width_step = ps->GetScoped<Float64>("lyafit.widthfitstep");
    m_opt_lya_fit_delta_min = ps->GetScoped<Float64>("lyafit.deltafitmin");
    m_opt_lya_fit_delta_max = ps->GetScoped<Float64>("lyafit.deltafitmax");
    m_opt_lya_fit_delta_step = ps->GetScoped<Float64>("lyafit.deltafitstep");

    SetSecondpassContinuumFitPrms();
  }
  if (isContinuumComponentTplfitxx()) {
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
  m_model =
      std::make_shared<CSpectrumModel>(CSpectrumModel(m_Elements, m_inputSpc));
  m_continuumManager = std::make_shared<CContinuumManager>(
      CContinuumManager(m_inputSpc, m_model, *(m_lambdaRange)));
  SetContinuumComponent(m_ContinuumComponent);

  // "New style" rules initialization:
  m_Regulament.CreateRulesFromJSONFiles();
  m_Regulament.EnableRulesAccordingToParameters(m_rulesoption);

  SetFittingMethod(m_fittingmethod);
  // Load the line catalog
  Log.LogDebug("About to load line catalog.");
  if (m_rigidity == "rules") {
    // load the regular catalog
    LoadCatalog(m_RestLineList);
    SetLSF();
  } else { //"tplratio" and "tplcorr"
    // load the tplratio catalog with only 1 element for all lines
    // LoadCatalogOneMultiline(restLineList);
    // load the tplratio catalog with 2 elements: 1 for the Em lines + 1 for the
    // Abs lines
    LoadCatalogTwoMultilinesAE(m_RestLineList);

    SetLSF();
    m_CatalogTplRatio = *(Context.GetTplRatioCatalog());
    initTplratioCatalogs(Context.GetParameterStore()->GetScoped<bool>(
        "linemodel.tplratio_ismfit"));
    SetTplratio_PriorHelper();
  }

  LogCatalogInfos();

  // TODO restore check the continuum flux axis for NaN
}
// hook
void CLineModelFitting::initTplratioCatalogs(Int32 opt_tplratio_ismFit) {
  // TODO: use the passed tplRatioCatalog
  // TODO: check if m_CatalogTplRatio changes between iterations

  m_CatalogTplRatio.Init(opt_tplratio_ismFit,
                         m_continuumManager->getIsmCorrectionFromTpl(),
                         m_NSigmaSupport);
  m_CatalogTplRatio.InitLineCorrespondingAmplitudes(m_Elements);
  SetMultilineNominalAmplitudesFast(0);
  // m_RestLineList = m_CatalogTplRatio.GetRestLinesList(0);
  // LoadCatalog(m_RestLineList);
  // LogCatalogInfos();

  // Resize tplratio buffers
  m_ChisquareTplratio.resize(m_CatalogTplRatio.GetCatalogsCount());
  m_ScaleMargCorrTplratio.resize(m_CatalogTplRatio.GetCatalogsCount());
  m_StrongELPresentTplratio.resize(m_CatalogTplRatio.GetCatalogsCount());
  m_StrongHalphaELPresentTplratio.resize(m_CatalogTplRatio.GetCatalogsCount());
  m_NLinesAboveSNRTplratio.resize(m_CatalogTplRatio.GetCatalogsCount());
  m_FittedAmpTplratio.resize(m_CatalogTplRatio.GetCatalogsCount());
  m_LyaAsymCoeffTplratio.resize(m_CatalogTplRatio.GetCatalogsCount());
  m_LyaWidthCoeffTplratio.resize(m_CatalogTplRatio.GetCatalogsCount());
  m_LyaDeltaCoeffTplratio.resize(m_CatalogTplRatio.GetCatalogsCount());
  m_LyaIgmIdxTplratio.resize(m_CatalogTplRatio.GetCatalogsCount());
  m_FittedErrorTplratio.resize(m_CatalogTplRatio.GetCatalogsCount());
  m_MtmTplratio.resize(m_CatalogTplRatio.GetCatalogsCount());
  m_DtmTplratio.resize(m_CatalogTplRatio.GetCatalogsCount());
  m_LinesLogPriorTplratio.resize(m_CatalogTplRatio.GetCatalogsCount());
  for (Int32 ktplratio = 0; ktplratio < m_CatalogTplRatio.GetCatalogsCount();
       ktplratio++) {
    m_FittedAmpTplratio[ktplratio].resize(m_Elements.size());
    m_FittedErrorTplratio[ktplratio].resize(m_Elements.size());
    m_MtmTplratio[ktplratio].resize(m_Elements.size());
    m_DtmTplratio[ktplratio].resize(m_Elements.size());
    m_LyaAsymCoeffTplratio[ktplratio].resize(m_Elements.size());
    m_LyaWidthCoeffTplratio[ktplratio].resize(m_Elements.size());
    m_LyaDeltaCoeffTplratio[ktplratio].resize(m_Elements.size());
    m_LyaIgmIdxTplratio[ktplratio].resize(m_Elements.size());
    m_LinesLogPriorTplratio[ktplratio].resize(m_Elements.size());
  }

  m_tplratioLeastSquareFast = false;
}

void CLineModelFitting::logParameters() {
  Log.LogInfo(Formatter() << "m_pass" << m_pass);
  Log.LogInfo(Formatter() << " m_enableAmplitudeOffsets"
                          << m_enableAmplitudeOffsets);
  Log.LogInfo(Formatter() << " m_LambdaOffsetMin" << m_LambdaOffsetMin);
  Log.LogInfo(Formatter() << " m_LambdaOffsetMax" << m_LambdaOffsetMax);
  Log.LogInfo(Formatter() << " m_LambdaOffsetStep" << m_LambdaOffsetStep);
  Log.LogInfo(Formatter() << " m_enableLambdaOffsetsFit"
                          << m_enableLambdaOffsetsFit);

  Log.LogInfo(Formatter() << " m_opt_lya_forcefit" << m_opt_lya_forcefit);
  Log.LogInfo(Formatter() << " m_opt_lya_forcedisablefit"
                          << m_opt_lya_forcedisablefit);
  Log.LogInfo(Formatter() << " m_opt_lya_fit_asym_min"
                          << m_opt_lya_fit_asym_min);
  Log.LogInfo(Formatter() << " m_opt_lya_fit_asym_max"
                          << m_opt_lya_fit_asym_max);
  Log.LogInfo(Formatter() << " m_opt_lya_fit_asym_step"
                          << m_opt_lya_fit_asym_step);
  Log.LogInfo(Formatter() << " m_opt_lya_fit_width_min"
                          << m_opt_lya_fit_width_min);
  Log.LogInfo(Formatter() << " m_opt_lya_fit_width_max"
                          << m_opt_lya_fit_width_max);
  Log.LogInfo(Formatter() << " m_opt_lya_fit_width_step"
                          << m_opt_lya_fit_width_step);
  Log.LogInfo(Formatter() << " m_opt_lya_fit_delta_min"
                          << m_opt_lya_fit_delta_min);
  Log.LogInfo(Formatter() << " m_opt_lya_fit_delta_max"
                          << m_opt_lya_fit_delta_max);
  Log.LogInfo(Formatter() << " m_opt_lya_fit_delta_step"
                          << m_opt_lya_fit_delta_step);

  Log.LogInfo(Formatter() << " m_opt_firstpass_forcedisableMultipleContinuumfit"
                          << m_opt_firstpass_forcedisableMultipleContinuumfit);
  Log.LogInfo(Formatter() << " m_opt_firstpass_forcedisableTplratioISMfit "
                          << m_opt_firstpass_forcedisableTplratioISMfit);
  Log.LogInfo(Formatter() << "m_opt_firstpass_fittingmethod "
                          << m_opt_firstpass_fittingmethod);
  Log.LogInfo(Formatter() << "m_opt_secondpass_fittingmethod"
                          << m_opt_secondpass_fittingmethod);

  Log.LogInfo(Formatter() << " m_opt_haprior" << m_opt_haprior);

  Log.LogInfo(Formatter() << "ContinuumComponent=" << m_ContinuumComponent);
  Log.LogInfo(Formatter() << "LineWidthType=" << m_LineWidthType);
  Log.LogInfo(Formatter() << "NSigmaSupport=" << m_NSigmaSupport);
  Log.LogInfo(Formatter() << "velocityEmission=" << m_velocityEmission);
  Log.LogInfo(Formatter() << "velocityAbsorption=" << m_velocityAbsorption);
  Log.LogInfo(Formatter() << "velocityEmissionInit=" << m_velocityEmissionInit);
  Log.LogInfo(Formatter() << "velocityAbsorptionInit="
                          << m_velocityAbsorptionInit);

  Log.LogInfo(Formatter() << "nominalWidthDefaultEmission="
                          << m_nominalWidthDefaultEmission);
  Log.LogInfo(Formatter() << "nominalWidthDefaultAbsorption="
                          << m_nominalWidthDefaultAbsorption);

  Log.LogInfo(Formatter() << "fittingmethod=" << m_fittingmethod);

  Log.LogInfo(Formatter() << "rulesoption=" << m_rulesoption);
  Log.LogInfo(Formatter() << "rigidity=" << m_rigidity);
  Log.LogInfo(Formatter() << "forcedisableTplratioISMfit="
                          << m_forcedisableTplratioISMfit);

  // Log.LogInfo(Formatter()<<"tplCatalog="<<m_tplCatalog);
  // Log.LogInfo(Formatter()<<"tplCategoryList="<<m_tplCategoryList);
  Log.LogInfo(Formatter() << "tplratioBestTplName=" << m_tplratioBestTplName);
  Log.LogInfo(Formatter() << "tplratioBestTplIsmCoeff="
                          << m_tplratioBestTplIsmCoeff);
  Log.LogInfo(Formatter() << "tplratioBestTplAmplitude="
                          << m_tplratioBestTplAmplitude);
  Log.LogInfo(Formatter() << "tplratioBestTplDtm=" << m_tplratioBestTplDtm);
  Log.LogInfo(Formatter() << "tplratioBestTplMtm=" << m_tplratioBestTplMtm);
  Log.LogInfo(Formatter()
              << "tplratioLeastSquareFast="
              << m_tplratioLeastSquareFast); // for rigidity=tplratio: switch to
                                             // use fast least square estimation

  Log.LogInfo(Formatter() << "secondpass_fitContinuum_dustfit="
                          << m_secondpass_fitContinuum_dustfit);
  Log.LogInfo(Formatter() << "secondpass_fitContinuum_igm="
                          << m_secondpass_fitContinuum_igm);

  //  Log.LogInfo(Formatter()<<"fitContinuum_tplFitPolyCoeffs="<<m_fitContinuum_tplFitPolyCoeffs);
  //  // only used with
  // m_fitContinuum_option==2 for now
  Log.LogInfo(Formatter() << "forcedisableMultipleContinuumfit="
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
    m_forceDisableLyaFitting = true;
    m_forcedisableTplratioISMfit = m_opt_firstpass_forcedisableTplratioISMfit;
    m_forcedisableMultipleContinuumfit =
        m_opt_firstpass_forcedisableMultipleContinuumfit;
    SetFittingMethod(m_opt_firstpass_fittingmethod);
    m_forceLyaFitting = false;
  }
  if (iPass == 2) {
    m_forceDisableLyaFitting = m_opt_lya_forcedisablefit;
    m_forcedisableTplratioISMfit = false;
    m_forcedisableMultipleContinuumfit = false;
    SetFittingMethod(m_opt_secondpass_fittingmethod);
    m_forceLyaFitting = m_opt_lya_forcefit;
    Log.LogDetail(
        "    model: set forceLyaFitting ASYMFIT for Tpl-ratio mode : %d",
        m_forceLyaFitting);
  }
  if (iPass == 3) {
    m_forceDisableLyaFitting = false;

    m_forcedisableMultipleContinuumfit = false;
    m_model->m_enableAmplitudeOffsets = true;
    m_enableAmplitudeOffsets = true;
    m_fitter->enableAmplitudeOffsets();
  }

  return true;
}
Int32 CLineModelFitting::GetPassNumber() const { return m_pass; }

void CLineModelFitting::SetForcedisableTplratioISMfit(bool opt) {
  m_forcedisableTplratioISMfit = opt;
}

/**
 * @brief CLineModelFitting::getFluxDirectIntegration
 * Integrates the flux (F-continuum) in a lbda range around the center observed
 * wavelength of the line. The wavelength range is defined by the instrument
 * resolution and a hardcoded nsigma factor
 * @param eIdx
 * @param subeIdx
 * @return
 */
void CLineModelFitting::getFluxDirectIntegration(const TInt32List &eIdx_list,
                                                 const TInt32List &subeIdx_list,
                                                 bool substract_abslinesmodel,
                                                 Float64 &fluxdi,
                                                 Float64 &snrdi) const {

  fluxdi = NAN;
  snrdi = NAN;

  Int32 nlines = eIdx_list.size();
  if (nlines != subeIdx_list.size())
    THROWG(INTERNAL_ERROR, " index sizes do not match");
  TInt32RangeList indexRangeList =
      getlambdaIndexesUnderLines(eIdx_list, subeIdx_list, N_SIGMA_SUPPORT_DI);

  if (!indexRangeList.size())
    THROWG(INTERNAL_ERROR, "empty indexRanges ");

  const CSpectrumSpectralAxis &spectralAxis = m_inputSpc->GetSpectralAxis();
  const auto &ContinuumFluxAxis = m_model->getContinuumFluxAxis();
  CSpectrumFluxAxis continuumFlux(ContinuumFluxAxis.GetSamplesCount());
  CSpectrumFluxAxis absLinesModelFlux;
  if (substract_abslinesmodel)
    absLinesModelFlux =
        getModel(CLine::nType_Absorption); // contains the continuum
                                           // and abs lines model;

  // polynome are shared between overlapping CElements
  //  find polynome belonging to this CElement
  // TODO: correct this for more than one CElement
  TPolynomCoeffs polynom_coeffs{0., 0., 0.};
  if (m_enableAmplitudeOffsets)
    polynom_coeffs = getPolynomCoeffs(eIdx_list[0]);

  // compute continuum
  for (const auto &r : indexRangeList)
    for (Int32 t = r.GetBegin(); t <= r.GetEnd(); t++) {
      continuumFlux[t] =
          substract_abslinesmodel ? absLinesModelFlux[t] : ContinuumFluxAxis[t];
      if (m_enableAmplitudeOffsets)
        continuumFlux[t] +=
            polynom_coeffs.x0 + polynom_coeffs.x1 * spectralAxis[t] +
            polynom_coeffs.x2 * spectralAxis[t] * spectralAxis[t];
    }
  // substarct continuum from spectrum flux
  const auto &SpcFluxAxis = m_model->getSpcFluxAxis();
  CSpectrumFluxAxis fluxMinusContinuum(SpcFluxAxis);
  for (const auto &r : indexRangeList)
    for (Int32 t = r.GetBegin(); t <= r.GetEnd(); t++)
      fluxMinusContinuum[t] -= continuumFlux[t];

  // estimate the integrated flux between obs. spectrum and continuum:
  // trapezoidal intg
  Float64 sumFlux = 0.0;
  Float64 sumErr = 0.0;
  integrateFluxes_usingTrapez(fluxMinusContinuum, indexRangeList, sumFlux,
                              sumErr);

  if (sumErr <= 0)
    return;

  fluxdi = sumFlux;
  snrdi = std::abs(fluxdi) / sqrt(sumErr);
  return;
}
/**
 * @brief
 *
 * @param eIdx_list
 * @param subeIdx_list
 * @param sigma_support
 * @return TInt32RangeList
 */
TInt32RangeList
CLineModelFitting::getlambdaIndexesUnderLines(const TInt32List &eIdx_list,
                                              const TInt32List &subeIdx_list,
                                              Float64 sigma_support) const {

  const CSpectrumSpectralAxis &spectralAxis = m_inputSpc->GetSpectralAxis();

  TInt32RangeList indexRangeList(eIdx_list.size());
  for (Int32 i = 0; i < eIdx_list.size(); i++) {
    Int32 eIdx = eIdx_list[i];
    Int32 subeIdx = subeIdx_list[i];

    Float64 mu = NAN;
    Float64 LineWidth = NAN;
    m_Elements[eIdx]->getObservedPositionAndLineWidth(subeIdx, m_Redshift, mu,
                                                      LineWidth);

    Float64 winsizeAngstrom = LineWidth * sigma_support;

    indexRangeList[i] = CLineModelElement::EstimateIndexRange(
        spectralAxis, mu, *(m_lambdaRange), winsizeAngstrom);
  }

  if (eIdx_list.size() == 1)
    return indexRangeList;
  TInt32RangeList nonOverlappingIndexRangeList =
      TInt32Range::joinIntersections(std::move(indexRangeList));

  return nonOverlappingIndexRangeList;
}

/**
 * @brief Construct a new clinemodelfitting::integratefluxes usingtrapez
 * object
 *
 * @param fluxMinusContinuum
 * @param indexRangeList
 * @param sumFlux
 * @param sumErr
 */
void CLineModelFitting::integrateFluxes_usingTrapez(
    const CSpectrumFluxAxis &fluxMinusContinuum,
    const TInt32RangeList &indexRangeList, Float64 &sumFlux,
    Float64 &sumErr) const {

  sumFlux = 0.0;
  sumErr = 0.0;

  const CSpectrumSpectralAxis &spectralAxis = m_inputSpc->GetSpectralAxis();
  if (spectralAxis.GetSamplesCount() < 2)
    THROWG(INTERNAL_ERROR, "Not enough samples in spectral axis");

  const auto &ErrorNoContinuum = m_inputSpc->GetFluxAxis().GetError();
  for (auto &r : indexRangeList) {
    for (Int32 t = r.GetBegin(); t < r.GetEnd(); t++) {
      Float64 trapweight = (spectralAxis[t + 1] - spectralAxis[t]) * 0.5;
      sumFlux +=
          trapweight * (fluxMinusContinuum[t + 1] + fluxMinusContinuum[t]);

      Float64 ea = ErrorNoContinuum[t] * ErrorNoContinuum[t];
      Float64 eb = ErrorNoContinuum[t + 1] * ErrorNoContinuum[t + 1];
      sumErr += trapweight * trapweight * (eb + ea);
    }
  }
  return;
}

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
  m_Redshift = redshift;

  if (reinterpolatedContinuum) {
    m_continuumManager->reinterpolateContinuum(redshift);
  }
}

const std::string &CLineModelFitting::getTplratio_bestTplName() const {
  return m_tplratioBestTplName;
}
Float64 CLineModelFitting::getTplratio_bestTplIsmCoeff() const {
  return m_tplratioBestTplIsmCoeff;
}

Float64 CLineModelFitting::getTplratio_bestAmplitude() const {
  return m_tplratioBestTplAmplitude;
}

Float64 CLineModelFitting::getTplratio_bestDtm() const {
  return m_tplratioBestTplDtm;
}

Float64 CLineModelFitting::getTplratio_bestMtm() const {
  return m_tplratioBestTplMtm;
}

Int32 CLineModelFitting::getTplratio_count() const {
  if (m_rigidity == "rules")
    return 0;

  return m_CatalogTplRatio.GetCatalogsCount();
}

const TFloat64List &CLineModelFitting::getTplratio_priors() {
  if (m_rigidity == "rules") {
    static TFloat64List dumb;
    return dumb;
  }
  return m_CatalogTplRatio.getCatalogsPriors();
}

const TFloat64List &CLineModelFitting::GetChisquareTplratio() const {
  return m_ChisquareTplratio;
}

/**
 * @brief CLineModelFitting::GetPriorLinesTplratio
 * WARNING: as stated in fit(), the prior is valid in this code structure only
 * for tplratio with only 1 element containing the EL component, hence using
 * idx=0 here
 * @return
 */
TFloat64List CLineModelFitting::GetPriorLinesTplratio() const {
  TFloat64List plinestplratio;
  Int32 eltIdx = 0;
  for (Int32 ktpl = 0; ktpl < m_LinesLogPriorTplratio.size(); ktpl++) {
    plinestplratio.push_back(m_LinesLogPriorTplratio[ktpl][eltIdx]);
  }
  return plinestplratio;
}

const TFloat64List &CLineModelFitting::GetScaleMargTplratio() const {
  return m_ScaleMargCorrTplratio;
}

const TBoolList &CLineModelFitting::GetStrongELPresentTplratio() const {
  return m_StrongELPresentTplratio;
}

const TBoolList &CLineModelFitting::getHaELPresentTplratio() const {
  return m_StrongHalphaELPresentTplratio;
}

const TInt32List &CLineModelFitting::GetNLinesAboveSNRTplratio() const {
  return m_NLinesAboveSNRTplratio;
}

bool CLineModelFitting::initModelAtZ(
    Float64 redshift, const CSpectrumSpectralAxis &spectralAxis) {
  m_Redshift = redshift;
  m_model->setRedshift(redshift);

  // prepare the elements support
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    m_Elements[iElts]->resetAsymfitParams();
    m_Elements[iElts]->prepareSupport(spectralAxis, redshift, *(m_lambdaRange));
  }

  return true;
}

bool CLineModelFitting::setTplratioModel(Int32 itplratio,
                                         bool enableSetVelocity) {
  SetMultilineNominalAmplitudesFast(itplratio);

  if (enableSetVelocity) {
    // Set the velocities from templates: todo auto switch when velfit is ON
    m_CatalogTplRatio.GetCatalogVelocities(itplratio, m_velocityEmission,
                                           m_velocityAbsorption);
  }

  Log.LogDebug("    model : setTplratioModel, loaded: %d = %s", itplratio,
               m_CatalogTplRatio.GetCatalogName(itplratio).c_str());
  return true;
}

void CLineModelFitting::SetTplratio_PriorHelper() {
  CAutoScope a = CAutoScope(Context.m_ScopeStack, "linemodel");
  CAutoScope b = CAutoScope(Context.m_ScopeStack, "tplratio");
  CAutoScope c = CAutoScope(Context.m_ScopeStack, "priors");
  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();

  m_tplratio_priorhelper = std::make_shared<CPriorHelper>();
  m_tplratio_priorhelper->Init(ps->GetScoped<std::string>("catalog_dirpath"),
                               0);
  m_tplratio_priorhelper->SetBetaA(ps->GetScoped<Float64>("betaA"));
  m_tplratio_priorhelper->SetBetaTE(ps->GetScoped<Float64>("betaTE"));
  m_tplratio_priorhelper->SetBetaZ(ps->GetScoped<Float64>("betaZ"));
}

bool CLineModelFitting::setTplratioAmplitude(const TFloat64List &ampsElts,
                                             const TFloat64List &errorsElts) {
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    m_Elements[iElts]->SetFittedAmplitude(ampsElts[iElts], errorsElts[iElts]);
  }
  return true;
}

bool CLineModelFitting::initDtd() {
  //  m_dTransposeDLambdaRange = TLambdaRange(*(m_lambdaRange));
  m_dTransposeDLambdaRange = *(m_lambdaRange);
  if (m_ContinuumComponent == "tplfit" ||
      m_ContinuumComponent == "tplfitauto") {
    m_dTransposeD = EstimateDTransposeD("raw");
  } else {
    m_dTransposeD = EstimateDTransposeD("nocontinuum");
  }
  m_likelihood_cstLog = EstimateLikelihoodCstLog();
  return true;
}

void CLineModelFitting::prepareAndLoadContinuum(Int32 k, Float64 redshift) {
  if (m_ContinuumComponent == "nocontinuum")
    return;

  m_model->PrepareContinuum();

  if (isContinuumComponentTplfitxx()) // the support has to be already computed
                                      // when LoadFitContinuum() is called
  {
    m_continuumManager->initObserveGridContinuumFlux(
        m_inputSpc->GetSampleCount());
    Int32 autoselect = m_ContinuumComponent == "tplfitauto";
    m_continuumManager->LoadFitContinuum(k, autoselect, redshift);
  }
}

void CLineModelFitting::computeSpectrumFluxWithoutContinuum() {
  m_model->initModelWithContinuum();
  m_model->substractContToFlux();
}

/**
 * @brief
 *
 * @param idx
 * @param _merit
 * @param _meritprior
 */
void CLineModelFitting::updateTplratioResults(
    Int32 idx, Float64 _merit, Float64 _meritprior,
    TFloat64List &bestTplratioMerit, TFloat64List &bestTplratioMeritPrior) {
  bestTplratioMerit[idx] = _merit;
  bestTplratioMeritPrior[idx] = _meritprior;
  m_ChisquareTplratio[idx] = _merit;
  m_ScaleMargCorrTplratio[idx] = getScaleMargCorrection();
  m_StrongELPresentTplratio[idx] = GetModelStrongEmissionLinePresent();
  // given that Ha is a strong emission line,
  if (m_opt_haprior > 0.) // check first that haprior is activated
    m_StrongHalphaELPresentTplratio[idx] = m_StrongELPresentTplratio[idx]
                                               ? GetModelHaStrongest()
                                               : false; // result per tplratio

  TStringList strongELSNRAboveCut; // = getLinesAboveSNR(3.5); //this
                                   // is costing a lot of processing
                                   // time, so deactivated for now.
  m_NLinesAboveSNRTplratio[idx] = strongELSNRAboveCut.size();

  // Saving the model A, errorA, and dtm, mtm, ... (for all tplratios,
  // needed ?) NB: this is only needed for the index=savedIdxFitted
  // ultimately
  for (Int32 iElt = 0; iElt < m_Elements.size(); iElt++) {
    bool savedAmp = false;
    bool allampzero = true;
    m_FittedAmpTplratio[idx][iElt] = NAN;
    m_FittedErrorTplratio[idx][iElt] = NAN;
    m_DtmTplratio[idx][iElt] = NAN;
    m_MtmTplratio[idx][iElt] = NAN;
    m_LyaAsymCoeffTplratio[idx][iElt] = NAN;
    m_LyaWidthCoeffTplratio[idx][iElt] = NAN;
    m_LyaDeltaCoeffTplratio[idx][iElt] = NAN;
    m_LyaIgmIdxTplratio[idx][iElt] = undefIdx;
    m_LinesLogPriorTplratio[idx][iElt] = _meritprior;
    Int32 nLines = m_Elements[iElt]->GetSize();

    for (Int32 j = 0; j < nLines; j++) {
      if (savedAmp) {
        break;
      }
      Float64 amp = m_Elements[iElt]->GetFittedAmplitude(j);
      if (amp > 0) {
        allampzero = false;
      }
      if (amp > 0 && !m_Elements[iElt]->IsOutsideLambdaRange(j)) {

        Float64 amp_error = m_Elements[iElt]->GetFittedAmplitudeErrorSigma(j);
        Float64 nominal_amp = m_Elements[iElt]->GetNominalAmplitude(j);
        m_FittedAmpTplratio[idx][iElt] = amp / nominal_amp;
        Log.LogDebug("    model : fit tplratio mode, tplratio_fittedamp: %e",
                     m_FittedAmpTplratio[idx][iElt]);

        m_FittedErrorTplratio[idx][iElt] = amp_error / nominal_amp;
        m_DtmTplratio[idx][iElt] = m_Elements[iElt]->GetSumCross();
        m_MtmTplratio[idx][iElt] = m_Elements[iElt]->GetSumGauss();

        TAsymParams params = m_Elements[iElt]->GetAsymfitParams(0);
        m_LyaAsymCoeffTplratio[idx][iElt] = params.alpha;
        m_LyaWidthCoeffTplratio[idx][iElt] = params.sigma;
        m_LyaDeltaCoeffTplratio[idx][iElt] = params.delta;

        TSymIgmParams params_igm = m_Elements[iElt]->GetSymIgmParams(0);
        m_LyaIgmIdxTplratio[idx][iElt] = params_igm.m_igmidx;

        savedAmp = true;
        break;
      }
    }
    // TODO: this case should be treated more
    // carefully, save dtm, mtm, and more...
    if (allampzero && !savedAmp)
      m_FittedAmpTplratio[idx][iElt] = 0.0;
  }
  return;
}

/**
 * @brief :copy the values for ebmv=ebmv_fixed (=0)
 *
 * @param idx
 */
void CLineModelFitting::duplicateTplratioResult(
    Int32 idx, TFloat64List &bestTplratioMerit,
    TFloat64List &bestTplratioMeritPrior) {
  bestTplratioMerit[idx] = bestTplratioMerit[idx - 1];
  bestTplratioMeritPrior[idx] = bestTplratioMeritPrior[idx - 1];
  m_ChisquareTplratio[idx] = m_ChisquareTplratio[idx - 1];
  m_ScaleMargCorrTplratio[idx] = m_ScaleMargCorrTplratio[idx - 1];
  m_StrongELPresentTplratio[idx] = m_StrongELPresentTplratio[idx - 1];
  m_StrongHalphaELPresentTplratio[idx] =
      m_StrongHalphaELPresentTplratio[idx - 1];
  m_NLinesAboveSNRTplratio[idx] = m_NLinesAboveSNRTplratio[idx - 1];

  for (Int32 iElt = 0; iElt < m_Elements.size(); iElt++) {
    m_FittedAmpTplratio[idx][iElt] = m_FittedAmpTplratio[idx - 1][iElt];
    m_FittedErrorTplratio[idx][iElt] = m_FittedErrorTplratio[idx - 1][iElt];
    m_DtmTplratio[idx][iElt] = m_DtmTplratio[idx - 1][iElt];
    m_MtmTplratio[idx][iElt] = m_MtmTplratio[idx - 1][iElt];
    m_LyaAsymCoeffTplratio[idx][iElt] = m_LyaAsymCoeffTplratio[idx - 1][iElt];
    m_LyaWidthCoeffTplratio[idx][iElt] = m_LyaWidthCoeffTplratio[idx - 1][iElt];
    m_LyaDeltaCoeffTplratio[idx][iElt] = m_LyaDeltaCoeffTplratio[idx - 1][iElt];
    m_LyaIgmIdxTplratio[idx][iElt] = m_LyaIgmIdxTplratio[idx - 1][iElt];
    m_LinesLogPriorTplratio[idx][iElt] = m_LinesLogPriorTplratio[idx - 1][iElt];
  }
  return;
}

Float64 CLineModelFitting::computelogLinePriorMerit(
    Int32 itratio,
    const std::vector<CPriorHelper::SPriorTZE> &logPriorDataTplRatio) {

  if (!logPriorDataTplRatio.size())
    return 0.;

  // lines prior
  Float64 _meritprior = -2. * logPriorDataTplRatio[itratio].betaTE *
                        logPriorDataTplRatio[itratio].logprior_precompTE;
  _meritprior += -2. * logPriorDataTplRatio[itratio].betaA *
                 logPriorDataTplRatio[itratio].logprior_precompA;
  _meritprior += -2. * logPriorDataTplRatio[itratio].betaZ *
                 logPriorDataTplRatio[itratio].logprior_precompZ;
  if (logPriorDataTplRatio[itratio].A_sigma <= 0.0)
    return _meritprior;

  Float64 ampl = 0.0;
  for (Int32 iElt = 0; iElt < m_Elements.size(); iElt++) {
    bool foundAmp = false;
    Int32 nLines = m_Elements[iElt]->GetSize();
    for (Int32 j = 0; j < nLines; j++) {
      Float64 amp = m_Elements[iElt]->GetFittedAmplitude(j);
      if (amp > 0 && !m_Elements[iElt]->IsOutsideLambdaRange(j)) {
        Float64 nominal_amp = m_Elements[iElt]->GetNominalAmplitude(j);
        ampl = amp / nominal_amp;
        foundAmp = true;
        break;
      }
    }
    /*if (foundAmp)
      break;
      //Didier: probably this is missing here??? I suppose we are
      // looking for the first non-null and valid amplitude?
     */
  }
  _meritprior += logPriorDataTplRatio[itratio].betaA *
                 (ampl - logPriorDataTplRatio[itratio].A_mean) *
                 (ampl - logPriorDataTplRatio[itratio].A_mean) /
                 (logPriorDataTplRatio[itratio].A_sigma *
                  logPriorDataTplRatio[itratio].A_sigma);

  return _meritprior;
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

  Int32 ntplratio = 1; // multiple fitting steps for rigidity=tplratio/tplratio
  std::vector<CPriorHelper::SPriorTZE> logPriorDataTplRatio;
  if (m_rigidity == "tplratio") {
    ntplratio = m_CatalogTplRatio.GetCatalogsCount();

    if (!m_tplratio_priorhelper->mInitFailed) {
      // prior initilization for tplratio EL only
      if (m_Elements.size() > 1)
        THROWG(INTERNAL_ERROR, "model: Unable to use tplratio line priors "
                               "with nElts>1 for now");
      // NB: this could be done if the EL element idx in searched (see later
      // in the itratio loop, UV Abs lines would be not affected by priors
      // then)

      for (Int32 itratio = 0; itratio < ntplratio; itratio++) {
        // prepare the lines prior data
        Int32 ebvfilter = m_CatalogTplRatio.GetIsmIndex(itratio);
        CPriorHelper::SPriorTZE logPriorData;
        std::string tplrationame = m_CatalogTplRatio.GetCatalogName(itratio);
        bool retGetPrior = m_tplratio_priorhelper->GetTZEPriorData(
            tplrationame, ebvfilter, redshift, logPriorData);
        if (retGetPrior == false)
          THROWG(INTERNAL_ERROR,
                 "model: Failed to get prior for chi2 solvecontinuum.");
        else
          logPriorDataTplRatio.push_back(logPriorData);
      }
    }
  }
  Int32 savedIdxFitted = -1; // for rigidity=tplratio

  Int32 nContinuum = 1;
  Int32 savedIdxContinuumFitted = -1; // for continuum tplfit
  if (isContinuumComponentTplfitxx() && !m_forcedisableMultipleContinuumfit)
    nContinuum = m_continuumManager->m_opt_fitcontinuum_maxCount;
  // 'on the fly' initialization
  Float64 bestMerit = INFINITY;
  Float64 bestMeritPrior = 0.0;
  TFloat64List bestTplratioMerit(ntplratio, INFINITY);
  TFloat64List bestTplratioMeritPrior(ntplratio, 0.0);

  for (Int32 k = 0; k < nContinuum; k++) {

    Float64 _merit = INFINITY;
    Float64 _meritprior = 0.; // only relevant for "tplratio"
    prepareAndLoadContinuum(k, redshift);
    if (m_ContinuumComponent != "nocontinuum")
      computeSpectrumFluxWithoutContinuum();

    if (m_enableAmplitudeOffsets)
      m_Elements.prepareAmplitudeOffset();

    for (Int32 itratio = 0; itratio < ntplratio; itratio++) {
      if (m_rigidity == "tplratio") {
        if (m_forcedisableTplratioISMfit && itratio > 0 &&
            m_CatalogTplRatio.GetIsmIndex(itratio) > 0) {
          duplicateTplratioResult(itratio, bestTplratioMerit,
                                  bestTplratioMeritPrior);
          continue;
        }
        setTplratioModel(itratio, false);
        // prepare the Lya width and asym coefficients if the asymfit profile
        // option is met INFO: tpl-shape are often ASYMFIXED in the tplratio
        // catalog files, for the lyaE profile, as of 2016-01-11 INFO:
        // tplratio can override the lyafitting, see m_opt_lya_forcefit
        setLyaProfile(m_Redshift,
                      m_CatalogTplRatio.GetCatalog(itratio).GetList(), true);
      } else {
        // prepare the Lya width and asym coefficients if the asymfit profile
        // option is met
        setLyaProfile(m_Redshift, m_RestLineList);
      }
      // generate random amplitudes
      m_fitter->fit(redshift);

      std::string bestTplratioName = "undefined";
      if (m_rigidity == "rules") {
        // Apply rules,
        applyRules(enableLogging);
        m_model->refreshModel();
        _merit = getLeastSquareMerit();
      }

      // correct lines amplitude with tplratioPrior (tpl-corr): Warning: Rules
      // must all be deactivated
      if (m_rigidity == "tplcorr") {
        m_model->refreshModel();
        // create spectrum model
        modelSolution = GetModelSolution(); // computed only to get
                                            // amplitude/err!
        continuumModelSolution =
            m_continuumManager->GetContinuumModelSolution(); // not used!
        TFloat64List correctedAmplitudes(modelSolution.Amplitudes.size());
        m_CatalogTplRatio.GetBestFit(m_RestLineList, modelSolution.Amplitudes,
                                     modelSolution.AmplitudesUncertainties,
                                     correctedAmplitudes, bestTplratioName);
        for (Int32 iRestLine = 0; iRestLine < m_RestLineList.size();
             iRestLine++) {
          Int32 subeIdx = undefIdx;
          Int32 eIdx = m_Elements.findElementIndex(iRestLine, subeIdx);
          if (eIdx == undefIdx || subeIdx == undefIdx)
            continue;

          Float64 er =
              modelSolution
                  .AmplitudesUncertainties[subeIdx]; // not modifying the
                                                     // fitting error for now
          Float64 nominalAmp = m_Elements[eIdx]->GetNominalAmplitude(subeIdx);
          m_Elements[eIdx]->SetFittedAmplitude(
              correctedAmplitudes[iRestLine] / nominalAmp, er);
        }
        m_model->refreshModel();
        _merit = getLeastSquareMerit();
      }

      if (m_rigidity == "tplratio") {
        if (!enableLogging && m_tplratioLeastSquareFast)
          _merit = getLeastSquareMeritFast();
        else {
          m_model->refreshModel();
          _merit = getLeastSquareMerit();
        }

        _meritprior = computelogLinePriorMerit(itratio, logPriorDataTplRatio);

        if (_merit + _meritprior <
            bestTplratioMerit[itratio] + bestTplratioMeritPrior[itratio]) {
          // update result variables
          updateTplratioResults(itratio, _merit, _meritprior, bestTplratioMerit,
                                bestTplratioMeritPrior);
        }
      }

      if (bestMerit + bestMeritPrior > _merit + _meritprior) {
        bestMerit = _merit;
        bestMeritPrior = _meritprior;
        savedIdxContinuumFitted = k;
        Int32 modelSolutionLevel =
            m_rigidity == "rules" ? Int32(enableLogging) : 0;
        modelSolution = GetModelSolution(modelSolutionLevel);
        continuumModelSolution =
            m_continuumManager->GetContinuumModelSolution();

        if (m_rigidity == "tplcorr")
          m_tplratioBestTplName = bestTplratioName;

        if (m_rigidity == "tplratio") {
          savedIdxFitted = itratio;
          m_tplratioBestTplName =
              m_CatalogTplRatio.GetCatalogName(savedIdxFitted);
          m_tplratioBestTplIsmCoeff =
              m_CatalogTplRatio.CalzettiInitFailed()
                  ? NAN
                  : m_CatalogTplRatio.GetIsmCoeff(savedIdxFitted);
          m_tplratioBestTplAmplitude =
              m_FittedAmpTplratio[savedIdxFitted][0]; // Should be only 1 elt
                                                      // in tpl ratio mode...
          m_tplratioBestTplDtm =
              m_DtmTplratio[savedIdxFitted]
                           [0]; // Should be only 1 elt in tpl ratio mode...
          m_tplratioBestTplMtm =
              m_MtmTplratio[savedIdxFitted]
                           [0]; // Should be only 1 elt in tpl ratio mode...
        }
      }
      if (m_ContinuumComponent == "nocontinuum")
        m_model->reinitModel();
    }
  }

  if (!enableLogging)
    return bestMerit;

  if (isContinuumComponentTplfitxx()) {
    if (m_fittingmethod != "svdlc" && nContinuum > 1) {
      Int32 autoselect = m_ContinuumComponent == "tplfitauto";
      // TODO savedIdxContinuumFitted=-1 if rigidity!=tplratio
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

  if (m_rigidity == "tplratio") {
    bool retSetMultiAmplFast =
        SetMultilineNominalAmplitudesFast(savedIdxFitted);
    if (!retSetMultiAmplFast) {
      Log.LogError("Linemodel: tplratio, Unable to set Multiline "
                   "NominalAmplitudes from Tplratio !");
    }

    // Set the velocities from templates: todo auto switch when velfit is ON
    // m_CatalogTplRatio.GetCatalogVelocities(savedIdxFitted,
    // m_velocityEmission, m_velocityAbsorption);
    for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
      Log.LogDetail("    model - Linemodel: tplratio = %d (%s, with "
                    "ebmv=%.3f), and A=%e",
                    savedIdxFitted, m_tplratioBestTplName.c_str(),
                    m_tplratioBestTplIsmCoeff,
                    m_FittedAmpTplratio[savedIdxFitted][iElts]);
      m_Elements[iElts]->SetFittedAmplitude(
          m_FittedAmpTplratio[savedIdxFitted][iElts],
          m_FittedErrorTplratio[savedIdxFitted][iElts]);
      m_Elements[iElts]->SetSumCross(m_DtmTplratio[savedIdxFitted][iElts]);
      m_Elements[iElts]->SetSumGauss(m_MtmTplratio[savedIdxFitted][iElts]);
    }

    // Lya
    for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++)
      m_Elements[iElts]->SetAsymfitParams(
          {m_LyaWidthCoeffTplratio[savedIdxFitted][iElts],
           m_LyaAsymCoeffTplratio[savedIdxFitted][iElts],
           m_LyaDeltaCoeffTplratio[savedIdxFitted][iElts]});

    for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++)
      m_Elements[iElts]->SetSymIgmParams(
          TSymIgmParams(m_LyaIgmIdxTplratio[savedIdxFitted][iElts], redshift));

    m_model->refreshModel();

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
  if (m_fittingmethod == "hybrid")
    m_fitter = std::make_shared<CHybridFitter>(
        CHybridFitter(m_Elements, m_inputSpc, m_lambdaRange, m_model));
  else if (m_fittingmethod == "svd")
    m_fitter = std::make_shared<CSvdFitter>(
        CSvdFitter(m_Elements, m_inputSpc, m_lambdaRange, m_model));
  else if (m_fittingmethod == "svdlc")
    m_fitter = std::make_shared<CSvdlcFitter>(CSvdlcFitter(
        m_Elements, m_inputSpc, m_lambdaRange, m_model, m_continuumManager));
  else if (m_fittingmethod == "svdlcp2")
    m_fitter = std::make_shared<CSvdlcFitter>(CSvdlcFitter(
        m_Elements, m_inputSpc, m_lambdaRange, m_model, m_continuumManager, 2));

  else if (m_fittingmethod == "ones")
    m_fitter = std::make_shared<COnesFitter>(
        COnesFitter(m_Elements, m_inputSpc, m_lambdaRange, m_model));
  else if (m_fittingmethod == "random")
    m_fitter = std::make_shared<CRandomFitter>(
        CRandomFitter(m_Elements, m_inputSpc, m_lambdaRange, m_model));
  else if (m_fittingmethod == "individual")
    m_fitter = std::make_shared<CIndividualFitter>(
        CIndividualFitter(m_Elements, m_inputSpc, m_lambdaRange, m_model));
  else
    THROWG(INTERNAL_ERROR, Formatter()
                               << "Unknown fitting method " << fitMethod);
}

void CLineModelFitting::SetLeastSquareFastEstimationEnabled(Int32 enabled) {
  m_tplratioLeastSquareFast = enabled;
}

void CLineModelFitting::SetAbsLinesLimit(Float64 limit) {
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    m_Elements[iElts]->SetAbsLinesLimit(limit);
  }
}

CSpectrumFluxAxis CLineModelFitting::getModel(Int32 lineTypeFilter) const {
  const CSpectrumSpectralAxis &spectralAxis = m_inputSpc->GetSpectralAxis();
  CSpectrumFluxAxis modelfluxAxis(spectralAxis.GetSamplesCount());

  Int32 nElements = m_Elements.size();
  for (Int32 iElts = 0; iElts < nElements; iElts++) {
    m_Elements[iElts]->initSpectrumModel(modelfluxAxis,
                                         m_model->getContinuumFluxAxis());
  }

  for (Int32 iElts = 0; iElts < nElements; iElts++) {
    Int32 lineType = m_Elements[iElts]->m_Lines[0].GetType();
    if (lineTypeFilter == -1 || lineTypeFilter == lineType) {
      m_Elements[iElts]->addToSpectrumModel(spectralAxis, modelfluxAxis,
                                            m_model->getContinuumFluxAxis(),
                                            m_Redshift);
    }
  }

  return modelfluxAxis;
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

/**
 * \brief Returns a sorted set of line indices present in the supports of the
 *argument. Create a vector named indexes. If the argument ind is an index to
 *m_Elements that IsOutSideLambdaRange, return indexes. For each entry in
 *m_Elements: If the entry has a different linetype than the line
 *corresponding to ind, go to the next entry. If the entry
 *IsOutsideLambdaRange, go to the enxt entry. For each subentry in the support
 *of entry: For each subsubentry in the support of ind: If the overlap in the
 *spectralAxis is smaller than -1 * overlapThres * winsize, add entry to
 *indexes. Sort indexes, remove duplicates from indexes, and return indexes.
 **/
TInt32List
CLineModelFitting::getOverlappingElementsBySupport(Int32 ind,
                                                   Float64 overlapThres) const {
  const CSpectrumSpectralAxis &spectralAxis = m_inputSpc->GetSpectralAxis();
  TInt32List indexes;

  if (m_Elements[ind]->IsOutsideLambdaRange()) {
    indexes.push_back(ind);
    return indexes;
  }
  TInt32RangeList refsupport = m_Elements[ind]->getSupport();
  const CLine &line = m_Elements[ind]->m_Lines[0];
  Int32 linetype = line.GetType();
  Float64 mu = line.GetPosition() * (1 + m_Redshift);
  Float64 c =
      m_Elements[ind]->GetLineWidth(mu, m_Redshift, line.GetIsEmission());
  Float64 winsize = line.GetProfile().GetNSigmaSupport() * c;
  Float64 overlapThresholdMin = winsize * overlapThres;
  // overlapThresholdMin = 0.0;

  Int32 x1 = 0;
  Int32 y1 = 0;
  Int32 x2 = 0;
  Int32 y2 = 0;
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    if (m_Elements[iElts]->m_Lines[0].GetType() != linetype) {
      continue;
    }

    if (m_Elements[iElts]->IsOutsideLambdaRange()) {
      continue;
    }
    TInt32RangeList s = m_Elements[iElts]->getSupport();
    for (Int32 iS = 0; iS < s.size(); iS++) {
      for (Int32 iRefS = 0; iRefS < refsupport.size(); iRefS++) {
        x1 = refsupport[iRefS].GetBegin();
        x2 = refsupport[iRefS].GetEnd();
        y1 = s[iS].GetBegin();
        y2 = s[iS].GetEnd();

        // Log.LogInfo( "hybrid fit: iRefS=%d - support=%d,%d", iRefS, x1,
        // x2); Log.LogInfo( "hybrid fit: iS=%d - support=%d,%d", iS, y1, y2);

        //                if( std::max(x1,y1) < std::min(x2,y2) ){
        //                    indexes.push_back(iElts);
        //                    break;
        //                }

        Float64 max = spectralAxis[std::max(x1, y1)];
        Float64 min = spectralAxis[std::min(x2, y2)];
        if (max - min < -overlapThresholdMin) {
          indexes.push_back(iElts);
          break;
        }
      }
    }
  }

  std::sort(indexes.begin(), indexes.end());
  indexes.erase(std::unique(indexes.begin(), indexes.end()), indexes.end());

  return indexes;
}

/**
 * @brief CLineModelFitting::setLyaProfile
 * If a Lya line is present with SYMIGM profile, fit igmIdx
 * If a Lya line is present with ASYMFIT profile, fit the width and asymmetry
 * parameters If a Lya line is present with ASYMFIXED profile, set the width
 * and asymmetry parameters according to profile parameters given in the
 * string
 * @param redshift, catalog (could be tplratio or linecatalog)
 * @return 1 if successfully fitted, 0 if error, 2 if Lya not present, 3 if Lya
 * not configured to be fitted in the catalog
 * note: this applies to all asym or symigm profiles (for symigm: all lines
 * below Lya)
 */
void CLineModelFitting::setLyaProfile(Float64 redshift,
                                      const CLineCatalog::TLineVector &catalog,
                                      bool tplratio) {

  auto idxLineIGM_ = m_Elements.getIgmLinesIndices();
  auto const idxEltIGM = std::move(idxLineIGM_.front());
  std::vector<TInt32List> idxLineIGM(
      std::make_move_iterator(idxLineIGM_.begin() + 1),
      std::make_move_iterator(idxLineIGM_.end()));

  if (idxEltIGM.empty())
    return;

  // assuming only one asymfit/fixed profile
  Int32 idxLyaE = idxEltIGM.front();
  Int32 idxLineLyaE = idxLineIGM.front().front();

  if (!m_Elements[idxLyaE]->IsOutsideLambdaRange(idxLineLyaE)) {
    const auto &profile =
        m_Elements[idxLyaE]->m_Lines[idxLineLyaE].GetProfile();
    if (profile.isAsym())
      setAsymProfile(idxLyaE, idxLineLyaE, redshift, catalog, tplratio);
  }

  for (Int32 i = 0; i < idxEltIGM.size(); ++i) {
    const auto &Elt = m_Elements[idxEltIGM[i]];
    if (!Elt->IsOutsideLambdaRange()) {
      TInt32List &idxLine = idxLineIGM[i];
      auto end =
          std::remove_if(idxLine.begin(), idxLine.end(), [Elt](Int32 idx) {
            return !Elt->m_Lines[idx].GetProfile().isSymIgm();
          });
      idxLine.erase(end, idxLine.end());
      if (!idxLine.empty())
        setSymIgmProfile(idxEltIGM[i], idxLine, redshift);
    }
  }
}

void CLineModelFitting::setAsymProfile(Int32 idxLyaE, Int32 idxLineLyaE,
                                       Float64 redshift,
                                       const CLineCatalog::TLineVector &catalog,
                                       bool tplratio) {
  Int32 lineIndex =
      getLineIndexInCatalog(idxLyaE, idxLineLyaE, catalog, tplratio);
  if (lineIndex == undefIdx)
    return;

  // finding or setting the correct profile
  CLineProfile_ptr profile;
  if (m_forceDisableLyaFitting && catalog[lineIndex].GetProfile().isAsymFit())
    // convert asymfit to asymfixed profile
    profile = dynamic_cast<const CLineProfileASYMFIT *>(
                  &catalog[lineIndex].GetProfile())
                  ->cloneToASYM();
  else if (m_forceLyaFitting && catalog[lineIndex].GetProfile().isAsymFixed())
    // convert asymfixed to asymfit
    profile =
        dynamic_cast<const CLineProfileASYM *>(&catalog[lineIndex].GetProfile())
            ->cloneToASYMFIT();
  else
    profile = catalog[lineIndex].GetProfile().Clone();

  bool doasymfit = profile->isAsymFit();

  m_Elements[idxLyaE]->m_Lines[idxLineLyaE].SetProfile(std::move(profile));

  if (!doasymfit)
    return;

  // find the best width and asym coeff. parameters
  TAsymParams bestfitParams = fitAsymParameters(redshift, idxLyaE, idxLineLyaE);

  // set the associated Lya members in the element definition
  m_Elements[idxLyaE]->SetAsymfitParams(bestfitParams);
}

void CLineModelFitting::setSymIgmProfile(Int32 iElts,
                                         const TInt32List &idxLineIGM,
                                         Float64 redshift) {

  bool fixedIGM = isContinuumComponentTplfitxx();

  // set to false when continuum is fitted to null
  fixedIGM &= m_continuumManager->isContFittedToNull();

  Int32 bestigmidx = fixedIGM
                         ? m_continuumManager->getFittedMeiksinIndex()
                         : fitAsymIGMCorrection(redshift, iElts, idxLineIGM);

  m_Elements[iElts]->SetSymIgmParams(TSymIgmParams(bestigmidx, redshift));
}

Int32 CLineModelFitting::getLineIndexInCatalog(
    Int32 iElts, Int32 idxLine, const CLineCatalog::TLineVector &catalog,
    bool tplratio) const {
  Int32 lineIndex = undefIdx;
  if (tplratio) {
    // get index of line inside tplratio catalog
    const std::string &strID = m_Elements[iElts]->m_Lines[idxLine].GetStrID();
    lineIndex = std::find_if(catalog.begin(), catalog.end(),
                             [strID](const CLine &line) {
                               return line.GetStrID() == strID;
                             }) -
                catalog.begin();
    if (lineIndex >= catalog.size())
      lineIndex = undefIdx;
  } else {
    lineIndex = m_Elements[iElts]->m_LineCatalogIndexes[idxLine];
    if (lineIndex < 0 || lineIndex >= catalog.size())
      THROWG(INTERNAL_ERROR, "Lya idx out-of-bound");
  }

  return lineIndex;
}

TAsymParams CLineModelFitting::fitAsymParameters(Float64 redshift,
                                                 Int32 idxLyaE,
                                                 const Int32 &idxLineLyaE) {

  const CSpectrumSpectralAxis &spectralAxis = m_inputSpc->GetSpectralAxis();
  // 3. find the best width and asym coeff. parameters
  Float64 widthCoeffStep = m_opt_lya_fit_width_step;
  Float64 widthCoeffMin = m_opt_lya_fit_width_min;
  Float64 widthCoeffMax = m_opt_lya_fit_width_max;
  Int32 nWidthSteps =
      int((widthCoeffMax - widthCoeffMin) / widthCoeffStep + 1.5);
  Float64 asymCoeffStep = m_opt_lya_fit_asym_step;
  Float64 asymCoeffMin = m_opt_lya_fit_asym_min;
  Float64 asymCoeffMax = m_opt_lya_fit_asym_max;
  Int32 nAsymSteps = int((asymCoeffMax - asymCoeffMin) / asymCoeffStep + 1.5);
  Float64 deltaStep = m_opt_lya_fit_delta_step;
  Float64 deltaMin = m_opt_lya_fit_delta_min;
  Float64 deltaMax = m_opt_lya_fit_delta_max;
  Int32 nDeltaSteps = int((deltaMax - deltaMin) / deltaStep + 1.5);

  TAsymParams bestparams = {widthCoeffMin, asymCoeffMin, deltaMin};
  Float64 meritMin = DBL_MAX;

  TInt32List filterEltsIdxLya(1, idxLyaE);

  for (Int32 iDelta = 0; iDelta < nDeltaSteps; iDelta++) {
    Float64 delta = deltaMin + deltaStep * iDelta;
    for (Int32 iWidth = 0; iWidth < nWidthSteps; iWidth++) {
      Float64 asymWidthCoeff = widthCoeffMin + widthCoeffStep * iWidth;
      for (Int32 iAsym = 0; iAsym < nAsymSteps; iAsym++) {
        Float64 asymAlphaCoeff = asymCoeffMin + asymCoeffStep * iAsym;
        m_Elements[idxLyaE]->SetAsymfitParams(
            {asymWidthCoeff, asymAlphaCoeff, delta});

        // idxLineLyaE = -1;
        m_Elements[idxLyaE]->fitAmplitude(
            spectralAxis, m_model->getSpcFluxAxisNoContinuum(),
            m_model->getContinuumFluxAxis(), redshift, idxLineLyaE);

        Float64 m = m_dTransposeD;
        if (1) {

          m_model->refreshModelUnderElements(filterEltsIdxLya, idxLineLyaE);
          m = m_model->getModelErrorUnderElement(idxLyaE);
        } else {
          m = getLeastSquareMeritFast(idxLyaE);
        }
        if (m < meritMin) {
          meritMin = m;
          bestparams = m_Elements[idxLyaE]->GetAsymfitParams(0);
        }

        Log.LogDebug("Fitting Lya Profile: width=%f, asym=%f, delta=%f",
                     asymWidthCoeff, asymAlphaCoeff, delta);
        Log.LogDebug("Fitting Lya Profile: merit=%e", m);
        Log.LogDebug("Fitting Lya Profile: idxLyaE=%d, idxLineLyaE=%d", idxLyaE,
                     idxLineLyaE);
      }
    }
  }
  Log.LogDebug("Lya Profile found: width=%f, asym=%f, delta=%f",
               bestparams.sigma, bestparams.alpha, bestparams.delta);
  return bestparams;
}

Int32 CLineModelFitting::fitAsymIGMCorrection(Float64 redshift, Int32 iElts,
                                              const TInt32List &idxLine) {

  const CSpectrumSpectralAxis &spectralAxis = m_inputSpc->GetSpectralAxis();
  if (spectralAxis[0] / (1 + redshift) > RESTLAMBDA_LYA)
    return -1;

  Float64 meritMin = DBL_MAX;
  Int32 bestIgmIdx = -1;

  Int32 igmCount =
      m_Elements[iElts]->getLineProfile(idxLine.front()).getIGMIdxCount();
  for (Int32 igmIdx = 0; igmIdx < igmCount; igmIdx++) {
    m_Elements[iElts]->SetSymIgmParams(TSymIgmParams(igmIdx, redshift));
    m_Elements[iElts]->fitAmplitude(spectralAxis,
                                    m_model->getSpcFluxAxisNoContinuum(),
                                    m_model->getContinuumFluxAxis(), redshift);

    m_model->refreshModelUnderElements(TInt32List(1, iElts));
    Float64 m = m_model->getModelErrorUnderElement(iElts);

    if (m < meritMin) {
      meritMin = m;
      bestIgmIdx = igmIdx;
    }
  }
  return bestIgmIdx;
}

/**
 * \brief Accumulates the squared differences between model and spectrum in
 *the argument lambdaRange and returns the sum.
 **/
Float64 CLineModelFitting::getLeastSquareMerit() const {
  const CSpectrumSpectralAxis &spcSpectralAxis = m_inputSpc->GetSpectralAxis();
  const CSpectrumFluxAxis &spcFluxAxis = m_model->getSpcFluxAxis();
  const CSpectrumFluxAxis &modelFluxAxis =
      m_model->GetModelSpectrum().GetFluxAxis();

  // Int32 numDevs = 0;
  Float64 fit = 0.0;
  const Float64 *Ymodel = modelFluxAxis.GetSamples();
  const Float64 *Yspc = spcFluxAxis.GetSamples();
  const auto &ErrorNoContinuum = m_inputSpc->GetFluxAxis().GetError();
  Float64 diff = 0.0;

  Int32 imin = spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange->GetBegin());
  Int32 imax = spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange->GetEnd());
  for (Int32 j = imin; j < imax; j++) {
    // numDevs++;
    diff = (Yspc[j] - Ymodel[j]);
    fit += (diff * diff) / (ErrorNoContinuum[j] * ErrorNoContinuum[j]);
    //        if ( 1E6 * diff < ErrorNoContinuum[j] )
    //        {
    //            Log.LogDebug( "Warning: noise is at least 6 orders greater
    //            than the residue!" ); Log.LogDebug(
    //            "CLineModelFitting::getLeastSquareMerit diff = %f", diff );
    //            Log.LogDebug( "CLineModelFitting::getLeastSquareMerit
    //            ErrorNoContinuum[%d] = %f", j, ErrorNoContinuum[j] );
    //        }
  }

  if (isContinuumComponentTplfitxx()) {
    fit += m_continuumManager->getFitSum(); // unconditionnal sum (if photometry
                                            // disabled, will sum 0.0)
  }

  Log.LogDebug("CLineModelFitting::getLeastSquareMerit fit = %f", fit);
  if (std::isnan(fit)) {
    Log.LogError("CLineModelFitting::getLeastSquareMerit: NaN value found on "
                 "the lambdarange = (%f, %f)",
                 m_lambdaRange->GetBegin(), m_lambdaRange->GetEnd());
    Log.LogError("CLineModelFitting::getLeastSquareMerit: NaN value found on "
                 "the true observed spectral axis lambdarange = (%f, %f)",
                 spcSpectralAxis[imin], spcSpectralAxis[imax]);
    for (Int32 j = imin; j < imax; j++) {
      if (std::isnan(Yspc[j])) {
        Log.LogError("CLineModelFitting::getLeastSquareMerit: NaN value found "
                     "for the observed spectrum at lambda=%f",
                     spcSpectralAxis[j]);
        break;
      }
      if (std::isnan(Ymodel[j])) {
        Log.LogError("CLineModelFitting::getLeastSquareMerit: NaN value found "
                     "for the model at lambda=%f",
                     spcSpectralAxis[j]);
        break;
      }

      if (std::isnan(ErrorNoContinuum[j])) {
        Log.LogError("CLineModelFitting::getLeastSquareMerit: NaN value found "
                     "for the sqrt(variance) at lambda=%f",
                     spcSpectralAxis[j]);
        break;
      }
      if (ErrorNoContinuum[j] == 0.0) {
        Log.LogError("CLineModelFitting::getLeastSquareMerit: 0 value found "
                     "for the sqrt(variance) at lambda=%f",
                     spcSpectralAxis[j]);
        break;
      }
    }
    THROWG(INTERNAL_ERROR, "NaN value found");
  }
  return fit;
}

Float64 CLineModelFitting::getLeastSquareContinuumMerit() const {
  const CSpectrumSpectralAxis &spcSpectralAxis = m_inputSpc->GetSpectralAxis();
  const CSpectrumFluxAxis &spcFluxAxis = m_model->getSpcFluxAxis();
  const auto &ErrorNoContinuum = m_inputSpc->GetFluxAxis().GetError();

  Int32 numDevs = 0;
  Float64 fit = 0.0;
  const Float64 *YCont = as_const(m_model->getContinuumFluxAxis()).GetSamples();
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

/**
 * \brief Get the squared difference by fast method proposed by D. Vibert
 **/
Float64 CLineModelFitting::getLeastSquareMeritFast(Int32 idxLine) const {
  Float64 fit = getLeastSquareContinuumMeritFast();

  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    if (idxLine != undefIdx && idxLine != iElts) {
      continue;
    }
    Float64 dtm = m_Elements[iElts]->GetSumCross();
    Float64 mtm = m_Elements[iElts]->GetSumGauss();
    Float64 a = m_Elements[iElts]->GetFitAmplitude();
    Float64 term1 = a * a * mtm;
    Float64 term2 = -2. * a * dtm;
    fit += term1 + term2;
  }

  Log.LogDebug("CLineModelFitting::getLeastSquareMerit fit fast = %f", fit);
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
 * \brief Get the scale marginalization correction
 *
 * WARNING: (todo-check) for lm-rules or lm-free, if hybrid method was used to
 *fit, mtm and dtm are not estimated for now...
 **/
Float64 CLineModelFitting::getScaleMargCorrection(Int32 idxLine) const {
  Float64 corr = 0.0;

  // scale marg for continuum
  // corr += getContinuumScaleMargCorrection();

  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    if (idxLine != undefIdx && idxLine != iElts) {
      continue;
    }
    if (m_Elements[iElts]->IsOutsideLambdaRange() == true) {
      continue;
    }
    // if(m_Elements[iElts]->GetElementAmplitude()<=0.0){
    //     continue;
    // }

    Float64 mtm = m_Elements[iElts]->GetSumGauss();
    if (mtm > 0.0) {
      corr += log(mtm);
    }
  }

  return corr;
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
 * @brief CLineModelFitting::getLinesAboveSNR
 * Only considering a list of strong emission lines for now
 * @param snrcut
 * @return
 */
TStringList CLineModelFitting::getLinesAboveSNR(Float64 snrcut) const {
  TInt32List eIdx_oii;
  TInt32List subeIdx_oii;
  TInt32List eIdx_ciii;
  TInt32List subeIdx_ciii;

  Float64 snr_ha = NAN;
  Float64 snr_oii = NAN;
  Float64 snr_oiiia = NAN;
  // snr_oiiib = NAN;
  Float64 snr_hb = NAN;
  Float64 snr_ciii = NAN;
  Float64 snr_lya = NAN;

  for (Int32 iRestLine = 0; iRestLine < m_RestLineList.size(); iRestLine++) {
    Int32 subeIdx = undefIdx;
    Int32 eIdx = m_Elements.findElementIndex(iRestLine, subeIdx);
    if (eIdx < 0 || subeIdx < 0 ||
        m_Elements[eIdx]->IsOutsideLambdaRange(subeIdx))
      continue;

    TPolynomCoeffs polynom_coeffs = getPolynomCoeffs(eIdx);
    Float64 cont = m_Elements[eIdx]->GetContinuumAtCenterProfile(
        subeIdx, m_inputSpc->GetSpectralAxis(), m_Redshift,
        m_model->getContinuumFluxAxis(), polynom_coeffs);
    Float64 mu = NAN;
    Float64 sigma = NAN;
    m_Elements[eIdx]->getObservedPositionAndLineWidth(subeIdx, m_Redshift, mu,
                                                      sigma, false);

    Float64 flux = NAN;
    Float64 fluxError = NAN;
    Float64 fluxDI = NAN;
    Float64 snrDI = NAN;
    TInt32List eIdx_line(1, eIdx);
    TInt32List subeIdx_line(1, subeIdx);
    bool isEmission =
        m_RestLineList[iRestLine].GetType() == CLine::nType_Emission;
    bool opt_cont_substract_abslinesmodel = isEmission;
    getFluxDirectIntegration(eIdx_line, subeIdx_line,
                             opt_cont_substract_abslinesmodel, fluxDI, snrDI);

    if (m_RestLineList[iRestLine].GetName() == linetags::halpha_em &&
        isEmission)
      snr_ha = snrDI;

    if (m_RestLineList[iRestLine].GetName() == linetags::oIIIa_em && isEmission)
      snr_oiiia = snrDI;

    if (m_RestLineList[iRestLine].GetName() == linetags::hbeta_em && isEmission)
      snr_hb = snrDI;

    if (m_RestLineList[iRestLine].GetName() == linetags::lya_em && isEmission)
      snr_lya = snrDI;

    if ((m_RestLineList[iRestLine].GetName() == linetags::oII3726_em ||
         m_RestLineList[iRestLine].GetName() == linetags::oII3729_em) &&
        isEmission) {
      // here we only cover the fluxDI case.
      eIdx_oii.push_back(eIdx);
      subeIdx_oii.push_back(subeIdx);
      fluxDI = NAN;
      snrDI = NAN;
      opt_cont_substract_abslinesmodel = false;
      getFluxDirectIntegration(eIdx_oii, subeIdx_oii,
                               opt_cont_substract_abslinesmodel, fluxDI, snrDI);

      snr_oii = snrDI;
    }
    if ((m_RestLineList[iRestLine].GetName() == linetags::cIII1907_em ||
         m_RestLineList[iRestLine].GetName() == linetags::cIII1909_em) &&
        isEmission) {
      // here we only cover the fluxDI case.
      eIdx_ciii.push_back(eIdx);
      subeIdx_ciii.push_back(subeIdx);
      fluxDI = NAN;
      snrDI = NAN;
      opt_cont_substract_abslinesmodel = false;
      getFluxDirectIntegration(eIdx_ciii, subeIdx_ciii,
                               opt_cont_substract_abslinesmodel, fluxDI, snrDI);

      snr_ciii = snrDI;
    }
  }

  // sum snr ht cut
  Int32 nhtcut = 0;
  TStringList str_above_cut;
  if (snr_ha > snrcut) {
    nhtcut++;
    str_above_cut.push_back("Ha");
  }
  if (snr_oiiia > snrcut) {
    nhtcut++;
    str_above_cut.push_back("OIIIa");
  }
  if (snr_hb > snrcut) {
    nhtcut++;
    str_above_cut.push_back("Hb");
  }
  if (snr_oii > snrcut) {
    nhtcut++;
    str_above_cut.push_back("OII");
  }
  if (snr_lya > snrcut) {
    nhtcut++;
    str_above_cut.push_back("Lya");
  }
  if (snr_ciii > snrcut) {
    nhtcut++;
    str_above_cut.push_back("CIII");
  }

  return str_above_cut;
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
 * @brief CLineModelFitting::GetModelStrongLinePresent
 * @return 1 if there is 1 strong emission line present
 */
bool CLineModelFitting::GetModelStrongEmissionLinePresent() const {
  if (!m_RestLineList.size())
    THROWG(INTERNAL_ERROR, "m_RestframeList is empty");

  bool isStrongPresent = false;

  TInt32List validEltsIdx = m_Elements.GetModelValidElementsIndexes();
  for (Int32 iValidElts = 0; iValidElts < validEltsIdx.size(); iValidElts++) {
    Int32 iElts = validEltsIdx[iValidElts];
    Int32 nlines = m_Elements[iElts]->GetLines().size();
    for (Int32 lineIdx = 0; lineIdx < nlines; lineIdx++) {
      if (!m_Elements[iElts]->m_Lines[lineIdx].GetIsEmission() ||
          !m_Elements[iElts]->m_Lines[lineIdx].GetIsStrong()) {
        continue;
      }

      Float64 amp = m_Elements[iElts]->GetFittedAmplitude(lineIdx);
      if (amp > 0.0) {
        isStrongPresent = true;
        Log.LogDebug("    model: GetModelStrongEmissionLinePresent - found "
                     "Strong EL: %s",
                     m_Elements[iElts]->m_Lines[lineIdx].GetName().c_str());
        break;
      }
    }
    if (isStrongPresent)
      break;
  }

  return isStrongPresent;
}

/**
 * @brief
 * @return 1 if ha em is the strongest line present
 */
bool CLineModelFitting::GetModelHaStrongest() const {

  Float64 ampMax = -DBL_MAX;
  std::string ampMaxLineTag = "";

  TInt32List validEltsIdx = m_Elements.GetModelValidElementsIndexes();
  for (Int32 iValidElts = 0; iValidElts < validEltsIdx.size(); iValidElts++) {
    Int32 iElts = validEltsIdx[iValidElts];
    Int32 nlines = m_Elements[iElts]->GetLines().size();
    for (Int32 lineIdx = 0; lineIdx < nlines; lineIdx++) {
      if (!m_Elements[iElts]->m_Lines[lineIdx].GetIsEmission()) {
        continue;
      }

      Float64 amp = m_Elements[iElts]->GetFittedAmplitude(lineIdx);
      if (amp > 0. && amp > ampMax) {
        ampMaxLineTag = m_Elements[iElts]->m_Lines[lineIdx].GetName().c_str();
        ampMax = amp;
      }
    }
  }

  bool isHaStrongest = (!std::isnan(ampMax) && ampMax > 0. &&
                        ampMaxLineTag == linetags::halpha_em);
  if (isHaStrongest) {
    Log.LogDebug("    model: GetModelHaStrongest - found to be true with "
                 "ampMax=%e (for line=Halpha)",
                 ampMax);
  }
  return isHaStrongest;
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

/**
 * /brief Calls the rules' methods depending on the JSON options.
 * If m_rulesoption is "no", do nothing.
 * If either "balmer" or "all" is in the rules string, call
 *ApplyBalmerRuleLinSolve. If "all" or "ratiorange" is in the rules string,
 *call ApplyAmplitudeRatioRangeRule parameterized for OII. If "all" or
 *"strongweak" is in the rules string, call ApplyStrongHigherWeakRule for
 *emission and then for absorption.
 **/
void CLineModelFitting::applyRules(bool enableLogs) {
  if (m_rulesoption == "no") {
    return;
  }

  m_Regulament.EnableLogs(enableLogs);
  m_Regulament.Apply(m_Elements);
}

const TStringList &CLineModelFitting::GetModelRulesLog() const {
  return m_Regulament.GetLogs();
}

/*
Reset all the model value to the previous solution found.
return 0 if every was ok; else -1
*/
void CLineModelFitting::LoadModelSolution(
    const CLineModelSolution &modelSolution) {

  setRedshift(modelSolution.Redshift, false);
  SetVelocityEmission(modelSolution.EmissionVelocity);
  SetVelocityAbsorption(modelSolution.AbsorptionVelocity);

  const CSpectrumSpectralAxis &spectralAxis = m_inputSpc->GetSpectralAxis();
  for (Int32 iRestLine = 0; iRestLine < m_RestLineList.size(); iRestLine++) {
    Int32 eIdx = modelSolution.ElementId[iRestLine];
    if (eIdx == undefIdx)
      continue;

    m_Elements[eIdx]->prepareSupport(spectralAxis, m_Redshift,
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
        m_Elements[eIdx]->GetObservedPosition(subeIdx, m_Redshift);
    modelSolution.Velocity[iRestLine] = m_Elements[eIdx]->GetVelocity();
    modelSolution.Offset[iRestLine] =
        m_Elements[eIdx]->m_Lines[subeIdx].GetOffset();

    if (opt_level) // brief, to save processing time, do not estimate fluxes
                   // and high level line properties
    {
      modelSolution.FittingError[iRestLine] =
          m_model->getModelErrorUnderElement(eIdx);
      TPolynomCoeffs polynom_coeffs = getPolynomCoeffs(eIdx);
      // save polynom info to output them in hdf5, mainly to recontruct
      // linemeas model
      modelSolution.continuum_pCoeff0[iRestLine] = polynom_coeffs.x0;
      modelSolution.continuum_pCoeff1[iRestLine] = polynom_coeffs.x1;
      modelSolution.continuum_pCoeff2[iRestLine] = polynom_coeffs.x2;

      Float64 cont = m_Elements[eIdx]->GetContinuumAtCenterProfile(
          subeIdx, m_inputSpc->GetSpectralAxis(), m_Redshift,
          m_model->getContinuumFluxAxis(), polynom_coeffs);
      modelSolution.CenterContinuumFlux[iRestLine] = cont;
      modelSolution.ContinuumError[iRestLine] =
          m_model->GetContinuumError(eIdx, subeIdx);
      Float64 mu = NAN;
      Float64 sigma = NAN;
      m_Elements[eIdx]->getObservedPositionAndLineWidth(
          subeIdx, m_Redshift, mu, sigma,
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
      getFluxDirectIntegration(eIdx_line, subeIdx_line,
                               opt_cont_substract_abslinesmodel, fluxDI, snrDI);
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
          getFluxDirectIntegration(eIdx_oii, subeIdx_oii,
                                   opt_cont_substract_abslinesmodel, fluxDI,
                                   snrDI);

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
  modelSolution.EmissionVelocity = m_velocityEmission;
  modelSolution.AbsorptionVelocity = m_velocityAbsorption;
  modelSolution.Redshift = m_Redshift;

  TStringList strongELSNRAboveCut = TStringList(); // getLinesAboveSNR(3.5);
  modelSolution.NLinesAboveSnrCut = strongELSNRAboveCut.size();

  return modelSolution;
}

/**
 * @brief Look for polynom coeffs corresponding to one specific Line
 * Here,we assume that we already fitted one line at a time, which is not
 * really the case
 * TODO: take into consideration lines fitted together for the element in
 * question
 *
 * @param eIdx
 * @return TPolynomCoeffs
 */
TPolynomCoeffs CLineModelFitting::getPolynomCoeffs(Int32 eIdx) const {
  TInt32List xInds = m_Elements.getSupportIndexes({Int32(eIdx)});
  TPolynomCoeffs polynom_coeffs = {0., 0., 0.};
  if (xInds.size() && m_Elements.m_ampOffsetsCoeffs.size()) {
    Int32 idxAmpOffset = m_Elements.getIndexAmpOffset(xInds[0]);
    polynom_coeffs = m_Elements.m_ampOffsetsCoeffs[idxAmpOffset];
  }
  return polynom_coeffs;
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

Float64 CLineModelFitting::GetRedshift() const { return m_Redshift; }

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

void CLineModelFitting::SetContinuumComponent(std::string component) {
  m_ContinuumComponent = component;
  m_model->SetContinuumComponent(component);
  m_continuumManager->setContinuumComponent(component);
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

/**
 * @brief CLineModelFitting::SetMultilineNominalAmplitudesFast
 * This method sets the linemodel unique elt nominal amplitudes to the
 * corresponding value of the iCatalog st catalog. INFO: fast method,
 * InitLineCorrespondence() should have been called previously with the same
 * LineModelElementList arg.
 * @param iCatalog
 * @return
 */
bool CLineModelFitting::SetMultilineNominalAmplitudesFast(Int32 iCatalog) {
  if (iCatalog < 0) {
    return false;
  }
  Float64 nominalAmp = 0.0;
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    Int32 nLines = m_Elements[iElts]->GetSize();
    for (Int32 j = 0; j < nLines; j++) {
      nominalAmp = m_CatalogTplRatio
                       .getNominalAmplitudeCorrespondance()[iElts][iCatalog][j];
      m_Elements[iElts]->SetNominalAmplitude(j, nominalAmp);
    }
  }
  return true;
}
