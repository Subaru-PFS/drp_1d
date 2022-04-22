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
#include "RedshiftLibrary/line/catalogsTplShape.h"
#include "RedshiftLibrary/line/line.h"
#include "RedshiftLibrary/line/regulament.h"
#include "RedshiftLibrary/linemodel/element.h"
#include "RedshiftLibrary/linemodel/lmfitfunctions.h"

#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/continuum/irregularsamplingmedian.h"
#include "RedshiftLibrary/extremum/extremum.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include <boost/chrono/thread_clock.hpp>
#include <boost/format.hpp>

#include <float.h>
#include <gsl/gsl_blas.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_spline.h>
#include <math.h>

#include <boost/numeric/conversion/bounds.hpp>
#include <gsl/gsl_multifit_nlin.h>
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
CLineModelFitting::CLineModelFitting(
    const CSpectrum &spectrum, const TFloat64Range &lambdaRange,
    const CTemplateCatalog &tplCatalog, const TStringList &tplCategoryList,
    const CLineCatalog::TLineVector &restLineList,
    const std::string &opt_fittingmethod,
    const std::string &opt_continuumcomponent,
    Float64 opt_continuum_neg_threshold, const std::string &widthType,
    Float64 nsigmasupport, Float64 velocityEmission, Float64 velocityAbsorption,
    const std::string &opt_rules, const std::string &opt_rigidity,
    Int32 amplitudeOffsetsDegree)
    : m_inputSpc(spectrum), m_lambdaRange(lambdaRange),
      m_tplCatalog(tplCatalog), m_tplCategoryList(tplCategoryList),
      m_RestLineList(restLineList), m_fittingmethod(opt_fittingmethod),
      m_opt_fitcontinuum_neg_threshold(opt_continuum_neg_threshold),
      m_LineWidthType(widthType), m_NSigmaSupport(nsigmasupport),
      m_velocityEmission(velocityEmission),
      m_velocityAbsorption(velocityAbsorption),
      m_velocityEmissionInit(m_velocityEmission),
      m_velocityAbsorptionInit(m_velocityAbsorption), m_rulesoption(opt_rules),
      m_rigidity(opt_rigidity), m_Regulament(),
      m_ErrorNoContinuum(m_spcFluxAxisNoContinuum.GetError()),
      m_templateFittingOperator(spectrum, lambdaRange),
      m_enableAmplitudeOffsets(false),
      m_AmplitudeOffsetsDegree(amplitudeOffsetsDegree) {

  // m_nominalWidthDefaultEmission = 1.15;// suited to new pfs simulations
  m_nominalWidthDefaultEmission = 13.4; // euclid 1 px
  m_nominalWidthDefaultAbsorption = m_nominalWidthDefaultEmission;

  m_enableLambdaOffsetsFit =
      true; // enable lambdaOffsetFit. Once enabled, the offset fixed value or
            // the fitting on/off switch is done through the offset calibration
            // file.

  m_SpectrumModel = CSpectrum(spectrum);
  m_SpcCorrectedUnderLines = CSpectrum(spectrum);
  const Int32 spectrumSampleCount = spectrum.GetSampleCount();
  m_SpcFluxAxis.SetSize(spectrumSampleCount);
  Log.LogDetail("    model: Continuum winsize found is %.2f A",
                m_inputSpc.GetMedianWinsize());

  m_ContinuumFluxAxis.SetSize(spectrumSampleCount);
  m_SpcFluxAxisModelDerivVelEmi.SetSize(spectrumSampleCount);
  m_SpcFluxAxisModelDerivVelAbs.SetSize(spectrumSampleCount);
  const CSpectrumFluxAxis &spectrumFluxAxis = spectrum.GetFluxAxis();

  m_spcFluxAxisNoContinuum.SetSize(spectrumSampleCount);
  m_ErrorNoContinuum = spectrumFluxAxis.GetError(); // sets the error vector

  SetContinuumComponent(opt_continuumcomponent);

  // NB: fitContinuum_option: this is the initialization (default value),
  // eventually overriden in SetFitContinuum_FitStore() when a fitStore gets
  // available
  m_fitContinuum_option =
      0; // 0=interactive fitting, 1=use precomputed fit store, 2=use fixed
         // values (typical use for second pass recompute)

  if (opt_fittingmethod == "lmfit") {
    m_lmfit_noContinuumTemplate = (m_ContinuumComponent == "fromspectrum");
    m_lmfit_bestTemplate = true;
    m_lmfit_fitContinuum = true;
    m_lmfit_fitEmissionVelocity = true;
    m_lmfit_fitAbsorptionVelocity = true;
    Log.LogInfo(
        "LMfit parameters : tplAlreadyFit %d ; dobestTemplate %d, fitContinuum "
        "%d, fitEmissionVelocity %d, fitAbsorptionVelocity %d",
        m_lmfit_noContinuumTemplate, m_lmfit_bestTemplate, m_lmfit_fitContinuum,
        m_lmfit_fitEmissionVelocity, m_lmfit_fitAbsorptionVelocity);
  }
  //*/

  // "New style" rules initialization:
  m_Regulament.CreateRulesFromJSONFiles();
  m_Regulament.EnableRulesAccordingToParameters(m_rulesoption);

  // Load the line catalog
  Log.LogDebug("About to load line catalog.");
  if (m_rigidity != "tplshape") {
    // load the regular catalog
    LoadCatalog(restLineList);
  } else {
    // load the tplshape catalog with only 1 element for all lines
    // LoadCatalogOneMultiline(restLineList);
    // load the tplshape catalog with 2 elements: 1 for the Em lines + 1 for the
    // Abs lines
    LoadCatalogTwoMultilinesAE(restLineList);
  }
  LogCatalogInfos();

  /*
  //check the continuum flux axis
  const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel.GetSpectralAxis();
  for(Int32 j=0; j<m_ContinuumFluxAxis.GetSamplesCount(); j++)
  {
      if(isnan(m_ContinuumFluxAxis[j]))
      {
          THROWG(INTERNAL_ERROR, "NaN value found for the ContinuumFluxAxis at
  lambda=%f", spectralAxis[j] );
      }
  }
  */

  SetLSF();
}

// hook
bool CLineModelFitting::initTplratioCatalogs(Int32 opt_tplratio_ismFit) {
  // TODO: use the passed tplRatioCatalog
  // TODO: check if m_CatalogTplShape changes between iterations
  bool ret =
      m_CatalogTplShape.Init(opt_tplratio_ismFit,
                             m_tplCatalog.GetTemplate(m_tplCategoryList[0], 0)
                                 ->m_ismCorrectionCalzetti,
                             m_NSigmaSupport);
  if (!ret) {
    Log.LogError(
        "Unable to initialize the the tpl-shape catalogs. aborting...");
    return false;
  }
  m_CatalogTplShape.InitLineCorrespondingAmplitudes(m_Elements);
  SetMultilineNominalAmplitudesFast(0);
  // m_CatalogTplShape.SetMultilineNominalAmplitudes( *this, 0 );
  // m_RestLineList = m_CatalogTplShape.GetRestLinesList(0);
  // LoadCatalog(m_RestLineList);
  // LogCatalogInfos();

  // Resize tplshape buffers
  m_ChisquareTplshape.resize(m_CatalogTplShape.GetCatalogsCount());
  m_ScaleMargCorrTplshape.resize(m_CatalogTplShape.GetCatalogsCount());
  m_StrongELPresentTplshape.resize(m_CatalogTplShape.GetCatalogsCount());
  m_StrongHalphaELPresentTplshape.resize(m_CatalogTplShape.GetCatalogsCount());
  m_NLinesAboveSNRTplshape.resize(m_CatalogTplShape.GetCatalogsCount());
  m_FittedAmpTplshape.resize(m_CatalogTplShape.GetCatalogsCount());
  m_LyaAsymCoeffTplshape.resize(m_CatalogTplShape.GetCatalogsCount());
  m_LyaWidthCoeffTplshape.resize(m_CatalogTplShape.GetCatalogsCount());
  m_LyaDeltaCoeffTplshape.resize(m_CatalogTplShape.GetCatalogsCount());
  m_LyaIgmIdxTplshape.resize(m_CatalogTplShape.GetCatalogsCount());
  m_FittedErrorTplshape.resize(m_CatalogTplShape.GetCatalogsCount());
  m_MtmTplshape.resize(m_CatalogTplShape.GetCatalogsCount());
  m_DtmTplshape.resize(m_CatalogTplShape.GetCatalogsCount());
  m_LinesLogPriorTplshape.resize(m_CatalogTplShape.GetCatalogsCount());
  for (Int32 ktplshape = 0; ktplshape < m_CatalogTplShape.GetCatalogsCount();
       ktplshape++) {
    m_FittedAmpTplshape[ktplshape].resize(m_Elements.size());
    m_FittedErrorTplshape[ktplshape].resize(m_Elements.size());
    m_MtmTplshape[ktplshape].resize(m_Elements.size());
    m_DtmTplshape[ktplshape].resize(m_Elements.size());
    m_LyaAsymCoeffTplshape[ktplshape].resize(m_Elements.size());
    m_LyaWidthCoeffTplshape[ktplshape].resize(m_Elements.size());
    m_LyaDeltaCoeffTplshape[ktplshape].resize(m_Elements.size());
    m_LyaIgmIdxTplshape[ktplshape].resize(m_Elements.size());
    m_LinesLogPriorTplshape[ktplshape].resize(m_Elements.size());
  }

  m_tplshapeLeastSquareFast = false;
  return true;
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
    m_fittingmethod = m_opt_firstpass_fittingmethod;
    m_forceLyaFitting = false;
  }
  if (iPass == 2) {
    m_forceDisableLyaFitting = m_opt_lya_forcedisablefit;
    m_forcedisableTplratioISMfit = false;
    m_forcedisableMultipleContinuumfit = false;

    m_fittingmethod = m_opt_secondpass_fittingmethod;
    m_forceLyaFitting = m_opt_lya_forcefit;
    Log.LogDetail(
        "    model: set forceLyaFitting ASYMFIT for Tpl-ratio mode : %d",
        m_forceLyaFitting);
  }
  if (iPass == 3) {
    m_forceDisableLyaFitting = false;

    m_forcedisableMultipleContinuumfit = false;
    m_enableAmplitudeOffsets = true;
  }

  return true;
}
Int32 CLineModelFitting::GetPassNumber() const { return m_pass; }

void CLineModelFitting::SetForcedisableTplratioISMfit(bool opt) {
  m_forcedisableTplratioISMfit = opt;
}

/**
 * \brief Returns a pointer to m_SpectrumModel.
 **/
const CSpectrum &CLineModelFitting::GetModelSpectrum() const {
  return m_SpectrumModel;
}

CSpectrum CLineModelFitting::GetSpectrumModelContinuum() const {
  CSpectrum spc(m_SpectrumModel);
  spc.SetFluxAxis(m_ContinuumFluxAxis);
  return spc;
}

/**
 * \brief Returns a pointer to a spectrum containing the observed spectrum with
 *the fitted lines subtracted
 **/
const CSpectrum &
CLineModelFitting::GetObservedSpectrumWithLinesRemoved(Int32 lineTypeFilter) {
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  const CSpectrumFluxAxis &fluxAxis = m_SpectrumModel.GetFluxAxis();

  //
  TFloat64List fluxAndContinuum(spectralAxis.GetSamplesCount(), 0.0);
  if (lineTypeFilter == CLine::nType_Emission) {
    refreshModel(CLine::nType_Absorption);
    fluxAndContinuum = fluxAxis.GetSamplesVector();
    refreshModel(CLine::nType_Emission);
  } else if (lineTypeFilter == CLine::nType_Absorption) {
    refreshModel(CLine::nType_Emission);
    fluxAndContinuum = fluxAxis.GetSamplesVector();
    refreshModel(CLine::nType_Absorption);
  }

  CSpectrumFluxAxis fluxAxisNothingUnderLines(
      m_SpcCorrectedUnderLines.GetSampleCount());
  TAxisSampleList &Y = fluxAxisNothingUnderLines.GetSamplesVector();
  const auto &SpcFluxAxis = m_SpcFluxAxis;
  const auto &ContinuumFluxAxis = m_ContinuumFluxAxis;
  for (Int32 t = 0; t < spectralAxis.GetSamplesCount(); t++) {
    Y[t] = SpcFluxAxis[t] - fluxAxis[t] + ContinuumFluxAxis[t];
  }

  bool enableSmoothBlendContinuUnderLines = true;
  if (enableSmoothBlendContinuUnderLines) {
    Float64 alphaMax = 0.9; // alpha blend = 0: only lineSubtractedFlux,
                            // alpha=1: only continuum

    TInt32List validEltsIdx = m_Elements.GetModelValidElementsIndexes();
    TInt32List nonZeroValidEltsIdx;
    for (Int32 j = 0; j < validEltsIdx.size(); j++) {
      Int32 lineType = m_Elements[validEltsIdx[j]]->m_Lines[0].GetType();
      if (!(lineTypeFilter == -1 || lineTypeFilter == lineType)) {
        continue;
      }

      if (m_Elements[validEltsIdx[j]]->GetElementAmplitude() > 0.0) {
        nonZeroValidEltsIdx.push_back(validEltsIdx[j]);
      }
    }
    if (nonZeroValidEltsIdx.size() > 0) {
      TInt32List supportIdxes =
          m_Elements.getSupportIndexes(nonZeroValidEltsIdx);
      if (supportIdxes.size() > 0) {
        for (Int32 i = 1; i < supportIdxes.size(); i++) {
          Float64 weighting = GetWeightingAnyLineCenterProximity(
              supportIdxes[i], nonZeroValidEltsIdx);
          Float64 alpha = alphaMax * weighting;
          Y[supportIdxes[i]] = (1. - alpha) * Y[supportIdxes[i]] +
                               alpha * (fluxAndContinuum[supportIdxes[i]]);
        }
      }
    }
  }

  m_SpcCorrectedUnderLines.SetFluxAxis(std::move(fluxAxisNothingUnderLines));

  return m_SpcCorrectedUnderLines;
}

/**
 * @brief GetContinuumError
 * Estimate the error on the continuum in a given window around the line center.
 * 1. calculate the observed spectrum flux with lines subtracted (fitted line
 * model)
 * 2. estimate the std in a window corresponding to the sliding median window
 * size (default/hardcoded=150A)
 * @param subeIdx
 * @param spectralAxis
 * @param spectrumfluxAxis
 * @param redshift
 * @param continuumfluxAxis
 * @return -1: zero samples found for error estimation, NAN: not enough samples
 * found for error estimation
 */
Float64 CLineModelFitting::GetContinuumError(Int32 eIdx, Int32 subeIdx) {
  Int32 nMinValueForErrorEstimation = 10;

  const CSpectrum &noLinesSpectrum = GetObservedSpectrumWithLinesRemoved();
  const CSpectrumFluxAxis &noLinesFluxAxis = noLinesSpectrum.GetFluxAxis();
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  const TFloat64Range lambdaRange =
      spectralAxis.GetLambdaRange(); // using the full wavelength range for this
                                     // error estimation
  const auto &ContinuumFluxAxis = m_ContinuumFluxAxis;
  Float64 winsizeAngstrom = 150.;

  Float64 mu = m_Elements[eIdx]->GetObservedPosition(subeIdx, m_Redshift);
  TInt32Range indexRange = CLineModelElement::EstimateIndexRange(
      spectralAxis, mu, lambdaRange, winsizeAngstrom);

  // estimate sum square error between continuum and nolines-spectrum
  Float64 sum = 0.0;
  Int32 nsum = 0;
  for (Int32 t = indexRange.GetBegin(); t < indexRange.GetEnd(); t++) {
    Float64 diff = noLinesFluxAxis[t] - ContinuumFluxAxis[t];
    sum += diff * diff;
    nsum++;
  }

  Float64 error = NAN;
  if (nsum >= nMinValueForErrorEstimation) {
    error = sqrt(sum / nsum);
  }

  return error;
}

/**
 * \brief Returns a pointer to the (re-)estimated continuum flux.
 **/
const CSpectrumFluxAxis &CLineModelFitting::GetModelContinuum() const {
  return m_ContinuumFluxAxis;
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
  TInt32List indexes =
      getlambdaIndexesUnderLines(eIdx_list, subeIdx_list, N_SIGMA_SUPPORT_DI);

  if (!indexes.size())
    THROWG(INTERNAL_ERROR, "empty indexes");

  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  const auto &ContinuumFluxAxis = m_ContinuumFluxAxis;
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
  for (const auto &t : indexes) {
    continuumFlux[t] =
        substract_abslinesmodel ? absLinesModelFlux[t] : ContinuumFluxAxis[t];
    if (m_enableAmplitudeOffsets)
      continuumFlux[t] += polynom_coeffs.x0 +
                          polynom_coeffs.x1 * spectralAxis[t] +
                          polynom_coeffs.x2 * spectralAxis[t] * spectralAxis[t];
  }
  // substarct continuum from spectrum flux
  const auto &SpcFluxAxis = m_SpcFluxAxis;
  CSpectrumFluxAxis fluxMinusContinuum(SpcFluxAxis);
  for (const auto &t : indexes) {
    fluxMinusContinuum[t] -= continuumFlux[t];
  }
  // estimate the integrated flux between obs. spectrum and continuum:
  // trapezoidal intg
  Float64 sumFlux = 0.0;
  Float64 sumErr = 0.0;
  integrateFluxes_usingTrapez(fluxMinusContinuum, indexes, sumFlux, sumErr);

  if (sumErr <= 0)
    return;

  fluxdi = sumFlux;
  snrdi = std::abs(fluxdi) / sqrt(sumErr);
  return;
}

TInt32List CLineModelFitting::getlambdaIndexesUnderLines(
    const TInt32List &eIdx_list, const TInt32List &subeIdx_list,
    const Float64 &sigma_support) const {

  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();

  TInt32List indexes;
  for (Int32 kl = 0; kl < eIdx_list.size(); kl++) {
    Int32 eIdx = eIdx_list[kl];
    Int32 subeIdx = subeIdx_list[kl];

    Float64 mu = NAN;
    Float64 LineWidth = NAN;
    m_Elements[eIdx]->getObservedPositionAndLineWidth(subeIdx, m_Redshift, mu,
                                                      LineWidth);

    Float64 winsizeAngstrom = LineWidth * sigma_support;

    TInt32Range indexRange = CLineModelElement::EstimateIndexRange(
        spectralAxis, mu, m_lambdaRange, winsizeAngstrom);

    auto newindices = indexRange.SpreadOver(1);
    indexes.insert(indexes.end(), newindices.begin(), newindices.end());
  }
  // remove duplicate indexes
  std::sort(indexes.begin(), indexes.end());
  indexes.erase(std::unique(indexes.begin(), indexes.end()), indexes.end());
  return indexes;
}

/**
 * @brief Construct a new clinemodelfitting::integratefluxes usingtrapez object
 *
 * @param fluxMinusContinuum
 * @param indexes
 * @param sumFlux
 * @param sumErr
 */
void CLineModelFitting::integrateFluxes_usingTrapez(
    const CSpectrumFluxAxis &fluxMinusContinuum, const TInt32List &indexes,
    Float64 &sumFlux, Float64 &sumErr) const {

  sumFlux = 0.0;
  sumErr = 0.0;

  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  const auto &ErrorNoContinuum = m_ErrorNoContinuum;
  for (const auto &t : indexes) {
    Float64 trapweight = (spectralAxis[t + 1] - spectralAxis[t]) * 0.5;
    sumFlux += trapweight * (fluxMinusContinuum[t + 1] + fluxMinusContinuum[t]);

    Float64 ea = ErrorNoContinuum[t] * ErrorNoContinuum[t];
    Float64 eb = ErrorNoContinuum[t + 1] * ErrorNoContinuum[t + 1];
    sumErr += trapweight * trapweight * (eb + ea);
  }
  return;
}

Float64 CLineModelFitting::getModelFluxVal(Int32 idx) const {
  const CSpectrumFluxAxis &modelFluxAxis = m_SpectrumModel.GetFluxAxis();
  if (idx < modelFluxAxis.GetSamplesCount()) {
    return modelFluxAxis[idx];
  }

  THROWG(INTERNAL_ERROR, "out-of-bound index");
}

Float64 CLineModelFitting::getModelFluxDerivEltVal(Int32 DerivEltIdx,
                                                   Int32 idx) const {

  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  Int32 iElts = DerivEltIdx;
  if (idx < spectralAxis.GetSamplesCount()) {
    Float64 derivateVal = m_Elements[iElts]->GetModelDerivAmplitudeAtLambda(
        spectralAxis[idx], m_Redshift, m_ContinuumFluxAxis[idx]);
    return derivateVal;
  }

  THROWG(INTERNAL_ERROR, "out-of-bound index");
}

Float64
CLineModelFitting::getModelFluxDerivContinuumAmpEltVal(Int32 DerivEltIdx,
                                                       Int32 idx) const {

  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  Int32 iElts = DerivEltIdx;
  if (idx < spectralAxis.GetSamplesCount()) {
    Float64 derivateVal = m_Elements[iElts]->GetModelDerivContinuumAmpAtLambda(
        spectralAxis[idx], m_Redshift, m_observeGridContinuumFlux[idx]);
    return derivateVal;
  }

  THROWG(INTERNAL_ERROR, "out-of-bound index");
}

Float64 CLineModelFitting::getModelFluxDerivZContinuumVal(Int32 idx) const {
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  if (idx < spectralAxis.GetSamplesCount()) {
    Float64 derivateVal = m_ContinuumFluxAxis[idx] * (-1) * spectralAxis[idx] /
                          (1 + m_Redshift) / (1 + m_Redshift);
    return derivateVal;
  }
  // if(idx< spectralAxis.GetSamplesCount()){
  //   return m_unscaleContinuumFluxAxisDerivZ[idx] *
  //   m_fitContinuum_tplFitAmplitude;
  // }
  THROWG(INTERNAL_ERROR, "out-of-bound index");
}

// void CLineModelFitting::calculateUnscaleContinuumDerivZ(){
//
//   CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel.GetSpectralAxis();
//
//   if(m_unscaleContinuumFluxAxisDerivZ == NULL){
//     m_unscaleContinuumFluxAxisDerivZ = new
//     Float64[spectralAxis.GetSamplesCount()]();
//     if(m_unscaleContinuumFluxAxisDerivZ == NULL){
//       Log.LogError("unable to allocate m_unscaleContinuumFluxAxisDerivZ");
//       return;
//     }
//   }
//   for (int idx = 0;idx<spectralAxis.GetSamplesCount(); idx ++ ){
//     m_unscaleContinuumFluxAxisDerivZ[idx] = m_observeGridContinuumFlux[idx] *
//     spectralAxis[idx] / (1+m_Redshift) / (1+m_Redshift);
//   }
//
// }

/*
 * Calculate the partial deriv of the flux by z, when the continuum is not a
 * variable of z
 */
Float64
CLineModelFitting::getModelFluxDerivZEltValNoContinuum(Int32 DerivEltIdx,
                                                       Int32 idx) const {
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  Int32 iElts = DerivEltIdx;
  if (idx < spectralAxis.GetSamplesCount()) {
    Float64 derivateVal = m_Elements[iElts]->GetModelDerivZAtLambdaNoContinuum(
        spectralAxis[idx], m_Redshift, m_ContinuumFluxAxis[idx]);
    return derivateVal;
  }

  THROWG(INTERNAL_ERROR, "out-of-bound index");
}

/*
 * Calculate the partial deriv of the flux by z, when the continuum is a
 * variable of z, the partial deriv of z of the continuum at this lambda is
 * given by continuumFluxDerivZ
 */
Float64
CLineModelFitting::getModelFluxDerivZEltVal(Int32 DerivEltIdx, Int32 idx,
                                            Float64 continuumFluxDerivZ) const {
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  Int32 iElts = DerivEltIdx;
  if (idx < spectralAxis.GetSamplesCount()) {
    Float64 derivateVal = m_Elements[iElts]->GetModelDerivZAtLambda(
        spectralAxis[idx], m_Redshift, m_ContinuumFluxAxis[idx],
        continuumFluxDerivZ);
    return derivateVal;
  }

  THROWG(INTERNAL_ERROR, "out-of-bound index");
}

/*
Change on this function, m_SpcFluxAxisModelDerivVel has been remove,
replace by m_SpcFluxAxisModelDerivVelEmi and m_SpcFluxAxisModelDerivVelAbs
*/
Float64 CLineModelFitting::getModelFluxDerivVelVal(Int32 idx) const {
  if (idx < m_SpcFluxAxisModelDerivVelEmi.GetSamplesCount()) {
    return m_SpcFluxAxisModelDerivVelEmi[idx] +
           m_SpcFluxAxisModelDerivVelAbs[idx];
  }

  THROWG(INTERNAL_ERROR, "out-of-bound index");
}

Float64 CLineModelFitting::getModelFluxDerivVelEmissionVal(Int32 idx) const {
  if (idx < m_SpcFluxAxisModelDerivVelEmi.GetSamplesCount()) {
    return m_SpcFluxAxisModelDerivVelEmi[idx];
  }

  THROWG(INTERNAL_ERROR, "out-of-bound index");
}

Float64 CLineModelFitting::getModelFluxDerivVelAbsorptionVal(Int32 idx) const {
  if (idx < m_SpcFluxAxisModelDerivVelAbs.GetSamplesCount()) {
    return m_SpcFluxAxisModelDerivVelAbs[idx];
  }

  THROWG(INTERNAL_ERROR, "out-of-bound index");
}

/**
 * \brief For each line in each group of the argument, finds the associated line
 *in the catalog and saves this information to m_Elements. Converts the argument
 *restLineList to a group list. For each entry in this list: For each line in
 *this entry: Finds the index in the catalog from the line name and type. Saves
 *the line, the catalog index and the nominal amplitude for the line thusly
 *associated to this line. If at least one line was found, save this result in
 *m_Elements.
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

void CLineModelFitting::LoadFitContinuumOneTemplate(
    const std::shared_ptr<const CTemplate> &tpl) {
  Float64 merit = INFINITY;

  Float64 fitContinuumAmplitude = NAN;
  Float64 fitContinuumAmplitudeError = NAN;
  Float64 fitContinuumAmplitudeSigma = NAN;
  Float64 FitEbmvCoeff = -1.0;
  Int32 fitMeiksinIdx = -1;
  Float64 fitDtM = NAN;
  Float64 fitMtM = NAN;
  Float64 fitLogprior = NAN;
  Float64 overlapThreshold = 1.0;

  bool ignoreLinesSupport = m_secondpass_fitContinuum_outsidelinesmask;
  std::vector<CMask> maskList;
  if (ignoreLinesSupport) {
    maskList.resize(1);
    maskList[0] = getOutsideLinesMask();
  }
  TFloat64List redshifts(1, m_Redshift);
  std::string opt_interp = "precomputedfinegrid"; //"lin";
  Int32 opt_extinction = m_secondpass_fitContinuum_igm;
  Int32 opt_dustFit = m_secondpass_fitContinuum_dustfit;

  if (m_observeGridContinuumFlux.empty()) {
    THROWG(INTERNAL_ERROR, "Elementlist, cannot SolveContinuum "
                           "without m_observeGridContinuumFlux");
  }
  bool ret =
      SolveContinuum(tpl, redshifts, overlapThreshold, maskList, opt_interp,
                     opt_extinction, opt_dustFit, merit, fitContinuumAmplitude,
                     fitContinuumAmplitudeError, fitContinuumAmplitudeSigma,
                     FitEbmvCoeff, fitMeiksinIdx, fitDtM, fitMtM, fitLogprior);
  /*//debug
  // export for debug
  FILE* fspc = fopen( "Continuum.txt", "w+" );
  for( Int32 t=0;t<m_ContinuumFluxAxis.GetSamplesCount();t++)
  {
      fprintf( fspc, "%f %f %f\n", spectralAxis[t],m_ContinuumFluxAxis[t]);
  }
  fclose( fspc );
  //*/
  Log.LogInfo("LMfit : Continuum amplitude set Init at %0.00f",
              fitContinuumAmplitude);
  Log.LogInfo("LMfit : Continuum amplitude error set Init at %0.00f",
              fitContinuumAmplitudeError);
  Log.LogInfo("LMfit : Continuum amplitude/error is set at %0.00f",
              fitContinuumAmplitudeSigma);
  Log.LogInfo("LMfit : Solve succes %s", (ret ? "true" : "false"));
  ApplyContinuumOnGrid(tpl, m_Redshift);
}

/**
 * \brief Generates a continuum from the fitting with a set of templates : uses
 * the templatefitting operator
 * TODO: LoadFitContinuum should be limited to reading continuum values from the
 * variable class, especially that we want that continuum fitting results are
 * saved in tplfitStore container outside CElementList and these stores will be
 * injected in the class whenever required !
 * TODO: study this possibility before doing the change
 */
void CLineModelFitting::LoadFitContinuum(Int32 icontinuum, Int32 autoSelect) {
  Log.LogDebug("Elementlist, m_fitContinuum_option=%d", m_fitContinuum_option);
  if (m_observeGridContinuumFlux.empty()) {
    THROWG(INTERNAL_ERROR,
           "Cannot loadfitcontinuum without precomputedGridTplFlux");
  }

  if (m_fitContinuum_option ==
      1) { // using precomputed fit store, i.e., fitValues
    CTemplatesFitStore::TemplateFitValues fitValues =
        m_fitContinuum_tplfitStore->GetFitValues(m_Redshift, icontinuum);
    if (fitValues.tplName.empty()) {
      THROWG(INTERNAL_ERROR, "Empty template name");
    }

    m_fitContinuum_tplName = fitValues.tplName;
    m_fitContinuum_tplFitAmplitude = fitValues.fitAmplitude;
    m_fitContinuum_tplFitAmplitudeError = fitValues.fitAmplitudeError;
    m_fitContinuum_tplFitMerit = fitValues.merit;
    m_fitContinuum_tplFitMerit_phot = fitValues.chiSquare_phot;
    m_fitContinuum_tplFitEbmvCoeff = fitValues.ismEbmvCoeff;
    m_fitContinuum_tplFitMeiksinIdx = fitValues.igmMeiksinIdx;
    m_fitContinuum_tplFitRedshift = m_Redshift;
    m_fitContinuum_tplFitDtM = fitValues.fitDtM;
    m_fitContinuum_tplFitMtM = fitValues.fitMtM;
    m_fitContinuum_tplFitLogprior = fitValues.logprior;
    m_fitContinuum_tplFitPolyCoeffs = {};

    m_fitContinuum_tplFitAmplitudeSigmaMAX =
        m_fitContinuum_tplfitStore->m_fitContinuum_fitAmplitudeSigmaMAX;
    m_fitContinuum_tplFitSNRMax =
        m_fitContinuum_tplfitStore->m_fitContinuum_tplFitSNRMax;
    m_opt_fitcontinuum_maxCount =
        m_fitContinuum_tplfitStore->m_opt_fitcontinuum_maxCount;

  } else if (m_fitContinuum_option == 2) {
    // values unmodified nothing to do
  } else {
    THROWG(INTERNAL_ERROR, "Cannot parse fitContinuum_option");
  }

  if (!m_fitContinuum_tplName.empty()) {
    // Retrieve the best template
    std::shared_ptr<const CTemplate> tpl = m_tplCatalog.GetTemplateByName(
        m_tplCategoryList, m_fitContinuum_tplName);
    if (tpl->GetName() == m_fitContinuum_tplName) {
      if (autoSelect) {
        // Float64 contsnr = getFitContinuum_snr();
        /*Float64 contsnr = m_fitContinuum_tplFitSNRMax;
        m_fitContinuum_tplFitAlpha = 1.0;
        if(contsnr>50.)
        {
            m_fitContinuum_tplFitAlpha=0.0;
        }*/
        m_fitContinuum_tplFitAlpha = 0.0;
        if (m_fitContinuum_tplFitAmplitudeSigmaMAX <
            m_opt_fitcontinuum_neg_threshold)
          m_fitContinuum_tplFitAlpha = 1.0; // switch to spectrum continuum
      }

      ApplyContinuumOnGrid(tpl, m_fitContinuum_tplFitRedshift);

      setFitContinuum_tplAmplitude(m_fitContinuum_tplFitAmplitude,
                                   m_fitContinuum_tplFitAmplitudeError,
                                   m_fitContinuum_tplFitPolyCoeffs);

      Log.LogDebug("    model : LoadFitContinuum, loaded: %s",
                   m_fitContinuum_tplName.c_str());
      Log.LogDebug("    model : LoadFitContinuum, loaded with A=%e",
                   m_fitContinuum_tplFitAmplitude);
      Log.LogDebug("    model : LoadFitContinuum, loaded with A_error=%e",
                   m_fitContinuum_tplFitAmplitudeError);
      Log.LogDebug("    model : LoadFitContinuum, loaded with DustCoeff=%e",
                   m_fitContinuum_tplFitEbmvCoeff);
      Log.LogDebug("    model : LoadFitContinuum, loaded with MeiksinIdx=%d",
                   m_fitContinuum_tplFitMeiksinIdx);
      Log.LogDebug("    model : LoadFitContinuum, loaded with dtm=%e",
                   m_fitContinuum_tplFitDtM);
      Log.LogDebug("    model : LoadFitContinuum, loaded with mtm=%e",
                   m_fitContinuum_tplFitMtM);
      Log.LogDebug("    model : LoadFitContinuum, loaded with logprior=%e",
                   m_fitContinuum_tplFitLogprior);
      Float64 tplfitsnr = NAN;
      if (m_fitContinuum_tplFitMtM > 0.0) {
        tplfitsnr =
            m_fitContinuum_tplFitDtM / std::sqrt(m_fitContinuum_tplFitMtM);
      }
      Log.LogDebug("    model : LoadFitContinuum, loaded with snr=%e",
                   tplfitsnr);
    } else {
      THROWG(
          INTERNAL_ERROR,
          Formatter()
              << "Failed to load-fit continuum. Failed to find best template="
              << m_fitContinuum_tplName.c_str());
    }
  } else {
    THROWG(INTERNAL_ERROR, Formatter()
                               << "Failed to load-fit continuum for cfitopt="
                               << m_fitContinuum_option);
  }
}

void CLineModelFitting::setFitContinuum_tplAmplitude(
    Float64 tplAmp, Float64 tplAmpErr, const TFloat64List &polyCoeffs) {
  const CSpectrumSpectralAxis &spcSpectralAxis =
      m_SpectrumModel.GetSpectralAxis();
  const auto &SpcFluxAxis = m_SpcFluxAxis;

  Float64 alpha =
      m_fitContinuum_tplFitAlpha; // alpha blend = 1: only
                                  // m_inputSpc->GetContinuumFluxAxis(),
                                  // alpha=0: only tplfit

  m_fitContinuum_tplFitAmplitude = tplAmp;
  m_fitContinuum_tplFitAmplitudeError = tplAmpErr;
  m_fitContinuum_tplFitPolyCoeffs = polyCoeffs;
  for (Int32 k = 0; k < m_ContinuumFluxAxis.GetSamplesCount(); k++) {
    if (alpha == 1.0)
      m_ContinuumFluxAxis[k] = 0.;
    else
      m_ContinuumFluxAxis[k] =
          (1. - alpha) * m_observeGridContinuumFlux[k] * tplAmp;
    if (alpha != 0.0) {
      m_ContinuumFluxAxis[k] += alpha * m_inputSpc.GetContinuumFluxAxis()[k];
    }
    if (alpha != 1.0) {
      Float64 lbdaTerm = 1.0;
      for (Int32 kCoeff = 0; kCoeff < polyCoeffs.size(); kCoeff++) {
        m_ContinuumFluxAxis[k] += (1. - alpha) * polyCoeffs[kCoeff] * lbdaTerm;
        lbdaTerm *= spcSpectralAxis[k];
      }
    }
    m_spcFluxAxisNoContinuum[k] =
        SpcFluxAxis[k] - as_const(m_ContinuumFluxAxis)[k];
  }
}

/*
Change the actual value of redshift.
the continuum can be reinterpolate.
*/
void CLineModelFitting::setRedshift(Float64 redshift,
                                    bool reinterpolatedContinuum) {
  m_Redshift = redshift;

  if (reinterpolatedContinuum) {
    std::shared_ptr<const CTemplate> tpl = m_tplCatalog.GetTemplateByName(
        m_tplCategoryList, m_fitContinuum_tplName);
    if (tpl->GetName() == m_fitContinuum_tplName) {
      ApplyContinuumOnGrid(tpl, redshift);
    }
  }
}

/**
 * Apply the template continuum by interpolating the grid as define in Init
 * Continuum
 */
Int32 CLineModelFitting::ApplyContinuumOnGrid(
    const std::shared_ptr<const CTemplate> &tpl, Float64 zcontinuum) {
  m_fitContinuum_tplName = tpl->GetName();
  Int32 n = tpl->GetSampleCount();

  Int32 idxDust = -1;
  if (m_fitContinuum_tplFitEbmvCoeff > 0.) {
    if (tpl->CalzettiInitFailed()) {
      THROWG(INTERNAL_ERROR, "  no calzetti calib. file in template");
    }
    idxDust = tpl->m_ismCorrectionCalzetti->GetEbmvIndex(
        m_fitContinuum_tplFitEbmvCoeff);
  }
  const CSpectrumSpectralAxis &tplSpectralAxis = tpl->GetSpectralAxis();
  TFloat64Range range(tplSpectralAxis[0], tplSpectralAxis[n - 1]);

  std::string inter_opt = "spline";
  Float64 overlapThreshold = 1., amplitude = 1.;
  std::shared_ptr<CModelSpectrumResult> spcmodel =
      m_templateFittingOperator.ComputeSpectrumModel(
          tpl, zcontinuum, m_fitContinuum_tplFitEbmvCoeff,
          m_fitContinuum_tplFitMeiksinIdx, amplitude, inter_opt,
          overlapThreshold);
  if (spcmodel == nullptr)
    THROWG(INTERNAL_ERROR, "Couldnt compute spectrum model");

  // m_observeGridContinuumFlux should be a CSpectrumFluxAxis not AxisSampleList
  m_observeGridContinuumFlux = std::move((*spcmodel).ModelFlux);

  return 0;
}

bool CLineModelFitting::SolveContinuum(
    const std::shared_ptr<const CTemplate> &tpl, const TFloat64List &redshifts,
    Float64 overlapThreshold, std::vector<CMask> maskList,
    std::string opt_interp, Int32 opt_extinction, Int32 opt_dustFit,
    Float64 &merit, Float64 &fitAmplitude, Float64 &fitAmplitudeError,
    Float64 &fitAmplitudeSigma, Float64 &FitEbmvCoeff, Int32 &fitMeiksinIdx,
    Float64 &fitDtM, Float64 &fitMtM, Float64 &fitLogprior) {
  CPriorHelper::TPriorZEList zePriorData;
  bool retGetPrior = m_fitContinuum_priorhelper->GetTplPriorData(
      tpl->GetName(), redshifts, zePriorData);
  if (retGetPrior == false) {
    THROWG(INTERNAL_ERROR,
           "    model: Failed to get prior for chi2 solvecontinuum.");
  }
  bool keepigmism = false;
  if (FitEbmvCoeff + fitMeiksinIdx != -2) {
    keepigmism = true;
  }

  // Compute merit function
  // Log.LogInfo("Solving continuum for %s at z=%.4e", tpl.GetName().c_str(),
  // redshifts[0]); CRef<CChisquareResult>  chisquareResult =
  // (CChisquareResult*)chiSquare.ExportChi2versusAZ( _spc, _tpl, lambdaRange,
  // redshifts, overlapThreshold );
  m_templateFittingOperator.SetRedshifts(redshifts);
  auto templateFittingResult =
      std::dynamic_pointer_cast<CTemplateFittingResult>(
          m_templateFittingOperator.Compute(
              tpl, overlapThreshold, maskList, opt_interp, opt_extinction,
              opt_dustFit, zePriorData, keepigmism, FitEbmvCoeff,
              fitMeiksinIdx));

  if (!templateFittingResult) {

    Log.LogError("Failed to compute chi square value");
    return false;
  } else {
    // Store results
    merit = templateFittingResult->ChiSquare[0];
    fitAmplitude = templateFittingResult->FitAmplitude[0];
    fitAmplitudeError = templateFittingResult->FitAmplitudeError[0];
    fitAmplitudeSigma = templateFittingResult->FitAmplitudeSigma[0];
    FitEbmvCoeff = templateFittingResult->FitEbmvCoeff[0];
    fitMeiksinIdx = templateFittingResult->FitMeiksinIdx[0];
    fitDtM = templateFittingResult->FitDtM[0];
    fitMtM = templateFittingResult->FitMtM[0];
    fitLogprior = templateFittingResult->LogPrior[0];
    return true;
  }
}

Int32 CLineModelFitting::LoadFitContaminantTemplate(const CTemplate &tpl) {
  Float64 fitAmplitude = NAN;
  Float64 fitAmplitudeError = NAN;
  Float64 fitAmplitudeSigma = NAN; // nan cause dependent on the the two above
  Float64 overlapThreshold = 1.0;

  bool ignoreLinesSupport = false;
  std::vector<CMask> maskList;
  if (ignoreLinesSupport) {
    maskList.resize(1);
    maskList[0] = getOutsideLinesMask();
  }
  std::string opt_interp = "lin";
  Int32 opt_extinction = 0;
  Int32 opt_dustFit = -1;

  // prepare observed spectrum without continuum
  CSpectrum::EType savetype = m_inputSpc.GetType();
  m_inputSpc.SetType(CSpectrum::nType_noContinuum);

  // prepare tpl contaminant
  const CSpectrumSpectralAxis &spcSpectralAxis =
      m_SpectrumModel.GetSpectralAxis();
  const std::string &category = "emission";
  m_tplContaminantSpcRebin =
      std::make_shared<CTemplate>("contaminantrebin", category);
  CSpectrumFluxAxis tplContaminantRebinFluxAxis(
      spcSpectralAxis.GetSamplesCount());
  CSpectrumSpectralAxis tplContaminantRebinSpcAxis(
      spcSpectralAxis.GetSamplesCount());

  // prepare INTERPOLATED contaminant
  Int32 n = tpl.GetSampleCount();
  CSpectrumFluxAxis tplFluxAxis = tpl.GetFluxAxis();
  const CSpectrumSpectralAxis &tplSpectralAxis = tpl.GetSpectralAxis();
  const Float64 *Ysrc = tplFluxAxis.GetSamples();
  const Float64 *Xsrc = tplSpectralAxis.GetSamples();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, n);
  gsl_spline_init(spline, Xsrc, Ysrc, n);
  gsl_interp_accel *accelerator = gsl_interp_accel_alloc();
  Int32 k = 0;
  Float64 x = 0.0;
  for (k = 0; k < spcSpectralAxis.GetSamplesCount(); k++) {
    x = spcSpectralAxis[k];
    tplContaminantRebinSpcAxis[k] = spcSpectralAxis[k];
    if (x < tplSpectralAxis[0] || x > tplSpectralAxis[n - 1]) {
      tplContaminantRebinFluxAxis[k] = 0.0;
    } else {
      tplContaminantRebinFluxAxis[k] = gsl_spline_eval(spline, x, accelerator);
    }
  }
  gsl_spline_free(spline);
  gsl_interp_accel_free(accelerator);

  m_tplContaminantSpcRebin->SetSpectralAndFluxAxes(
      std::move(tplContaminantRebinSpcAxis),
      std::move(tplContaminantRebinFluxAxis));

  //*
  // Fit contaminant template AMPLITUDE
  m_templateFittingOperator.SetRedshifts(
      TFloat64List(1, 0.0)); // fitting an already redshifted model
  auto templateFittingResult =
      std::dynamic_pointer_cast<CTemplateFittingResult>(
          m_templateFittingOperator.Compute(
              m_tplContaminantSpcRebin, overlapThreshold, maskList, opt_interp,
              opt_extinction, opt_dustFit));

  if (!templateFittingResult) {

    Log.LogError("    model: contaminant - Failed to compute chi square value");
    return false;
  } else {
    // Extract amplitude and amplitude error
    fitAmplitude = templateFittingResult->FitAmplitude[0];
    fitAmplitudeError = templateFittingResult->FitAmplitudeError[0];
    fitAmplitudeSigma = templateFittingResult->FitAmplitudeSigma[0];
  }
  Log.LogInfo("    model: contaminant raw fit amplitude is %f.", fitAmplitude);
  //*/
  if (fitAmplitude > 1.0) {
    Log.LogInfo(
        "    model: contaminant raw fit amplitude was %f. NOW SET TO 1.0",
        fitAmplitude);
    fitAmplitude = 1.0;
  }
  //*/
  // fitAmplitude = 1.0;
  Log.LogInfo("    model: contaminant fit amplitude = %f", fitAmplitude);
  for (k = 0; k < spcSpectralAxis.GetSamplesCount(); k++) {
    tplContaminantRebinFluxAxis[k] =
        fitAmplitude * tplContaminantRebinFluxAxis[k];
  }

  /*//debug:
  FILE* f = fopen( "contaminantFittedInterpolated.txt", "w+" );
  for(k=0; k<m_SpcFluxAxis.GetSamplesCount(); k++)
  {
      fprintf( f, "%f\t%e\n", spcSpectralAxis[k],
  tplContaminantRebinFluxAxis[k]);
  }
  fclose( f );
  //*/

  m_inputSpc.SetType(savetype);
  if (m_ContinuumComponent == "nocontinuum" ||
      m_ContinuumComponent == "fromspectrum" ||
      m_ContinuumComponent == "tplfit" ||
      m_ContinuumComponent == "tplfitauto") {
    for (k = 0; k < m_SpcFluxAxis.GetSamplesCount(); k++) {
      m_SpcFluxAxis[k] =
          m_inputSpc.GetFluxAxis()[k] - tplContaminantRebinFluxAxis[k];
    }
  }
  return 1;
}

std::shared_ptr<CModelSpectrumResult>
CLineModelFitting::GetContaminantSpectrumResult() const {
  std::shared_ptr<CModelSpectrumResult> resultcont =
      std::make_shared<CModelSpectrumResult>(*m_tplContaminantSpcRebin);
  return resultcont;
}

const std::string &CLineModelFitting::getFitContinuum_tplName() const {
  return m_fitContinuum_tplName;
}

Float64 CLineModelFitting::getFitContinuum_tplAmplitude() const {
  return m_fitContinuum_tplFitAmplitude;
}

Float64 CLineModelFitting::getFitContinuum_tplAmplitudeError() const {
  return m_fitContinuum_tplFitAmplitudeError;
}

// This SNR estimate maybe needs to use observed spectrum with lines removed ?
Float64 CLineModelFitting::getFitContinuum_snr() const {
  Float64 snr = -1.; // shoudnt be NAN here?
  if (m_fitContinuum_tplFitMtM > 0.) {
    snr = m_fitContinuum_tplFitDtM / std::sqrt(m_fitContinuum_tplFitMtM);
  }
  return snr;
}

Float64 CLineModelFitting::getFitContinuum_tplMerit() const {
  return m_fitContinuum_tplFitMerit;
}

Float64 CLineModelFitting::getFitContinuum_tplMeritPhot() const {
  return m_fitContinuum_tplFitMerit_phot;
}

Float64 CLineModelFitting::getFitContinuum_tplIsmEbmvCoeff() const {
  return m_fitContinuum_tplFitEbmvCoeff;
}

Float64 CLineModelFitting::getFitContinuum_tplIgmMeiksinIdx() const {
  return m_fitContinuum_tplFitMeiksinIdx;
}

void CLineModelFitting::SetContinuumComponent(std::string component) {
  m_ContinuumComponent = component;

  const Int32 spectrumSampleCount = m_inputSpc.GetSampleCount();

  if (m_ContinuumComponent == "nocontinuum") {
    m_fitContinuum_tplName = "nocontinuum"; // to keep track in resultstore
    // the continuum is set to zero and the observed spectrum is the spectrum
    // without continuum
    m_spcFluxAxisNoContinuum =
        m_inputSpc
            .GetRawFluxAxis(); // m_inputSpc.GetWithoutContinuumFluxAxis();
    m_SpcFluxAxis = m_spcFluxAxisNoContinuum;
    m_SpectrumModel.SetFluxAxis(CSpectrumFluxAxis(spectrumSampleCount));
  }
  if (m_ContinuumComponent == "fromspectrum") {
    m_fitContinuum_tplName = "fromspectrum"; // to keep track in resultstore
    // the continuum is set to the spectrum continuum and the observed spectrum
    // is the raw spectrum
    m_spcFluxAxisNoContinuum = m_inputSpc.GetWithoutContinuumFluxAxis();
    m_SpcFluxAxis = m_inputSpc.GetRawFluxAxis();
    m_ContinuumFluxAxis = m_inputSpc.GetContinuumFluxAxis();
    m_SpectrumModel.SetFluxAxis(m_ContinuumFluxAxis);
  }
  if (m_ContinuumComponent == "tplfit" ||
      m_ContinuumComponent == "tplfitauto") {
    // the continuum is set to zero and the observed spectrum is the raw
    // spectrum
    m_SpcFluxAxis = m_inputSpc.GetRawFluxAxis();
    m_SpectrumModel.SetFluxAxis(CSpectrumFluxAxis(spectrumSampleCount));
  }
}

Int32 CLineModelFitting::SetFitContinuum_FitStore(
    const std::shared_ptr<const CTemplatesFitStore> &fitStore) {
  m_fitContinuum_option = 1; // enable use of the fit store
  Log.LogDetail("Elementlist: enabling fitContinuum store.");
  m_fitContinuum_tplfitStore = fitStore;
  return 1;
}

const std::shared_ptr<const CTemplatesFitStore> &
CLineModelFitting::GetFitContinuum_FitStore() const {
  return m_fitContinuum_tplfitStore;
}

Int32 CLineModelFitting::SetFitContinuum_PriorHelper(
    const std::shared_ptr<const CPriorHelper> &priorhelper) {
  m_fitContinuum_priorhelper = priorhelper;
  return 1;
}

void CLineModelFitting::SetFitContinuum_Option(Int32 opt) {
  m_fitContinuum_option = opt;
}

void CLineModelFitting::SetFitContinuum_SNRMax(Float64 snr_max) {
  m_fitContinuum_tplFitSNRMax = snr_max;
}

Int32 CLineModelFitting::GetFitContinuum_Option() const {
  return m_fitContinuum_option;
}

void CLineModelFitting::SetFitContinuum_FitValues(
    std::string tplfit_name, Float64 tplfit_amp, Float64 tplfit_amperr,
    Float64 tplfit_chi2, Float64 tplfit_chi2_phot, Float64 tplfit_ebmv,
    Int32 tplfit_meiksinidx, Float64 tplfit_continuumredshift,
    Float64 tplfit_dtm, Float64 tplfit_mtm, Float64 tplfit_logprior,
    const TFloat64List &polyCoeffs) {
  m_fitContinuum_tplName = tplfit_name;
  m_fitContinuum_tplFitAmplitude = tplfit_amp;
  m_fitContinuum_tplFitAmplitudeError = tplfit_amperr;
  m_fitContinuum_tplFitMerit = tplfit_chi2;
  m_fitContinuum_tplFitMerit_phot = tplfit_chi2_phot;
  m_fitContinuum_tplFitEbmvCoeff = tplfit_ebmv;
  m_fitContinuum_tplFitMeiksinIdx = tplfit_meiksinidx;
  m_fitContinuum_tplFitRedshift = tplfit_continuumredshift;

  m_fitContinuum_tplFitDtM = tplfit_dtm;
  m_fitContinuum_tplFitMtM = tplfit_mtm;
  m_fitContinuum_tplFitLogprior = tplfit_logprior;
  m_fitContinuum_tplFitPolyCoeffs = polyCoeffs;
}

/**
 * \brief This function prepares the continuum for use in the fit with the line
 *elements. Rebin with PFG buffer Find and apply amplitude factor from
 *previously fitted tpl
 **/
void CLineModelFitting::PrepareContinuum() {
  const CSpectrumSpectralAxis &targetSpectralAxis =
      m_SpectrumModel.GetSpectralAxis();
  TAxisSampleList &Yrebin = m_ContinuumFluxAxis.GetSamplesVector();

  if (m_ContinuumComponent == "nocontinuum")
    return;

  if (m_ContinuumComponent == "tplfit" ||
      m_ContinuumComponent == "tplfitauto" || m_fittingmethod == "lmfit") {
    m_observeGridContinuumFlux.resize(m_SpectrumModel.GetSampleCount());
  } else {
    Yrebin = m_inputSpc.GetContinuumFluxAxis().GetSamplesVector();
  }
}

const std::string &CLineModelFitting::getTplshape_bestTplName() const {
  return m_tplshapeBestTplName;
}
Float64 CLineModelFitting::getTplshape_bestTplIsmCoeff() const {
  return m_tplshapeBestTplIsmCoeff;
}

Float64 CLineModelFitting::getTplshape_bestAmplitude() const {
  return m_tplshapeBestTplAmplitude;
}

Float64 CLineModelFitting::getTplshape_bestDtm() const {
  return m_tplshapeBestTplDtm;
}

Float64 CLineModelFitting::getTplshape_bestMtm() const {
  return m_tplshapeBestTplMtm;
}

Int32 CLineModelFitting::getTplshape_count() const {
  if (m_rigidity != "tplshape") {
    return 0;
  }
  return m_CatalogTplShape.GetCatalogsCount();
}

const TFloat64List &CLineModelFitting::getTplshape_priors() {
  if (m_rigidity != "tplshape") {
    static TFloat64List dumb;
    return dumb;
  }
  return m_CatalogTplShape.getCatalogsPriors();
}

const TFloat64List &CLineModelFitting::GetChisquareTplshape() const {
  return m_ChisquareTplshape;
}

/**
 * @brief CLineModelFitting::GetPriorLinesTplshape
 * WARNING: as stated in fit(), the prior is valid in this code structure only
 * for tplratio with only 1 element containing the EL component, hence using
 * idx=0 here
 * @return
 */
TFloat64List CLineModelFitting::GetPriorLinesTplshape() const {
  TFloat64List plinestplshape;
  Int32 eltIdx = 0;
  for (Int32 ktpl = 0; ktpl < m_LinesLogPriorTplshape.size(); ktpl++) {
    plinestplshape.push_back(m_LinesLogPriorTplshape[ktpl][eltIdx]);
  }
  return plinestplshape;
}

const TFloat64List &CLineModelFitting::GetScaleMargTplshape() const {
  return m_ScaleMargCorrTplshape;
}

const TBoolList &CLineModelFitting::GetStrongELPresentTplshape() const {
  return m_StrongELPresentTplshape;
}

const TBoolList &CLineModelFitting::getHaELPresentTplshape() const {
  return m_StrongHalphaELPresentTplshape;
}

const TInt32List &CLineModelFitting::GetNLinesAboveSNRTplshape() const {
  return m_NLinesAboveSNRTplshape;
}

bool CLineModelFitting::initModelAtZ(
    Float64 redshift, const CSpectrumSpectralAxis &spectralAxis) {
  m_Redshift = redshift;

  // prepare the elements support
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    m_Elements[iElts]->resetAsymfitParams();
    m_Elements[iElts]->prepareSupport(spectralAxis, redshift, m_lambdaRange);
  }

  return true;
}

bool CLineModelFitting::setTplshapeModel(Int32 itplshape,
                                         bool enableSetVelocity) {
  SetMultilineNominalAmplitudesFast(itplshape);

  if (enableSetVelocity) {
    // Set the velocities from templates: todo auto switch when velfit is ON
    m_CatalogTplShape.GetCatalogVelocities(itplshape, m_velocityEmission,
                                           m_velocityAbsorption);
  }

  Log.LogDebug("    model : setTplshapeModel, loaded: %d = %s", itplshape,
               m_CatalogTplShape.GetCatalogName(itplshape).c_str());
  return true;
}

Int32 CLineModelFitting::SetTplshape_PriorHelper(
    const std::shared_ptr<const CPriorHelper> &priorhelper) {
  m_tplshape_priorhelper = priorhelper;
  return 1;
}

bool CLineModelFitting::setTplshapeAmplitude(const TFloat64List &ampsElts,
                                             const TFloat64List &errorsElts) {
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    m_Elements[iElts]->SetFittedAmplitude(ampsElts[iElts], errorsElts[iElts]);
  }
  return true;
}

bool CLineModelFitting::initDtd() {
  m_dTransposeDLambdaRange = m_lambdaRange;
  if (m_ContinuumComponent == "tplfit" ||
      m_ContinuumComponent == "tplfitauto") {
    m_dTransposeD = EstimateDTransposeD("raw");
  } else {
    m_dTransposeD = EstimateDTransposeD("nocontinuum");
  }
  m_likelihood_cstLog = EstimateLikelihoodCstLog();
  return true;
}

/**
 * \brief Prepares the context and fits the Linemodel to the spectrum, returning
 *the merit of the fit. Prepare the continuum. Initialize the model spectrum.
 * Prepare the elements.
 * Fit the amplitudes of each element independently.
 * Fit the amplitude of all elements together with iterative solver: Nelder Mead
 *Simplex. Fit the amplitude of all elements together with linear solver:
 *gsl_multifit_wlinear. Fit the amplitudes of each element independently, unless
 *there is overlap. Apply a continuum iterative re-estimation with lines removed
 *from the initial spectrum. Apply rules. Create spectrum model. Return merit.
 **/
Float64 CLineModelFitting::fit(Float64 redshift,
                               CLineModelSolution &modelSolution,
                               CContinuumModelSolution &continuumModelSolution,
                               Int32 contreest_iterations, bool enableLogging) {
  // initialize the model spectrum
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();

  if (m_fittingmethod == "lmfit") {
    SetVelocityEmission(m_velocityEmissionInit);
    SetVelocityAbsorption(m_velocityAbsorptionInit);
  }

  initModelAtZ(redshift, spectralAxis);

  if (!(m_dTransposeDLambdaRange.GetBegin() == m_lambdaRange.GetBegin() &&
        m_dTransposeDLambdaRange.GetEnd() == m_lambdaRange.GetEnd())) {
    initDtd();
  }

  Int32 nfitting = 1; // multiple fitting steps for rigidity=tplshape/tplratio
  std::vector<CPriorHelper::SPriorTZE> logPriorDataTplShape;
  if (m_rigidity == "tplshape") {
    nfitting = m_CatalogTplShape.GetCatalogsCount();

    if (!m_tplshape_priorhelper->mInitFailed) {
      // prior initilization for tplshape EL only
      if (m_Elements.size() > 1) {
        THROWG(INTERNAL_ERROR, "model: Unable to use tplshape line priors "
                               "with nElts>1 for now");
        // NB: this could be done if the EL element idx in searched (see later
        // in the ifitting loop, UV Abs lines would be not affected by priors
        // then)
      }
      for (Int32 ifitting = 0; ifitting < nfitting; ifitting++) {
        // prepare the lines prior data
        Int32 ebvfilter = m_CatalogTplShape.GetIsmIndex(ifitting);
        CPriorHelper::SPriorTZE logPriorData;
        std::string tplrationame = m_CatalogTplShape.GetCatalogName(ifitting);
        bool retGetPrior = m_tplshape_priorhelper->GetTZEPriorData(
            tplrationame, ebvfilter, redshift, logPriorData);
        if (retGetPrior == false) {
          THROWG(INTERNAL_ERROR,
                 "model: Failed to get prior for chi2 solvecontinuum.");
        } else {
          logPriorDataTplShape.push_back(logPriorData);
        }
      }
    }
  }
  Int32 savedIdxFitted = -1; // for rigidity=tplshape

  Int32 ncontinuumfitting = 1;
  Int32 savedIdxContinuumFitted = -1; // for continuum tplfit
  if ((m_ContinuumComponent == "tplfit" ||
       m_ContinuumComponent == "tplfitauto") &&
      !m_forcedisableMultipleContinuumfit) {
    ncontinuumfitting = m_opt_fitcontinuum_maxCount;
  }

  Float64 merit = INFINITY; // initializing 'on the fly' best-merit
  Float64 meritprior = 0.0; // initializing 'on the fly' best-merit-prior
  TFloat64List meritTplratio(
      nfitting, INFINITY); // initializing 'on the fly' best-merit per tplratio
  TFloat64List meritPriorTplratio(
      nfitting, 0.0); // initializing 'on the fly' best-merit-prior per tplratio
  for (Int32 icontfitting = 0; icontfitting < ncontinuumfitting;
       icontfitting++) {
    if (m_ContinuumComponent != "nocontinuum") {
      PrepareContinuum();
    }
    if (m_ContinuumComponent == "tplfit" ||
        m_ContinuumComponent ==
            "tplfitauto") // the support has to be already computed when
                          // LoadFitContinuum() is called
    {
      Int32 autoselect = (m_ContinuumComponent == "tplfitauto");
      LoadFitContinuum(icontfitting, autoselect);
    }
    if (m_ContinuumComponent != "nocontinuum") {
      CSpectrumFluxAxis modelFluxAxis = m_ContinuumFluxAxis;
      const auto &SpcFluxAxis = m_SpcFluxAxis;
      const auto &ContinuumFluxAxis = m_ContinuumFluxAxis;
      m_SpectrumModel.SetFluxAxis(std::move(modelFluxAxis));
      for (Int32 i = 0; i < m_SpectrumModel.GetSampleCount(); i++) {
        m_spcFluxAxisNoContinuum[i] = SpcFluxAxis[i] - ContinuumFluxAxis[i];
      }
    }

    if (m_enableAmplitudeOffsets) {
      m_Elements.prepareAmplitudeOffset();
    }

    for (Int32 ifitting = 0; ifitting < nfitting; ifitting++) {
      if (m_rigidity == "tplshape") {
        if (m_forcedisableTplratioISMfit && ifitting > 0 &&
            m_CatalogTplShape.GetIsmCoeff(ifitting) > 0) {
          // copy the values for ebmv=ebmv_fixed (=0) here
          m_ChisquareTplshape[ifitting] = m_ChisquareTplshape[ifitting - 1];
          m_ScaleMargCorrTplshape[ifitting] =
              m_ScaleMargCorrTplshape[ifitting - 1];
          m_StrongELPresentTplshape[ifitting] =
              m_StrongELPresentTplshape[ifitting - 1];
          m_StrongHalphaELPresentTplshape[ifitting] =
              m_StrongHalphaELPresentTplshape[ifitting - 1];
          m_NLinesAboveSNRTplshape[ifitting] =
              m_NLinesAboveSNRTplshape[ifitting - 1];
          meritTplratio[ifitting] = meritTplratio[ifitting - 1];
          meritPriorTplratio[ifitting] = meritPriorTplratio[ifitting - 1];
          for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
            m_FittedAmpTplshape[ifitting][iElts] =
                m_FittedAmpTplshape[ifitting - 1][iElts];
            m_FittedErrorTplshape[ifitting][iElts] =
                m_FittedErrorTplshape[ifitting - 1][iElts];
            m_DtmTplshape[ifitting][iElts] = m_DtmTplshape[ifitting - 1][iElts];
            m_MtmTplshape[ifitting][iElts] = m_MtmTplshape[ifitting - 1][iElts];

            m_LyaAsymCoeffTplshape[ifitting][iElts] =
                m_LyaAsymCoeffTplshape[ifitting - 1][iElts];
            m_LyaWidthCoeffTplshape[ifitting][iElts] =
                m_LyaWidthCoeffTplshape[ifitting - 1][iElts];
            m_LyaDeltaCoeffTplshape[ifitting][iElts] =
                m_LyaDeltaCoeffTplshape[ifitting - 1][iElts];
            m_LyaIgmIdxTplshape[ifitting][iElts] =
                m_LyaIgmIdxTplshape[ifitting - 1][iElts];
            m_LinesLogPriorTplshape[ifitting][iElts] =
                m_LinesLogPriorTplshape[ifitting - 1][iElts];
          }
          continue;
        }
        setTplshapeModel(ifitting, false);
        // prepare the Lya width and asym coefficients if the asymfit profile
        // option is met INFO: tpl-shape are often ASYMFIXED in the tplshape
        // catalog files, for the lyaE profile, as of 2016-01-11 INFO: tplshape
        // can override the lyafitting, see m_opt_lya_forcefit
        setLyaProfile(m_Redshift,
                      m_CatalogTplShape.GetCatalog(ifitting).GetList(), true);
      } else {
        if (!m_forceDisableLyaFitting) // for asymFit, disable fitting
        {
          // prepare the Lya width and asym coefficients if the asymfit profile
          // option is met
          setLyaProfile(redshift, m_RestLineList);
        }
      }
      // generate random amplitudes
      if (m_fittingmethod == "random") {
        srand(time(0));
        Float64 randNumFloat = (Float64)rand() / (Float64)(RAND_MAX);

        Float64 coeffAmpEmission = pow(10.0, randNumFloat * 3.0 - 1.0);
        randNumFloat = (Float64)rand() / (Float64)(RAND_MAX);
        Float64 coeffAmpAbsorption = pow(10.0, randNumFloat * 1.0 - 1.0);
        Log.LogInfo("\nLineModel simulation: coeffAmpEmission = %.2f",
                    coeffAmpEmission);
        Log.LogInfo("LineModel simulation: coeffAmpAbsorption = %.2f",
                    coeffAmpAbsorption);
        // fit the model amplitudes individually
        for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
          Float64 meanContinuum = getContinuumMeanUnderElement(iElts);
          Float64 err = 1e-22;
          Float64 amax = meanContinuum;
          if (m_Elements[iElts]->m_Lines[0].GetType() ==
              CLine::nType_Absorption) {
            amax = meanContinuum * 0.5 * coeffAmpAbsorption;
          } else {
            amax = meanContinuum * coeffAmpEmission;
          }
          randNumFloat = (Float64)rand() / (Float64)(RAND_MAX);
          Float64 a = randNumFloat * amax;
          if (a < 0.0) {
            a = 0.0;
          }
          // get the max nominal amplitude
          Int32 nLines = m_Elements[iElts]->GetSize();
          Float64 maxNominalAmp = -1.0;
          for (Int32 j = 0; j < nLines; j++) {
            if (maxNominalAmp < m_Elements[iElts]->GetNominalAmplitude(j)) {
              maxNominalAmp = m_Elements[iElts]->GetNominalAmplitude(j);
            }
          }

          m_Elements.SetElementAmplitude(iElts, a / maxNominalAmp, err);
        }
      }

      // fit the amplitudes of each element independently
      if (m_fittingmethod == "individual") {
        // fit the model amplitudes individually
        for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
          m_Elements[iElts]->fitAmplitude(spectralAxis,
                                          m_spcFluxAxisNoContinuum,
                                          m_ContinuumFluxAxis, redshift);

          // optionnally log some fitting details
          if (0) {
            Log.LogDebug("    model: elt #%d individual fit", iElts);
            Float64 dbg_amp = m_Elements[iElts]->GetElementAmplitude();
            Log.LogDebug("    model:     fitted elt amp = %f", dbg_amp);
            Int32 nSubE = m_Elements[iElts]->GetSize();
            for (Int32 iSubElts = 0; iSubElts < nSubE; iSubElts++) {
              Log.LogDebug("    model:     sub #%d - fitted amp = %f", iSubElts,
                           m_Elements[iElts]->GetFittedAmplitude(iSubElts));
              Log.LogDebug("    model:     sub #%d - outside range = %d",
                           iSubElts,
                           m_Elements[iElts]->IsOutsideLambdaRange(iSubElts));
            }
            Log.LogDebug("    model:     dtm = %e",
                         m_Elements[iElts]->GetSumCross());
            Log.LogDebug("    model:     dtd = %e",
                         m_Elements[iElts]->GetSumGauss());
          }
        }
      }

      // fit the amplitude of all elements together (but Emission or Absorption
      // separately) with iterative   solver: lmfit
      if (m_fittingmethod == "lmfit") {
        // Log.LogInfo( "Linemodel: Fitting method lmfit, z= %f", m_Redshift);
        SetVelocityEmission(m_velocityEmissionInit);
        SetVelocityAbsorption(m_velocityAbsorptionInit);
        m_Redshift = redshift;
        Float64 bestMerit = INFINITY;
        std::shared_ptr<CLmfitController>
            bestController; // at each run of lmfit we create a controller, we
                            // keep the best evaluation
        std::vector<CLmfitController *> controllers = createLmfitControllers();
        if (controllers.size() == 0) {
          Log.LogError("LMfit : No Controller created");
        }
        for (Int32 i = 0; i < controllers.size(); i++) {
          CLmfitController *controller = controllers[i];
          // Log.LogInfo("Continuum Template use : %s",
          // controller->getTemplate()->GetName().c_str());
          if (!controller->isNoContinuum() &&
              !controller->isContinuumLoaded()) {
            LoadFitContinuumOneTemplate(controller->getTemplate());
          }
          // adding element base on configuration
          TInt32List validEltsIdx = m_Elements.GetModelValidElementsIndexes();
          for (Int32 iElt = 0; iElt < validEltsIdx.size(); iElt++) {
            // if(controller->isLineTypeVelocityFitted(m_RestLineList[m_Elements[validEltsIdx[iElt]]->m_LineCatalogIndexes[0]].GetType()))
            // {
            controller->addElement(validEltsIdx[iElt]);
            m_Elements[validEltsIdx[iElt]]->fitAmplitude(
                spectralAxis, m_spcFluxAxisNoContinuum, m_ContinuumFluxAxis,
                redshift);
            // }
          }

          controller->calculatedIndices();
          bool NegVal = true;
          while (NegVal) {
            Int32 retVal = fitAmplitudesLmfit(m_SpcFluxAxis, controller);
            // Log.LogInfo("LineModel LMfit: retVal = %d", retVal);
            NegVal = controller->removeNegAmpLine();
            // Log.LogInfo("LineModel LMfit: NegVal = %s", (NegVal?"t":"f") );
            // Log.LogInfo("LineModel LMfit: merit = %.0f ; bestMerit : %.0f ",
            // controller->getMerit(), bestMerit );
            if (retVal == 0 && !NegVal && controller->getMerit() < bestMerit) {

              bestMerit = controller->getMerit();
              bestController = std::shared_ptr<CLmfitController>(controller);
              // Log.LogInfo("Pointer best Controller setted");
            }
            // NegVal = false;//no while
          }
        }

        // We set the best model found
        if (bestController) {
          if (bestController->isContinuumFitted()) {
            if (bestController->isRedshiftFitted()) {
              m_Redshift = bestController->getRedshift();
            }
            std::shared_ptr<const CTemplate> tpl =
                bestController->getTemplate();

            ApplyContinuumOnGrid(tpl, m_Redshift);
            TFloat64List polyCoeffs;
            setFitContinuum_tplAmplitude(bestController->getContinuumAmp(),
                                         bestController->getContinuumAmpErr(),
                                         polyCoeffs);
            // setFitContinuum_tplAmplitude( bestController->getContinuumAmp());
          } else if (bestController->isRedshiftFitted()) {
            setRedshift(bestController->getRedshift(),
                        bestController->isContinuumLoaded());
          }
          TInt32List bestFilteredEltsIdx = bestController->getFilteredIdx();

          for (Int32 iElt = 0; iElt < bestFilteredEltsIdx.size(); iElt++) {
            Float64 amp = bestController->getLineAmp(iElt);
            Float64 ampErr = bestController->getLineAmpErr(iElt);
            // Log.LogInfo("LMfit : amp [%d] = %.1f, z = %f", iElt, amp,
            // m_Redshift);
            m_Elements.SetElementAmplitude(bestFilteredEltsIdx[iElt], amp,
                                           ampErr);
          }
          if (bestController->isEmissionVelocityFitted()) {
            SetVelocityEmission(bestController->getEmissionVelocity());
            // Log.LogInfo("LMfit : VelEmission = %.1f, z = %f",
            // bestController->getEmissionVelocity(), m_Redshift);
          }
          if (bestController->isAbsorptionVelocityFitted()) {
            SetVelocityAbsorption(bestController->getAbsorptionVelocity());
            // Log.LogInfo("LMfit : VelAbsorption= %.1f, z = %f",
            // bestController->getAbsorptionVelocity(), m_Redshift);
          }
          // prepare the elements support
          for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
            m_Elements[iElts]->prepareSupport(spectralAxis, redshift,
                                              m_lambdaRange);
          }
          refreshModelInitAllGrid();
          // modelSolution = GetModelSolution();
          while (!controllers.empty()) {
            CLmfitController *controller = controllers.back();
            controllers.pop_back();
            delete controller;
          }

        } else {
          Log.LogError("LineModel LMfit: not able to fit values at z %f",
                       m_Redshift);
          // continue;
          while (!controllers.empty()) {
            CLmfitController *controller = controllers.back();
            controllers.pop_back();
            delete controller;
          }
          return INFINITY;
        }
      }

      // fit the amplitude of all elements together with linear solver:
      // gsl_multifit_wlinear
      if (m_fittingmethod == "svd") {
        TInt32List validEltsIdx = m_Elements.GetModelValidElementsIndexes();
        TFloat64List ampsfitted;
        TFloat64List errorsfitted;
        fitAmplitudesLinSolveAndLambdaOffset(
            validEltsIdx, spectralAxis, m_spcFluxAxisNoContinuum,
            m_ContinuumFluxAxis, ampsfitted, errorsfitted,
            m_enableLambdaOffsetsFit);
      }

      // fit the amplitude of all elements AND continuum amplitude together with
      // linear solver: gsl_multifit_wlinear
      if (m_fittingmethod == "svdlc" || m_fittingmethod == "svdlcp2") {
        // 1. fit only the current continuum
        // prepare continuum on the observed grid
        Log.LogDebug("    model: fitting svdlc, with continuum-tpl=%s",
                     m_fitContinuum_tplName.c_str());

        // re-interpolate the continuum on the grid
        std::shared_ptr<const CTemplate> tpl = m_tplCatalog.GetTemplateByName(
            m_tplCategoryList, m_fitContinuum_tplName);
        ApplyContinuumOnGrid(tpl, m_fitContinuum_tplFitRedshift);

        m_fitContinuum_tplFitAmplitude = 1.0;
        m_fitContinuum_tplFitAmplitudeError = 1.0;
        TFloat64List polyCoeffs_unused;
        setFitContinuum_tplAmplitude(m_fitContinuum_tplFitAmplitude,
                                     m_fitContinuum_tplFitAmplitudeError,
                                     polyCoeffs_unused);

        TInt32List validEltsIdx = m_Elements.GetModelValidElementsIndexes();
        TFloat64List ampsfitted;
        TFloat64List errorsfitted;
        Float64 chi2_cl = INFINITY;
        Int32 fitc_polyOrder = -1;
        if (m_fittingmethod == "svdlcp2") {
          fitc_polyOrder = 2;
        }
        fitAmplitudesLinesAndContinuumLinSolve(
            validEltsIdx, spectralAxis, m_SpcFluxAxis, m_ContinuumFluxAxis,
            ampsfitted, errorsfitted, chi2_cl, fitc_polyOrder);
        Log.LogDebug("    model: fitting svdlc done");

        m_fitContinuum_tplFitAmplitude = ampsfitted[validEltsIdx.size()];
        TFloat64List polyCoeffs;
        for (Int32 kpoly = validEltsIdx.size() + 1; kpoly < ampsfitted.size();
             kpoly++) {
          polyCoeffs.push_back(ampsfitted[kpoly]);
        }
        setFitContinuum_tplAmplitude(m_fitContinuum_tplFitAmplitude,
                                     m_fitContinuum_tplFitAmplitudeError,
                                     polyCoeffs);
        CSpectrumFluxAxis modelFluxAxis = m_ContinuumFluxAxis;
        const auto &SpcFluxAxis = m_SpcFluxAxis;
        const auto &ContinuumFluxAxis = m_ContinuumFluxAxis;
        m_SpectrumModel.SetFluxAxis(std::move(modelFluxAxis));
        for (Int32 i = 0; i < m_SpectrumModel.GetSampleCount(); i++) {
          m_spcFluxAxisNoContinuum[i] = SpcFluxAxis[i] - ContinuumFluxAxis[i];
        }

        // 2. loop on the continuum templates from catalog
        /*
        bool stop_loop_tpl = false;
        Float64 bestChi2 = INFINITY;
        Float64 bestAmp = 1.0;
        std::string bestTplName = "undefined";
        for( Int32 i=0; i<m_tplCategoryList.size(); i++ )
        {
            if(stop_loop_tpl)
            {
                break;
            }
            std::string category = m_tplCategoryList[i];

            for( Int32 j=0; j<m_tplCatalog.GetTemplateCount( category ); j++ )
            {
                if(stop_loop_tpl)
                {
                    break;
                }
                const CTemplate& tpl = m_tplCatalog.GetTemplate( category, j );

                //prepare continuum on the observed grid
                m_fitContinuum_tplFitEbmvCoeff = 0;
                m_fitContinuum_tplFitMeiksinIdx = 0;
                m_fitContinuum_tplName = tpl.GetName();
                Log.LogDebug( "    model: fitting svdlc, with continuum-tpl=%s",
        m_fitContinuum_tplName.c_str()); setRedshift(m_Redshift, true);
                m_fitContinuum_tplFitAmplitude = 1.0;
                setFitContinuum_tplAmplitude(m_fitContinuum_tplFitAmplitude);

                TInt32List validEltsIdx = GetModelValidElementsIndexes();
                TFloat64List ampsfitted;
                TFloat64List errorsfitted;
                Float64 chi2_cl = INFINITY;
                fitAmplitudesLinesAndContinuumLinSolve(validEltsIdx,
                                                       lambdaRange,
                                                       spectralAxis,
                                                       m_SpcFluxAxis,
                                                       m_ContinuumFluxAxis,
                                                       ampsfitted,
                                                       errorsfitted,
                                                       chi2_cl);
                Log.LogDebug( "    model: fitting svdlc done");

                if(chi2_cl<bestChi2)
                {
                    m_fitContinuum_tplFitAmplitude =
        ampsfitted[ampsfitted.size()-1]; if(1) //not needed here ?
                    {
                        Log.LogDebug( "    model: fitting svdlc, amp found=%e",
        m_fitContinuum_tplFitAmplitude); Log.LogDebug( "    model: fitting
        svdlc, chi2 found=%e", chi2_cl); m_fitContinuum_tplFitMerit = -1;
                        //m_fitContinuum_tplFitEbmvCoeff = -1;
                        //m_fitContinuum_tplFitMeiksinIdx = -1;
                        m_fitContinuum_tplFitDtM = -1;
                        m_fitContinuum_tplFitMtM = -1;
                        setFitContinuum_tplAmplitude(m_fitContinuum_tplFitAmplitude);
                        //PrepareContinuum(m_Redshift);

                        for(Int32 i=0; i<modelFluxAxis.GetSamplesCount(); i++)
                        {
                            modelFluxAxis[i] = m_ContinuumFluxAxis[i];
                            m_spcFluxAxisNoContinuum[i] =
        m_SpcFluxAxis[i]-m_ContinuumFluxAxis[i];
                        }
                        //refreshModel();
                    }
                    bestChi2 = chi2_cl;
                    bestAmp = m_fitContinuum_tplFitAmplitude;
                    bestTplName = m_fitContinuum_tplName;
                }
                //stop_loop_tpl = true; //only process one template
            }
        }

        m_fitContinuum_tplName = bestTplName;
        setRedshift(m_Redshift, true);
        m_fitContinuum_tplFitAmplitude = bestAmp;
        setFitContinuum_tplAmplitude(m_fitContinuum_tplFitAmplitude);
        //PrepareContinuum(m_Redshift);
        for(Int32 i=0; i<modelFluxAxis.GetSamplesCount(); i++)
        {
            modelFluxAxis[i] = m_ContinuumFluxAxis[i];
            m_spcFluxAxisNoContinuum[i] =
        m_SpcFluxAxis[i]-m_ContinuumFluxAxis[i];
        }
        //*/
      }

      // fit the amplitudes of each element independently, unless there is
      // overlap
      if (m_fittingmethod == "hybrid") {
        fitAmplitudesHybrid(spectralAxis, m_spcFluxAxisNoContinuum,
                            m_ContinuumFluxAxis, redshift);

        // apply a continuum iterative re-estimation with lines removed from the
        // initial spectrum
        Int32 nIt = contreest_iterations;
        Int32 it = 0;
        while (it < nIt) {
          applyRules();

          //*
          // iterative continuum estimation :: RAW SLOW METHOD
          refreshModel();
          Float64 enhanceLines = 0;
          //*
          if (nIt > 2 * it && nIt > 3.0 && it <= 3) {
            enhanceLines = 2.0 - ((Float64)it * 0.33);
          }

          //*/
          /*
          if(it==0 && nIt>1){
              enhanceLines = 1.5;
          }
          */
          EstimateSpectrumContinuum(enhanceLines);
          CSpectrumFluxAxis modelFluxAxis = m_ContinuumFluxAxis;
          const auto &SpcFluxAxis = m_SpcFluxAxis;
          const auto &ContinuumFluxAxis = m_ContinuumFluxAxis;
          m_SpectrumModel.SetFluxAxis(std::move(modelFluxAxis));
          for (Int32 i = 0; i < m_SpectrumModel.GetSampleCount(); i++) {
            m_spcFluxAxisNoContinuum[i] = SpcFluxAxis[i] - ContinuumFluxAxis[i];
          }
          //*/

          /*
      //iterative continuum estimation approx: APPROX. METHOD
      TInt32List validEltsIdx = GetModelValidElementsIndexes();
      TInt32List refreshIdxs = ReestimateContinuumApprox(validEltsIdx);
      refreshModelAfterContReestimation(refreshIdxs, modelFluxAxis,
      spcFluxAxisNoContinuum);
      //*/

          /*
      //iterative continuum estimation approx: FAST METHOD
      TInt32List validEltsIdx = GetModelValidElementsIndexes();
      TInt32List highSNRvalidEltsIdx;
      Float64 snrthres = 5.0;
      for( Int32 i=0; i<validEltsIdx.size(); i++ )
      {
          Int32 eltIdx = validEltsIdx[i];
          bool isSnrHigh = false;
          Int32 nlines = m_Elements[eltIdx]->GetSize();
          for(Int32 iline=0; iline<nlines; iline++)
          {
              Float64 A = m_Elements[eltIdx]->GetFittedAmplitude(iline);
              Float64 Sigma =
      m_Elements[eltIdx]->GetFittedAmplitudeErrorSigma(iline);

              Float64 snr = A/Sigma;
              if(snr > snrthres){
                  isSnrHigh = true;
                  break;
              }
          }
          if(isSnrHigh){
              highSNRvalidEltsIdx.push_back(eltIdx);
          }
      }
      TInt32List refreshIdxs =
      ReestimateContinuumUnderLines(highSNRvalidEltsIdx);
      refreshModelAfterContReestimation(refreshIdxs, modelFluxAxis,
      spcFluxAxisNoContinuum);
      //*/

          fitAmplitudesHybrid(spectralAxis, m_spcFluxAxisNoContinuum,
                              m_ContinuumFluxAxis, redshift);
          it++;
        }
      }

      // set all the amplitudes to 1.0
      if (m_fittingmethod == "ones") {
        for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
          m_Elements[iElts]->SetFittedAmplitude(1.0, 1.0);
        }
      }

      if (m_rigidity == "rules") {
        // Apply rules,
        applyRules(enableLogging);

        refreshModel();
        // create spectrum model
        Int32 modelSolutionLevel = Int32(enableLogging);
        modelSolution = GetModelSolution(modelSolutionLevel);
        continuumModelSolution = GetContinuumModelSolution();

        merit = getLeastSquareMerit();
      }

      // correct lines amplitude with tplshapePrior (tpl-corr): Warning: Rules
      // must all be deactivated
      if (m_rigidity == "tplcorr") {
        refreshModel();
        // create spectrum model
        modelSolution = GetModelSolution();
        continuumModelSolution = GetContinuumModelSolution();
        // Log.LogInfo( "LineModel Infos: TPLCORR");
        TFloat64List correctedAmplitudes;
        correctedAmplitudes.resize(modelSolution.Amplitudes.size());
        std::string bestTplName = "";

        m_CatalogTplShape.GetBestFit(modelSolution.Lines,
                                     modelSolution.Amplitudes,
                                     modelSolution.AmplitudesUncertainties,
                                     correctedAmplitudes, bestTplName);
        for (Int32 iRestLine = 0; iRestLine < m_RestLineList.size();
             iRestLine++) {
          Int32 subeIdx = undefIdx;
          Int32 eIdx = m_Elements.findElementIndex(iRestLine, subeIdx);
          if (eIdx == undefIdx || subeIdx == undefIdx) {
            Log.LogInfo("indexes are  -1");
            continue;
          }
          Float64 er = m_Elements[eIdx]->GetFittedAmplitudeErrorSigma(
              subeIdx); // not modifying the fitting error for now
          Float64 nominalAmp = m_Elements[eIdx]->GetNominalAmplitude(subeIdx);
          m_Elements[eIdx]->SetFittedAmplitude(
              correctedAmplitudes[iRestLine] / nominalAmp, er);
        }
        refreshModel();
        modelSolution = GetModelSolution();
        continuumModelSolution = GetContinuumModelSolution();
        m_tplshapeBestTplName = bestTplName;

        merit = getLeastSquareMerit();
      }

      if (m_rigidity == "tplshape") {
        Float64 _merit = INFINITY;
        // if(enableLogging || m_ContinuumComponent == "tplfit")
        if (enableLogging) {
          refreshModel();
          _merit = getLeastSquareMerit();
          //_merit = getLeastSquareMeritFast();
        } else {
          if (m_tplshapeLeastSquareFast) {
            _merit = getLeastSquareMeritFast();
          } else {
            refreshModel();
            _merit = getLeastSquareMerit();
          }
        }

        // lines prior
        Float64 _meritprior = 0.0;
        if (logPriorDataTplShape.size() > 0) {
          _meritprior += -2. * logPriorDataTplShape[ifitting].betaTE *
                         logPriorDataTplShape[ifitting].logprior_precompTE;
          _meritprior += -2. * logPriorDataTplShape[ifitting].betaA *
                         logPriorDataTplShape[ifitting].logprior_precompA;
          _meritprior += -2. * logPriorDataTplShape[ifitting].betaZ *
                         logPriorDataTplShape[ifitting].logprior_precompZ;
          if (logPriorDataTplShape[ifitting].A_sigma > 0.0) {
            //
            Float64 ampl = 0.0;
            for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
              bool foundAmp = false;
              Int32 nLines = m_Elements[iElts]->GetSize();
              for (Int32 j = 0; j < nLines; j++) {
                if (foundAmp) {
                  break;
                }
                Float64 amp = m_Elements[iElts]->GetFittedAmplitude(j);
                if (amp > 0 && !m_Elements[iElts]->IsOutsideLambdaRange(j)) {
                  Float64 nominal_amp =
                      m_Elements[iElts]->GetNominalAmplitude(j);
                  ampl = amp / nominal_amp;
                  foundAmp = true;
                  break;
                }
              }
            }
            //
            Float64 logPa = logPriorDataTplShape[ifitting].betaA *
                            (ampl - logPriorDataTplShape[ifitting].A_mean) *
                            (ampl - logPriorDataTplShape[ifitting].A_mean) /
                            (logPriorDataTplShape[ifitting].A_sigma *
                             logPriorDataTplShape[ifitting].A_sigma);
            _meritprior += logPa;
          }
        }

        if ((_merit + _meritprior) <
            meritTplratio[ifitting] + meritPriorTplratio[ifitting]) {
          meritTplratio[ifitting] = _merit;
          meritPriorTplratio[ifitting] = _meritprior;
          m_ChisquareTplshape[ifitting] = _merit;
          m_ScaleMargCorrTplshape[ifitting] = getScaleMargCorrection();
          m_StrongELPresentTplshape[ifitting] =
              GetModelStrongEmissionLinePresent();
          // given that Ha is a strong emission line,
          if (m_opt_haprior > 0.) // check first that haprior is activated
            m_StrongHalphaELPresentTplshape[ifitting] =
                m_StrongELPresentTplshape[ifitting]
                    ? GetModelHaStrongest()
                    : false; // result per tplratio

          TStringList strongELSNRAboveCut; // = getLinesAboveSNR(3.5); //this is
                                           // costing a lot of processing time,
                                           // so deactivated for now.
          m_NLinesAboveSNRTplshape[ifitting] = strongELSNRAboveCut.size();

          // Saving the model A, errorA, and dtm, mtm, ... (for all tplratios,
          // needed ?) NB: this is only needed for the index=savedIdxFitted
          // ultimately
          for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
            m_LinesLogPriorTplshape[ifitting][iElts] = _meritprior;

            bool savedAmp = false;
            // init
            //*
            m_FittedAmpTplshape[ifitting][iElts] = NAN;
            m_FittedErrorTplshape[ifitting][iElts] = NAN;
            m_DtmTplshape[ifitting][iElts] = NAN;
            m_MtmTplshape[ifitting][iElts] = NAN;
            m_LyaAsymCoeffTplshape[ifitting][iElts] = NAN;
            m_LyaWidthCoeffTplshape[ifitting][iElts] = NAN;
            m_LyaDeltaCoeffTplshape[ifitting][iElts] = NAN;
            m_LyaIgmIdxTplshape[ifitting][iElts] = undefIdx;
            bool allampzero = true;
            //*/

            Int32 nLines = m_Elements[iElts]->GetSize();
            for (Int32 j = 0; j < nLines; j++) {
              if (savedAmp) {
                break;
              }
              Float64 amp = m_Elements[iElts]->GetFittedAmplitude(j);
              if (amp > 0) {
                allampzero = false;
              }
              if (amp > 0 && !m_Elements[iElts]->IsOutsideLambdaRange(j)) {

                Float64 amp_error =
                    m_Elements[iElts]->GetFittedAmplitudeErrorSigma(j);
                Float64 nominal_amp = m_Elements[iElts]->GetNominalAmplitude(j);
                m_FittedAmpTplshape[ifitting][iElts] = amp / nominal_amp;
                Log.LogDebug(
                    "    model : fit tplratio mode, tplratio_fittedamp: %e",
                    m_FittedAmpTplshape[ifitting][iElts]);

                m_FittedErrorTplshape[ifitting][iElts] =
                    amp_error / nominal_amp;
                m_DtmTplshape[ifitting][iElts] =
                    m_Elements[iElts]->GetSumCross();
                m_MtmTplshape[ifitting][iElts] =
                    m_Elements[iElts]->GetSumGauss();

                TAsymParams params = m_Elements[iElts]->GetAsymfitParams(0);
                m_LyaAsymCoeffTplshape[ifitting][iElts] = params.alpha;
                m_LyaWidthCoeffTplshape[ifitting][iElts] = params.sigma;
                m_LyaDeltaCoeffTplshape[ifitting][iElts] = params.delta;

                TSymIgmParams params_igm =
                    m_Elements[iElts]->GetSymIgmParams(0);
                m_LyaIgmIdxTplshape[ifitting][iElts] = params_igm.m_igmidx;

                savedAmp = true;
                break;
              }
            }

            if (allampzero &&
                !savedAmp) // TODO: this case should be treated more carefully,
                           // save dtm, mtm, and more...
            {
              m_FittedAmpTplshape[ifitting][iElts] = 0.0;
            }
          }
        }

        if (merit + meritprior > _merit + _meritprior) {
          merit = _merit;
          meritprior = _meritprior;
          savedIdxContinuumFitted = icontfitting;
          savedIdxFitted = ifitting;

          modelSolution = GetModelSolution();
          continuumModelSolution = GetContinuumModelSolution();
          m_tplshapeBestTplName =
              m_CatalogTplShape.GetCatalogName(savedIdxFitted);
          m_tplshapeBestTplIsmCoeff =
              m_CatalogTplShape.GetIsmCoeff(savedIdxFitted);
          m_tplshapeBestTplAmplitude =
              m_FittedAmpTplshape[savedIdxFitted][0]; // Should be only 1 elt in
                                                      // tpl ratio mode...
          m_tplshapeBestTplDtm =
              m_DtmTplshape[savedIdxFitted]
                           [0]; // Should be only 1 elt in tpl ratio mode...
          m_tplshapeBestTplMtm =
              m_MtmTplshape[savedIdxFitted]
                           [0]; // Should be only 1 elt in tpl ratio mode...
        }
      } else
        savedIdxContinuumFitted = 0;

      // if(m_rigidity=="tplcorr")
      //{
      //   Float64 tplshapePriorCoeff = m_CatalogTplShape.GetBestFit(
      //   modelSolution.Lines, modelSolution.Amplitudes,  ); Log.LogDebug(
      //   "Linemodel: tplshapePriorCoeff = %f", tplshapePriorCoeff); merit =
      //   sqrt(merit*merit/2.0-log(tplshapePriorCoeff));
      // }

      if (m_ContinuumComponent == "nocontinuum") {
        reinitModel();
      }
    }
  }

  if (enableLogging) {
    if (m_ContinuumComponent == "tplfit" ||
        m_ContinuumComponent == "tplfitauto") {
      if (ncontinuumfitting > 1 && m_fittingmethod != "svdlc") {
        Int32 autoselect = 0;
        if (m_ContinuumComponent == "tplfit") {
          autoselect = 0;
        } else if (m_ContinuumComponent == "tplfitauto") {
          autoselect = 1;
        }
        // TODO savedIdxContinuumFitted=-1 if rigidity!=tplshape
        LoadFitContinuum(savedIdxContinuumFitted, autoselect);
      }
      if (m_fittingmethod == "svdlc") {
        //
      }
      Log.LogDetail("    model - Linemodel: fitcontinuum = %d (%s, with "
                    "ebmv=%.3f), and A=%e",
                    savedIdxContinuumFitted, m_fitContinuum_tplName.c_str(),
                    m_fitContinuum_tplFitEbmvCoeff,
                    m_fitContinuum_tplFitAmplitude);
    }

    if (m_rigidity == "tplshape") {
      // m_CatalogTplShape.SetMultilineNominalAmplitudes( *this, savedIdxFitted
      // );
      bool retSetMultiAmplFast =
          SetMultilineNominalAmplitudesFast(savedIdxFitted);
      if (!retSetMultiAmplFast) {
        Log.LogError("Linemodel: tplshape, Unable to set Multiline "
                     "NominalAmplitudes from Tplshape !");
      }

      // Set the velocities from templates: todo auto switch when velfit is ON
      // m_CatalogTplShape.GetCatalogVelocities(savedIdxFitted,
      // m_velocityEmission, m_velocityAbsorption);
      for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
        Log.LogDetail("    model - Linemodel: tplratio = %d (%s, with "
                      "ebmv=%.3f), and A=%e",
                      savedIdxFitted, m_tplshapeBestTplName.c_str(),
                      m_tplshapeBestTplIsmCoeff,
                      m_FittedAmpTplshape[savedIdxFitted][iElts]);
        m_Elements[iElts]->SetFittedAmplitude(
            m_FittedAmpTplshape[savedIdxFitted][iElts],
            m_FittedErrorTplshape[savedIdxFitted][iElts]);
        m_Elements[iElts]->SetSumCross(m_DtmTplshape[savedIdxFitted][iElts]);
        m_Elements[iElts]->SetSumGauss(m_MtmTplshape[savedIdxFitted][iElts]);
      }

      // Lya
      for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
        m_Elements[iElts]->SetAsymfitParams(
            {m_LyaWidthCoeffTplshape[savedIdxFitted][iElts],
             m_LyaAsymCoeffTplshape[savedIdxFitted][iElts],
             m_LyaDeltaCoeffTplshape[savedIdxFitted][iElts]});
      }

      for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
        m_Elements[iElts]->SetSymIgmParams(TSymIgmParams(
            m_LyaIgmIdxTplshape[savedIdxFitted][iElts], redshift));
      }

      refreshModel();

      Int32 modelSolutionLevel = Int32(enableLogging);
      modelSolution = GetModelSolution(modelSolutionLevel);
      continuumModelSolution = GetContinuumModelSolution();
    }
  }
  return merit;
}

std::vector<CLmfitController *> CLineModelFitting::createLmfitControllers() {
  std::vector<CLmfitController *> useLmfitControllers;
  if (m_lmfit_noContinuumTemplate) {
    useLmfitControllers.push_back(new CLmfitController(
        m_lmfit_fitEmissionVelocity, m_lmfit_fitAbsorptionVelocity));
  } else {
    if (m_lmfit_bestTemplate) {
      LoadFitContinuum(-1, 0);
      std::shared_ptr<const CTemplate> tpl = m_tplCatalog.GetTemplateByName(
          m_tplCategoryList, m_fitContinuum_tplName);
      if (m_fitContinuum_tplName == tpl->GetName()) {
        bool continumLoaded = true;
        useLmfitControllers.push_back(new CLmfitController(
            tpl, continumLoaded, m_lmfit_fitContinuum,
            m_lmfit_fitEmissionVelocity, m_lmfit_fitAbsorptionVelocity));
      }
    } else {
      bool continumLoaded = false;
      for (Int32 i = 0; i < m_tplCategoryList.size(); i++) {
        std::string category = m_tplCategoryList[i];

        for (Int32 j = 0; j < m_tplCatalog.GetTemplateCount(category); j++) {
          std::shared_ptr<const CTemplate> tpl =
              m_tplCatalog.GetTemplate(category, j);
          useLmfitControllers.push_back(new CLmfitController(
              tpl, continumLoaded, m_lmfit_fitContinuum,
              m_lmfit_fitEmissionVelocity, m_lmfit_fitAbsorptionVelocity));
        }
      }
    }
  }
  return useLmfitControllers;
}

void CLineModelFitting::SetSecondpassContinuumFitPrms(Int32 dustfit,
                                                      Int32 meiksinfit,
                                                      Int32 outsidelinemask,
                                                      Int32 observedFrame) {
  m_secondpass_fitContinuum_dustfit = dustfit;
  m_secondpass_fitContinuum_igm = meiksinfit;
  m_secondpass_fitContinuum_outsidelinesmask =
      outsidelinemask; // hardcoded deactivated because ortho templates are used
  m_secondpass_fitContinuum_observedFrame = observedFrame;

  if (1) {
    Log.LogDetail(
        "Elementlist: SetSecondpassContinuumFitPrms fitContinuum_dustfit = %d",
        m_secondpass_fitContinuum_dustfit);
    Log.LogDetail(
        "Elementlist: SetSecondpassContinuumFitPrms fitContinuum_igm = %d",
        m_secondpass_fitContinuum_igm);
    Log.LogDetail("Elementlist: SetSecondpassContinuumFitPrms "
                  "fitContinuum_outsidelinemask = %d",
                  m_secondpass_fitContinuum_outsidelinesmask);
    Log.LogDetail("Elementlist: SetSecondpassContinuumFitPrms "
                  "fitContinuum_observedFrame = %d",
                  m_secondpass_fitContinuum_observedFrame);
  }
}

void CLineModelFitting::SetFittingMethod(const std::string &fitMethod) {
  m_fittingmethod = fitMethod;
}

void CLineModelFitting::SetLeastSquareFastEstimationEnabled(Int32 enabled) {
  m_tplshapeLeastSquareFast = enabled;
}

void CLineModelFitting::SetAbsLinesLimit(Float64 limit) {
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    m_Elements[iElts]->SetAbsLinesLimit(limit);
  }
}

/**
 * \brief Init the whole spectrum model with continuum.
 **/
void CLineModelFitting::reinitModel() {
  CSpectrumFluxAxis modelFluxAxis = m_SpectrumModel.GetFluxAxis();

  /*
  //check the continuum flux axis
  const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel.GetSpectralAxis();
  for(Int32 j=0; j<m_ContinuumFluxAxis.GetSamplesCount(); j++)
  {
      if(isnan(m_ContinuumFluxAxis[j]))
      {
          THROWG(INTERNAL_ERROR,
  "NaN value found for the ContinuumFluxAxis at
  lambda=%f", spectralAxis[j] ); break;
      }
  }*/

  // init spectrum model with continuum
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    m_Elements[iElts]->initSpectrumModel(modelFluxAxis, m_ContinuumFluxAxis);
  }

  m_SpectrumModel.SetFluxAxis(std::move(modelFluxAxis));
}

/**
 * \brief Init the argument elements from the spectrum model with continuum.
 **/
void CLineModelFitting::reinitModelUnderElements(
    const TInt32List &filterEltsIdx, Int32 lineIdx) {
  Int32 iElts;
  CSpectrumFluxAxis modelFluxAxis = m_SpectrumModel.GetFluxAxis();
  // init spectrum model with continuum
  for (Int32 i = 0; i < filterEltsIdx.size(); i++) {
    iElts = filterEltsIdx[i];
    m_Elements[iElts]->initSpectrumModel(modelFluxAxis, m_ContinuumFluxAxis,
                                         lineIdx);
  }
  m_SpectrumModel.SetFluxAxis(std::move(modelFluxAxis));
}

CSpectrumFluxAxis CLineModelFitting::getModel(Int32 lineTypeFilter) const {
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  CSpectrumFluxAxis modelfluxAxis(spectralAxis.GetSamplesCount());

  Int32 nElements = m_Elements.size();
  for (Int32 iElts = 0; iElts < nElements; iElts++) {
    m_Elements[iElts]->initSpectrumModel(modelfluxAxis, m_ContinuumFluxAxis);
  }

  for (Int32 iElts = 0; iElts < nElements; iElts++) {
    Int32 lineType = m_Elements[iElts]->m_Lines[0].GetType();
    if (lineTypeFilter == -1 || lineTypeFilter == lineType) {
      m_Elements[iElts]->addToSpectrumModel(spectralAxis, modelfluxAxis,
                                            m_ContinuumFluxAxis, m_Redshift);
    }
  }

  return modelfluxAxis;
}

/**
 * \brief Adds a new model to each m_Elements entry.
 * Calls reinitModel.
 * For each entry in m_Elements, addToSpectrumModel using the reinitModel output
 *as arguments.
 **/
void CLineModelFitting::refreshModel(Int32 lineTypeFilter) {
  reinitModel();

  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  CSpectrumFluxAxis modelFluxAxis = m_SpectrumModel.GetFluxAxis();

  bool enableMinMaxLog = false;
  if (enableMinMaxLog) {
    // check the model min/max values
    // Int32 imin = spectralAxis.GetIndexAtWaveLength(12500.);
    // Int32 imax = spectralAxis.GetIndexAtWaveLength(18500.);
    Int32 imin = 0;
    Int32 imax = modelFluxAxis.GetSamplesCount();
    Float64 fmin = DBL_MAX;
    Float64 fmax = -DBL_MAX;
    for (Int32 j = imin; j < imax; j++) {
      if (as_const(modelFluxAxis)[j] < fmin) {
        fmin = as_const(modelFluxAxis)[j];
      }
      if (modelFluxAxis[j] > fmax) {
        fmax = as_const(modelFluxAxis)[j];
      }
    }
    Log.LogDebug("CLineModelFitting::refreshModel AFTER REINIT: model min=%e "
                 "and model max=%e",
                 fmin, fmax);
  }

  /*
  //check the model reinited for nan values
  Int32 imin = spectralAxis.GetIndexAtWaveLength(12500.);
  Int32 imax = spectralAxis.GetIndexAtWaveLength(18500.);
  //Int32 imin = 0;
  //Int32 imax = modelFluxAxis.GetSamplesCount();

  for(Int32 j=imin; j<imax; j++)
  {
      if(isnan(modelFluxAxis[j]))
      {
          THROWG(INTERNAL_ERROR,
  "CLineModelFitting::refreshModel: NaN value found for the reinited model
  spectrum at lambda=%f", spectralAxis[j] ); break;
      }
  }
  //*/

  /*
  CSpectrumFluxAxis contFluxAxisWithAmpOffset =
  CSpectrumFluxAxis(m_ContinuumFluxAxis.GetSamplesCount()); for( Int32 i=0;
  i<m_ContinuumFluxAxis.GetSamplesCount(); i++ )
  {
      contFluxAxisWithAmpOffset[i] = m_ContinuumFluxAxis[i];
  }
  //*/
  if (m_enableAmplitudeOffsets) {
    // add amplitude offsets
    m_Elements.addToSpectrumAmplitudeOffset(m_SpectrumModel.GetSpectralAxis(),
                                            modelFluxAxis);
    // addToSpectrumAmplitudeOffset(contFluxAxisWithAmpOffset);
  }

  // create spectrum model
  Int32 nElements = m_Elements.size();
  for (Int32 iElts = 0; iElts < nElements; iElts++) {
    Int32 lineType = m_Elements[iElts]->m_Lines[0].GetType();
    if (lineTypeFilter == -1 || lineTypeFilter == lineType) {
      m_Elements[iElts]->addToSpectrumModel(spectralAxis, modelFluxAxis,
                                            m_ContinuumFluxAxis, m_Redshift);
    }

    /*
    //check the model for nan values
    Int32 imin = spectralAxis.GetIndexAtWaveLength(12500.);
    Int32 imax = spectralAxis.GetIndexAtWaveLength(18500.);
    //Int32 imin = 0;
    //Int32 imax = modelFluxAxis.GetSamplesCount();

    for(Int32 j=imin; j<imax; j++)
    {
        if(isnan(modelFluxAxis[j]))
        {
            THROWG(INTERNAL_ERROR,
    "CLineModelFitting::refreshModel: NaN value found for the model spectrum at
    lambda=%f", spectralAxis[j] ); break;
        }
    }
    //*/
  }

  if (enableMinMaxLog) {
    // check the model min/max values
    Int32 imin = spectralAxis.GetIndexAtWaveLength(12500.);
    Int32 imax = spectralAxis.GetIndexAtWaveLength(18500.);
    // Int32 imin = 0;
    // Int32 imax = modelFluxAxis.GetSamplesCount();
    Float64 fmin = DBL_MAX;
    Float64 fmax = -DBL_MAX;
    for (Int32 j = imin; j < imax; j++) {
      if (as_const(modelFluxAxis)[j] < fmin) {
        fmin = as_const(modelFluxAxis)[j];
      }
      if (as_const(modelFluxAxis)[j] > fmax) {
        fmax = as_const(modelFluxAxis)[j];
      }
    }
    Log.LogDebug("CLineModelFitting::refreshModel AFTER addToSpectrumModel: "
                 "model min=%e and model max=%e",
                 fmin, fmax);
  }

  m_SpectrumModel.SetFluxAxis(std::move(modelFluxAxis));
}

/**
 * \brief refreshing all the grid and Adds a new model to each m_Elements entry
 *. Calls iterate on model flux to set it equal to continuum. For each entry in
 *m_Elements, addToSpectrumModel using the reinitModel output as arguments.
 **/
void CLineModelFitting::refreshModelInitAllGrid() {

  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  CSpectrumFluxAxis modelFluxAxis = m_ContinuumFluxAxis;

  // create spectrum model
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    m_Elements[iElts]->addToSpectrumModel(spectralAxis, modelFluxAxis,
                                          m_ContinuumFluxAxis, m_Redshift);
  }

  m_SpectrumModel.SetFluxAxis(std::move(modelFluxAxis));
}

/**
 * \brief Adds a new model to each m_Elements entry specified on the argument.
 * Works as refreshModel.
 **/
void CLineModelFitting::refreshModelUnderElements(
    const TInt32List &filterEltsIdx, Int32 lineIdx) {
  reinitModelUnderElements(filterEltsIdx, lineIdx);
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  CSpectrumFluxAxis modelFluxAxis = m_SpectrumModel.GetFluxAxis();
  // create spectrum model
  Int32 iElts;
  for (Int32 i = 0; i < filterEltsIdx.size(); i++) {
    iElts = filterEltsIdx[i];
    m_Elements[iElts]->addToSpectrumModel(
        spectralAxis, modelFluxAxis, m_ContinuumFluxAxis, m_Redshift, lineIdx);
  }
  m_SpectrumModel.SetFluxAxis(std::move(modelFluxAxis));
}

void CLineModelFitting::refreshModelDerivVelEmissionUnderElements(
    const TInt32List &filterEltsIdx) {
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  TInt32List supportIdxes = m_Elements.getSupportIndexes(filterEltsIdx);

  for (Int32 i = 0; i < supportIdxes.size(); i++) {
    m_SpcFluxAxisModelDerivVelEmi[supportIdxes[i]] = 0.0;
  }

  // create spectrum model partial derivate vs sigma
  Int32 iElts;
  for (Int32 i = 0; i < filterEltsIdx.size(); i++) {
    iElts = filterEltsIdx[i];
    m_Elements[iElts]->addToSpectrumModelDerivVel(
        spectralAxis, m_SpcFluxAxisModelDerivVelEmi, m_ContinuumFluxAxis,
        m_Redshift, true);
  }
}

void CLineModelFitting::refreshModelDerivVelAbsorptionUnderElements(
    const TInt32List &filterEltsIdx) {
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  TInt32List supportIdxes = m_Elements.getSupportIndexes(filterEltsIdx);

  for (Int32 i = 0; i < supportIdxes.size(); i++) {
    m_SpcFluxAxisModelDerivVelAbs[supportIdxes[i]] = 0.0;
  }

  // create spectrum model partial derivate vs sigma
  Int32 iElts;
  for (Int32 i = 0; i < filterEltsIdx.size(); i++) {
    iElts = filterEltsIdx[i];
    m_Elements[iElts]->addToSpectrumModelDerivVel(
        spectralAxis, m_SpcFluxAxisModelDerivVelAbs, m_ContinuumFluxAxis,
        m_Redshift, false);
  }
}

void CLineModelFitting::refreshModelDerivVelUnderElements(
    const TInt32List &filterEltsIdx) {
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  TInt32List supportIdxes = m_Elements.getSupportIndexes(filterEltsIdx);

  for (Int32 i = 0; i < supportIdxes.size(); i++) {
    m_SpcFluxAxisModelDerivVelEmi[supportIdxes[i]] = 0.0;
    m_SpcFluxAxisModelDerivVelAbs[supportIdxes[i]] = 0.0;
  }

  // create spectrum model partial derivate vs sigma
  Int32 iElts;
  for (Int32 i = 0; i < filterEltsIdx.size(); i++) {
    iElts = filterEltsIdx[i];
    m_Elements[iElts]->addToSpectrumModelDerivVel(
        spectralAxis, m_SpcFluxAxisModelDerivVelEmi, m_ContinuumFluxAxis,
        m_Redshift, true);
    m_Elements[iElts]->addToSpectrumModelDerivVel(
        spectralAxis, m_SpcFluxAxisModelDerivVelAbs, m_ContinuumFluxAxis,
        m_Redshift, false);
  }
}

void CLineModelFitting::setModelSpcObservedOnSupportZeroOutside() {
  m_Redshift = 0.0;

  // initialize the model spectrum
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  CSpectrumFluxAxis modelFluxAxis(
      m_SpectrumModel.GetSampleCount()); // set to zeros

  // prepare the elements
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    m_Elements[iElts]->prepareSupport(spectralAxis, m_Redshift, m_lambdaRange);
  }

  TInt32List validEltsIdx = m_Elements.GetModelValidElementsIndexes();
  TInt32List supportIdxes = m_Elements.getSupportIndexes(validEltsIdx);
  const auto &SpcFluxAxis = m_SpcFluxAxis;
  for (Int32 i = 0; i < supportIdxes.size(); i++) {
    modelFluxAxis[supportIdxes[i]] = SpcFluxAxis[supportIdxes[i]];
  }

  m_SpectrumModel.SetFluxAxis(std::move(modelFluxAxis));
}

/**
 * \brief Creates and returns a Mask with 0 in the lines support, 1 under the
 *lines
 **/
CMask CLineModelFitting::getOutsideLinesMask() const {
  CMask _mask;
  // initialize the model spectrum
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
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
 * input: which = 1: uses the spectrum flux continuum subtracted to compute STD
 * input: which = 2: uses the spectrum error to compute STD
 **/
Float64 CLineModelFitting::getOutsideLinesSTD(Int32 which) const {
  if (which != 1 && which != 2) {
    Log.LogError("    model: getOutsideLinesSTD - Failed to parse input "
                 "argument, which");
    return -1;
  }

  CMask _mask = getOutsideLinesMask();

  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  Float64 sum2 = 0.0;
  Int32 nsum = 0;
  Int32 imin = spectralAxis.GetIndexAtWaveLength(m_lambdaRange.GetBegin());
  Int32 imax = spectralAxis.GetIndexAtWaveLength(m_lambdaRange.GetEnd());
  const auto &spcFluxAxisNoContinuum = m_spcFluxAxisNoContinuum;
  const auto &ErrorNoContinuum = m_ErrorNoContinuum;
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
 * \brief Tries to fit subelements considering their overlap.
 * For each entry in GetModelValidElementsIndexes:
 *   If subelement in the entry already fitted, go for the next entry.
 *   getOverlappingElements for the fitted subelements.
 *   If the overlap is smaller than 2, call fitAmplitude on the entry.
 *   If the overlap is greater than or equal to 2:
 *     Call fitAmplitudeLinSolve with the subelements as argument.
 *     Store all non-negative fits.
 *     Set to 0.0 all negative fits.
 *     If the size of non-negative fits is 1, call the entry's fitAmplitude.
 *     If the size of non-negative fits is not 1:
 *       If the size of non-negative fits is greater than 1:
 *         Call fitAmplitudesLinSolve with the indexes of the non-negative
 *subelements. If the above call return is different than 1: For each
 *non-negative subelement, if the amplitude fitted is greater than 0, call
 *fitAmplitude on its entry. Else, SetElementAmplitude to 0. Update the index of
 *already-fitted subelements.
 **/
Int32 CLineModelFitting::fitAmplitudesHybrid(
    const CSpectrumSpectralAxis &spectralAxis,
    const CSpectrumFluxAxis &spcFluxAxisNoContinuum,
    const CSpectrumFluxAxis &continuumfluxAxis, Float64 redshift) {
  const bool verbose = true;
  TInt32List validEltsIdx = m_Elements.GetModelValidElementsIndexes();
  TInt32List indexesFitted;
  for (Int32 iValidElts = 0; iValidElts < validEltsIdx.size(); iValidElts++) {
    Int32 iElts = validEltsIdx[iValidElts];
    // skip if already fitted
    bool alreadyfitted = false;
    for (Int32 i = 0; i < indexesFitted.size(); i++) {
      if (iElts == indexesFitted[i]) {
        alreadyfitted = true;
        break;
      }
    }
    if (alreadyfitted) {
      continue;
    }
    // do the fit on the ovelapping elements
    Float64 overlapThres =
        0.15; // 15% seemed necessary for Ha/SII complex when lines are very
              // wide (either because of PSF or source size)
    TInt32List overlappingInds = m_Elements.getOverlappingElements(
        iElts, indexesFitted, m_Redshift, overlapThres);

    // setting the fitting group info
    for (Int32 ifit = 0; ifit < overlappingInds.size(); ifit++) {
      std::string fitGroupTag = boost::str(boost::format("hy%d") % iValidElts);
      m_Elements[overlappingInds[ifit]]->m_fittingGroupInfo = fitGroupTag;
    }

    if (verbose) {
      // Log.LogDebug( "Redshift: %f", m_Redshift);
      Log.LogDebug("    model: hybrid fit: #%d - N overlapping=%d", iValidElts,
                   overlappingInds.size());
      for (Int32 ifit = 0; ifit < overlappingInds.size(); ifit++) {
        Log.LogDebug("    model: hybrid fit:     overlapping #%d - eltIdx=%d",
                     ifit, overlappingInds[ifit]);
      }
    }
    if (!m_enableAmplitudeOffsets && overlappingInds.size() < 2) {
      if (verbose) {
        Log.LogDebug("    model: hybrid fit:     Individual fit");
      }
      m_Elements[iElts]->fitAmplitudeAndLambdaOffset(
          spectralAxis, spcFluxAxisNoContinuum, continuumfluxAxis, redshift,
          undefIdx, m_enableLambdaOffsetsFit, m_LambdaOffsetStep,
          m_LambdaOffsetMin, m_LambdaOffsetMax);
    } else {
      if (verbose) {
        Log.LogDebug("    model: hybrid fit:     SVD fit");
      }
      // fit individually: mainly for mtm and dtm estimation
      for (Int32 ifit = 0; ifit < overlappingInds.size(); ifit++) {
        m_Elements[overlappingInds[ifit]]->fitAmplitudeAndLambdaOffset(
            spectralAxis, spcFluxAxisNoContinuum, continuumfluxAxis, redshift,
            undefIdx, m_enableLambdaOffsetsFit, m_LambdaOffsetStep,
            m_LambdaOffsetMin, m_LambdaOffsetMax);
      }
      TFloat64List ampsfitted;
      TFloat64List errorsfitted;
      Int32 retVal = fitAmplitudesLinSolveAndLambdaOffset(
          overlappingInds, spectralAxis, spcFluxAxisNoContinuum,
          continuumfluxAxis, ampsfitted, errorsfitted,
          m_enableLambdaOffsetsFit);
      // if all the amplitudes fitted don't have the same sign, do it separately
      TInt32List overlappingIndsSameSign;
      if (retVal != 1 && ampsfitted.size() > 0) {
        for (Int32 ifit = 0; ifit < overlappingInds.size(); ifit++) {
          if (ampsfitted[ifit] > 0) {
            overlappingIndsSameSign.push_back(overlappingInds[ifit]);
            // m_Elements[overlappingInds[ifit]]->fitAmplitude(spectralAxis,
            // spcFluxAxisNoContinuum, redshift);
          } else {
            m_Elements.SetElementAmplitude(overlappingInds[ifit], 0.0,
                                           errorsfitted[ifit]);
          }
        }
        // fit the rest of the overlapping elements (same sign) together
        if (!m_enableAmplitudeOffsets && overlappingIndsSameSign.size() == 1) {
          m_Elements[overlappingIndsSameSign[0]]->fitAmplitudeAndLambdaOffset(
              spectralAxis, spcFluxAxisNoContinuum, continuumfluxAxis, redshift,
              undefIdx, m_enableLambdaOffsetsFit, m_LambdaOffsetStep,
              m_LambdaOffsetMin, m_LambdaOffsetMax);
        } else if (overlappingIndsSameSign.size() > 0) {
          Int32 retVal2 = fitAmplitudesLinSolveAndLambdaOffset(
              overlappingIndsSameSign, spectralAxis, spcFluxAxisNoContinuum,
              continuumfluxAxis, ampsfitted, errorsfitted,
              m_enableLambdaOffsetsFit);

          if (retVal2 != 1) {
            for (Int32 ifit = 0; ifit < overlappingIndsSameSign.size();
                 ifit++) {
              if (ampsfitted[ifit] > 0) {
                m_Elements[overlappingIndsSameSign[ifit]]
                    ->fitAmplitudeAndLambdaOffset(
                        spectralAxis, spcFluxAxisNoContinuum, continuumfluxAxis,
                        redshift, undefIdx, m_enableLambdaOffsetsFit,
                        m_LambdaOffsetStep, m_LambdaOffsetMin,
                        m_LambdaOffsetMax);
              } else {
                m_Elements.SetElementAmplitude(overlappingIndsSameSign[ifit],
                                               0.0, errorsfitted[ifit]);
              }
            }
          }
        }
      }
    }

    // update the already fitted list
    for (Int32 i = 0; i < overlappingInds.size(); i++) {
      indexesFitted.push_back(overlappingInds[i]);
    }
  }

  if (m_opt_enable_improveBalmerFit) {
    improveBalmerFit();
  }

  return 0;
}

/**
 * @brief CLineModelFitting::estimateMeanSqFluxAndGradient
 * @param varPack: the variables of the model being fitted
 * @param normFactor
 * @param filteredEltsIdx: index of the elements included in the fit, also sets
 * the size of varPack (n=nElts+1)
 * @param xInds: indexes of the samples where the data of the model is fitted
 * @param lineType: E or A
 * @param fluxdata: data to be fitted (already reshaped, no need to use xInds
 * for this vector)
 * @param msqBuffer: buffer for fast computing of the meansquare
 * @param f: output meansquare residual
 * @param g: output meansquare gradient residual
 * @return
 */
Int32 CLineModelFitting::estimateMeanSqFluxAndGradient(
    const Float64 *varPack, const Float64 normFactor,
    const TInt32List &filteredEltsIdx, const TInt32List &xInds, Int32 lineType,
    Float64 *fluxdata, Float64 *msqBuffer, Float64 &f, Float64 *g) {
  // update the linemodel amplitudes
  for (Int32 iElt = 0; iElt < filteredEltsIdx.size(); iElt++) {
    Float64 amp = varPack[iElt] / normFactor;
    m_Elements.SetElementAmplitude(filteredEltsIdx[iElt], amp, 0.0);
  }
  // update the linemodel velocity/linewidth
  Int32 idxVelocity = filteredEltsIdx.size();
  Float64 velocity = varPack[idxVelocity];
  if (lineType == CLine::nType_Emission) {
    SetVelocityEmission(velocity);
  } else {
    SetVelocityAbsorption(velocity);
  }

  Int32 nsamples = xInds.size();
  // retrieve the model
  refreshModelUnderElements(filteredEltsIdx);
  f = 0.0;
  for (Int32 i = 0; i < nsamples; i++) {
    Float64 Yi = getModelFluxVal(xInds[i]) * normFactor;
    msqBuffer[i] = Yi - fluxdata[i];
    f += msqBuffer[i] * msqBuffer[i];
  }

  refreshModelDerivVelUnderElements(filteredEltsIdx);
  for (Int32 iElt = 0; iElt < filteredEltsIdx.size(); iElt++) {
    g[iElt] = 0.0;
  }
  for (Int32 i = 0; i < nsamples; i++) {
    for (Int32 iElt = 0; iElt < filteredEltsIdx.size(); iElt++) {
      Float64 dm = getModelFluxDerivEltVal(filteredEltsIdx[iElt], xInds[i]);
      Float64 grad = 2 * dm * msqBuffer[i];
      g[iElt] += grad;
    }
    //*
    Int32 iElt = filteredEltsIdx.size();
    // getModelFluxDerivVelVal change for LMFit, emssion and aborption have
    // there own deriv sigma value
    Float64 dm = getModelFluxDerivVelVal(xInds[i]) * normFactor;
    Float64 grad = 2 * dm * msqBuffer[i];
    g[iElt] += grad;
    //*/
  }

  return 0;
}

// temporary stuff, dev of the lmfit method: int

void print_state(size_t iter, gsl_multifit_fdfsolver *s) {
  printf("iter: %3lu x = % 15.8f % 15.8f % 15.8f "
         "|f(x)| = %g\n",
         iter, gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1),
         gsl_vector_get(s->x, 2), gsl_blas_dnrm2(s->f));
}

Int32 CLineModelFitting::fitAmplitudesLmfit(const CSpectrumFluxAxis &fluxAxis,
                                            CLmfitController *controller)
// TFloat64List& ampsfitted, TFloat64List& ampsErrfitted, Float64&
// velocityFitted, Float64& continuumAmpFitted, Float64& merit, Int32 lineType
{
  // http://www.gnu.org/software/gsl/manual/html_node/Example-programs-for-Nonlinear-Least_002dSquares-Fitting.html
  bool verbose = true;

  if (verbose) {
    // Log.LogInfo("fitAmplitudesLmfit");
  }
  TInt32List filteredEltsIdx = controller->getFilteredIdx();
  Int32 nddl = filteredEltsIdx.size();
  if (nddl < 1) {
    Log.LogError("Linemodel LMfit: No line amplitude to fit");
    return -1;
  }

  TInt32List xInds;
  if (!controller->isContinuumFitted()) {
    xInds = m_Elements.getSupportIndexes(filteredEltsIdx);
  } else {
    xInds = TInt32List(fluxAxis.GetSamplesCount());
    // boost::push_back(xInds, fluxAxis.GetSamplesCount());
    std::iota(std::begin(xInds), std::end(xInds), 0);
  }

  const Float64 *flux = fluxAxis.GetSamples();
  Float64 amplitudeContinuumInit = m_fitContinuum_tplFitAmplitude;

  // create gsl solver object
  const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
  gsl_multifit_fdfsolver *s;
  int status, info;
  size_t i;
  size_t n; // n samples on the support, /* number of data points to fit */
  size_t p =
      controller
          ->getNumberParameters(); // DOF = n amplitudes to fit (1 for each
                                   // element) + 1 (EL velocity) + 1 Absorption
                                   // VEl + 1 Amplitude Template

  n = xInds.size();
  if (n < nddl) {
    controller->resizeAmpsLine();
    Log.LogError("LineModel LMfit: not enough smaples on support");
    for (Int32 iddl = 0; iddl < nddl; iddl++) {
      controller->setAmpLine(iddl, 0.0, 0.0);
    }
    return -1;
  }

  // set the norm factor
  Int32 idx = 0;
  Float64 maxabsval = DBL_MIN;
  Float64 normFactor;
  // if(controller->needCalculateNormFactor()){
  for (i = 0; i < n; i++) {
    idx = xInds[i];

    if (maxabsval < std::abs(flux[idx])) {
      maxabsval = std::abs(flux[idx]);
    }
  }
  normFactor = 1.0 / maxabsval;
  controller->setNormFactor(normFactor);
  // }else{
  //   // the norm factor is the same as the last call
  //   normFactor=controller->getNormFactor();
  // }

  if (verbose) {
    // Log.LogDetail("normFactor = '%.3e'", normFactor);
  }

  gsl_matrix *J = gsl_matrix_alloc(n, p);
  gsl_matrix *covar = gsl_matrix_alloc(p, p);
  double y[n], weights[n];
  struct lmfitdata d = {
      n, y, this, xInds, m_observeGridContinuumFlux.data(), controller};
  gsl_multifit_function_fdf f;

  Float64 *x_init = (Float64 *)calloc(p, sizeof(Float64));
  // initialize lmfit with previously estimated individual/hybrid fit method
  Float64 bestAmpLine = -1.0;
  for (Int32 kp = 0; kp < nddl; kp++) {
    Float64 amp = m_Elements[filteredEltsIdx[kp]]->GetElementAmplitude();
    if (amp > bestAmpLine) {
      bestAmpLine = amp;
    }
  }
  controller->setNormAmpLine(1.0);
  //    if(bestAmpLine>0.){
  //        controller->setNormAmpLine(1/bestAmpLine);//setNormAmpLine(controller->lineAmp_LmToModel(bestAmpLine));
  //        if(verbose){
  //            Log.LogDetail("LineModel LMfit: normAmpLine : %f",
  //            controller->getNormAmpLine());
  //        }
  //    }
  for (Int32 kp = 0; kp < nddl; kp++) {
    Float64 ampInitGuess =
        m_Elements[filteredEltsIdx[kp]]
            ->GetElementAmplitude(); // std::max(m_Elements[filteredEltsIdx[kp]]->GetElementAmplitude()
                                     // *1.5, bestAmpLine*0.001) ;
    if (ampInitGuess < 0.0) {
      Log.LogDetail("Amp negative %d , set to 0", ampInitGuess);
      ampInitGuess = 0.0;
    }
    // Float64 ampInitGuess = 0.0;

    x_init[kp] = controller->lineAmp_ModelToLm(ampInitGuess); //*normFactor;
    if (verbose) {
      Log.LogDetail(
          "LineModel LMfit: set init guess amp [%d] Model = %f lmfit: %f",
          filteredEltsIdx[kp], ampInitGuess, x_init[kp]);
    }
    x_init[kp] = ampInitGuess * normFactor;
  }

  if (controller->isEmissionVelocityFitted()) {
    Float64 normEmiFactor = 1.0; // bestAmpLine/GetVelocityEmission();//1.0;//
    controller->setNormEmiFactor(normEmiFactor);
    x_init[controller->getIndEmissionVel()] =
        controller->emiVel_ModelToLm(GetVelocityEmission()); // bestAmpLine;
    if (verbose) {
      Log.LogDetail("LineModel LMfit: normEmiFactor = %f", normEmiFactor);
      Log.LogDetail("LineModel LMfit: set init guess Emission velocity = %f ",
                    GetVelocityEmission());
    }
  }
  if (controller->isAbsorptionVelocityFitted()) {
    Float64 normAbsFactor =
        1.0; // bestAmpLine/GetVelocityAbsorption();//1.0f;//

    controller->setNormAbsFactor(normAbsFactor);
    x_init[controller->getIndAbsorptionVel()] =
        controller->absVel_ModelToLm(GetVelocityAbsorption()); // bestAmpLine;
    if (verbose) {
      Log.LogDetail("LineModel LMfit: normAbsFactor = %f", normAbsFactor);
      Log.LogDetail("LineModel LMfit: set init guess Absorption velocity = %f ",
                    GetVelocityAbsorption());
    }
  }
  if (controller->isContinuumFitted()) {
    x_init[controller->getIndContinuumAmp()] =
        controller->continuumAmp_ModelToLm(
            getFitContinuum_tplAmplitude()); //*normFactor;
    if (verbose) {
      Log.LogDetail("LineModel LMfit: set init guess continuum amp = %f ",
                    getFitContinuum_tplAmplitude());
    }
  }

  if (controller->isRedshiftFitted()) {
    x_init[controller->getIndRedshift()] = m_Redshift;
    if (verbose) {
      Log.LogDetail("LineModel LMfit: set init guess redshift = %f ",
                    m_Redshift);
    }
  }

  gsl_vector_view x = gsl_vector_view_array(x_init, p);
  gsl_vector_view w = gsl_vector_view_array(weights, n);
  // const gsl_rng_type * type;
  // gsl_rng * r;
  gsl_vector *res_f;
  double chi, chi0;

  const double xtol = 1e-8;
  const double gtol = 1e-8;
  const double ftol = 1e-8;
  Int32 maxIterations = 250;

  // gsl_rng_env_setup();

  // type = gsl_rng_default;
  // r = gsl_rng_alloc (type);

  f.f = &lmfit_f;
  f.df = &lmfit_df;
  f.n = n;
  f.p = p;
  f.params = &d;

  // This is the data to be fitted
  Float64 ei;
  const auto &ErrorNoContinuum = m_ErrorNoContinuum;
  for (i = 0; i < n; i++) {
    idx = xInds[i];
    ei = ErrorNoContinuum[idx] * normFactor;
    weights[i] = 1.0 / (ei * ei);
    y[i] = flux[idx] * normFactor;
  }

  s = gsl_multifit_fdfsolver_alloc(T, n, p);

  /* initialize solver with starting point and weights */
  gsl_multifit_fdfsolver_wset(s, &f, &x.vector, &w.vector);

  /* compute initial residual norm */
  res_f = gsl_multifit_fdfsolver_residual(s);
  chi0 = gsl_blas_dnrm2(res_f);

  /* solve the system with a maximum of maxIterations iterations */
  status =
      gsl_multifit_fdfsolver_driver(s, maxIterations, xtol, gtol, ftol, &info);

  gsl_multifit_fdfsolver_jac(s, J);
  gsl_multifit_covar(J, 0.0, covar);

  /* compute final residual norm */
  chi = gsl_blas_dnrm2(res_f);

  if (verbose) {
#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar, i, i))

    Log.LogDebug("summary from method '%s'", gsl_multifit_fdfsolver_name(s));
    Log.LogDebug("number of iterations: %zu", gsl_multifit_fdfsolver_niter(s));
    Log.LogDebug("function evaluations: %zu", f.nevalf);
    Log.LogDebug("Jacobian evaluations: %zu", f.nevaldf);
    Log.LogDebug("reason for stopping: %s", (info == 1) ? "small step size"
                                            : (info == 2)
                                                ? "small gradient"
                                                : "small change in f");
    Log.LogDebug("initial |f(x)| = %g", chi0);
    Log.LogDebug("final   |f(x)| = %g", chi);

    {
      double dof = n - p;
      // double c = GSL_MAX_DBL(1, chi / sqrt(dof));

      Log.LogDebug("chisq/dof = %g", pow(chi, 2.0) / dof);

      for (Int32 k = 0; k < p; k++) {
        if (FIT(k) < 1e-3) {
          // Log.LogDebug("A %d     = %.3e +/- %.8f", k, FIT(k), c*ERR(k));
        } else {
          // Log.LogDebug("A %d     = %.5f +/- %.8f", k, FIT(k), c*ERR(k));
        }
      }
    }
    Log.LogDebug("status = %s (%d)", gsl_strerror(status), status);
  }

  //=======================Result treament ====================
  Int32 ampPos = 1;
  controller->resizeAmpsLine();
  Float64 dof = n - p;
  Float64 c = GSL_MAX_DBL(1, chi / sqrt(dof));
  for (Int32 iddl = 0; iddl < nddl; iddl++) {
    Float64 a = controller->lineAmp_LmToModel(
        gsl_vector_get(s->x, iddl)); /// normFactor;
    Float64 sigma = controller->lineAmp_LmToModel(
        c * sqrt(gsl_matrix_get(covar, iddl, iddl))); /// normFactor;
    if (a < 0.0) {
      ampPos = 0;
    }
    controller->setAmpLine(iddl, a, sigma);
    if (verbose) {
      Log.LogDetail("LineModel LMfit: set final amp [%] = %f with err = %f",
                    iddl, a, sigma);
    }
  }
  //
  // // Analyse if result is good or not
  // if(controller->isContinuumFitted()){
  //   Float64 continuumTplAmp =
  //   sqrt(gsl_vector_get(s->x,controller->getIndContinuumAmp()));///normFactor;
  //   if(continuumTplAmp<0){
  //     ampPos =0;
  //   }
  // }
  //
  // if(controller->isEmissionVelocityFitted()){
  //   Float64 vel =
  //   gsl_vector_get(s->x,controller->getIndEmissionVel())/controller->getNormEmiFactor();
  //   if(vel<0.0){
  //     ampPos =0;
  //   }
  // }
  // if(controller->isAbsorptionVelocityFitted()){
  //   Float64 vel =
  //   gsl_vector_get(s->x,controller->getIndAbsorptionVel())/controller->getNormAbsFactor();
  //   if(vel<0.0){
  //     ampPos =0;
  //   }
  // }

  Int32 outputStatus = 0;

  // finally populate the fitting results to the output
  if (ampPos && status == 0) {
    controller->setMerit(chi / normFactor / normFactor);
    if (verbose) {
      Log.LogDetail("LineModel LMfit: Result accepted");
    }
    if (controller->isContinuumFitted()) {
      Int32 continuumId = controller->getIndContinuumAmp();
      Float64 continuumTplAmp = controller->continuumAmp_LmToModel(
          gsl_vector_get(s->x, continuumId)); /// normFactor;
      Float64 continuumTplAmpErr = controller->continuumAmp_LmToModel(
          c * sqrt(gsl_matrix_get(covar, continuumId,
                                  continuumId))); /// normFactor;
      controller->setContinummAmp(continuumTplAmp, continuumTplAmpErr);
      if (verbose) {
        Log.LogDetail(
            "LineModel LMfit: continuumTplAmp calculated =%.3f err =%.3f",
            continuumTplAmp, continuumTplAmpErr);
      }
    }

    if (controller->isEmissionVelocityFitted()) {
      Int32 id = controller->getIndEmissionVel();
      Float64 vel = controller->emiVel_LmToModel(gsl_vector_get(s->x, id));
      Float64 errVel =
          controller->emiVel_LmToModel(c * sqrt(gsl_matrix_get(covar, id, id)));
      controller->setVelocityEmission(vel, errVel);
      if (verbose) {
        Log.LogDetail(
            "LineModel LMfit: Emission velocity found = %.3f with err = %.3f",
            vel, errVel);
      }
    }
    // TODO ajouter des condition sur les err?
    if (controller->isAbsorptionVelocityFitted()) {
      Int32 id = controller->getIndAbsorptionVel();
      Float64 vel = controller->absVel_LmToModel(gsl_vector_get(s->x, id));
      Float64 errVel =
          controller->absVel_LmToModel(c * sqrt(gsl_matrix_get(covar, id, id)));
      controller->setVelocityAbsorption(vel, errVel);
      if (verbose) {
        Log.LogDetail(
            "LineModel LMfit: Absorption velocity found = %.3f with err = %.3f",
            vel, errVel);
      }
    }

    if (controller->isRedshiftFitted()) {
      Int32 id = controller->getIndRedshift();
      Float64 vel = gsl_vector_get(s->x, id);
      Float64 errVel = c * sqrt(gsl_matrix_get(covar, id, id));
      controller->setRedshift(vel, errVel);
      if (verbose) {
        Log.LogDetail("LineModel LMfit: Redshift found = %f with err = %f", vel,
                      errVel);
      }
    }

    outputStatus = 0;
  } else {
    Log.LogDetail("LineModel LMfit: Result dump");
    controller->setMerit(INFINITY);

    if (controller->isEmissionVelocityFitted()) {
      controller->setVelocityEmission(m_velocityEmissionInit, 0);
    }

    if (controller->isAbsorptionVelocityFitted()) {
      controller->setVelocityEmission(m_velocityAbsorptionInit, 0);
    }
    if (controller->isContinuumFitted()) {
      controller->setContinummAmp(amplitudeContinuumInit, 0);
    }
    outputStatus = 1;
  }

  gsl_multifit_fdfsolver_free(s);
  gsl_matrix_free(covar);
  gsl_matrix_free(J);
  // gsl_rng_free (r);

  return outputStatus;
}

Float64 CLineModelFitting::GetWeightingAnyLineCenterProximity(
    Int32 sampleIndex, const TInt32List &EltsIdx) const {
  Float64 maxWeight = 0.0;
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  Float64 currentLbda = spectralAxis[sampleIndex];

  for (Int32 i = 0; i < EltsIdx.size(); i++) {
    Int32 iElts = EltsIdx[i];

    // getTheoreticalSupport reads from the 2 class variables
    TInt32RangeList s = m_Elements[iElts]->getTheoreticalSupport();
    for (Int32 iS = 0; iS < s.size(); iS++) {
      Float64 weight = 0;
      TInt32Range range = s[iS];
      if (sampleIndex > range.GetBegin() && sampleIndex < range.GetEnd()) {
        Float64 minRangeLbda = spectralAxis[range.GetBegin()];
        Float64 maxRangeLbda = spectralAxis[range.GetEnd()];

        Float64 fullInterval = maxRangeLbda - minRangeLbda;
        if (fullInterval > 0.0) {
          Float64 leftInterval = currentLbda - minRangeLbda;
          Float64 rightInterval = maxRangeLbda - currentLbda;
          Float64 coeff = std::abs(leftInterval - rightInterval);

          weight = 1. - (Float64)coeff / (Float64)fullInterval;
        }
      }

      if (maxWeight < weight) {
        maxWeight = weight;
      }
    }
  }

  return maxWeight;
}

/**
 * \brief Returns a sorted set of line indices present in the supports of the
 *argument. Create a vector named indexes. If the argument ind is an index to
 *m_Elements that IsOutSideLambdaRange, return indexes. For each entry in
 *m_Elements: If the entry has a different linetype than the line corresponding
 *to ind, go to the next entry. If the entry IsOutsideLambdaRange, go to the
 *enxt entry. For each subentry in the support of entry: For each subsubentry in
 *the support of ind: If the overlap in the spectralAxis is smaller than -1 *
 *overlapThres * winsize, add entry to indexes. Sort indexes, remove duplicates
 *from indexes, and return indexes.
 **/
TInt32List
CLineModelFitting::getOverlappingElementsBySupport(Int32 ind,
                                                   Float64 overlapThres) const {
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
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

        // Log.LogInfo( "hybrid fit: iRefS=%d - support=%d,%d", iRefS, x1, x2);
        // Log.LogInfo( "hybrid fit: iS=%d - support=%d,%d", iS, y1, y2);

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

Int32 CLineModelFitting::fitAmplitudesLinSolveAndLambdaOffset(
    TInt32List EltsIdx, const CSpectrumSpectralAxis &spectralAxis,
    const CSpectrumFluxAxis &fluxAxis,
    const CSpectrumFluxAxis &continuumfluxAxis, TFloat64List &ampsfitted,
    TFloat64List &errorsfitted, bool enableOffsetFitting) {
  Int32 ret = -1;
  Int32 nSteps =
      int((m_LambdaOffsetMax - m_LambdaOffsetMin) / m_LambdaOffsetStep + 0.5);

  bool atLeastOneOffsetToFit = false;
  if (enableOffsetFitting) {
    for (Int32 iE = 0; iE < EltsIdx.size(); iE++) {
      Float64 nLines = m_Elements[iE]->m_Lines.size();
      for (Int32 iR = 0; iR < nLines; iR++) {
        // check if the line is to be fitted
        if (m_Elements[iE]->m_Lines[iR].GetOffsetFitEnabled()) {
          atLeastOneOffsetToFit = true;
          break;
        }
      }
    }
  }

  if (!atLeastOneOffsetToFit) {
    nSteps = 1;
  }

  Float64 bestMerit = DBL_MAX;
  Int32 idxBestMerit = -1;
  for (Int32 iO = 0; iO < nSteps; iO++) {
    // set offset value
    if (atLeastOneOffsetToFit) {
      Float64 offset = m_LambdaOffsetMin + m_LambdaOffsetStep * iO;
      for (Int32 iE = 0; iE < EltsIdx.size(); iE++) {
        Float64 nLines = m_Elements[iE]->m_Lines.size();

        for (Int32 iR = 0; iR < nLines; iR++) {
          if (m_Elements[iE]->m_Lines[iR].GetOffsetFitEnabled()) {
            m_Elements[iE]->m_Lines[iR].SetOffset(offset);
          }
        }
      }
    }

    // fit for this offset
    ret = fitAmplitudesLinSolve(EltsIdx, spectralAxis, fluxAxis,
                                continuumfluxAxis, ampsfitted, errorsfitted);

    // check fitting
    if (atLeastOneOffsetToFit) {
      Float64 sumFit = 0.0;
      refreshModelUnderElements(EltsIdx);

      // todo: replace lambdarange using elements limits for speed
      for (Int32 iE = 0; iE < EltsIdx.size(); iE++) {
        Float64 _fit = m_Elements.getModelErrorUnderElement(
            iE, fluxAxis, m_SpectrumModel.GetFluxAxis());
        sumFit += _fit;
      }
      if (sumFit < bestMerit) {
        bestMerit = sumFit;
        idxBestMerit = iO;
      }
    }
  }

  if (idxBestMerit >= 0 && atLeastOneOffsetToFit) {
    // set offset value
    if (atLeastOneOffsetToFit) {
      Float64 offset = m_LambdaOffsetMin + m_LambdaOffsetStep * idxBestMerit;
      for (Int32 iE = 0; iE < EltsIdx.size(); iE++) {
        Int32 nLines = m_Elements[iE]->m_Lines.size();
        for (Int32 iR = 0; iR < nLines; iR++) {
          if (m_Elements[iE]->m_Lines[iR].GetOffsetFitEnabled()) {
            m_Elements[iE]->m_Lines[iR].SetOffset(offset);
          }
        }
      }
    }
    // fit again for this offset
    ret = fitAmplitudesLinSolve(EltsIdx, spectralAxis, fluxAxis,
                                continuumfluxAxis, ampsfitted, errorsfitted);
  }

  return ret;
}

/**
 * \brief Use GSL to fit linearly the elements listed in argument EltsIdx.
 * If size of argument EltsIdx is less than 1 return -1.
 **/
Int32 CLineModelFitting::fitAmplitudesLinSolve(
    const TInt32List &EltsIdx, const CSpectrumSpectralAxis &spectralAxis,
    const CSpectrumFluxAxis &fluxAxis,
    const CSpectrumFluxAxis &continuumfluxAxis, TFloat64List &ampsfitted,
    TFloat64List &errorsfitted) {
  // boost::chrono::thread_clock::time_point start_prep =
  // boost::chrono::thread_clock::now();

  const bool verbose = false;
  bool useAmpOffset = m_enableAmplitudeOffsets;
  Int32 idxAmpOffset = -1;

  Int32 idx = 0;

  Int32 nddl = EltsIdx.size();
  if (nddl < 1) {
    return -1;
  }
  TInt32List xInds = m_Elements.getSupportIndexes(EltsIdx);
  if (xInds.size() < 1) {
    return -1;
  }

  if (useAmpOffset) {
    nddl += m_AmplitudeOffsetsDegree + 1;
    // find the amplitudeOffset Support that corresponds to these elts
    idxAmpOffset = m_Elements.getIndexAmpOffset(xInds[0]);
  }

  for (Int32 iddl = 0; iddl < EltsIdx.size(); iddl++) {
    m_Elements.SetElementAmplitude(EltsIdx[iddl], 1.0, 0.0);
  }

  const Float64 *spectral = spectralAxis.GetSamples();
  const Float64 *flux = fluxAxis.GetSamples();
  const auto &ErrorNoContinuum = m_ErrorNoContinuum;

  // Linear fit
  int i, n;
  Float64 fval;
  double chisq;
  gsl_matrix *X, *cov;
  gsl_vector *y, *w, *c;

  n = xInds.size();
  if (n < nddl) {
    ampsfitted.resize(nddl);
    errorsfitted.resize(nddl);
    for (Int32 iddl = 0; iddl < nddl; iddl++) {
      ampsfitted[iddl] = 0.0;
      errorsfitted[iddl] = 1e12; // some high number
    }
    return -1;
  }

  X = gsl_matrix_alloc(n, nddl);
  y = gsl_vector_alloc(n);
  w = gsl_vector_alloc(n);
  c = gsl_vector_alloc(nddl);
  cov = gsl_matrix_alloc(nddl, nddl);

  // Normalize
  Float64 maxabsval = DBL_MIN;
  for (i = 0; i < n; i++) {
    idx = xInds[i];
    if (maxabsval < std::abs(flux[idx])) {
      maxabsval = std::abs(flux[idx]);
    }
  }
  Float64 normFactor = 1.0 / maxabsval;
  if (verbose) {
    Log.LogDetail("normFactor = '%.3e'\n", normFactor);
  }

  // Prepare the fit data
  for (i = 0; i < n; i++) {
    double xi, yi, ei;
    idx = xInds[i];
    xi = spectral[idx];
    yi = flux[idx] * normFactor;
    ei = ErrorNoContinuum[idx] * normFactor;

    gsl_vector_set(y, i, yi);
    gsl_vector_set(w, i, 1.0 / (ei * ei));

    for (Int32 iddl = 0; iddl < EltsIdx.size(); iddl++) {
      fval = m_Elements[EltsIdx[iddl]]->getModelAtLambda(
          xi, m_Redshift, continuumfluxAxis[idx]);
      gsl_matrix_set(X, i, iddl, fval);

      if (verbose) {
        Log.LogDebug("fval = '%.3e'", fval);
      }
    }

    if (useAmpOffset) {
      gsl_matrix_set(X, i, EltsIdx.size(), 1.0);
      if (m_AmplitudeOffsetsDegree == 0)
        continue;
      gsl_matrix_set(X, i, EltsIdx.size() + 1, xi);
      if (m_AmplitudeOffsetsDegree == 1)
        continue;
      gsl_matrix_set(X, i, EltsIdx.size() + 2, xi * xi);
    }
  }

  {
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n, nddl);
    gsl_multifit_wlinear(X, w, y, c, cov, &chisq, work);
    gsl_multifit_linear_free(work);
  }

  if (verbose) {
#define C(i) (gsl_vector_get(c, (i)))
#define COV(i, j) (gsl_matrix_get(cov, (i), (j)))
    if (1) {
      Log.LogDebug("# best fit: Y = %g X1 + %g X2 ...", C(0), C(1));
      Log.LogDebug("# covariance matrix:");
      Log.LogDebug("[");
      Log.LogDebug("  %+.5e, %+.5e", COV(0, 0), COV(0, 1));
      Log.LogDebug("  %+.5e, %+.5e", COV(1, 0), COV(1, 1));

      //        Log.LogDebug("[ %+.5e, %+.5e, %+.5e  \n", COV(0,0), COV(0,1),
      //        COV(0,2)); Log.LogDebug("  %+.5e, %+.5e, %+.5e  \n", COV(1,0),
      //        COV(1,1), COV(1,2)); Log.LogDebug("  %+.5e, %+.5e, %+.5e ]\n",
      //        COV(2,0), COV(2,1), COV(2,2));

      Log.LogDebug("]");
      Log.LogDebug("# chisq/n = %g", chisq / n);
    }

    for (Int32 iddl = 0; iddl < nddl; iddl++) {
      Float64 a = gsl_vector_get(c, iddl) / normFactor;
      if (verbose) {
        Log.LogDetail("# Found amplitude %d: %+.5e", iddl, a);
      }
    }
  }

  Int32 sameSign = 1;
  Float64 a0 = gsl_vector_get(c, 0) / normFactor;
  for (Int32 iddl = 1; iddl < EltsIdx.size(); iddl++) {
    Float64 a = gsl_vector_get(c, iddl) / normFactor;
    Float64 product = a0 * a;
    if (product < 0) {
      sameSign = 0;
    }
  }

  if (verbose) {
    Log.LogDetail("# Found amplitudes with sameSign=%d", sameSign);
  }
  if (sameSign) {
    for (Int32 iddl = 0; iddl < EltsIdx.size(); iddl++) {
      Float64 a = gsl_vector_get(c, iddl) / normFactor;
      Float64 cova = gsl_matrix_get(cov, iddl, iddl);
      Float64 sigma = sqrt(cova) / normFactor;
      m_Elements.SetElementAmplitude(EltsIdx[iddl], a, sigma);
    }
    // refreshModel();
  } else {
    ampsfitted.resize(EltsIdx.size());
    errorsfitted.resize(EltsIdx.size());
    for (Int32 iddl = 0; iddl < EltsIdx.size(); iddl++) {
      Float64 a = gsl_vector_get(c, iddl) / normFactor;
      Float64 cova = gsl_matrix_get(cov, iddl, iddl);
      Float64 sigma = sqrt(cova) / normFactor;
      m_Elements.SetElementAmplitude(EltsIdx[iddl], a, sigma);
      ampsfitted[iddl] = (a);
      errorsfitted[iddl] = (sigma);
    }
  }

  if (useAmpOffset) {
    Float64 x0 = gsl_vector_get(c, EltsIdx.size()) / normFactor;
    Float64 x1 = 0.0;
    Float64 x2 = 0.0;
    if (m_AmplitudeOffsetsDegree > 0)
      x1 = gsl_vector_get(c, EltsIdx.size() + 1) / normFactor;
    if (m_AmplitudeOffsetsDegree > 1)
      x2 = gsl_vector_get(c, EltsIdx.size() + 2) / normFactor;
    m_Elements.setAmplitudeOffsetsCoeffsAt(idxAmpOffset, {x0, x1, x2});
  }

  gsl_matrix_free(X);
  gsl_vector_free(y);
  gsl_vector_free(w);
  gsl_vector_free(c);
  gsl_matrix_free(cov);

  return sameSign;
}

/**
 * @brief CLineModelFitting::fitAmplitudesLinesAndContinuumLinSolve
 * @param EltsIdx:  elements to be fitted
 * @param lambdaRange
 * @param spectralAxis
 * @param fluxAxis: must contain the observed spectrum
 * @param continuumfluxAxis: must contain the continuum to be amp-fitted
 * @param ampsfitted: output
 * @param errorsfitted: output
 * @param polyOrder: order of the polynom to be fitted along with the continuum
 * (-1 = disabled)
 * @return
 *
 * WARNING: not sure about fitting abs. lines with this method...
 */
Int32 CLineModelFitting::fitAmplitudesLinesAndContinuumLinSolve(
    const TInt32List &EltsIdx, const CSpectrumSpectralAxis &spectralAxis,
    const CSpectrumFluxAxis &fluxAxis,
    const CSpectrumFluxAxis &continuumfluxAxis, TFloat64List &ampsfitted,
    TFloat64List &errorsfitted, Float64 &chisquare, Int32 polyOrder) {
  // boost::chrono::thread_clock::time_point start_prep =
  // boost::chrono::thread_clock::now();

  bool verbose = false;
  Int32 idx = 0;

  Int32 nddl =
      EltsIdx.size() + 1; // number of param to be fitted=nlines+continuum
  if (nddl < 2) {
    return -1;
  }

  if (polyOrder >= 0) {
    nddl += polyOrder + 1;
  }

  TInt32List xInds;
  for (Int32 i = 0; i < spectralAxis.GetSamplesCount(); i++) {
    if (spectralAxis[i] >= m_lambdaRange.GetBegin() &&
        spectralAxis[i] <= m_lambdaRange.GetEnd()) {
      xInds.push_back(i);
    }
  }
  if (xInds.size() < 1) {
    return -1;
  }

  for (Int32 iddl = 0; iddl < EltsIdx.size(); iddl++) {
    m_Elements.SetElementAmplitude(EltsIdx[iddl], 1.0, 0.0);
  }

  const Float64 *spectral = spectralAxis.GetSamples();
  const Float64 *flux = fluxAxis.GetSamples();
  const auto &ErrorNoContinuum = m_ErrorNoContinuum;

  // Linear fit
  int i, n;
  Float64 fval;
  double chisq;
  gsl_matrix *X, *cov;
  gsl_vector *y, *w, *c;

  n = xInds.size();
  if (n < nddl) {
    ampsfitted.resize(nddl);
    errorsfitted.resize(nddl);
    for (Int32 iddl = 0; iddl < nddl; iddl++) {
      ampsfitted[iddl] = 0.0;
      errorsfitted[iddl] = 1e12; // some high number
    }
    return -1;
  }

  X = gsl_matrix_alloc(n, nddl);
  y = gsl_vector_alloc(n);
  w = gsl_vector_alloc(n);
  c = gsl_vector_alloc(nddl);
  cov = gsl_matrix_alloc(nddl, nddl);

  // Normalize
  Float64 maxabsval = DBL_MIN;
  for (i = 0; i < n; i++) {
    idx = xInds[i];
    if (maxabsval < std::abs(flux[idx])) {
      maxabsval = std::abs(flux[idx]);
    }
  }
  Float64 normFactor = 1.0 / maxabsval;
  if (verbose) {
    Log.LogDetail("normFactor = '%.3e'\n", normFactor);
  }

  // Prepare the fit data
  for (i = 0; i < n; i++) {
    double xi, yi, ei, ci;
    idx = xInds[i];
    xi = spectral[idx];
    yi = flux[idx] * normFactor;
    ei = ErrorNoContinuum[idx] * normFactor;
    ci = continuumfluxAxis[idx];

    for (Int32 iddl = 0; iddl < EltsIdx.size(); iddl++) {
      fval = m_Elements[EltsIdx[iddl]]->getModelAtLambda(xi, m_Redshift, ci);
      gsl_matrix_set(X, i, iddl, fval);

      if (false && verbose) {
        Log.LogDebug("fval = '%.3e'", fval);
      }
    }

    gsl_matrix_set(X, i, EltsIdx.size(), ci);

    if (polyOrder >= 0) {
      for (Int32 kCoeff = 0; kCoeff < polyOrder + 1; kCoeff++) {
        Float64 vect = 1.0;
        for (Int32 kt = 0; kt < kCoeff; kt++) {
          vect *= xi;
        }
        gsl_matrix_set(X, i, EltsIdx.size() + 1 + kCoeff, vect);
      }
    }

    gsl_vector_set(y, i, yi);
    gsl_vector_set(w, i, 1.0 / (ei * ei));
  }

  //
  // boost::chrono::thread_clock::time_point stop_prep =
  // boost::chrono::thread_clock::now();
  // Float64 duration_prep =
  // boost::chrono::duration_cast<boost::chrono::microseconds>(stop_prep -
  // start_prep).count();
  // boost::chrono::thread_clock::time_point start_fit =
  // boost::chrono::thread_clock::now();

  gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n, nddl);
  gsl_multifit_wlinear(X, w, y, c, cov, &chisq, work);
  gsl_multifit_linear_free(work);

  //
  // boost::chrono::thread_clock::time_point stop_fit =
  // boost::chrono::thread_clock::now(); Float64 duration_fit =
  // boost::chrono::duration_cast<boost::chrono::microseconds>(stop_fit -
  // start_fit).count(); Log.LogInfo("LineModel linear fit: prep = %.3f - fit =
  // %.3f", duration_prep, duration_fit);

  if (verbose) {
#define C(i) (gsl_vector_get(c, (i)))
#define COV(i, j) (gsl_matrix_get(cov, (i), (j)))
    if (1) {
      Log.LogDebug("# best fit: Y = %g X1 + %g X2 ...", C(0), C(1));
      Log.LogDebug("# covariance matrix:");
      Log.LogDebug("[");
      Log.LogDebug("  %+.5e, %+.5e", COV(0, 0), COV(0, 1));
      Log.LogDebug("  %+.5e, %+.5e", COV(1, 0), COV(1, 1));

      //        Log.LogDebug("[ %+.5e, %+.5e, %+.5e  \n", COV(0,0), COV(0,1),
      //        COV(0,2)); Log.LogDebug("  %+.5e, %+.5e, %+.5e  \n", COV(1,0),
      //        COV(1,1), COV(1,2)); Log.LogDebug("  %+.5e, %+.5e, %+.5e ]\n",
      //        COV(2,0), COV(2,1), COV(2,2));

      Log.LogDebug("]");
      Log.LogDebug("# chisq = %g", chisq);
      Log.LogDebug("# chisq/n = %g", chisq / n);
    }

    for (Int32 iddl = 0; iddl < nddl; iddl++) {
      Float64 a = gsl_vector_get(c, iddl) / normFactor;
      if (verbose) {
        Log.LogDetail("# Found amplitude %d: %+.5e", iddl, a);
      }
    }
  }

  Int32 sameSign = 1;
  Float64 a0 = gsl_vector_get(c, 0) / normFactor;
  for (Int32 iddl = 1; iddl < EltsIdx.size(); iddl++) {
    Float64 a = gsl_vector_get(c, iddl) / normFactor;
    Float64 product = a0 * a;
    if (product < 0) {
      sameSign = 0;
    }
  }

  if (verbose) {
    Log.LogDetail("# Found n=%d amplitudes with sameSign=%d", EltsIdx.size(),
                  sameSign);
  }

  ampsfitted.resize(nddl);
  errorsfitted.resize(nddl);
  for (Int32 iddl = 0; iddl < nddl; iddl++) {
    Float64 a = gsl_vector_get(c, iddl) / normFactor;
    Float64 cova = gsl_matrix_get(cov, iddl, iddl);
    Float64 sigma = sqrt(cova) / normFactor;
    if (iddl < EltsIdx.size()) {
      m_Elements.SetElementAmplitude(EltsIdx[iddl], a, sigma);
    }
    ampsfitted[iddl] = (a);
    errorsfitted[iddl] = (sigma);
  }

  if (verbose && polyOrder >= 0) {
    for (Int32 kCoeff = 0; kCoeff < polyOrder + 1; kCoeff++) {
      Float64 p = gsl_vector_get(c, EltsIdx.size() + 1 + kCoeff) / normFactor;
      Log.LogDetail("# Found p%d poly amplitude = %+.5e", kCoeff, p);
    }
  }

  if (verbose) {
    Log.LogDetail("# Returning (L+C) n=%d amplitudes", ampsfitted.size());
  }

  gsl_matrix_free(X);
  gsl_vector_free(y);
  gsl_vector_free(w);
  gsl_vector_free(c);
  gsl_matrix_free(cov);

  chisquare = chisq;
  return sameSign;
}

/**
 * @brief CLineModelFitting::setLyaProfile
 * If a Lya line is present with SYMIGM profile, fit igmIdx
 * If a Lya line is present with ASYMFIT profile, fit the width and asymmetry
 * parameters If a Lya line is present with ASYMFIXED profile, set the width and
 * asymmetry parameters according to profile parameters given in the string
 * @param redshift, catalog (could be tplshape or linecatalog)
 * @return 1 if successfully fitted, 0 if error, 2 if Lya not present, 3 if Lya
 * not configured to be fitted in the catalog
 */
Int32 CLineModelFitting::setLyaProfile(Float64 redshift,
                                       const CLineCatalog::TLineVector &catalog,
                                       bool tplratio)

{

  // 1. retrieve the Lya index
  std::string lyaTag = linetags::lya_em;

  Int32 idxLineLyaE = undefIdx;
  Int32 idxLyaE = m_Elements.findElementIndex(lyaTag, -1, idxLineLyaE);
  if (idxLyaE < 0 || idxLineLyaE < 0) {
    return 2; // Lya alpha not found
  }

  Int32 lineIndex = undefIdx;
  // get index of lya inside tplshape catalog
  if (tplratio) {
    lineIndex = std::find_if(catalog.begin(), catalog.end(),
                             [lyaTag](const CLine &line) {
                               return line.GetName() == lyaTag;
                             }) -
                catalog.begin();
    if (lineIndex == undefIdx || lineIndex > catalog.size() - 1)
      return 2; // Lya alpha not found in tplshape catalog
  } else {
    lineIndex = m_Elements[idxLyaE]->m_LineCatalogIndexes[idxLineLyaE];
  }

  if (lineIndex < 0 || lineIndex > catalog.size() - 1)
    THROWG(INTERNAL_ERROR, "Lye idx out-of-bound");

  // finding or setting the correct profile
  CLineProfile_ptr profile;
  if (m_forceLyaFitting && catalog[lineIndex].GetProfile().isAsymFixed())
    profile = std::unique_ptr<CLineProfileASYMFIT>(new CLineProfileASYMFIT(
        m_NSigmaSupport, catalog[lineIndex].GetAsymParams(), "mean"));
  else
    profile = catalog[lineIndex].GetProfile().Clone();

  if (m_Elements[idxLyaE]->IsOutsideLambdaRange())
    return 3;

  bool doasymfit = profile->isAsymFit();
  bool fixedIGM =
      m_ContinuumComponent == "tplfit" || m_ContinuumComponent == "tplfitauto";
  bool doIGMfit = profile->GetName() == SYMIGM && !fixedIGM;

  if (profile->GetName() == SYMIGM && fixedIGM)
    profile->SetSymIgmParams(
        TSymIgmParams(m_fitContinuum_tplFitMeiksinIdx, redshift));

  m_Elements[idxLyaE]->m_Lines[idxLineLyaE].SetProfile(std::move(profile));

  if (!doasymfit && !doIGMfit)
    return 3;

  TInt32List filterEltsIdxLya;
  filterEltsIdxLya.push_back(idxLyaE);
  if (doasymfit) {
    // 3. find the best width and asym coeff. parameters
    TAsymParams bestfitParams =
        fitAsymParameters(redshift, idxLyaE, filterEltsIdxLya, idxLineLyaE);
    // 4. set the associated Lya members in the element definition
    m_Elements[idxLyaE]->SetAsymfitParams(bestfitParams);
  }

  if (doIGMfit) {
    Int32 bestigmidx =
        fitAsymIGMCorrection(redshift, idxLyaE, filterEltsIdxLya, idxLineLyaE);
    m_Elements[idxLyaE]->SetSymIgmParams(TSymIgmParams(bestigmidx, redshift));
  }
  return 1;
}

Int32 CLineModelFitting::fitAsymIGMCorrection(
    const Float64 &redshift, const Int32 &idxLyaE,
    const TInt32List &filterEltsIdxLya, const Int32 &idxLineLyaE) {

  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  // TODO: to check with Didier
  if (spectralAxis[0] / (1 + redshift) > RESTLAMBDA_LYA)
    return -1;

  Float64 meritMin = DBL_MAX;
  Int32 bestIgmIdx = -1;

  Int32 igmCount =
      m_Elements[idxLyaE]->getLineProfile(idxLineLyaE).getIGMIdxCount();
  for (Int32 igmIdx = 0; igmIdx < igmCount; igmIdx++) {
    m_Elements[idxLyaE]->SetSymIgmParams(TSymIgmParams(igmIdx, redshift));
    m_Elements[idxLyaE]->fitAmplitude(spectralAxis, m_spcFluxAxisNoContinuum,
                                      m_ContinuumFluxAxis, redshift,
                                      idxLineLyaE);

    refreshModelUnderElements(filterEltsIdxLya, idxLineLyaE);
    Float64 m = m_Elements.getModelErrorUnderElement(
        idxLyaE, m_SpcFluxAxis, m_SpectrumModel.GetFluxAxis());

    if (m < meritMin) {
      meritMin = m;
      bestIgmIdx = igmIdx;
    }
  }
  return bestIgmIdx;
}

TAsymParams CLineModelFitting::fitAsymParameters(
    const Float64 &redshift, const Int32 &idxLyaE,
    const TInt32List &filterEltsIdxLya, const Int32 &idxLineLyaE) {

  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
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

  bool verbose = false;
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
            spectralAxis, m_spcFluxAxisNoContinuum, m_ContinuumFluxAxis,
            redshift, idxLineLyaE);

        Float64 m = m_dTransposeD;
        if (1) {
          refreshModelUnderElements(filterEltsIdxLya, idxLineLyaE);
          m = m_Elements.getModelErrorUnderElement(
              idxLyaE, m_SpcFluxAxis, m_SpectrumModel.GetFluxAxis());
        } else {
          m = getLeastSquareMeritFast(idxLyaE);
        }
        if (m < meritMin) {
          meritMin = m;
          bestparams = m_Elements[idxLyaE]->GetAsymfitParams(0);
        }

        if (verbose) {
          Log.LogInfo("Fitting Lya Profile: width=%f, asym=%f, delta=%f",
                      asymWidthCoeff, asymAlphaCoeff, delta);
          Log.LogInfo("Fitting Lya Profile: merit=%e", m);
          Log.LogInfo("Fitting Lya Profile: idxLyaE=%d, idxLineLyaE=%d",
                      idxLyaE, idxLineLyaE);
        }
      }
    }
  }
  if (verbose) {
    Log.LogInfo("Lya Profile found: width=%f, asym=%f, delta=%f",
                bestparams.sigma, bestparams.alpha, bestparams.delta);
  }
  return bestparams;
}

/**
 * @brief CLineModelFitting::ReestimateContinuumUnderLines
 * For each line, reestimate the continuum using the original median routines on
 * a sub segment of the original spectrum
 * - the subsegment is a contiguous segment around the support of all the
 * elements
 * @param EltsIdx
 * @return
 */
TInt32List
CLineModelFitting::ReestimateContinuumUnderLines(const TInt32List &EltsIdx) {
  if (EltsIdx.size() < 1) {
    TInt32List empty;
    return empty;
  }

  // smoothing factor in continuum median filter
  Float64 smoof = 150;

  // modify m_ContinuumFluxAxis
  CSpectrumFluxAxis &fluxAxisModified = m_ContinuumFluxAxis;
  Float64 *Ycont = fluxAxisModified.GetSamples();
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();

  TInt32List xInds = m_Elements.getSupportIndexes(EltsIdx);
  Int32 minInd = xInds[0];
  Int32 maxInd = xInds[xInds.size() - 1];
  //
  Int32 iminMerge =
      spectralAxis.GetIndexAtWaveLength(spectralAxis[minInd] - 2 * smoof);
  if (iminMerge == 0) {
    iminMerge = 1;
  }
  Int32 imaxMerge =
      spectralAxis.GetIndexAtWaveLength(spectralAxis[maxInd] + 2 * smoof);
  //
  //
  Int32 iminMerge2 =
      spectralAxis.GetIndexAtWaveLength(spectralAxis[minInd] - smoof);
  if (iminMerge2 == 0) {
    iminMerge2 = 1;
  }
  Int32 imaxMerge2 =
      spectralAxis.GetIndexAtWaveLength(spectralAxis[maxInd] + smoof);
  //
  Int32 imin =
      spectralAxis.GetIndexAtWaveLength(spectralAxis[minInd] - 3 * smoof);
  if (imin == 0) {
    imin = 1;
  }
  Int32 imax =
      spectralAxis.GetIndexAtWaveLength(spectralAxis[maxInd] + 3 * smoof);
  Int32 spcSize = imax - imin + 1;

  // prepare the line model with no continuum;
  CSpectrumFluxAxis modelFluxAxisTmp = m_SpectrumModel.GetFluxAxis();
  for (Int32 i = 0; i < spcSize; i++) {
    modelFluxAxisTmp[i + imin] = 0.0;
  }
  for (Int32 idx = 0; idx < EltsIdx.size(); idx++) {
    Int32 eltIdx = EltsIdx[idx];
    m_Elements[eltIdx]->addToSpectrumModel(spectralAxis, modelFluxAxisTmp,
                                           m_ContinuumFluxAxis, m_Redshift);
  }

  // gather the modified indexes
  TInt32List modifiedIdxs;

  // create the spcBuffer with spectrum minus the lines
  CSpectrum spcBuffer;
  CSpectrumSpectralAxis _SpectralAxis = CSpectrumSpectralAxis(spcSize, false);
  CSpectrumFluxAxis _FluxAxis = CSpectrumFluxAxis(spcSize);
  const auto &SpcFluxAxis = m_SpcFluxAxis;

  for (Int32 i = 0; i < spcSize; i++) {
    (_SpectralAxis)[i] = spectralAxis[i + imin];
    (_FluxAxis)[i] =
        SpcFluxAxis[i + imin] - as_const(modelFluxAxisTmp)[i + imin];
    //                if( error!= NULL ){
    //                    (*_FluxAxis).GetError()[i] = tmpError[i];
    //                }
  }
  spcBuffer.SetSpectralAndFluxAxes(std::move(_SpectralAxis),
                                   std::move(_FluxAxis));

  /*
  // export for debug
  FILE* fspc = fopen( "ReestimateContinuumUnderLines_correctedSpc_dbg.txt", "w+"
  ); Float64 coeffSaveSpc = 1e16; for( Int32
  t=0;t<spcBuffer.GetSampleCount();t++)
  {
      fprintf( fspc, "%f %f %f\n", t, spcBuffer.GetSpectralAxis()[t],
  (m_SpcFluxAxis[t+imin])*coeffSaveSpc,
  (spcBuffer.GetFluxAxis()[t])*coeffSaveSpc);
  }
  fclose( fspc );
  //*/

  // apply continuum routine on this spcbuffer
  CContinuumIrregularSamplingMedian continuum;
  CSpectrumFluxAxis fluxAxisWithoutContinuumCalc;
  continuum.RemoveContinuum(spcBuffer, fluxAxisWithoutContinuumCalc);

  /*
  // export for debug
  FILE* f = fopen( "continuum_reestimated_underlines_dbg.txt", "w+" );
  Float64 coeffSave = 1e16;
  for( Int32 t=iminMerge;t<imaxMerge;t++)
  {
      fprintf( f, "%f %f %f\n", t, spectralAxis[t],
  (m_inputSpc->GetContinuumFluxAxis()[t])*coeffSave,
  (spcBuffer.GetFluxAxis()[t-imin] -
  fluxAxisWithoutContinuumCalc[t-imin])*coeffSave);
  }
  fclose( f );
  //*/

  Float64 modified = 0.0;
  Float64 coeff = 0.0;
  // merge raw continuum free with the newly calculated cont. under the line,
  // (todo: with cross-fade on the borders)
  for (Int32 i = iminMerge; i < imaxMerge; i++) {
    modified = spcBuffer.GetFluxAxis()[i - imin] -
               as_const(fluxAxisWithoutContinuumCalc)[i - imin];
    coeff = 1.0;

    if (i <= iminMerge2) {
      coeff = (Float64(i - iminMerge) / Float64(iminMerge2 - iminMerge));
    } else if (i >= imaxMerge2) {
      coeff = 1.0 - (Float64(i - imaxMerge2) / Float64(imaxMerge - imaxMerge2));
    }

    Ycont[i] = coeff * modified + (1 - coeff) * as_const(Ycont)[i];
    modifiedIdxs.push_back(i);
  }

  std::sort(modifiedIdxs.begin(), modifiedIdxs.end());
  modifiedIdxs.erase(std::unique(modifiedIdxs.begin(), modifiedIdxs.end()),
                     modifiedIdxs.end());
  return modifiedIdxs;
}

/**
 * \brief Modifies the fluxAxis where SNR is good enough, and returns the list
 *of modified indexes of fluxAxis.
 **/
TInt32List
CLineModelFitting::ReestimateContinuumApprox(const TInt32List &EltsIdx) {
  // smoothing factor in continuum median filter
  Float64 smoof = 150;

  //
  TInt32List modifiedIdxs;

  // modify m_ContinuumFluxAxis
  CSpectrumFluxAxis &fluxAxisModified = m_ContinuumFluxAxis;
  Float64 *Ycont = fluxAxisModified.GetSamples();
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();

  for (Int32 idx = 0; idx < EltsIdx.size(); idx++) {
    Int32 eltIdx = EltsIdx[idx];

    Int32 nlines = m_Elements[eltIdx]->GetSize();
    for (Int32 iline = 0; iline < nlines; iline++) {
      TInt32Range s = m_Elements[eltIdx]->getSupportSubElt(iline);
      Int32 smin = s.GetBegin();
      Int32 smax = s.GetEnd();
      if (smax - smin < 2) {
        continue;
      }
      Int32 imin =
          spectralAxis.GetIndexAtWaveLength(spectralAxis[smin] - smoof);
      if (imin == 0) {
        imin = 1;
      }
      Int32 imax =
          spectralAxis.GetIndexAtWaveLength(spectralAxis[smax] + smoof);
      Int32 sSize = imax - imin + 1;

      Float64 A = m_Elements[eltIdx]->GetFittedAmplitude(iline);
      Float64 Sigma = m_Elements[eltIdx]->GetFittedAmplitudeErrorSigma(iline);

      if (A <= 0 ||
          std::abs(Sigma) > std::abs(A)) { // todo: check this error sigma rule,
                                           // should we add a sigma thres. ?
        continue;
      }
      A *= m_Elements[eltIdx]->GetSignFactor(iline);

      Float64 mu = NAN;
      Float64 sigma = NAN;
      m_Elements[eltIdx]->getObservedPositionAndLineWidth(iline, m_Redshift, mu,
                                                          sigma, false);
      Float64 integratedA = A * sigma * sqrt(2 * M_PI);
      Float64 coeffA = integratedA / (Float64)sSize;

      Float64 term = 0.0;
      for (Int32 i = imin; i < imax; i++) {
        Float64 dx = spectralAxis[imin] - spectralAxis[imin - 1];
        if (i > smin && i < smax) {
          term = coeffA / dx;
        } else if (i > imin && i <= smin) {
          term = ((Float64(i - imin) / Float64(smin - imin))) * coeffA / dx;
        } else if (i >= smax && i < imax) {
          term =
              (1.0 - (Float64(i - smax) / Float64(imax - smax))) * coeffA / dx;
        } else {
          term = 0.0;
        }
        Ycont[i] -= term;

        modifiedIdxs.push_back(i);
      }
    }
  }

  std::sort(modifiedIdxs.begin(), modifiedIdxs.end());
  modifiedIdxs.erase(std::unique(modifiedIdxs.begin(), modifiedIdxs.end()),
                     modifiedIdxs.end());
  return modifiedIdxs;
}

/**
 * \brief Copies the continuum flux to the model flux, and the no continuum flux
 *receives the value of spectrum flux minus the continuum flux.
 **/
void CLineModelFitting::refreshModelAfterContReestimation(
    const TInt32List &xInds, CSpectrumFluxAxis &modelFluxAxis,
    CSpectrumFluxAxis &spcFluxAxisNoContinuum) const {
  Int32 n = xInds.size();
  const auto &ContinuumFluxAxis = m_ContinuumFluxAxis;
  const auto &SpcFluxAxis = m_SpcFluxAxis;
  Int32 idx = 0;
  for (Int32 i = 0; i < n; i++) {
    idx = xInds[i];

    modelFluxAxis[idx] = ContinuumFluxAxis[idx];
    spcFluxAxisNoContinuum[idx] = SpcFluxAxis[idx] - ContinuumFluxAxis[idx];
  }
}

/**
 * \brief Accumulates the squared differences between model and spectrum in the
 *argument lambdaRange and returns the sum.
 **/
Float64 CLineModelFitting::getLeastSquareMerit() const {
  const CSpectrumSpectralAxis &spcSpectralAxis =
      m_SpectrumModel.GetSpectralAxis();
  const CSpectrumFluxAxis &spcFluxAxis = m_SpcFluxAxis;
  const CSpectrumFluxAxis &modelFluxAxis = m_SpectrumModel.GetFluxAxis();

  // Int32 numDevs = 0;
  Float64 fit = 0.0;
  const Float64 *Ymodel = modelFluxAxis.GetSamples();
  const Float64 *Yspc = spcFluxAxis.GetSamples();
  const auto &ErrorNoContinuum = m_ErrorNoContinuum;
  Float64 diff = 0.0;

  Int32 imin = spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange.GetBegin());
  Int32 imax = spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange.GetEnd());
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

  if (m_ContinuumComponent == "tplfit" ||
      m_ContinuumComponent == "tplfitauto") {
    fit += m_fitContinuum_tplFitMerit_phot; // unconditionnal sum (if photometry
                                            // disabled, will sum 0.0)
    fit += m_fitContinuum_tplFitLogprior;
  }

  Log.LogDebug("CLineModelFitting::getLeastSquareMerit fit = %f", fit);
  if (std::isnan(fit)) {
    Log.LogError("CLineModelFitting::getLeastSquareMerit: NaN value found on "
                 "the lambdarange = (%f, %f)",
                 m_lambdaRange.GetBegin(), m_lambdaRange.GetEnd());
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
  const CSpectrumSpectralAxis &spcSpectralAxis =
      m_SpectrumModel.GetSpectralAxis();
  const CSpectrumFluxAxis &spcFluxAxis = m_SpcFluxAxis;
  const auto &ErrorNoContinuum = m_ErrorNoContinuum;

  Int32 numDevs = 0;
  Float64 fit = 0.0;
  const Float64 *YCont = as_const(m_ContinuumFluxAxis).GetSamples();
  const Float64 *Yspc = spcFluxAxis.GetSamples();
  Float64 diff = 0.0;

  Float64 imin = spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange.GetBegin());
  Float64 imax = spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange.GetEnd());
  for (Int32 j = imin; j < imax; j++) {
    numDevs++;
    diff = (Yspc[j] - YCont[j]);
    fit += (diff * diff) / (ErrorNoContinuum[j] * ErrorNoContinuum[j]);
  }

  if (m_ContinuumComponent == "tplfit" ||
      m_ContinuumComponent == "tplfitauto") {
    fit += m_fitContinuum_tplFitLogprior;
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

  if (m_ContinuumComponent == "tplfit" ||
      m_ContinuumComponent == "tplfitauto") {
    Float64 term1 = m_fitContinuum_tplFitAmplitude *
                    m_fitContinuum_tplFitAmplitude * m_fitContinuum_tplFitMtM;
    Float64 term2 =
        -2. * m_fitContinuum_tplFitAmplitude * m_fitContinuum_tplFitDtM;
    fit += term1 + term2;
    fit += m_fitContinuum_tplFitLogprior;
  }

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

Float64 CLineModelFitting::getContinuumScaleMargCorrection() const {
  Float64 corr = 0.0;

  // scale marg for continuum
  if (m_ContinuumComponent == "tplfit" ||
      m_ContinuumComponent ==
          "tplfitauto") // the support has to be already computed when
                        // LoadFitContinuum() is called
  {
    corr += log(m_fitContinuum_tplFitMtM);
  }

  return corr;
}

/**
 * \brief Returns the number of spectral samples between lambdaRange.
 **/
// TODO rename this ! not a simple getter
Int32 CLineModelFitting::getSpcNSamples() const {
  const CSpectrumSpectralAxis &spcSpectralAxis =
      m_SpectrumModel.GetSpectralAxis();

  Int32 numDevs = 0;

  Float64 imin = spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange.GetBegin());
  Float64 imax = spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange.GetEnd());

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
  const CSpectrumFluxAxis &spcFluxAxis = m_SpcFluxAxis;
  const CSpectrumFluxAxis &modelFluxAxis = m_SpectrumModel.GetFluxAxis();
  const CSpectrumNoiseAxis &ErrorNoContinuum = m_ErrorNoContinuum;

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
        subeIdx, m_SpectrumModel.GetSpectralAxis(), m_Redshift,
        m_ContinuumFluxAxis, polynom_coeffs);
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
    TInt32Range support =
        TInt32Range(supportList[k].GetBegin(), supportList[k].GetEnd());

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

    // relevant if contributes to the 'thresRatio'*100 percent (ex. 80%) of the
    // SumSNR value
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

  const CSpectrumFluxAxis &modelFluxAxis = m_SpectrumModel.GetFluxAxis();
  const Float64 *Ymodel = modelFluxAxis.GetSamples();
  const auto &ErrorNoContinuum = m_ErrorNoContinuum;
  const auto &ContinuumFluxAxis = m_ContinuumFluxAxis;

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

Float64 CLineModelFitting::getContinuumMeanUnderElement(Int32 eltId) const {
  Int32 n = 0;
  Float64 m = 0.0;
  // Float64 sumErr=0.0;

  TInt32RangeList support;
  Int32 iElts = eltId;
  {
    if (m_Elements[iElts]->IsOutsideLambdaRange()) {
      return 0.0;
    }
    TInt32RangeList s = m_Elements[iElts]->getSupport();
    for (Int32 iS = 0; iS < s.size(); iS++) {
      support.push_back(s[iS]);
    }
  }

  const auto &ContinuumFluxAxis = m_ContinuumFluxAxis;
  // const auto & ErrorNoContinuum = m_ErrorNoContinuum;
  for (Int32 iS = 0; iS < support.size(); iS++) {
    for (Int32 j = support[iS].GetBegin(); j < support[iS].GetEnd(); j++) {
      n++;
      // w = 1.0 / ErrorNoContinuum[j];
      // sumErr += w;
      // m += ContinuumFluxAxis[j] * w;
      m += ContinuumFluxAxis[j];
    }
  }

  return m / Float64(n);
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
 *ApplyBalmerRuleLinSolve. If "all" or "ratiorange" is in the rules string, call
 *ApplyAmplitudeRatioRangeRule parameterized for OII. If "all" or "strongweak"
 *is in the rules string, call ApplyStrongHigherWeakRule for emission and then
 *for absorption.
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
  if (m_RestLineList.size() != modelSolution.Lines.size())
    THROWG(INTERNAL_ERROR, "m_restLineList "
                           "and modelSolution.m_Lines have different size");

  setRedshift(modelSolution.Redshift, false);
  SetVelocityEmission(modelSolution.EmissionVelocity);
  SetVelocityAbsorption(modelSolution.AbsorptionVelocity);

  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  for (Int32 iRestLine = 0; iRestLine < m_RestLineList.size(); iRestLine++) {
    Int32 eIdx = modelSolution.ElementId[iRestLine];
    if (eIdx == undefIdx)
      continue;

    m_Elements[eIdx]->prepareSupport(spectralAxis, m_Redshift, m_lambdaRange);
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

    if (m_RestLineList[iRestLine] != modelSolution.Lines[iRestLine])
      THROWG(INTERNAL_ERROR, "m_restLineList and "
                             "modelSolution.m_Lines dont correspond");

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

  if (!std::isnan(modelSolution.LyaIgm)) {
    std::string lyaTag = ltags.lya_em;
    Int32 idxLyaE = m_Elements.findElementIndex(lyaTag);
    if (idxLyaE != undefIdx)
      m_Elements[idxLyaE]->SetSymIgmParams(
          {modelSolution.LyaIgm, modelSolution.Redshift});
  }
  return;
}

// return error: 1=can't find element index, 2=Abs_width not high enough
// compared to Em_width
Int32 CLineModelFitting::improveBalmerFit() {

  // Emission Balmer lines
  TStringList linetagsE;
  linetagsE.push_back(linetags::halpha_em);
  linetagsE.push_back(linetags::hbeta_em);
  linetagsE.push_back(linetags::hgamma_em);
  linetagsE.push_back(linetags::hdelta_em);
  // Absorption Balmer lines
  TStringList linetagsA;
  linetagsA.push_back(linetags::halpha_abs);
  linetagsA.push_back(linetags::hbeta_abs);
  linetagsA.push_back(linetags::hgamma_abs);
  linetagsA.push_back(linetags::hdelta_abs);
  // Additional lines to be fitted with the Balmer lines, WARNING: only EMISSION
  // for now !!
  TStringList linetagsNII;
  linetagsNII.push_back(linetags::niia_em);
  linetagsNII.push_back(linetags::niib_em);
  TStringList linetagsVoid;
  std::vector<TStringList> linetagsMore;
  linetagsMore.push_back(linetagsNII);
  linetagsMore.push_back(linetagsVoid);
  linetagsMore.push_back(linetagsVoid);
  linetagsMore.push_back(linetagsVoid);

  if (linetagsE.size() != linetagsA.size() ||
      linetagsE.size() != linetagsMore.size()) {
    return -1;
  }

  for (Int32 itag = 0; itag < linetagsE.size(); itag++) {
    std::string tagE = linetagsE[itag];
    std::string tagA = linetagsA[itag];

    Int32 ilineE = m_Elements.findElementIndex(tagE, CLine::nType_Emission);
    Int32 ilineA = m_Elements.findElementIndex(tagA, CLine::nType_Absorption);
    // Were the lines indexes found ?
    if (ilineE < 0 || ilineA < 0) {
      continue;
    }
    // for now only allow this process if Em and Abs line are single lines
    if (m_Elements[ilineE]->m_Lines.size() > 1 ||
        m_Elements[ilineA]->m_Lines.size() > 1) {
      continue;
    }
    Int32 subeIdxE = 0;
    Int32 subeIdxA = 0;

    // find the linesMore unique elements indexes
    TInt32List ilinesMore;
    for (Int32 imore = 0; imore < linetagsMore[itag].size(); imore++) {
      std::string tagMore = linetagsMore[itag][imore];
      Int32 ilineMore =
          m_Elements.findElementIndex(tagMore, CLine::nType_Emission);
      if (ilineMore < 0) {
        continue;
      }
      ilinesMore.push_back(ilineMore);
    }
    std::sort(ilinesMore.begin(), ilinesMore.end());
    ilinesMore.erase(std::unique(ilinesMore.begin(), ilinesMore.end()),
                     ilinesMore.end());
    for (Int32 imore = 0; imore < ilinesMore.size(); imore++) {
      Log.LogDebug("    model: balmerImprove more tags = %d",
                   ilinesMore[imore]);
    }

    // try if the width is significantly different: abs > em
    Float64 AbsVSEmWidthCoeffThreshold = 2.0;
    Float64 muE = NAN;
    Float64 muA = NAN;
    Float64 sigmaE = NAN;
    Float64 sigmaA = NAN;
    m_Elements[ilineE]->getObservedPositionAndLineWidth(
        subeIdxE, m_Redshift, muE, sigmaE,
        false); // do not apply Lya asym offset
    m_Elements[ilineA]->getObservedPositionAndLineWidth(
        subeIdxA, m_Redshift, muA, sigmaA,
        false); // do not apply Lya asym offset
    if (sigmaA < AbsVSEmWidthCoeffThreshold * sigmaE) {
      continue;
    }

    // simulatneous fit with linsolve
    Float64 modelErr_init = m_Elements.getModelErrorUnderElement(
        ilineA, m_SpcFluxAxis, m_SpectrumModel.GetFluxAxis());
    Float64 ampA = m_Elements[ilineA]->GetFittedAmplitude(0);
    Float64 amp_errorA = m_Elements[ilineA]->GetFittedAmplitudeErrorSigma(0);
    Float64 ampE = m_Elements[ilineE]->GetFittedAmplitude(0);
    Float64 amp_errorE = m_Elements[ilineE]->GetFittedAmplitudeErrorSigma(0);
    TFloat64List ampsMore;
    TFloat64List ampErrorsMore;
    for (Int32 imore = 0; imore < ilinesMore.size(); imore++) {
      Float64 amp = m_Elements[ilinesMore[imore]]->GetFittedAmplitude(0);
      Float64 ampErr =
          m_Elements[ilinesMore[imore]]->GetFittedAmplitudeErrorSigma(0);
      ampsMore.push_back(amp);
      ampErrorsMore.push_back(ampErr);
    }

    TInt32List eltsIdx;
    eltsIdx.push_back(ilineA);
    eltsIdx.push_back(ilineE);
    for (Int32 imore = 0; imore < ilinesMore.size(); imore++) {
      eltsIdx.push_back(ilinesMore[imore]);
    }
    TFloat64List ampsfitted;
    TFloat64List errorsfitted;
    fitAmplitudesLinSolve(eltsIdx, m_SpectrumModel.GetSpectralAxis(),
                          m_spcFluxAxisNoContinuum, m_ContinuumFluxAxis,
                          ampsfitted, errorsfitted);

    // decide if the fit is better than previous amps
    TInt32List elts;
    elts.push_back(ilineA);
    elts.push_back(ilineE);
    for (Int32 imore = 0; imore < ilinesMore.size(); imore++) {
      elts.push_back(ilinesMore[imore]);
    }
    refreshModelUnderElements(elts);
    Float64 modelErr_withfit = m_Elements.getModelErrorUnderElement(
        ilineA, m_SpcFluxAxis, m_SpectrumModel.GetFluxAxis());
    if (modelErr_withfit > modelErr_init) {
      Float64 nominal_ampA = m_Elements[ilineA]->GetNominalAmplitude(0);
      Float64 nominal_ampE = m_Elements[ilineE]->GetNominalAmplitude(0);
      m_Elements[ilineA]->SetFittedAmplitude(ampA / nominal_ampA,
                                             amp_errorA / nominal_ampA);
      m_Elements[ilineE]->SetFittedAmplitude(ampE / nominal_ampE,
                                             amp_errorE / nominal_ampE);
      for (Int32 imore = 0; imore < ilinesMore.size(); imore++) {
        Float64 nominal_ampMore =
            m_Elements[ilinesMore[imore]]->GetNominalAmplitude(0);
        m_Elements[ilinesMore[imore]]->SetFittedAmplitude(
            ampsMore[imore] / nominal_ampMore,
            ampErrorsMore[imore] / nominal_ampMore);
      }
    }
  }

  return 0;
}

/**
 * \brief Returns a CLineModelSolution object populated with the current
 *solutions.
 **/
CLineModelSolution CLineModelFitting::GetModelSolution(Int32 opt_level) {
  Int32 s = m_RestLineList.size();
  CLineModelSolution modelSolution;
  modelSolution.nDDL = m_Elements.GetModelNonZeroElementsNDdl();
  modelSolution.Lines = m_RestLineList; // copying the catalog?
  modelSolution.ElementId.resize(s);
  modelSolution.snrHa = NAN;
  modelSolution.lfHa = NAN;
  modelSolution.snrOII = NAN;
  modelSolution.lfOII = NAN;
  modelSolution.LyaWidthCoeff = NAN;
  modelSolution.LyaAlpha = NAN;
  modelSolution.LyaDelta = NAN;
  modelSolution.LyaIgm = undefIdx;
  modelSolution.Amplitudes = TFloat64List(s, NAN);
  modelSolution.AmplitudesUncertainties = TFloat64List(s, NAN);
  modelSolution.FittingError = TFloat64List(s, NAN);
  modelSolution.LambdaObs = TFloat64List(s, NAN);
  modelSolution.Offset = TFloat64List(s, NAN);
  modelSolution.Velocity = TFloat64List(s, NAN);
  modelSolution.CenterContinuumFlux = TFloat64List(s, NAN);
  modelSolution.ContinuumError = TFloat64List(s, NAN);
  modelSolution.Sigmas = TFloat64List(s, NAN);
  modelSolution.Fluxs = TFloat64List(s, NAN);
  modelSolution.FluxErrors = TFloat64List(s, NAN);
  modelSolution.FluxDirectIntegration = TFloat64List(s, NAN);
  modelSolution.FluxDirectIntegrationError = TFloat64List(s, NAN);
  modelSolution.OutsideLambdaRange = TBoolList(s, true);
  modelSolution.fittingGroupInfo = TStringList(s, "undefined");
  modelSolution.continuum_pCoeff0 = TFloat64List(s, NAN);
  modelSolution.continuum_pCoeff1 = TFloat64List(s, NAN);
  modelSolution.continuum_pCoeff2 = TFloat64List(s, NAN);

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

    if (opt_level) // brief, to save processing time, do not estimate fluxes and
                   // high level line properties
    {
      modelSolution.FittingError[iRestLine] =
          m_Elements.getModelErrorUnderElement(eIdx, m_SpcFluxAxis,
                                               m_SpectrumModel.GetFluxAxis());
      TPolynomCoeffs polynom_coeffs = getPolynomCoeffs(eIdx);
      // save polynom info to output them in hdf5, mainly to recontruct linemeas
      // model
      modelSolution.continuum_pCoeff0[iRestLine] = polynom_coeffs.x0;
      modelSolution.continuum_pCoeff1[iRestLine] = polynom_coeffs.x1;
      modelSolution.continuum_pCoeff2[iRestLine] = polynom_coeffs.x2;

      Float64 cont = m_Elements[eIdx]->GetContinuumAtCenterProfile(
          subeIdx, m_SpectrumModel.GetSpectralAxis(), m_Redshift,
          m_ContinuumFluxAxis, polynom_coeffs);
      modelSolution.CenterContinuumFlux[iRestLine] = cont;
      modelSolution.ContinuumError[iRestLine] =
          GetContinuumError(eIdx, subeIdx);
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
 * Here,we assume that we already fitted one line at a time, which is not really
 * the case
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

CContinuumModelSolution CLineModelFitting::GetContinuumModelSolution() const {
  CContinuumModelSolution continuumModelSolution;

  continuumModelSolution.tplName = m_fitContinuum_tplName;
  continuumModelSolution.tplEbmvCoeff = m_fitContinuum_tplFitEbmvCoeff;
  continuumModelSolution.tplMeiksinIdx = m_fitContinuum_tplFitMeiksinIdx;
  continuumModelSolution.tplAmplitude = m_fitContinuum_tplFitAmplitude;
  continuumModelSolution.tplAmplitudeError =
      m_fitContinuum_tplFitAmplitudeError;
  continuumModelSolution.tplMerit = m_fitContinuum_tplFitMerit;
  continuumModelSolution.tplMeritPhot = m_fitContinuum_tplFitMerit_phot;
  continuumModelSolution.tplDtm = m_fitContinuum_tplFitDtM;
  continuumModelSolution.tplMtm = m_fitContinuum_tplFitMtM;
  continuumModelSolution.tplLogPrior = m_fitContinuum_tplFitLogprior;
  continuumModelSolution.tplRedshift = m_fitContinuum_tplFitRedshift;
  continuumModelSolution.pCoeffs = m_fitContinuum_tplFitPolyCoeffs;

  return continuumModelSolution;
}

/**
 * \brief Returns the size of m_Elements.
 **/
Int32 CLineModelFitting::GetNElements() const {
  Int32 nddl = m_Elements.size();
  return nddl;
}

void CLineModelFitting::SetLSF() {
  const std::shared_ptr<const CLSF> &lsf = m_inputSpc.GetLSF();

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
 * \brief this function estimates the continuum after removal(interpolation) of
 *the flux samples under the lines for a given redshift
 * //todo: use lambdaRange in order to speed up the continuum estimation
 *process. (rmq: use lambdarange with some margin in order to not detetriorate
 *cont. estimation close to the borders)
 **/
void CLineModelFitting::EstimateSpectrumContinuum(Float64 opt_enhance_lines) {
  TInt32List validEltsIdx = m_Elements.GetModelValidElementsIndexes();
  // TInt32List xInds = getSupportIndexes( validEltsIdx );
  const CSpectrumSpectralAxis &spectralAxis = m_SpectrumModel.GetSpectralAxis();
  const auto &ContinuumFluxAxis = m_ContinuumFluxAxis;

  // create new spectrum, which is corrected under the lines

  /*
  //1. interp from previous values
  for( Int32 t=0;t<spectralAxis.GetSamplesCount();t++)
  {
      Y[t] = m_SpcFluxAxis[t];
  }
  Int32 idx = 0;
  Float64 valf = m_SpcFluxAxis[0];
  Float64 corrRatio=0.8;
  for (int i = 0; i < xInds.size(); i++)
  {
      idx = xInds[i];
      if(idx>0){
          if ( std::find(xInds.begin(), xInds.end(), idx-1) == xInds.end() )
          {
              valf=m_SpcFluxAxis[idx-1];
          }
      }
      Y[idx]= corrRatio*valf + (1.0-corrRatio)*m_SpcFluxAxis[idx];
  }
  //*/

  // 2. subtract lines from model
  // model for subtraction
  CSpectrumFluxAxis spcmodel4linefittingFluxAxis =
      m_SpectrumModel.GetFluxAxis();
  // CSpectrum spcmodel4linefitting = GetModelSpectrum();
  for (Int32 i = 0; i < spcmodel4linefittingFluxAxis.GetSamplesCount(); i++) {
    spcmodel4linefittingFluxAxis[i] -= ContinuumFluxAxis[i];
  }
  // optionnaly enhance the abs model component
  if (opt_enhance_lines > 0.0) {

    bool filterOnlyAbs = false;
    for (Int32 t = 0; t < spectralAxis.GetSamplesCount(); t++) {
      if (filterOnlyAbs) {
        if (spcmodel4linefittingFluxAxis[t] < 0.0) {
          spcmodel4linefittingFluxAxis[t] *= opt_enhance_lines;
        }
      } else {
        spcmodel4linefittingFluxAxis[t] *= opt_enhance_lines;
      }
    }
  }

  // subtract the lines component
  CSpectrumFluxAxis fluxAxisNothingUnderLines = m_SpcFluxAxis;
  for (Int32 t = 0; t < spectralAxis.GetSamplesCount(); t++) {
    fluxAxisNothingUnderLines[t] -= as_const(spcmodel4linefittingFluxAxis)[t];
  }

  // evaluate contiuum
  CSpectrum spcCorrectedUnderLines(m_inputSpc);
  spcCorrectedUnderLines.SetFluxAxis(std::move(fluxAxisNothingUnderLines));
  m_ContinuumFluxAxis = spcCorrectedUnderLines.GetContinuumFluxAxis();
}

/**
 * \brief this function returns the dtd value withing the wavelength range for a
 *given spcComponent
 *
 **/
// TODO rename this ! not a simple getter
Float64 CLineModelFitting::getDTransposeD() {
  if (!(m_dTransposeDLambdaRange.GetBegin() == m_lambdaRange.GetBegin() &&
        m_dTransposeDLambdaRange.GetEnd() == m_lambdaRange.GetEnd())) {
    initDtd();
  }

  return m_dTransposeD;
}

/**
 * \brief this function returns the dtd value withing the wavelength range for a
 *given spcComponent
 *
 **/
// TODO rename this ! not a simple getter
Float64 CLineModelFitting::getLikelihood_cstLog() {
  if (!(m_dTransposeDLambdaRange.GetBegin() == m_lambdaRange.GetBegin() &&
        m_dTransposeDLambdaRange.GetEnd() == m_lambdaRange.GetEnd())) {
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
      m_SpectrumModel.GetSpectralAxis();
  const CSpectrumFluxAxis &spcFluxAxis = m_SpcFluxAxis;
  const CSpectrumFluxAxis &spcFluxAxisNoContinuum = m_spcFluxAxisNoContinuum;
  const auto &ErrorNoContinuum = m_ErrorNoContinuum;

  Int32 numDevs = 0;
  Float64 dtd = 0.0;
  const Float64 *Yspc = spcFluxAxis.GetSamples();
  const Float64 *YspcNoContinuum = spcFluxAxisNoContinuum.GetSamples();
  Float64 flux = 0.0;

  Float64 imin = spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange.GetBegin());
  Float64 imax = spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange.GetEnd());
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
  const CSpectrumSpectralAxis &spcSpectralAxis =
      m_SpectrumModel.GetSpectralAxis();
  const CSpectrumFluxAxis &spcFluxAxis = m_SpectrumModel.GetFluxAxis();
  const auto &ErrorNoContinuum = m_ErrorNoContinuum;

  Int32 numDevs = 0;
  Float64 mtm = 0.0;
  const Float64 *Yspc = spcFluxAxis.GetSamples();
  Float64 diff = 0.0;

  Float64 imin = spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange.GetBegin());
  Float64 imax = spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange.GetEnd());
  for (Int32 j = imin; j < imax; j++) {
    numDevs++;
    { diff = Yspc[j]; }
    mtm += (diff * diff) / (ErrorNoContinuum[j] * ErrorNoContinuum[j]);
  }
  // Log.LogDebug( "CLineModelFitting::EstimateMTransposeM val = %f", mtm );

  return mtm;
}

/**
 * \brief this function estimates the likelihood_cstLog term withing the
 *wavelength range
 **/
Float64 CLineModelFitting::EstimateLikelihoodCstLog() const {
  const CSpectrumSpectralAxis &spcSpectralAxis =
      m_SpectrumModel.GetSpectralAxis();
  const auto &ErrorNoContinuum = m_ErrorNoContinuum;

  Int32 numDevs = 0;
  Float64 cstLog = 0.0;
  Float64 sumLogNoise = 0.0;

  Float64 imin = spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange.GetBegin());
  Float64 imax = spcSpectralAxis.GetIndexAtWaveLength(m_lambdaRange.GetEnd());
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
      nominalAmp = m_CatalogTplShape
                       .getNominalAmplitudeCorrespondance()[iElts][iCatalog][j];
      m_Elements[iElts]->SetNominalAmplitude(j, nominalAmp);
    }
  }
  return true;
}

/**
 * @brief CLineModelFitting::SetMultilineNominalAmplitudes
 * This method sets the linemodel unique elt nominal amplitudes to the
 * corresponding value of the iCatalog st catalog. INFO: slow method
 * @param iCatalog
 * @return
 */
bool CLineModelFitting::SetMultilineNominalAmplitudes(Int32 iCatalog) {
  // first set all amplitudes to 0.0
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    // get the max nominal amplitude
    Int32 nLines = m_Elements[iElts]->GetSize();
    for (Int32 j = 0; j < nLines; j++) {
      m_Elements[iElts]->SetNominalAmplitude(j, 0.0);
    }
  }
  // loop the amplitudes in the iLine_st catalog
  CLineCatalog::TLineVector currentCatalogLineList =
      m_CatalogTplShape.GetCatalog(iCatalog).GetList();
  Int32 nLines = currentCatalogLineList.size();
  for (Int32 kL = 0; kL < nLines; kL++) {
    Float64 nominalAmp = currentCatalogLineList[kL].GetNominalAmplitude();
    Float64 restLambda = currentCatalogLineList[kL].GetPosition();
    // here below we retrieve ism coeff and apply it to nominal Amp to then set
    // lines Nominal amplitudes
    Float64 dustCoeff =
        m_tplCatalog.GetTemplate(m_tplCategoryList[0], 0)
            ->m_ismCorrectionCalzetti->GetDustCoeff(
                m_CatalogTplShape.GetIsmIndex(iCatalog), restLambda);
    nominalAmp *= dustCoeff;
    // find line in the elementList
    for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
      // get the max nominal amplitude
      Int32 nLines = m_Elements[iElts]->GetSize();
      for (Int32 j = 0; j < nLines; j++) {

        if (m_Elements[iElts]->m_Lines[j].GetName() ==
            currentCatalogLineList[kL].GetName()) {
          m_Elements[iElts]->SetNominalAmplitude(j, nominalAmp);
        }
      }
    }
  }
  return true;
}
