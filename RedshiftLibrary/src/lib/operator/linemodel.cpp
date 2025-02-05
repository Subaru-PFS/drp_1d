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
#include <boost/chrono/thread_clock.hpp>
#include <boost/format.hpp>
#include <boost/numeric/conversion/bounds.hpp>

#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/common/flag.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/common/indexing.h"
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/common/vectorOperations.h"
#include "RedshiftLibrary/common/zgridparam.h"
#include "RedshiftLibrary/extremum/extremum.h"
#include "RedshiftLibrary/linemodel/lineratiomanager.h"
#include "RedshiftLibrary/linemodel/outsideLineMaskBuilder.h"
#include "RedshiftLibrary/linemodel/powerlawstore.h"
#include "RedshiftLibrary/linemodel/rulesmanager.h"
#include "RedshiftLibrary/linemodel/templatesfitstore.h"
#include "RedshiftLibrary/linemodel/templatesortho.h"
#include "RedshiftLibrary/linemodel/tplratiomanager.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/operator/linemodel.h"
#include "RedshiftLibrary/operator/modelphotvalueresult.h"
#include "RedshiftLibrary/operator/powerlaw.h"
#include "RedshiftLibrary/operator/powerlawresult.h"
#include "RedshiftLibrary/operator/templatefitting.h"
#include "RedshiftLibrary/operator/templatefittinglog.h"
#include "RedshiftLibrary/operator/templatefittingresult.h"
#include "RedshiftLibrary/operator/templatefittingwithphot.h"
#include "RedshiftLibrary/processflow/autoscope.h"
#include "RedshiftLibrary/processflow/context.h"
#include "RedshiftLibrary/processflow/inputcontext.h"
#include "RedshiftLibrary/processflow/parameterstore.h"
#include "RedshiftLibrary/spectrum/axis.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "RedshiftLibrary/statistics/deltaz.h"
#include "RedshiftLibrary/statistics/priorhelper.h"

using namespace NSEpic;
using namespace std;

/**
 * @brief COperatorLineModel::ComputeFirstPass
 * @return 0=no errors, -1=error
 */
std::shared_ptr<const CLineModelResult> const
COperatorLineModel::ComputeFirstPass() {
  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();
  m_opt_continuumcomponent =
      ps->GetScoped<std::string>("lineModel.continuumComponent");

  std::shared_ptr<const CTemplateCatalog> tplCatalog;
  tplCatalog = Context.GetTemplateCatalog();

  makeContinuumFittingOperator(m_redshifts);
  if (m_continuumFittingOperator->IsFFTProcessing()) { // create a default
    const TFloat64List &redshifts = m_redshifts;
    m_fittingManager = std::make_shared<CLineModelFitting>(
        std::make_shared<COperatorTemplateFitting>(redshifts));
  } else {
    m_fittingManager =
        std::make_shared<CLineModelFitting>(m_continuumFittingOperator);
  }

  Int32 nfitcontinuum = 0;
  if (m_opt_continuumcomponent.isTplFitXXX())
    nfitcontinuum = tplCatalog->GetTemplateCount(m_tplCategory);
  else if (m_opt_continuumcomponent.isPowerLawXXX())
    nfitcontinuum = 1;
  m_result->Init(m_redshifts, Context.getCLineMap(), nfitcontinuum,
                 m_fittingManager->getTplratio_count(),
                 m_fittingManager->getTplratio_priors());

  Log.LogInfo("  Operator-Linemodel: initialized");

  // commom between firstpass and secondpass processes
  // TODO not pretty, maybe move fitcontinuum_prior help building to operator
  m_phelperContinuum =
      m_fittingManager->getContinuumManager()->SetFitContinuum_PriorHelper();

  Log.LogInfo("  Operator-Linemodel: start processing");

  // Set model parameters to FIRST-PASS
  m_fittingManager->setPassMode(1);
  Log.LogInfo(
      "  Operator-Linemodel: ---------- ---------- ---------- ----------");
  Log.LogInfo("  Operator-Linemodel: now computing first-pass");
  Log.LogInfo(
      "  Operator-Linemodel: ---------- ---------- ---------- ----------");

  // fit continuum
  ////////////////////
  if (m_opt_continuumcomponent.isContinuumFit()) {

    Log.LogInfo(Formatter() << "Precompuute continuum fit");
    m_tplfitStore_firstpass = PrecomputeContinuumFit(m_redshifts);
  }

  m_result->nSpcSamples = m_fittingManager->computeSpcNSamples();
  m_result->dTransposeD = m_fittingManager->getDTransposeD();
  m_result->cstLog = m_fittingManager->getLikelihood_cstLog();

  Int32 contreest_iterations =
      ps->GetScoped<std::string>("lineModel.continuumReestimation") == "always"
          ? 1
          : 0;

  // Set model parameter: abs lines limit
  Float64 absLinesLimit = 1.0; //-1 to disable, 1.0 is typical
  m_fittingManager->SetAbsLinesLimit(absLinesLimit);
  Log.LogInfo(Formatter() << "  Operator-Linemodel: set abs lines limit to "
                          << absLinesLimit
                          << " (ex: -1 means "
                             "disabled)");

  // Set model parameter: continuum least-square estimation fast
  // note: this fast method requires continuum templates and linemodels to be
  // orthogonal. The velfit option turns this trickier...
  m_estimateLeastSquareFast = 0;
  // TODO should be called only with lineRatioType=tplratio
  m_fittingManager->m_lineRatioManager->SetLeastSquareFastEstimationEnabled(
      m_estimateLeastSquareFast);
  Log.LogInfo(
      Formatter() << "  Operator-Linemodel: set estimateLeastSquareFast to "
                  << m_estimateLeastSquareFast
                  << " (ex: "
                     "0 means disabled)");

  TBoolList allAmplitudesZero;

  boost::chrono::thread_clock::time_point start_mainloop =
      boost::chrono::thread_clock::now();

  //#pragma omp parallel for
  m_fittingManager->logParameters();
  for (Int32 i = 0; i < m_result->Redshifts.size(); i++) {

    m_result->ChiSquare[i] = m_fittingManager->fit(
        m_result->Redshifts[i], m_result->LineModelSolutions[i],
        m_result->ContinuumModelSolutions[i], contreest_iterations,
        m_opt_skipsecondpass);

    m_result->ScaleMargCorrection[i] =
        m_fittingManager->getScaleMargCorrection();
    // TODO create a method to make this or comparison ("mustComputeContinuum"
    // or something similar)
    if (m_opt_continuumcomponent.isContinuumFit())
      m_result->SetChisquareContinuumResult(i, m_tplfitStore_firstpass);

    if (m_fittingManager->getLineRatioType() == "tplRatio")
      m_result->SetChisquareTplratioResult(
          i, dynamic_pointer_cast<CTplratioManager>(
                 m_fittingManager->m_lineRatioManager));

    m_result->ChiSquareContinuum[i] =
        m_estimateLeastSquareFast
            ? m_fittingManager->getLeastSquareContinuumMeritFast()
            : m_fittingManager->getLeastSquareContinuumMerit();

    m_result->ScaleMargCorrectionContinuum[i] =
        m_fittingManager->getContinuumManager()
            ->getContinuumScaleMargCorrection();
    Log.LogDebug(Formatter() << "  Operator-Linemodel: Z interval "
                             << m_result->Redshifts[i] << ": and zvalue = " << i
                             << " Chi2 = " << m_result->ChiSquare[i]);

    // Flags on continuum and model amplitudes
    bool continuumAmplitudeZero = true;
    if (m_opt_continuumcomponent.isTplFitXXX() &&
        m_result->ContinuumModelSolutions[i].tplAmplitude > 0.0) // false if NAN
      continuumAmplitudeZero = false;
    if (m_opt_continuumcomponent.isPowerLawXXX())
      if (m_result->ContinuumModelSolutions[i].a1 > 0.0 ||
          m_result->ContinuumModelSolutions[i].a2 > 0.0)
        continuumAmplitudeZero = false;
    bool modelAmplitudesZero = true;
    auto it = std::find_if(m_result->LineModelSolutions[i].Amplitudes.cbegin(),
                           m_result->LineModelSolutions[i].Amplitudes.cend(),
                           [](Float64 amp) { return amp > 0.0; });
    if (it != m_result->LineModelSolutions[i].Amplitudes.end())
      modelAmplitudesZero = false;
    allAmplitudesZero.push_back(modelAmplitudesZero && continuumAmplitudeZero);
  }
  // Check if all amplitudes are zero for all z
  bool checkAllAmplitudes = AllAmplitudesAreZero(allAmplitudesZero);
  if (checkAllAmplitudes == true)
    THROWG(ErrorCode::NULL_MODEL,
           "Null amplitudes (continuum & model) at all z");

  boost::chrono::thread_clock::time_point stop_mainloop =
      boost::chrono::thread_clock::now();
  Float64 duration_mainloop =
      boost::chrono::duration_cast<boost::chrono::microseconds>(stop_mainloop -
                                                                start_mainloop)
          .count();
  Float64 duration_firstpass_seconds = duration_mainloop / 1e6;
  Log.LogInfo(Formatter() << "  Operator-Linemodel: first-pass done in "
                          << duration_firstpass_seconds << " sec");
  Log.LogInfo(Formatter() << "<proc-lm-firstpass><"
                          << (Int32)duration_firstpass_seconds << ">");
  return m_result;
}

bool COperatorLineModel::AllAmplitudesAreZero(const TBoolList &amplitudesZero) {
  auto it = std::find_if(amplitudesZero.cbegin(), amplitudesZero.cend(),
                         [](bool ampIsZero) { return !ampIsZero; });
  if (it != amplitudesZero.end())
    return false;
  return true;
}

bool COperatorLineModel::isfftprocessingActive(Int32 redshiftsTplFitCount) {
  if (m_opt_continuumcomponent.isPowerLawXXX()) {
    return false;
  }

  bool fftprocessing =
      (m_fittingManager == nullptr) || m_fittingManager->GetPassNumber() == 1
          ? m_opt_tplfit_fftprocessing
          : m_opt_tplfit_fftprocessing_secondpass;
  Log.LogDebug(Formatter()
               << "COperatorLineModel::isfftprocessingActive: redshtplfitsize "
               << redshiftsTplFitCount);
  // heuristic value corresponding to a threshold below wich fftprocessing is
  // considered slow nb of redshift samples
  const Int32 fftprocessing_min_sample_nb = 100;
  if ((redshiftsTplFitCount < fftprocessing_min_sample_nb) && fftprocessing) {
    fftprocessing = false;
    Log.LogInfo("COperatorLineModel::isfftprocessingActive: auto deselect fft "
                "processing (faster when only few redshifts calc. points)");
  }
  return fftprocessing;
}

void COperatorLineModel::fitContinuumTemplates(
    Int32 candidateIdx, const TFloat64List &redshiftsContinuumFit,
    std::vector<std::shared_ptr<const COperatorResult>> &chisquareResultsAllTpl,
    TStringList &chisquareResultsTplName) {
  std::shared_ptr<const CTemplateCatalog> tplCatalog =
      Context.GetTemplateCatalog();
  Float64 overlapThreshold = 1.0;
  std::string opt_interp = "preComputedFineGrid"; //"lin";
  Log.LogDetail(Formatter()
                << "COperatorLineModel::PrecomputeContinuumFit: fitContinuum "
                   "opt_interp = "
                << opt_interp);
  TInt32List meiksinIndices;
  TInt32List ebmvIndices;
  TTemplateConstRefList tplList;
  bool fftprocessing = isfftprocessingActive(redshiftsContinuumFit.size());
  if (m_fittingManager->GetPassNumber() == 2 &&
      m_continnuum_fit_option == EContinuumFit::reFitFirstPass) {
    // case where we only want to refit around the m_opt_fitcontinuum_maxN
    // best continuum from firstpass
    getContinuumInfoFromFirstpassFitStore(candidateIdx, meiksinIndices,
                                          ebmvIndices, tplList, fftprocessing);
  } else {
    if (m_opt_continuumcomponent.isPowerLawXXX())
      tplList.push_back({});
    else
      tplList = tplCatalog->GetOrthoTemplateList(TStringList{m_tplCategory},
                                                 fftprocessing);
    ebmvIndices.assign(tplList.size(), undefIdx);
    meiksinIndices.assign(tplList.size(), undefIdx);
    Log.LogDebug(Formatter()
                 << "Processing " << tplList.size() << " templates");
  }
  for (Int32 i = 0; i < tplList.size(); i++) {
    std::shared_ptr<const COperatorResult> templatefittingResult;
    std::string tplname;
    if (m_opt_continuumcomponent.isPowerLawXXX()) {
      auto castedOperator = std::dynamic_pointer_cast<COperatorPowerLaw>(
          m_continuumFittingOperator);
      tplname = "powerLaw";
      castedOperator->SetRedshifts(redshiftsContinuumFit);
      templatefittingResult =
          castedOperator->Compute(m_opt_tplfit_extinction, m_opt_tplfit_dustFit,
                                  m_opt_continuum_null_amp_threshold, "full",
                                  ebmvIndices[i], meiksinIndices[i]);
    } else {
      CPriorHelper::TPriorZEList zePriorData;
      tplname = tplList[i]->GetName();
      m_phelperContinuum->GetTplPriorData(tplname, redshiftsContinuumFit,
                                          zePriorData);
      Log.LogDebug(Formatter() << "Processing tpl " << tplname);
      tplList[i]->setRebinInterpMethod(opt_interp);

      auto castedOperator =
          std::dynamic_pointer_cast<COperatorTemplateFittingBase>(
              m_continuumFittingOperator);
      castedOperator->SetRedshifts(redshiftsContinuumFit);
      templatefittingResult = castedOperator->Compute(
          tplList[i], overlapThreshold, opt_interp, m_opt_tplfit_extinction,
          m_opt_tplfit_dustFit, m_opt_continuum_null_amp_threshold, zePriorData,
          ebmvIndices[i], meiksinIndices[i]);
    }

    if (!templatefittingResult) {
      THROWG(ErrorCode::INTERNAL_ERROR,
             Formatter() << "Failed "
                            "to compute chisquare value for tpl=%"
                         << tplname);
    } else {
      chisquareResultsAllTpl.push_back(templatefittingResult);
      chisquareResultsTplName.push_back(tplname);
    }
  }
  return;
}

// get name, Meiksin and ISM coeff for all continuum
// returns vectors of these entities
void COperatorLineModel::getContinuumInfoFromFirstpassFitStore(
    Int32 candidateIdx, TInt32List &meiksinIndices, TInt32List &ebmvIndices,
    TTemplateConstRefList &tplList, bool fft) const {

  meiksinIndices.assign(m_opt_fitcontinuum_maxN, undefIdx);
  ebmvIndices.assign(m_opt_fitcontinuum_maxN, undefIdx);

  if (m_fittingManager->GetPassNumber() != 2 ||
      m_continnuum_fit_option != EContinuumFit::reFitFirstPass)
    return;

  if (candidateIdx < 0 || candidateIdx >= m_firstpass_extremaResult->size())
    THROWG(ErrorCode::INTERNAL_ERROR, "Candidate index is out of range");

  std::shared_ptr<const CTemplateCatalog> tplCatalog =
      Context.GetTemplateCatalog();

  for (Int32 icontinuum = 0; icontinuum < m_opt_fitcontinuum_maxN;
       icontinuum++) {
    // get the closest lower or equal redshift in coarse grid
    Int32 coarseIdx = m_tplfitStore_firstpass->getClosestLowerRedshiftIndex(
        m_firstpass_extremaResult->Redshift(candidateIdx));

    CContinuumModelSolution fitValue =
        m_tplfitStore_firstpass->GetFitValues(coarseIdx, icontinuum);

    if (m_opt_continuumcomponent.isPowerLawXXX())
      tplList.push_back({});
    else
      tplList.push_back(tplCatalog->GetTemplateByName(
          TStringList{m_tplCategory}, fitValue.name, true, fft));

    if (m_opt_tplfit_extinction)
      meiksinIndices[icontinuum] = fitValue.meiksinIdx;
    if (m_opt_tplfit_dustFit)
      ebmvIndices[icontinuum] =
          Context.getFluxCorrectionCalzetti()->GetEbmvIndex(fitValue.ebmvCoef);
  }
  return;
}

void COperatorLineModel::makeContinuumFittingOperator(
    const TFloat64List &redshifts) {
  // Initializes template fitting operator

  bool fftprocessing = isfftprocessingActive(redshifts.size());

  if (m_opt_continuumcomponent.isPowerLawXXX()) {
    m_continuumFittingOperator = std::make_shared<COperatorPowerLaw>(redshifts);
    return;
  }

  if (fftprocessing) {
    if (m_continuumFittingOperator == nullptr ||
        !m_continuumFittingOperator
             ->IsFFTProcessing()) // else reuse the shared
                                  // pointer for secondpass
      m_continuumFittingOperator =
          std::make_shared<COperatorTemplateFittingLog>(redshifts);
    return;
  }

  if (m_continuumFittingOperator == nullptr ||
      m_continuumFittingOperator->IsFFTProcessing()) { // else reuse the shared
                                                       // pointer for secondpass
    if (m_opt_tplfit_use_photometry) {
      const std::shared_ptr<const CPhotBandCatalog> &photBandCat =
          Context.GetPhotBandCatalog();
      std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();
      m_continuumFittingOperator =
          std::make_shared<COperatorTemplateFittingPhot>(
              photBandCat, redshifts,
              ps->GetScoped<Float64>("lineModel.photometry.weight"));
    } else
      m_continuumFittingOperator =
          std::make_shared<COperatorTemplateFitting>(redshifts);
  }
}
/**
 * Estimate once for all the continuum amplitude which is only dependent
 * from the name, ism and igm indexes. This is useful when the model
 * fitting option corresponds to fitting separately the continuum and the
 * lines. In such case, playing with (fit) Lines parameters (velocity, line
 * offsets, lines amps, etc.) do not affect continuum amplitudes.. thus we
 * can save time  by fitting once-for-all the continuum amplitudes, prior to
 * fitting the lines.
 * @candidateIdx@ is also an indicator of pass mode
 * */
std::shared_ptr<CContinuumFitStore>
COperatorLineModel::PrecomputeContinuumFit(const TFloat64List &redshifts,
                                           Int32 candidateIdx) {
  std::shared_ptr<CContinuumFitStore> continuumFitStore;
  if (m_opt_continuumcomponent.isPowerLawXXX()) {
    continuumFitStore = make_shared<CPowerLawStore>(redshifts);
  } else {
    continuumFitStore = make_shared<CTemplatesFitStore>(redshifts);
  }
  const TFloat64List &redshiftsContinuumFit =
      continuumFitStore->GetRedshiftList();
  Log.LogInfo(
      Formatter() << "COperatorLineModel::PrecomputeContinuumFit: continuum "
                     "redshift list n="
                  << redshiftsContinuumFit.size());

  bool fftprocessing = isfftprocessingActive(redshiftsContinuumFit.size());

  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();

  bool ignoreLinesSupport =
      ps->GetScoped<bool>("lineModel.continuumFit.ignoreLineSupport");
  boost::chrono::thread_clock::time_point start_tplfitprecompute =
      boost::chrono::thread_clock::now();
  Log.LogInfo(Formatter()
              << "COperatorLineModel::PrecomputeContinuumFit: continuum tpl "
                 "fitting: min="
              << redshifts.front() << ", max=" << redshifts.back());

  Log.LogInfo(Formatter()
              << "COperatorLineModel::PrecomputeContinuumFit: fftprocessing = "
              << fftprocessing);
  Log.LogDetail(
      Formatter()
      << "COperatorLineModel::PrecomputeContinuumFit: fitContinuum_dustfit = "
      << m_opt_tplfit_dustFit);
  Log.LogDetail(
      Formatter()
      << "COperatorLineModel::PrecomputeContinuumFit: fitContinuum_igm = "
      << m_opt_tplfit_extinction);

  makeContinuumFittingOperator(redshifts);

  if (fftprocessing && ignoreLinesSupport == true) {
    ignoreLinesSupport = false;
    Flag.warning(WarningCode::FORCED_IGNORELINESUPPORT_TO_FALSE,
                 Formatter() << "  COperatorLineModel::" << __func__
                             << ": unable to ignoreLinesSupport if "
                                "fftProcessing. ignoreLinesSupport disabled");
  }
  if (ignoreLinesSupport) {
    m_fittingManager->getSpectraIndex()
        .setAtBegining(); // TODO multiobs, dummy implementation
    m_continuumFittingOperator->setMaskBuilder(
        std::make_shared<COutsideLineMaskBuilder>(
            m_fittingManager->getElementListVector(),
            m_fittingManager->getSpectraIndex()));
  }
  std::vector<std::shared_ptr<const COperatorResult>> chisquareResultsAllTpl;
  TStringList chisquareResultsTplName;
  fitContinuumTemplates(candidateIdx, redshiftsContinuumFit,
                        chisquareResultsAllTpl, chisquareResultsTplName);

  // fill the fit store with fitted values: only the best fitted values FOR
  // EACH TEMPLATE are used
  if (m_opt_continuumcomponent.isPowerLawXXX()) {
    const auto &chisquareResult =
        std::dynamic_pointer_cast<const CPowerLawResult>(
            chisquareResultsAllTpl[0]);
    // size 1
    Float64 bestFitSNR = 0.0;
    Int32 nredshiftsTplFitResults = redshiftsContinuumFit.size();
    for (Int32 i = 0; i < nredshiftsTplFitResults; i++) {
      Float64 redshift = redshiftsContinuumFit[i];
      std::dynamic_pointer_cast<CPowerLawStore>(continuumFitStore)
          ->Add(chisquareResult->FitEbmvCoeff[i],
                chisquareResult->FitMeiksinIdx[i], redshift,
                chisquareResult->ChiSquare[i],
                chisquareResult->ReducedChiSquare[i], chisquareResult->coefs[i],
                chisquareResult->SNR[i]);

      if (chisquareResult->SNR[i] > bestFitSNR)
        bestFitSNR = chisquareResult->SNR[i];
    }
  } else {
    Int32 nredshiftsTplFitResults = redshiftsContinuumFit.size();
    for (Int32 i = 0; i < nredshiftsTplFitResults; i++) { // zIDX
      Float64 redshift = redshiftsContinuumFit[i];

      for (Int32 j = 0; j < chisquareResultsAllTpl.size(); j++) { // n templates
        const auto &chisquareResult =
            std::dynamic_pointer_cast<const CTemplateFittingResult>(
                chisquareResultsAllTpl[j]);

        std::dynamic_pointer_cast<CTemplatesFitStore>(continuumFitStore)
            ->Add(chisquareResultsTplName[j], chisquareResult->FitEbmvCoeff[i],
                  chisquareResult->FitMeiksinIdx[i], redshift,
                  chisquareResult->ChiSquare[i],
                  chisquareResult->ReducedChiSquare[i],
                  chisquareResult->ChiSquarePhot[i],
                  chisquareResult->FitAmplitude[i],
                  chisquareResult->FitAmplitudeError[i],
                  chisquareResult->FitAmplitudeSigma[i],
                  chisquareResult->FitDtM[i], chisquareResult->FitMtM[i],
                  chisquareResult->LogPrior[i], chisquareResult->SNR[i]);
      }
    }
    Log.LogDetail(
        Formatter()
        << "COperatorLineModel::PrecomputeContinuumFit: continuumcount set to "
        << continuumFitStore->getContinuumCount());

    if (continuumFitStore->getContinuumCount() < m_opt_fitcontinuum_maxN) {
      THROWG(ErrorCode::INTERNAL_ERROR,
             "Couldn't compute the required continuum count");
    }
  }
  m_fittingManager->getContinuumManager()->SetFitContinuum_FitStore(
      continuumFitStore);

  boost::chrono::thread_clock::time_point stop_tplfitprecompute =
      boost::chrono::thread_clock::now();
  Float64 duration_tplfitprecompute =
      boost::chrono::duration_cast<boost::chrono::microseconds>(
          stop_tplfitprecompute - start_tplfitprecompute)
          .count();
  Float64 duration_tplfit_seconds = duration_tplfitprecompute / 1e6;
  Log.LogInfo(
      Formatter()
      << "COperatorLineModel::PrecomputeContinuumFit: tplfit-precompute "
         "done in "
      << duration_tplfit_seconds << " sec");
  Log.LogDetail(Formatter() << "<proc-lm-tplfit><"
                            << (Int32)duration_tplfit_seconds << ">");

  // raise exception or switch to fromspectrum if bad or negative continuum
  // switch to NoContinuum if null-compatible continuum
  if (m_fittingManager->GetPassNumber() == 1)
    evaluateAndUpdateContinuumComponent(continuumFitStore);

  return continuumFitStore;
}

bool COperatorLineModel::updateContinuumComponentIfBadChi2(
    CContinuumModelSolution const &fitValues) {

  // check if min_z reducedChisquare is not to big
  auto const min_chi2_r = fitValues.reducedChi2;
  bool const good_cont = min_chi2_r <= m_opt_continuum_bad_chi2_threshold;

  if (good_cont)
    return false;

  if (!m_opt_continuumcomponent.isAuto())
    THROWG(ErrorCode::BAD_CONTINUUMFIT,
           Formatter() << "Bad continuum fit found at z=" << fitValues.redshift
                       << ": best continuum tpl " << fitValues.name
                       << ", reduced chi2 = " << min_chi2_r);

  Flag.warning(
      WarningCode::FORCED_CONTINUUM_COMPONENT_TO_FROMSPECTRUM,
      Formatter() << "Switching to continuum filtered from spectrum since "
                     "bad continuum fit found at z="
                  << fitValues.redshift << ": best continuum tpl "
                  << fitValues.name << ", reduced chi2 = " << min_chi2_r);
  m_opt_continuumcomponent.set("fromSpectrum");
  m_fittingManager->setContinuumComponent(m_opt_continuumcomponent);
  return true;
}

bool COperatorLineModel::updateContinuumComponentIfNegative(
    Float64 max_fitamplitudeSigma, CContinuumModelSolution const &fitValues) {

  bool const negative_cont =
      max_fitamplitudeSigma < m_opt_continuum_neg_amp_threshold;

  if (!negative_cont)
    return false;

  if (!m_opt_continuumcomponent.isAuto())
    THROWG(ErrorCode::NEGATIVE_CONTINUUMFIT,
           Formatter() << "Negative "
                          "continuum amplitude found at z="
                       << fitValues.redshift << ": best continuum tpl "
                       << fitValues.name
                       << ", amplitude/error = " << max_fitamplitudeSigma);

  Flag.warning(
      WarningCode::FORCED_CONTINUUM_COMPONENT_TO_FROMSPECTRUM,
      Formatter() << ": Switching to spectrum continuum since Negative "
                     "continuum amplitude found at z="
                  << fitValues.redshift << ": best continuum tpl "
                  << fitValues.name
                  << ", amplitude/error = " << max_fitamplitudeSigma
                  << " & error = " << fitValues.tplAmplitudeError);
  m_opt_continuumcomponent.set("fromSpectrum");
  m_fittingManager->setContinuumComponent(m_opt_continuumcomponent);

  return true;
}

bool COperatorLineModel::updateContinuumComponentIfNotSignificant(
    Float64 max_fitamplitudeSigma, CContinuumModelSolution const &fitValues) {

  // check if continuum is too weak comparing to the preset threshold, or
  // falls within [thres_neg; thresh_null], at all z

  bool significant =
      max_fitamplitudeSigma >= m_opt_continuum_null_amp_threshold;
  if (significant)
    return false;

  Flag.warning(
      WarningCode::FORCED_CONTINUUM_TO_NOCONTINUUM,
      Formatter() << ": Switching to nocontinuum since close"
                     "to null or not enough negative continuum amplitude "
                     "found at z="
                  << fitValues.redshift << ": best continuum tpl "
                  << fitValues.name
                  << ", amplitude/error = " << max_fitamplitudeSigma);
  m_opt_continuumcomponent.set("noContinuum");
  m_fittingManager->setContinuumComponent(m_opt_continuumcomponent);

  return true;
}

void COperatorLineModel::evaluateAndUpdateContinuumComponent(
    std::shared_ptr<CContinuumFitStore> const &continuumFitStore) {

  // get smallest reduced chi2 at all z
  auto const &fitValues_of_min_chi2_r = continuumFitStore->FindMinReducedChi2();
  if (updateContinuumComponentIfBadChi2(fitValues_of_min_chi2_r))
    return;

  // get greatest continuum amplitude at all z
  auto const &[max_fitamplitudeSigma, fitValues_of_max_amplitude] =
      continuumFitStore->FindMaxAmplitudeSigma();

  if (updateContinuumComponentIfNegative(max_fitamplitudeSigma,
                                         fitValues_of_max_amplitude))
    return;

  updateContinuumComponentIfNotSignificant(max_fitamplitudeSigma,
                                           fitValues_of_max_amplitude);
}

void COperatorLineModel::SetFirstPassCandidates(
    const TCandidateZbyRank &zCandidates) {

  m_firstpass_extremaResult =
      std::make_shared<LineModelExtremaResult>(zCandidates);

  // now preparing the candidates extrema results
  for (Int32 i = 0; i < m_firstpass_extremaResult->size(); i++) {
    // find the index in the zaxis results
    Int32 idx = CIndexing<Float64>::getIndex(m_result->Redshifts,
                                             zCandidates[i].second->Redshift);

    // save basic fitting info from first pass
    m_firstpass_extremaResult->m_savedModelFittingResults[i] =
        std::make_shared<CLineModelSolution>(m_result->LineModelSolutions[idx]);
    // save the continuum fitting parameters from first pass
    // std::shared_ptr<const CContinuumModelSolution> csolution =
    //     std::make_shared<CContinuumModelSolution>(
    //         m_result->ContinuumModelSolutions[idx]);
    auto &candidate = m_firstpass_extremaResult->getRankedCandidatePtr(i);

    candidate->updateFromContinuumModelSolution(
        m_result->ContinuumModelSolutions[idx]);
    candidate->updateFromLineModelSolution(m_result->LineModelSolutions[idx]);
    //... TODO: more first pass results can be saved here if needed
  }
}

std::shared_ptr<const CLineModelResult> const
COperatorLineModel::ComputeSecondPass() {

  std::shared_ptr<const CTemplateCatalog> tplCatalog =
      Context.GetTemplateCatalog();
  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();

  boost::chrono::thread_clock::time_point start_secondpass =
      boost::chrono::thread_clock::now();
  // Set model parameters to SECOND-PASS
  m_fittingManager->setPassMode(2);
  CContinuumManager::EFitType savedFitContinuumOption =
      m_fittingManager->getContinuumManager()
          ->GetFitContinuum_Option(); // the first time was set in
                                      // precomputeContinuumFit
  Log.LogInfo(
      "  Operator-Linemodel: ---------- ---------- ---------- ----------");
  Log.LogInfo("  Operator-Linemodel: now computing second-pass");
  Log.LogInfo(
      "  Operator-Linemodel: ---------- ---------- ---------- ----------");

  std::string opt_continuumfit_method =
      ps->GetScoped<std::string>("lineModel.secondPass.continuumFit");
  std::string opt_continuumreest =
      ps->GetScoped<std::string>("lineModel.continuumReestimation");
  std::string opt_fittingmethod =
      ps->GetScoped<std::string>("lineModel.fittingMethod");
  m_continnuum_fit_option =
      CTwoPassSolve::str2ContinuumFit.at(opt_continuumfit_method);

  // extend z around the extrema
  buildExtendedRedshifts();

  // insert extendedRedshifts into m_redshifts
  updateRedshiftGridAndResults(m_result);

  // Deal with continuum, either recompute it or keep from first pass
  if (m_opt_continuumcomponent.isContinuumFit()) {
    // precompute only whenever required and whenever the result can be a
    // tplfitStore
    if (mustReFit(m_continnuum_fit_option)) {
      m_tplfitStore_secondpass.resize(m_firstpass_extremaResult->size());
      for (Int32 i = 0; i < m_firstpass_extremaResult->size(); i++) {
        m_tplfitStore_secondpass[i] =
            PrecomputeContinuumFit(m_extendedRedshifts[i], i);
        if (m_opt_continuumcomponent.isFromSpectrum())
          break; // when set to "fromSpectrum" by PrecomputeContinuumFit
                 // because negative continuum with tplfitauto
      }
    } else {
      // since precompute is not called all the time, secondpass candidates do
      // not have systematically a tplfitstore_secondpass copy the firstpass
      // tplfitstore into the secondpass tplfitstore
      m_tplfitStore_secondpass.resize(1);
      if (m_continnuum_fit_option == EContinuumFit::fromFirstPass)
        m_tplfitStore_secondpass[0] = m_tplfitStore_firstpass;
      m_fittingManager->getContinuumManager()->SetFitContinuum_FitStore(
          m_tplfitStore_firstpass);
      // duplicated: double make sure that these info are present in the
      // modelElement
    }
  }

  // now that we recomputed what should be recomputed, we define once for all
  // the secondpass
  //  estimate second pass parameters (mainly elv, alv...)
  EstimateSecondPassParameters();

  // recompute the fine grid results around the extrema
  RecomputeAroundCandidates(opt_continuumreest, m_continnuum_fit_option);

  // additional fitting with fittingmethod=svdlcp2
  if (m_opt_secondpasslcfittingmethod == "svdlc" ||
      m_opt_secondpasslcfittingmethod == "svdlcp2") {
    Log.LogInfo(Formatter()
                << "  Operator-Linemodel: now computing second-pass "
                << m_opt_secondpasslcfittingmethod
                << " on each secondPass candidate (n="
                << m_firstpass_extremaResult->size() << ")");
    auto saved_fitter = std::move(m_fittingManager->m_fitter);
    m_fittingManager->SetFittingMethod(m_opt_secondpasslcfittingmethod);
    RecomputeAroundCandidates(opt_continuumreest, EContinuumFit::fromFirstPass,
                              true);
    m_fittingManager->m_fitter = std::move(saved_fitter);

    Log.LogInfo("  Operator-Linemodel: now re-computing the final chi2 for "
                "each candidate");
    RecomputeAroundCandidates(opt_continuumreest, EContinuumFit::fromFirstPass);
  }

  boost::chrono::thread_clock::time_point stop_secondpass =
      boost::chrono::thread_clock::now();
  Float64 duration_secondpass =
      boost::chrono::duration_cast<boost::chrono::microseconds>(
          stop_secondpass - start_secondpass)
          .count();
  Float64 duration_secondpass_seconds = duration_secondpass / 1e6;
  Log.LogInfo(Formatter() << "  Operator-Linemodel: second pass done in "
                          << duration_secondpass_seconds << " sec");
  Log.LogInfo(Formatter() << "<proc-lm-secondpass><"
                          << (Int32)duration_secondpass_seconds << ">");

  m_fittingManager->getContinuumManager()->SetFitContinuum_Option(
      savedFitContinuumOption);

  return m_result;
}

std::shared_ptr<LineModelExtremaResult>
COperatorLineModel::buildExtremaResults(const TCandidateZbyRank &zCandidates,
                                        const std::string &opt_continuumreest) {
  CContinuumManager::EFitType savedFitContinuumOption =
      m_fittingManager->getContinuumManager()->GetFitContinuum_Option();
  Log.LogInfo("  Operator-Linemodel: Now storing extrema results");

  Int32 extremumCount = zCandidates.size();
  if (extremumCount > m_maxModelSaveCount) {
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "ExtremumCount " << extremumCount
                       << " is greater than maxModelSaveCount "
                       << m_maxModelSaveCount);
  }

  std::shared_ptr<LineModelExtremaResult> ExtremaResult =
      make_shared<LineModelExtremaResult>(zCandidates);

  Int32 savedModels = 0;

  Log.LogDetail(Formatter()
                << "  Operator-Linemodel: N extrema results will be saved : "
                << extremumCount);
  for (Int32 i = 0; i < extremumCount; i++) {
    std::string Id = zCandidates[i].first;
    std::string parentId = zCandidates[i].second->ParentId; // retrieve parentID
    Float64 z = zCandidates[i].second->Redshift;
    //  find the index in the zaxis results
    Int32 idx = CIndexing<Float64>::getIndex(m_result->Redshifts, z);
    Float64 m = m_result->ChiSquare[idx];

    Int32 i_1pass = undefIdx;
    if (parentId != "") { // second pass processed
      for (Int32 j = 0; j != m_firstpass_extremaResult->size(); ++j)
        if (m_firstpass_extremaResult->ID(j) == parentId)
          i_1pass = j;
      if (i_1pass == -1) {
        THROWG(ErrorCode::INTERNAL_ERROR,
               Formatter() << "Impossible to find the first pass extrema "
                              "id corresponding to 2nd pass extrema "
                           << Id);
      }
    }

    Log.LogInfo("");
    Log.LogInfo(Formatter()
                << "  Operator-Linemodel: Saving candidate #" << i << ", idx="
                << idx << ", z=" << m_result->Redshifts[idx] << ", m=" << m);

    m_fittingManager->getContinuumManager()->SetFitContinuum_FitValues(
        m_result->ContinuumModelSolutions[idx]);

    m_fittingManager->getContinuumManager()->SetFitContinuum_Option(
        CContinuumManager::EFitType::fixedValues);

    // reestimate the model (eventually with continuum reestimation) on the
    // extrema selected
    Int32 contreest_iterations = 0;
    if (opt_continuumreest == "always")
      contreest_iterations = 1;
    else if (opt_continuumreest == "onlyextrema") {
      contreest_iterations = 8; // 4
      if (Context.GetParameterStore()->GetScoped<bool>(
              "lineModel.skipSecondPass")) {
        contreest_iterations = 0;
        Flag.warning(WarningCode::FORCED_CONTINUUM_REESTIMATION_TO_NO,
                     "onlyextrema value for ContinuumReestimation is "
                     "unavailable since secondpass is skipped");
      }
    }

    /*    if (m_enableWidthFitByGroups && i_1pass != undefIdx) {
      // Bug: this option requires that prepareSupport is called prior to
      // calling GetModelVelFitGroups, which is not done here. TBSolved in #6623
      // absorption
      std::vector<TInt32List> idxVelfitGroups =
          m_fittingManager->getElementList().GetModelVelfitGroups(
              CLine::EType::nType_Absorption);
      std::string alv_list_str = "";
      for (Int32 kgroup = 0; kgroup < idxVelfitGroups.size(); kgroup++) {
        m_fittingManager->setVelocityAbsorptionByGroup(
            m_firstpass_extremaResult->GroupsALv[i_1pass][kgroup],
            idxVelfitGroups[kgroup]);
        alv_list_str.append(
            boost::str(boost::format("%.2f, ") %
                       m_firstpass_extremaResult->GroupsALv[i_1pass][kgroup]));
      }
      Log.LogInfo(Formatter()
                  << "    Operator-Linemodel: saveResults with groups alv="
                  << alv_list_str);
      // emission
      idxVelfitGroups = m_fittingManager->getElementList().GetModelVelfitGroups(
          CLine::EType::nType_Emission);
      std::string elv_list_str = "";
      for (Int32 kgroup = 0; kgroup < idxVelfitGroups.size(); kgroup++) {
        m_fittingManager->setVelocityEmissionByGroup(
            m_firstpass_extremaResult->GroupsELv[i_1pass][kgroup],
            idxVelfitGroups[kgroup]);
        elv_list_str.append(
            boost::str(boost::format("%.2f") %
                       m_firstpass_extremaResult->GroupsELv[i_1pass][kgroup]));
      }
      Log.LogInfo(Formatter()
                  << "    Operator-Linemodel: saveResults with groups elv="
                  << elv_list_str);
    } else {*/
    // m_fittingManager->SetVelocityEmission(m_firstpass_extremaResult->Elv[i_1pass]);
    // m_fittingManager->SetVelocityAbsorption(m_firstpass_extremaResult->Alv[i_1pass]);
    m_fittingManager->SetVelocityEmission(
        m_result->LineModelSolutions[idx].EmissionVelocity);
    m_fittingManager->SetVelocityAbsorption(
        m_result->LineModelSolutions[idx].AbsorptionVelocity);
    //}

    m_result->ChiSquare[idx] = m_fittingManager->fit(
        m_result->Redshifts[idx], m_result->LineModelSolutions[idx],
        m_result->ContinuumModelSolutions[idx], contreest_iterations, true);
    m_result->ScaleMargCorrection[idx] =
        m_fittingManager->getScaleMargCorrection();
    if (m_fittingManager->getLineRatioType() == "tplRatio")
      m_result->SetChisquareTplratioResult(
          idx, std::dynamic_pointer_cast<CTplratioManager>(
                   m_fittingManager->m_lineRatioManager));
    if (!m_estimateLeastSquareFast)
      m_result->ChiSquareContinuum[idx] =
          m_fittingManager->getLeastSquareContinuumMerit();
    else
      m_result->ChiSquareContinuum[idx] =
          m_fittingManager->getLeastSquareContinuumMeritFast();

    m_result->ScaleMargCorrectionContinuum[idx] =
        m_fittingManager->getContinuumManager()
            ->getContinuumScaleMargCorrection();
    if (m != m_result->ChiSquare[idx] &&
        m_fittingManager->getFittingMethod() != "random")
      THROWG(ErrorCode::INTERNAL_ERROR,
             Formatter() << "COperatorLineModel::" << __func__ << ": m (" << m
                         << " for idx=" << idx << ") !=chi2 ("
                         << m_result->ChiSquare[idx] << ")");
    m = m_result->ChiSquare[idx];
    // save the model result
    // WARNING: saving results TODO: this is currently wrong !! the model
    // saved corresponds to the bestchi2 model. PDFs should be combined
    // prior to exporting the best model for each extrema...
    Int32 maxModelSave = std::min(m_maxModelSaveCount, extremumCount);
    Int32 maxSaveNLinemodelContinua = maxModelSave;
    if (savedModels < maxModelSave) {
      // CModelSpectrumResult
      std::shared_ptr<CModelSpectrumResult> resultspcmodel =
          std::make_shared<CModelSpectrumResult>();
      Int32 overrideModelSavedType = 0;
      // 0=save model, (DEFAULT)
      // 1=save model with lines removed,
      // 2=save model with only Em. lines removed.
      for (auto &obs : m_fittingManager->getSpectraIndex()) {

        if (overrideModelSavedType == 0) {
          resultspcmodel->addModel(
              m_fittingManager->getSpectrumModel().GetModelSpectrum(),
              m_fittingManager->getSpectrum().getObsID());
        } else if (overrideModelSavedType == 1 || overrideModelSavedType == 2) {
          auto lineTypeFilter = CLine::EType::nType_All;
          if (overrideModelSavedType == 2)
            lineTypeFilter = CLine::EType::nType_Emission;

          resultspcmodel->addModel(
              m_fittingManager->getSpectrumModel()
                  .GetObservedSpectrumWithLinesRemoved(lineTypeFilter),
              m_fittingManager->getSpectrum().getObsID());
        }
      }
      ExtremaResult->m_savedModelSpectrumResults[i] = resultspcmodel;

      // below spectrumModel doesnt include identified lines
      auto &cont = m_result->ContinuumModelSolutions[idx];
      TPhotVal phot_values;
      if (!m_opt_continuumcomponent.isTplFitXXX() ||
          cont.name == "noContinuum") { // no photometry
        Log.LogDetail("photometry cannot be applied outside tplFit");
      } else {
        m_fittingManager->getSpectraIndex().setAtBegining();
        phot_values = m_fittingManager->getSpectrumModel().getPhotValues();
      }
      ExtremaResult->m_modelPhotValues[i] =
          std::make_shared<CModelPhotValueResult>(std::move(phot_values));
      ExtremaResult->m_savedModelFittingResults[i] =
          std::make_shared<CLineModelSolution>(
              m_result->LineModelSolutions[idx]);

      // CModelRulesResult
      if (m_fittingManager->getLineRatioType() == "rules") {
        ExtremaResult->m_savedModelRulesResults[i] =
            std::make_shared<CModelRulesResult>(
                std::dynamic_pointer_cast<CRulesManager>(
                    m_fittingManager->m_lineRatioManager)
                    ->GetModelRulesLog());
      }

      std::shared_ptr<CModelSpectrumResult> baselineResult =
          std::make_shared<CModelSpectrumResult>();
      for (auto &obs : m_fittingManager->getSpectraIndex()) {

        // Save the reestimated continuum, only the first
        // n=maxSaveNLinemodelContinua extrema
        baselineResult->addModel(m_fittingManager->getSpectrum()
                                     .GetSpectralAxis()
                                     .GetSamplesVector(),
                                 m_fittingManager->getSpectrumModel()
                                     .GetModelContinuum()
                                     .GetSamplesVector(),
                                 m_fittingManager->getSpectrum().getObsID());
      }
      ExtremaResult->m_savedModelContinuumSpectrumResults[i] = baselineResult;
      savedModels++;
    }
    auto candidate = ExtremaResult->getRankedCandidatePtr(i);
    // code here has been moved to TLineModelResult::updateFromModel
    candidate->updateFromModel(m_fittingManager, m_result,
                               m_estimateLeastSquareFast, idx);

    // save the continuum tpl fitting results

    candidate->updateFromContinuumModelSolution(
        *m_fittingManager->getContinuumFitValues());

    if (m_fittingManager->getLineRatioType() == "tplRatio")
      candidate->updateTplRatioFromModel(
          std::dynamic_pointer_cast<CTplratioManager>(
              m_fittingManager->m_lineRatioManager));
    // save the tplcorr/tplratio results
  }

  // ComputeArea2(ExtremaResult);

  m_fittingManager->getContinuumManager()->SetFitContinuum_Option(
      savedFitContinuumOption);

  return ExtremaResult;
}

/**
 * @brief COperatorLineModel::estimateSecondPassParameters
 * - Estimates best parameters: elv and alv
 * - Store parameters for further use
 *
 *
 * @return
 */
void COperatorLineModel::EstimateSecondPassParameters() {
  // setup velocity fitting

  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();
  const bool enableVelocityFitting =
      ps->GetScoped<bool>("lineModel.velocityFit");
  const std::string &opt_continuumreest =
      ps->GetScoped<std::string>("lineModel.continuumReestimation");

  m_fittingManager->logParameters();
  m_velocitySolutions =
      std::vector<TVelocityFitSolution>(m_firstpass_extremaResult->size());
  for (Int32 i = 0; i < m_firstpass_extremaResult->size(); i++) {
    Log.LogInfo("");
    Log.LogInfo(
        Formatter()
        << "  Operator-Linemodel: Second pass - estimate parameters for "
           "candidate #"
        << i);
    Float64 z = m_firstpass_extremaResult->Redshift(i);
    Float64 m = m_firstpass_extremaResult->ValProba(i);
    Log.LogInfo(Formatter()
                << "  Operator-Linemodel: redshift=" << z << " proba=%" << m);

    if (m_opt_continuumcomponent.isContinuumFit()) {
      // inject continuumFitValues of current candidate
      if (mustReFit(m_continnuum_fit_option))
        m_fittingManager->getContinuumManager()->SetFitContinuum_FitStore(
            m_tplfitStore_secondpass[i]);
      else {
        m_fittingManager->getContinuumManager()->SetFitContinuum_FitStore(
            nullptr);
        m_fittingManager->getContinuumManager()->SetFitContinuum_FitValues(
            m_firstpass_extremaResult->getRankedCandidateCPtr(i)
                ->fittedContinuum);
        m_fittingManager->getContinuumManager()->SetFitContinuum_Option(
            CContinuumManager::EFitType::fixedValues);
      }
    }
    // find the index in the zaxis results
    Int32 idx = CIndexing<Float64>::getIndex(m_result->Redshifts, z);

    // reestimate the model (eventually with continuum reestimation) on
    // the extrema selected
    Int32 contreest_iterations = (opt_continuumreest == "always") ? 1 : 0;

    // model.LoadModelSolution(m_result->LineModelSolutions[idx]);

    m_fittingManager->fit(
        m_result->Redshifts[idx], m_result->LineModelSolutions[idx],
        m_result->ContinuumModelSolutions[idx], contreest_iterations, false);

    if (!enableVelocityFitting) {
      // TODO move these copies to another structure
      //  Emission
      //       Int32 s = m_firstpass_extremaResult->GroupsELv[i].size();
      Float64 emVel = m_result->LineModelSolutions[idx].EmissionVelocity;

      // m_firstpass_extremaResult->GroupsELv[i].assign(s, emVel);

      // Absorption
      //      s = m_firstpass_extremaResult->GroupsALv[i].size();
      Float64 absVel = m_result->LineModelSolutions[idx].AbsorptionVelocity;

      m_velocitySolutions[i] = TVelocityFitSolution{emVel, absVel};
      // m_firstpass_extremaResult->GroupsALv[i].assign(s, absVel);

      continue; // move to the following candidate
    }
    // enableVelocityFitting == true
    fitVelocity(idx, i, contreest_iterations);
  }
}
/**
 * @brief
 *
 * @param Zidx : refers to redshift index on the zGrid
 * @param candidateIdx
 * @param contreest_iterations
 */
void COperatorLineModel::fitVelocity(Int32 Zidx, Int32 candidateIdx,
                                     Int32 contreest_iterations) {
  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();
  const Float64 velfitMinE =
      ps->GetScoped<Float64>("lineModel.emVelocityFitMin");
  const Float64 velfitMaxE =
      ps->GetScoped<Float64>("lineModel.emVelocityFitMax");
  const Float64 velfitStepE =
      ps->GetScoped<Float64>("lineModel.emVelocityFitStep");
  const Float64 velfitMinA =
      ps->GetScoped<Float64>("lineModel.absVelocityFitMin");
  const Float64 velfitMaxA =
      ps->GetScoped<Float64>("lineModel.absVelocityFitMax");
  const Float64 velfitStepA =
      ps->GetScoped<Float64>("lineModel.absVelocityFitStep");
  const std::string opt_lineRatioType =
      ps->GetScoped<std::string>("lineModel.lineRatioType");
  const std::string opt_fittingmethod =
      ps->GetScoped<std::string>("lineModel.fittingMethod");

  // once for all get indices of secondpass interval
  const Int32 half_nb_zsteps = 6;
  const Float64 z_front = m_extendedRedshifts[candidateIdx].front();
  const Float64 z_back = m_extendedRedshifts[candidateIdx].back();
  const Int32 idx_begin =
      CIndexing<Float64>::getIndex(m_result->Redshifts, z_front);
  const Int32 idx_end =
      CIndexing<Float64>::getIndex(m_result->Redshifts, z_back);
  const Int32 lowerzIdx = std::max(idx_begin, Zidx - half_nb_zsteps);
  const Int32 higherzIdx = std::min(idx_end, Zidx + half_nb_zsteps);

  Log.LogInfo(Formatter() << "  Operator-Linemodel: "
                             "velocity fitting bounds for Emission: min="
                          << velfitMinE << " - max=" << velfitMaxE
                          << " - step=" << velfitStepE);
  Log.LogInfo(Formatter() << "  Operator-Linemodel: "
                             "velocity fitting bounds for Absorption: min="
                          << velfitMinA << " - max=" << velfitMaxA
                          << " - "
                             "step="
                          << velfitStepA);

  auto saved_fitter = std::move(m_fittingManager->m_fitter);

  // fit the emission and absorption width by minimizing the
  // linemodel merit with linemodel "hybrid" fitting method
  m_fittingManager->SetFittingMethod("hybrid");
  if (opt_lineRatioType == "tplRatio") {
    m_fittingManager->SetFittingMethod("individual");
    std::dynamic_pointer_cast<CTplratioManager>(
        m_fittingManager->m_lineRatioManager)
        ->SetForcedisableTplratioISMfit(
            std::dynamic_pointer_cast<CTplratioManager>(
                m_fittingManager->m_lineRatioManager)
                ->m_opt_firstpass_forcedisableTplratioISMfit); // TODO: add
  }

  std::vector<TInt32List> idxVelfitGroups;
  Float64 vInfLim;
  Float64 vSupLim;
  Float64 vStep;
  std::string lineTypeStr = undefStr;
  auto lineTypeInt = CLine::EType::nType_All;
  for (Int32 iLineType = 0; iLineType < 2; iLineType++) {
    if (iLineType == 0) {
      vInfLim = velfitMinA;
      vSupLim = velfitMaxA;
      vStep = velfitStepA;
      lineTypeStr = "ABSORPTION";
      lineTypeInt = CLine::EType::nType_Absorption;
    } else {
      vInfLim = velfitMinE;
      vSupLim = velfitMaxE;
      vStep = velfitStepE;
      lineTypeStr = "EMISSION";
      lineTypeInt = CLine::EType::nType_Emission;
    }
    Log.LogDetail(Formatter()
                  << "  Operator-Linemodel: manualStep velocity fit "
                  << lineTypeStr << ", for z = " << m_result->Redshifts[Zidx]);
    if (m_enableWidthFitByGroups) {
      // Note: this option requires that prepareSupport is called prior to
      // calling GetModelVelFitGroups,which is the case here (through a prior
      // call to ::fit)
      idxVelfitGroups =
          m_fittingManager->getElementList().GetModelVelfitGroups(lineTypeInt);
      Log.LogDetail(Formatter()
                    << "  Operator-Linemodel: VelfitGroups " << lineTypeStr
                    << " - n = " << idxVelfitGroups.size());
      if (m_firstpass_extremaResult->size() > 1 && idxVelfitGroups.size() > 1)
        THROWG(ErrorCode::INTERNAL_ERROR,
               "  Operator-Linemodel: not allowed to use more than 1 "
               "group per E/A for "
               "more than 1 extremum (see .json "
               "lineModel.extremaCount)");
    } else {
      // create a dumb vector, condition enough to reach the fitting loop
      // idxVelfitGroups wont be used when m_enableWidthFitByGroups = false
      idxVelfitGroups.assign(1, TInt32List{});
    }
    // Prepare velocity grid to be checked
    const TFloat64List velfitlist = makeVelFitBins(vInfLim, vSupLim, vStep);

    for (Int32 kgroup = 0; kgroup < idxVelfitGroups.size(); kgroup++) {
      Log.LogDetail(Formatter()
                    << "  Operator-Linemodel: manualStep fitting group="
                    << kgroup);

      Float64 meritMin = DBL_MAX;
      Float64 vOptim = NAN;
      Float64 z_vOptim = NAN;

      for (Int32 idzTest = lowerzIdx; idzTest <= higherzIdx; ++idzTest) {
        Float64 zTest = m_result->Redshifts[idzTest];
        for (const auto vTest : velfitlist) {
          if (m_enableWidthFitByGroups) {
            m_fittingManager->setVelocityByGroup(vTest,
                                                 idxVelfitGroups[kgroup]);
          } else {
            if (iLineType == 0) {
              m_fittingManager->SetVelocityAbsorption(vTest);
            } else {
              m_fittingManager->SetVelocityEmission(vTest);
            }
          }
          Float64 meritv = m_fittingManager->fit(
              zTest, m_result->LineModelSolutions[Zidx],
              m_result->ContinuumModelSolutions[Zidx], contreest_iterations,
              false); // maybe m_result members
                      // should be replaced here
                      // by an unused variable

          Log.LogDebug(Formatter() << "  Operator-Linemodel: testing velocity: "
                                      "merit="
                                   << meritv << " for velocity = " << vTest);
          if (meritMin > meritv) {
            meritMin = meritv;
            vOptim = (iLineType == 0)
                         ? m_result->LineModelSolutions[Zidx].AbsorptionVelocity
                         : m_result->LineModelSolutions[Zidx].EmissionVelocity;
            z_vOptim = zTest;
          }
        }
      }
      if (std::isnan(vOptim))
        continue;

      Log.LogDetail(Formatter()
                    << "  Operator-Linemodel: best Velocity found = " << vOptim
                    << " for kgroup: " << kgroup);
      m_result->ChiSquare[Zidx] = meritMin;
      if (iLineType == 0) {
        if (m_enableWidthFitByGroups) {
          // m_fittingManager->setVelocityAbsorptionByGroup(
          //     vOptim, idxVelfitGroups[kgroup]);

          // m_firstpass_extremaResult->GroupsALv[candidateIdx][kgroup] =
          //     vOptim; // not sure yet how to deal with this
        } else
          m_fittingManager->SetVelocityAbsorption(vOptim);
        m_velocitySolutions[candidateIdx].Alv = vOptim;
        Log.LogDebug(
            Formatter()
            << "    Operator-Linemodel: secondpass_parameters extrema #"
            << candidateIdx << " set: alv= " << vOptim);
      } else {
        if (m_enableWidthFitByGroups) {
          // m_fittingManager->setVelocityEmissionByGroup(vOptim,
          //                                              idxVelfitGroups[kgroup]);

          // m_firstpass_extremaResult->GroupsELv[candidateIdx][kgroup] =
          // vOptim;
        } else
          m_fittingManager->SetVelocityEmission(vOptim);
        m_velocitySolutions[candidateIdx].Elv = vOptim;

        Log.LogDebug(Formatter()
                     << "    Operator-Linemodel: secondpass_parameters "
                        "extrema #"
                     << candidateIdx << " set: elv=" << vOptim
                     << " (for z-optim=" << z_vOptim);
      }
    }
  }
  // restore some params
  m_fittingManager->m_fitter = std::move(saved_fitter);
  if (m_fittingManager->getLineRatioType() == "tplRatio")
    std::dynamic_pointer_cast<CTplratioManager>(
        m_fittingManager->m_lineRatioManager)
        ->SetForcedisableTplratioISMfit(
            false); // TODO: coordinate with SetPassMode() ?
  return;
}

void COperatorLineModel::RecomputeAroundCandidates(
    const std::string &opt_continuumreest, const EContinuumFit tplfit_option,
    const bool overrideRecomputeOnlyOnTheCandidate) {

  if (m_firstpass_extremaResult->size() < 1) {
    THROWG(ErrorCode::INTERNAL_ERROR, "ExtremaResult is empty");
  }

  Log.LogInfo("");
  Log.LogInfo(
      Formatter() << "  Operator-Linemodel: Second pass - recomputing around n="
                  << m_firstpass_extremaResult->size() << "candidates");

  for (Int32 i = 0; i < m_firstpass_extremaResult->size(); i++) {
    Log.LogInfo("");
    Log.LogInfo(
        Formatter()
        << "  Operator-Linemodel: Second pass - recompute around Candidate #"
        << i);
    Log.LogInfo(Formatter()
                << "  Operator-Linemodel: ---------- /\\ ---------- ---------- "
                << "---------- Candidate #" << i);

    Float64 Z = m_firstpass_extremaResult->Redshift(i);

    if (m_enableWidthFitByGroups) {
      std::vector<TInt32List> idxVelfitGroups;
      // absorption
      idxVelfitGroups = m_fittingManager->getElementList().GetModelVelfitGroups(
          CLine::EType::nType_Absorption);
      std::string alv_list_str = "";
      for (Int32 kgroup = 0; kgroup < idxVelfitGroups.size(); kgroup++) {
        // m_fittingManager->setVelocityAbsorptionByGroup(
        //     m_firstpass_extremaResult->GroupsALv[i][kgroup],
        //     idxVelfitGroups[kgroup]);
        // alv_list_str.append(
        //     boost::str(boost::format("%.2f, ") %
        //                m_firstpass_extremaResult->GroupsALv[i][kgroup]));
      }
      Log.LogInfo(Formatter()
                  << "    Operator-Linemodel: recompute with groups alv="
                  << alv_list_str);
      // emission
      idxVelfitGroups = m_fittingManager->getElementList().GetModelVelfitGroups(
          CLine::EType::nType_Emission);
      std::string elv_list_str = "";
      for (Int32 kgroup = 0; kgroup < idxVelfitGroups.size(); kgroup++) {
        // m_fittingManager->setVelocityEmissionByGroup(
        // 					     m_firstpass_extremaResult->getRankedCandidate(i)->GroupsELv[kgroup],
        //     idxVelfitGroups[kgroup]);
        // elv_list_str.append(
        //     boost::str(boost::format("%.2f") %
        //                m_firstpass_extremaResult->getRankedCandidate(i)->GroupsELv[kgroup]));
      }
      Log.LogInfo(Formatter()
                  << "    Operator-Linemodel: recompute with groups elv="
                  << elv_list_str);

    } else {
      m_fittingManager->SetVelocityEmission(m_velocitySolutions[i].Elv);
      m_fittingManager->SetVelocityAbsorption(m_velocitySolutions[i].Alv);
      Log.LogInfo(Formatter() << "    Operator-Linemodel: recompute with elv="
                              << m_velocitySolutions[i].Elv
                              << ", alv=" << m_velocitySolutions[i].Alv);
    }

    if (m_opt_continuumcomponent.isContinuumFit()) {
      // fix some fitcontinuum values for this extremum
      if (tplfit_option == EContinuumFit::fromFirstPass) {
        m_fittingManager->getContinuumManager()->SetFitContinuum_FitStore(
            nullptr);
        m_fittingManager->getContinuumManager()->SetFitContinuum_FitValues(
            m_firstpass_extremaResult->getRankedCandidateCPtr(i)
                ->fittedContinuum);
        m_fittingManager->getContinuumManager()->SetFitContinuum_Option(
            CContinuumManager::EFitType::fixedValues);
      } else if (tplfit_option == EContinuumFit::retryAll ||
                 tplfit_option ==
                     EContinuumFit::reFitFirstPass) // for these cases we called
                                                    // precompute in secondpass,
                                                    // so we have new fitstore
        m_fittingManager->getContinuumManager()->SetFitContinuum_FitStore(
            m_tplfitStore_secondpass[i]);
    }

    Log.LogInfo(Formatter()
                << "    Operator-Linemodel: recompute with tplfit_option="
                << static_cast<Int32>(tplfit_option));

    // find the index in the zaxis results
    const Int32 idx = CIndexing<Float64>::getIndex(m_result->Redshifts, Z);

    // reestimate the model (eventually with continuum reestimation) on
    // the extrema selected
    Int32 contreest_iterations = 0;
    if (opt_continuumreest == "always")
      contreest_iterations = 1;
    else if (opt_continuumreest == "onlyextrema")
      contreest_iterations = 8;

    // finally compute the redshifts on the z-range around the extremum
    // m_fittingManager->SetFittingMethod("nofit");
    Int32 n_progresssteps = m_extendedRedshifts[i].size();
    Log.LogInfo(Formatter()
                << "    Operator-Linemodel: Fit n=" << n_progresssteps
                << " values for z in [" << m_extendedRedshifts[i].front()
                << "; " << m_extendedRedshifts[i].back() << "]");
    m_fittingManager->logParameters();
    for (const Float64 z : m_extendedRedshifts[i]) {
      const Int32 iz = CIndexing<Float64>::getIndex(m_result->Redshifts, z);
      Log.LogDetail(Formatter()
                    << "Fit for Extended redshift " << iz << ", z = " << z);

      m_result->ChiSquare[iz] = m_fittingManager->fit(
          m_result->Redshifts[iz], m_result->LineModelSolutions[iz],
          m_result->ContinuumModelSolutions[iz], contreest_iterations, false);
      m_result->ScaleMargCorrection[iz] =
          m_fittingManager->getScaleMargCorrection();
      if (m_opt_continuumcomponent.isContinuumFit()) {
        if (mustReFit(tplfit_option))
          m_result->SetChisquareContinuumResult(
              iz, m_fittingManager->getContinuumManager()
                      ->GetFitContinuum_FitStore());
        // nothing to do when fromfirstpass: keep
        // m_result->ChiSquareTplContinuum from first pass
      }
      if (m_fittingManager->getLineRatioType() == "tplRatio")
        m_result->SetChisquareTplratioResult(
            iz, std::dynamic_pointer_cast<CTplratioManager>(
                    m_fittingManager->m_lineRatioManager));
      if (!m_estimateLeastSquareFast) {
        m_result->ChiSquareContinuum[iz] =
            m_fittingManager->getLeastSquareContinuumMerit();
      } else {
        m_result->ChiSquareContinuum[iz] =
            m_fittingManager->getLeastSquareContinuumMeritFast();
      }
      m_result->ScaleMargCorrectionContinuum[iz] =
          m_fittingManager->getContinuumManager()
              ->getContinuumScaleMargCorrection();
    }
  }
}

void COperatorLineModel::Init(const TFloat64List &redshifts, Float64 zStep,
                              const bool zLogSampling) {

  m_tplCategory = Context.GetCurrentCategory();
  // initialize empty results so that it can be returned anyway in case of an
  // error
  m_result = std::make_shared<CLineModelResult>();

  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();
  CAutoScope autoscope(Context.m_ScopeStack, "lineModel");

  m_opt_continuumcomponent = ps->GetScoped<std::string>("continuumComponent");

  // below should be part of constructor
  m_redshifts = redshifts;
  // init relevant elements to generate secondpass intervals
  m_fineStep = zStep;
  m_zLogSampling = zLogSampling;

  if (Context.GetCurrentMethod() == "lineModelSolve") {
    m_secondPass_halfwindowsize =
        ps->GetScoped<Float64>("secondPass.halfWindowSize");

    m_opt_firstpass_fittingmethod =
        ps->GetScoped<std::string>("firstPass.fittingMethod");
  }
  //
  if (m_opt_continuumcomponent.isContinuumFit()) {

    m_opt_fitcontinuum_maxN = ps->GetScoped<Int32>("continuumFit.count");
    if (m_opt_continuumcomponent.isPowerLawXXX()) {
      m_opt_fitcontinuum_maxN = 1;
    }
    Log.LogDetail(Formatter()
                  << "  method Linemodel wit tplfit: fitcontinuum_maxN set to "
                  << m_opt_fitcontinuum_maxN);

    m_opt_tplfit_fftprocessing =
        ps->GetScoped<bool>("continuumFit.fftProcessing");
    m_opt_tplfit_fftprocessing_secondpass =
        m_opt_tplfit_fftprocessing; // TODO add a real parameter or remove
                                    // this member
    if (ps->HasScoped<bool>("enablePhotometry"))
      m_opt_tplfit_use_photometry = ps->GetScoped<bool>("enablePhotometry");
    m_opt_tplfit_dustFit = ps->GetScoped<bool>("continuumFit.ismFit");
    m_opt_tplfit_extinction = ps->GetScoped<bool>("continuumFit.igmFit");

    m_opt_tplfit_ignoreLinesSupport =
        ps->GetScoped<bool>("continuumFit.ignoreLineSupport");
    if (Context.GetCurrentMethod() == "lineModelSolve") {

      m_opt_secondpasslcfittingmethod =
          ps->GetScoped<std::string>("secondPassLcFittingMethod");
    }
    m_opt_continuum_neg_amp_threshold =
        ps->GetScoped<Float64>("continuumFit.negativeThreshold");

    m_opt_continuum_null_amp_threshold =
        ps->GetScoped<Float64>("continuumFit.nullThreshold");

    m_opt_skipsecondpass = ps->GetScoped<bool>("skipSecondPass");

    m_opt_continuum_bad_chi2_threshold =
        ps->GetScoped<Float64>("continuumFit.badChi2Threshold");
  }
}

/**
 * \brief Returns a non-negative value for the width that yields the least
 *squared difference between the flux and a exponentially decayed maximum
 *amplitude. Find the maximum flux amplitude. If this not greater than zero,
 *return zero. For each value of c within the range: Sum the squared
 *difference between the flux and the maximum amplitude with a exponential
 *decay parameterized by c. Save the minimal result. If the result is not
 *greater than zero, return zero. Return the result.
 **/
Float64
COperatorLineModel::FitBayesWidth(const CSpectrumSpectralAxis &spectralAxis,
                                  const CSpectrumFluxAxis &fluxAxis, Float64 z,
                                  Int32 start, Int32 end) {

  Float64 A = *std::max_element(fluxAxis.GetSamplesVector().begin() + start,
                                fluxAxis.GetSamplesVector().begin() + end);

  if (A <= 0)
    return 0.0;

  // c fitting iteration loop
  Float64 mu = z;
  Float64 c = 0.0001;
  Float64 cmax = 0.05;
  Int32 maxIteration = 500;
  Float64 cstepup = (cmax - c) / ((Float64)(maxIteration + 1));
  Float64 sum2 = boost::numeric::bounds<float>::highest();
  Float64 minsum2 = boost::numeric::bounds<float>::highest();
  Float64 minc = c;
  Int32 icmpt = 0;
  while (icmpt < maxIteration) {
    sum2 = 0.0;
    for (Int32 i = start; i < end; i++) {
      Float64 x = spectralAxis[i];
      Float64 Yi = A * exp(-1. * (x - mu) * (x - mu) / (2 * c * c));
      sum2 += pow(Yi - fluxAxis[i], 2.0);
      // sum2 += pow( Yi - flux[i] , 2.0 ) / pow( error[i], 2.0 );
    }
    if (sum2 < minsum2) {
      minc = c;
      minsum2 = sum2;
    }
    icmpt++;
    c = c + cstepup;
  }

  if (minc < 0)
    minc = 0;

  return minc;
}

CLineModelSolution COperatorLineModel::fitWidthByGroups(
    std::shared_ptr<const CInputContext> context, Float64 redshift) {
  /*  CDataStore &datastore = context.GetDataStore();
  const TFloat64Range &lambdaRange = context.GetLambdaRange();
  Float64 redshift_min =
  datastore.GetFloat64Range("redshiftRange").GetBegin(); Float64 redshift_max
  = datastore.GetFloat64Range("redshiftRange").GetEnd(); CLineModelSolution
  modelSolution; CLineModelSolution bestModelSolution;

  CContinuumModelSolution continuumModelSolution;

  m_fittingManager->SetFittingMethod("hybrid");

    if (opt_lineRatioType == "tplRatio")
    {
    m_fittingManager->SetFittingMethod("individual");
    }

  m_fittingManager->SetForcedisableTplratioISMfit(true);//m_fittingManager->m_opt_firstpass_forcedisableTplratioISMfit);
  //todo, add new param for this ?

  //TODO these params should be in a dedicated scope "velocityFit"
  datastore.PushScope("lineModel");
  const Float64&velfitMinE =
  datastore.GetScopedFloat64Param("emVelocityFitMin"); const
  Float64&velfitMaxE = datastore.GetScopedFloat64Param("emVelocityFitMax");
  const Float64&velfitStepE =
  datastore.GetScopedFloat64Param("emVelocityFitStep"); const
  Float64&velfitMinA = datastore.GetScopedFloat64Param("absVelocityFitMin");
  const Float64&velfitMaxA =
  datastore.GetScopedFloat64Param("absVelocityFitMax"); const
  Float64&velfitStepA = datastore.GetScopedFloat64Param("absVelocityFitStep");
  const Float64&opt_manvelfit_dzmin =
  datastore.GetScopedFloat64Param("manvelocityfitdzmin"); const
  Float64&opt_manvelfit_dzmax =
  datastore.GetScopedFloat64Param("manvelocityfitdzmax"); const
  Float64&opt_manvelfit_dzstep =
  datastore.GetScopedFloat64Param("manvelocityfitdzstep");
  datastore.PopScope();

  //TODO there should be a mapping between this variable and
  opt_continuumreest Int32 contreest_iterations = 0;
  // contreest_iterations = 1;


  //output
  TFloat64List GroupsElv;
  TFloat64List GroupsAlv;
  Float64 Elv,Alv;

  Float64 dzInfLim = roundf(opt_manvelfit_dzmin*10000)/10000;//set precision
  to 10^4 Float64 dzStep = opt_manvelfit_dzstep; Float64 dzSupLim =
  roundf(10000*opt_manvelfit_dzmax)/10000;

  TFloat64List velFitEList =
  TFloat64Range(velfitMinE,velfitMaxE).SpreadOverEpsilon(velfitStepE);
  TFloat64List velFitAList =
  TFloat64Range(velfitMinA,velfitMaxA).SpreadOverEpsilon(velfitStepA);

  // Prepare velocity grid to be checked
  // TODO log grid ?

  if (redshift + dzInfLim < redshift)
    {
      dzInfLim = redshift_min - redshift;
    }
  if (redshift + dzSupLim > redshift_max)
    {
      dzSupLim = redshift_max - redshift;
    }
  TFloat64List zList =
  TFloat64Range(redshift-opt_manvelfit_dzmin,redshift+opt_manvelfit_dzmax).SpreadOverEpsilon(opt_manvelfit_dzstep);

  fitVelocityByGroups(velFitEList,zList,CLine::EType::nType_Emission);
  fitVelocityByGroups(velFitAList,zList,CLine::EType::nType_Absorption);

*/
  CLineModelSolution clms;
  return clms;
}

CLineModelSolution COperatorLineModel::computeForLineMeas(
    std::shared_ptr<const CInputContext> inputContext,
    const TFloat64List &redshiftsGrid, Float64 &bestz) {
  std::shared_ptr<const CParameterStore> params =
      inputContext->GetParameterStore();
  if (params->GetScoped<bool>("lineModel.velocityFit") &&
      params->GetScoped<std::string>("lineModel.fittingMethod") != "lbfgsb")
    THROWG(ErrorCode::INVALID_PARAMETER,
           "velocityFit implemented only for lbfgsb ftting method");

  Int32 amplitudeOffsetsDegree =
      params->GetScoped<Int32>("lineModel.polynomialDegree");
  if (amplitudeOffsetsDegree < 0 || amplitudeOffsetsDegree > 2)
    THROWG(ErrorCode::INVALID_PARAMETER, "the polynomial degree "
                                         "parameter should be between 0 and 2");

  makeContinuumFittingOperator(m_redshifts);

  m_fittingManager = std::make_shared<CLineModelFitting>(
      m_continuumFittingOperator, ElementComposition::OneLine);

  // TODO handle igm coeff

  // does m_fittingManager->m_enableAmplitudeOffsets = true;
  m_fittingManager->setPassMode(3);

  m_estimateLeastSquareFast = 0;
  m_fittingManager->m_lineRatioManager->SetLeastSquareFastEstimationEnabled(
      m_estimateLeastSquareFast);

  m_fittingManager->initDtd();

  CLineModelSolution modelSolution;
  CContinuumModelSolution continuumModelSolution;
  CLineModelSolution bestModelSolution;

  Float64 bestScore = DBL_MAX;
  bestz = NAN;
  for (const Float64 &z : redshiftsGrid) {
    Log.LogDebug(Formatter() << "test with z=" << z);

    Float64 score = m_fittingManager->fit(z, modelSolution,
                                          continuumModelSolution, 0, true);

    if (score < bestScore) {
      bestScore = score;
      bestModelSolution = std::move(modelSolution);
      bestz = z;
    }
  }

  Log.LogInfo(Formatter() << "best z=" << bestz);

  return bestModelSolution;
}

/**
 * @brief Compute spectrum model.
 * TODO: currently it only works for linemeas since we do not the continuum
 * @param z : best redshift
 * @param bestModelSolution : linemodel solution corresponding to the best Z
 * @return std::shared_ptr<const CModelSpectrumResult>
 */
const CSpectrum &COperatorLineModel::getFittedModelWithoutcontinuum(
    const CLineModelSolution &bestModelSolution) {
  // make sure polynom info are correctly set. it s up to refresh model to use
  // these coeffs
  m_fittingManager->LoadModelSolution(bestModelSolution);
  m_fittingManager->refreshAllModels();
  m_fittingManager->getSpectraIndex()
      .setAtBegining(); // TODO dummy implementation, should return all models
  return m_fittingManager->getSpectrumModel().GetModelSpectrum();
}

std::shared_ptr<CContinuumFitStore const> const &
COperatorLineModel::getContinuumFitStoreFirstPass() const {
  // 1 tplfitstore per extrema result
  return m_tplfitStore_firstpass;
}

TFloat64List COperatorLineModel::makeVelFitBins(Float64 vInfLim,
                                                Float64 vSupLim,
                                                Float64 vStep) const {
  TFloat64List velFitBins;
  Int32 nVelSteps = 0;
  nVelSteps = std::floor((vSupLim - vInfLim) / vStep + 1);
  for (Int32 i = 0; i < nVelSteps; ++i)
    velFitBins.push_back(vInfLim + i * vStep);
  return velFitBins;
}