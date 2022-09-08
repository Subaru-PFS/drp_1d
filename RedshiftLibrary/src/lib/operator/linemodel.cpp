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
#include "RedshiftLibrary/operator/linemodel.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/extremum/extremum.h"
#include "RedshiftLibrary/linemodel/lineratiomanager.h"
#include "RedshiftLibrary/linemodel/rulesmanager.h"
#include "RedshiftLibrary/linemodel/templatesfitstore.h"
#include "RedshiftLibrary/linemodel/templatesortho.h"
#include "RedshiftLibrary/linemodel/tplratiomanager.h"
#include "RedshiftLibrary/operator/spectraFluxResult.h"
#include "RedshiftLibrary/operator/templatefitting.h"
#include "RedshiftLibrary/operator/templatefittinglog.h"
#include "RedshiftLibrary/operator/templatefittingresult.h"
#include "RedshiftLibrary/operator/templatefittingwithphot.h"
#include "RedshiftLibrary/spectrum/axis.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "RedshiftLibrary/statistics/deltaz.h"
#include "RedshiftLibrary/statistics/priorhelper.h"

#include "RedshiftLibrary/common/flag.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/common/indexing.h"
#include "RedshiftLibrary/log/log.h"

#include <boost/chrono/thread_clock.hpp>
#include <boost/format.hpp>

#include "RedshiftLibrary/processflow/autoscope.h"
#include "RedshiftLibrary/processflow/context.h"
#include "RedshiftLibrary/processflow/inputcontext.h"
#include "RedshiftLibrary/processflow/parameterstore.h"

#include <algorithm> //std::sort
#include <assert.h>
#include <boost/numeric/conversion/bounds.hpp>
#include <gsl/gsl_interp.h>

using namespace NSEpic;
using namespace std;

void COperatorLineModel::CreateRedshiftLargeGrid(
    Int32 ratio, TFloat64List &largeGridRedshifts) {
  for (Int32 i = 0; i < m_sortedRedshifts.size(); i += ratio) {
    largeGridRedshifts.push_back(m_sortedRedshifts[i]);
  }
  /*for(Int32 i = m_sortedRedshifts.size()-1;i>=0; i-=ratio){
      largeGridRedshifts.push_back(m_sortedRedshifts[i]);
  }
  std::reverse(largeGridRedshifts.begin(), largeGridRedshifts.end());//*/
  if (largeGridRedshifts.empty()) {
    m_enableFastFitLargeGrid = 0;
    Log.LogInfo("  Operator-Linemodel: FastFitLargeGrid auto disabled: "
                "raw %d redshifts will be calculated",
                m_sortedRedshifts.size());
  } else {
    Log.LogInfo("  Operator-Linemodel: FastFitLargeGrid enabled: %d redshifts "
                "will be calculated on the large grid (%d initially)",
                largeGridRedshifts.size(), m_sortedRedshifts.size());
  }
}
/**
 * @brief COperatorLineModel::ComputeFirstPass
 * @return 0=no errors, -1=error
 */
Int32 COperatorLineModel::ComputeFirstPass() {
  const CSpectrum &spectrum = *(Context.GetSpectrum());
  const CSpectrum &logSampledSpectrum =
      *(Context.GetRebinnedSpectrum()); // this is temporary
  std::shared_ptr<const CTemplateCatalog> tplCatalog =
      Context.GetTemplateCatalog();
  const CLineCatalogsTplRatio &tplRatioCatalog =
      *(Context.GetTplRatioCatalog());
  const std::shared_ptr<const CPhotBandCatalog> &photBandCat =
      Context.GetPhotBandCatalog();
  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();

  const Int32 &opt_twosteplargegridstep_ratio =
      ps->GetScoped<Int32>("linemodel.firstpass.largegridstepratio");
  m_opt_continuumcomponent =
      ps->GetScoped<std::string>("linemodel.continuumcomponent");
  TFloat64List largeGridRedshifts;
  // redefine redshift grid
  m_enableFastFitLargeGrid = 0;
  if (opt_twosteplargegridstep_ratio > 1) {
    m_enableFastFitLargeGrid = 1;
    CreateRedshiftLargeGrid(opt_twosteplargegridstep_ratio, largeGridRedshifts);
  } else {
    largeGridRedshifts = m_sortedRedshifts;
  }

  m_fittingManager = std::make_shared<CLineModelFitting>();

  // TODO: check option tplfit
  Int32 nfitcontinuum = 0;
  if (m_opt_continuumcomponent == "tplfit" ||
      m_opt_continuumcomponent == "tplfitauto")
    for (const auto &category : m_tplCategoryList)
      nfitcontinuum += tplCatalog->GetTemplateCount(category);
  m_result->Init(m_sortedRedshifts, Context.getLineVector(), nfitcontinuum,
                 m_fittingManager->getTplratio_count(),
                 m_fittingManager->getTplratio_priors());

  Log.LogInfo("  Operator-Linemodel: initialized");

  // commom between firstpass and secondpass processes
  // TODO not pretty, maybe move fitcontinuum_prior help building to operator
  m_phelperContinuum =
      m_fittingManager->m_continuumManager->SetFitContinuum_PriorHelper();

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
  if (m_opt_continuumcomponent == "tplfit" ||
      m_opt_continuumcomponent == "tplfitauto") {
    tplCatalog->m_orthogonal = 1;
    Log.LogInfo(Formatter() << "Precompuute continuum fit ortho="
                            << tplCatalog->m_orthogonal);
    m_tplfitStore_firstpass = PrecomputeContinuumFit(largeGridRedshifts);
    tplCatalog->m_orthogonal = 0;
  }

  m_result->nSpcSamples = m_fittingManager->getSpcNSamples();
  m_result->dTransposeD = m_fittingManager->getDTransposeD();
  m_result->cstLog = m_fittingManager->getLikelihood_cstLog();

  Int32 contreest_iterations = 0;
  if (ps->GetScoped<std::string>("linemodel.continuumreestimation") ==
      "always") {
    contreest_iterations = 1;
  }

  // Set model parameter: abs lines limit
  Float64 absLinesLimit = 1.0; //-1 to disable, 1.0 is typical
  m_fittingManager->SetAbsLinesLimit(absLinesLimit);
  Log.LogInfo("  Operator-Linemodel: set abs lines limit to %f (ex: -1 means "
              "disabled)",
              absLinesLimit);

  // Set model parameter: continuum least-square estimation fast
  // note: this fast method requires continuum templates and linemodels to be
  // orthogonal. The velfit option turns this trickier...
  m_estimateLeastSquareFast = 0;
  // TODO should be called only with lineRatioType=tplratio
  m_fittingManager->m_lineRatioManager->SetLeastSquareFastEstimationEnabled(
      m_estimateLeastSquareFast);
  Log.LogInfo("  Operator-Linemodel: set estimateLeastSquareFast to %d (ex: "
              "0 means disabled)",
              m_estimateLeastSquareFast);

  //
  TBoolList allAmplitudesZero;
  Int32 indexLargeGrid = 0;
  TFloat64List calculatedLargeGridRedshifts;
  TFloat64List calculatedLargeGridMerits;
  std::vector<TFloat64List> calculatedChiSquareTplratios(
      m_result->ChiSquareTplratios.size());
  std::vector<TFloat64List> calculatedChisquareTplContinuum(
      m_result->ChiSquareTplContinuum.size());

  boost::chrono::thread_clock::time_point start_mainloop =
      boost::chrono::thread_clock::now();

  //#pragma omp parallel for
  m_fittingManager->logParameters();
  for (Int32 i = 0; i < m_result->Redshifts.size(); i++) {
    if (m_enableFastFitLargeGrid == 0 ||
        m_result->Redshifts[i] == largeGridRedshifts[indexLargeGrid]) {
      m_result->ChiSquare[i] = m_fittingManager->fit(
          m_result->Redshifts[i], m_result->LineModelSolutions[i],
          m_result->ContinuumModelSolutions[i], contreest_iterations, false);
      calculatedLargeGridRedshifts.push_back(m_result->Redshifts[i]);
      calculatedLargeGridMerits.push_back(m_result->ChiSquare[i]);
      m_result->ScaleMargCorrection[i] =
          m_fittingManager->getScaleMargCorrection();
      if (m_opt_continuumcomponent == "tplfit" ||
          m_opt_continuumcomponent == "tplfitauto") {
        m_result->SetChisquareTplContinuumResult(i, m_tplfitStore_firstpass);
        for (Int32 k = 0, s = m_result->ChiSquareTplContinuum.size(); k < s;
             k++) {
          calculatedChisquareTplContinuum[k].push_back(
              m_result->ChiSquareTplContinuum[k][i]);
        }
      }
      if (m_fittingManager->getLineRatioType() == "tplratio") {
        m_result->SetChisquareTplratioResult(
            i, dynamic_pointer_cast<CTplratioManager>(
                   m_fittingManager->m_lineRatioManager));
      }
      for (Int32 k = 0, s = m_result->ChiSquareTplratios.size(); k < s; k++) {
        calculatedChiSquareTplratios[k].push_back(
            m_result->ChiSquareTplratios[k][i]);
      }
      if (!m_estimateLeastSquareFast) {
        m_result->ChiSquareContinuum[i] =
            m_fittingManager->getLeastSquareContinuumMerit();
      } else {
        m_result->ChiSquareContinuum[i] =
            m_fittingManager->getLeastSquareContinuumMeritFast();
      }
      m_result->ScaleMargCorrectionContinuum[i] =
          m_fittingManager->m_continuumManager
              ->getContinuumScaleMargCorrection();
      Log.LogDebug("  Operator-Linemodel: Z interval %d: Chi2 = %f", i,
                   m_result->ChiSquare[i]);
      indexLargeGrid++;
      // Log.LogInfo( "\nLineModel Infos: large grid step %d", i);
    } else {
      m_result->ChiSquare[i] = m_result->ChiSquare[i - 1] +
                               1e-2; // these values will be replaced by the
                                     // fine grid interpolation below...
      m_result->ScaleMargCorrection[i] = m_result->ScaleMargCorrection[i - 1];
      m_result->LineModelSolutions[i] = m_result->LineModelSolutions[i - 1];
      m_result->ContinuumModelSolutions[i] =
          m_result->ContinuumModelSolutions[i - 1];
      m_result->SetChisquareTplContinuumResultFromPrevious(i);
      m_result->SetChisquareTplratioResultFromPrevious(i);
      if (!m_estimateLeastSquareFast) {
        m_result->ChiSquareContinuum[i] =
            m_fittingManager->getLeastSquareContinuumMerit();
      } else {
        m_result->ChiSquareContinuum[i] =
            m_fittingManager->getLeastSquareContinuumMeritFast();
      }
      m_result->ScaleMargCorrectionContinuum[i] =
          m_result->ScaleMargCorrectionContinuum[i - 1];
    }
    // Flags on continuum and model amplitudes
    Int32 nbLines = m_result->LineModelSolutions[i].Amplitudes.size();
    bool continuumAmplitudeZero =
        (m_result->ContinuumModelSolutions[i].tplAmplitude <= 0.0);
    bool modelAmplitudesZero = true;
    for (Int32 l = 0; l < nbLines; l++) {
      modelAmplitudesZero =
          (modelAmplitudesZero &&
           m_result->LineModelSolutions[i].Amplitudes[l] <= 0.0);
    }
    allAmplitudesZero.push_back(modelAmplitudesZero && continuumAmplitudeZero);
  }
  // Check if all amplitudes are zero for all z
  bool checkAllAmplitudes =
      AllAmplitudesAreZero(allAmplitudesZero, m_result->Redshifts.size());
  if (checkAllAmplitudes == true) {
    THROWG(INTERNAL_ERROR, "Null amplitudes (continuum & model) at all z");
  }

  // now interpolate large grid merit results onto the fine grid
  if (m_result->Redshifts.size() > calculatedLargeGridMerits.size() &&
      calculatedLargeGridMerits.size() > 1) {
    interpolateLargeGridOnFineGrid(
        calculatedLargeGridRedshifts, m_result->Redshifts,
        calculatedLargeGridMerits, m_result->ChiSquare);
  }

  for (Int32 kt = 0, s = m_result->ChiSquareTplContinuum.size(); kt < s; kt++) {
    if (m_result->Redshifts.size() >
            calculatedChisquareTplContinuum[kt].size() &&
        calculatedChisquareTplContinuum[kt].size() > 1) {
      if (m_opt_continuumcomponent == "tplfit" ||
          m_opt_continuumcomponent == "tplfitauto")
        interpolateLargeGridOnFineGrid(calculatedLargeGridRedshifts,
                                       m_result->Redshifts,
                                       calculatedChisquareTplContinuum[kt],
                                       m_result->ChiSquareTplContinuum[kt]);
    }
  }

  for (Int32 kts = 0, s = m_result->ChiSquareTplratios.size(); kts < s; kts++) {
    if (m_result->Redshifts.size() > calculatedChiSquareTplratios[kts].size() &&
        calculatedChiSquareTplratios[kts].size() > 1) {
      interpolateLargeGridOnFineGrid(
          calculatedLargeGridRedshifts, m_result->Redshifts,
          calculatedChiSquareTplratios[kts], m_result->ChiSquareTplratios[kts]);
    }
  }

  // WARNING: HACK, first pass with continuum from spectrum.
  // model.SetContinuumComponent(opt_continuumcomponent);
  // model.InitFitContinuum();
  //
  boost::chrono::thread_clock::time_point stop_mainloop =
      boost::chrono::thread_clock::now();
  Float64 duration_mainloop =
      boost::chrono::duration_cast<boost::chrono::microseconds>(stop_mainloop -
                                                                start_mainloop)
          .count();
  Float64 duration_firstpass_seconds = duration_mainloop / 1e6;
  Log.LogInfo("  Operator-Linemodel: first-pass done in %.4e sec",
              duration_firstpass_seconds);
  Log.LogInfo("<proc-lm-firstpass><%d>", (Int32)duration_firstpass_seconds);

  return 0;
}

bool COperatorLineModel::AllAmplitudesAreZero(const TBoolList &amplitudesZero,
                                              Int32 nbZ) {
  bool areZero = true;
  for (Int32 iZ = 0; iZ < nbZ; iZ++) {
    areZero = (areZero && amplitudesZero[iZ]);
  }
  return areZero;
}

bool COperatorLineModel::isfftprocessingActive(Int32 redshiftsTplFitCount) {
  bool fftprocessing = m_fittingManager->GetPassNumber() == 1
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
    Int32 candidateIdx, const TFloat64List &redshiftsTplFit,
    const std::vector<CMask> &maskList,
    std::vector<std::shared_ptr<CTemplateFittingResult>>
        &chisquareResultsAllTpl,
    TStringList &chisquareResultsTplName) {
  std::shared_ptr<const CTemplateCatalog> tplCatalog =
      Context.GetTemplateCatalog();
  Float64 overlapThreshold = 1.0;
  std::string opt_interp = "precomputedfinegrid"; //"lin";
  Log.LogDetail("COperatorLineModel::PrecomputeContinuumFit: fitContinuum "
                "opt_interp = %s",
                opt_interp.c_str());
  TInt32List meiksinIndices;
  TInt32List ebmvIndices;
  TFloat64List EbmvCoeffs;
  TTemplateConstRefList tplList;

  bool keepismigm = false;
  if (m_fittingManager->GetPassNumber() == 2 && m_continnuum_fit_option == 3) {
    // case where we only want to refit around the m_opt_fitcontinuum_maxN
    // best continuum from firstpass
    keepismigm = getContinuumInfoFromFirstpassFitStore(
        candidateIdx, meiksinIndices, EbmvCoeffs, ebmvIndices, tplList);
  } else {
    tplList = tplCatalog->GetTemplateList(m_tplCategoryList);
    meiksinIndices.resize(tplList.size(), -1);
    ebmvIndices.resize(tplList.size(), m_opt_tplfit_dustFit ? -10 : -1);
    EbmvCoeffs.resize(tplList.size(), NAN);
  }
  Log.LogDebug(Formatter() << "Processing " << tplList.size() << " templates");

  for (Int32 i = 0; i < tplList.size(); i++) {
    std::string tplname = tplList[i]->GetName();
    Log.LogDebug(Formatter() << "Processing tpl " << tplname);

    CPriorHelper::TPriorZEList zePriorData;
    bool retGetPrior = m_phelperContinuum->GetTplPriorData(
        tplname, redshiftsTplFit, zePriorData);
    if (retGetPrior == false)
      THROWG(INTERNAL_ERROR, "Failed to get prior "
                             "for chi2 continuum precomp fit");

    m_templateFittingOperator->SetRedshifts(redshiftsTplFit);
    auto templatefittingResult =
        std::dynamic_pointer_cast<CTemplateFittingResult>(
            m_templateFittingOperator->Compute(
                tplList[i], overlapThreshold, maskList, opt_interp,
                m_opt_tplfit_extinction, ebmvIndices[i], zePriorData,
                keepismigm, EbmvCoeffs[i], meiksinIndices[i]));

    if (!templatefittingResult) {
      THROWG(INTERNAL_ERROR, Formatter()
                                 << "Failed "
                                    "to compute chisquare value for tpl=%"
                                 << tplname);
    } else {
      chisquareResultsAllTpl.push_back(templatefittingResult);
      chisquareResultsTplName.push_back(tplname);
    }
  }
  return;
}

// get tplName, Meiksin and ISM coeff for all continuum
// returns vectors of these entities
bool COperatorLineModel::getContinuumInfoFromFirstpassFitStore(
    Int32 candidateIdx, TInt32List &meiksinIndices, TFloat64List &EbmvCoeffs,
    TInt32List &ebmvIndices, TTemplateConstRefList &tplList) const {

  if (m_fittingManager->GetPassNumber() != 2 ||
      m_continnuum_fit_option != 3) // not secondpass or not refitfirstpass
    return false;

  if (candidateIdx < 0 ||
      candidateIdx > m_firstpass_extremaResult->size() - 1) {
    THROWG(INTERNAL_ERROR, "Candidate index is out of range");
  }
  std::shared_ptr<const CTemplateCatalog> tplCatalog =
      Context.GetTemplateCatalog();
  bool keepismigm = m_opt_tplfit_dustFit || m_opt_tplfit_extinction;
  // telling that we want to keep the ism and igm indexes

  meiksinIndices.resize(m_opt_fitcontinuum_maxN, -1);
  EbmvCoeffs.resize(m_opt_fitcontinuum_maxN, NAN);
  ebmvIndices.resize(m_opt_fitcontinuum_maxN, -1);
  for (Int32 icontinuum = 0; icontinuum < m_opt_fitcontinuum_maxN;
       icontinuum++) {
    // get the closest lower or equal redshift in coarse grid
    Int32 coarseIdx = m_tplfitStore_firstpass->getClosestLowerRedshiftIndex(
        m_firstpass_extremaResult->m_ranked_candidates[candidateIdx]
            .second->Redshift);

    CTemplatesFitStore::TemplateFitValues fitValue =
        m_tplfitStore_firstpass->GetFitValues(coarseIdx, icontinuum);

    tplList.push_back(
        tplCatalog->GetTemplateByName(m_tplCategoryList, fitValue.tplName));

    if (!keepismigm)
      continue;

    meiksinIndices[icontinuum] = fitValue.igmMeiksinIdx;
    EbmvCoeffs[icontinuum] = fitValue.ismEbmvCoeff;
    // access any template and retrieve the ismcorrection object
    if (m_opt_tplfit_dustFit)
      ebmvIndices[icontinuum] =
          tplCatalog->GetTemplate(m_tplCategoryList[0], 0)
              ->m_ismCorrectionCalzetti->GetEbmvIndex(EbmvCoeffs[icontinuum]);
  }
  return keepismigm;
}

/**
 * Estimate once for all the continuum amplitude which is only dependent
 * from the tplName, ism and igm indexes. This is useful when the model
 * fitting option corresponds to fitting separately the continuum and the
 * lines. In such case, playing with (fit) Lines parameters (velocity, line
 * offsets, lines amps, etc.) do not affect continuum amplitudes.. thus we
 * can save time  by fitting once-for-all the continuum amplitudes, prior to
 * fitting the lines.
 * @candidateIdx@ is also an indicator of pass mode
 * */
std::shared_ptr<CTemplatesFitStore>
COperatorLineModel::PrecomputeContinuumFit(const TFloat64List &redshifts,
                                           Int32 candidateIdx) {
  const CSpectrum &spectrum = *(Context.GetSpectrum());
  const CSpectrum &logSampledSpectrum =
      *(Context.GetRebinnedSpectrum()); // this is temporary
  const std::shared_ptr<const CPhotBandCatalog> &photBandCat =
      Context.GetPhotBandCatalog();
  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();

  std::shared_ptr<const CTemplateCatalog> tplCatalog =
      Context.GetTemplateCatalog();

  bool ignoreLinesSupport =
      ps->GetScoped<bool>("linemodel.continuumfit.ignorelinesupport");
  boost::chrono::thread_clock::time_point start_tplfitprecompute =
      boost::chrono::thread_clock::now();
  Log.LogInfo("COperatorLineModel::PrecomputeContinuumFit: continuum tpl "
              "fitting: min=%.5e, max=%.5e",
              redshifts.front(), redshifts.back());

  std::shared_ptr<CTemplatesFitStore> tplfitStore =
      make_shared<CTemplatesFitStore>(redshifts);
  const TFloat64List &redshiftsTplFit = tplfitStore->GetRedshiftList();
  Log.LogInfo("COperatorLineModel::PrecomputeContinuumFit: continuum tpl "
              "redshift list n=%d",
              redshiftsTplFit.size());

  Int32 n_tplfit = std::min(Int32(redshiftsTplFit.size()), 10);
  for (Int32 i = 0; i < n_tplfit; i++)
    Log.LogDebug("COperatorLineModel::PrecomputeContinuumFit: continuum tpl "
                 "redshift list[%d] = %f",
                 i, redshiftsTplFit[i]);

  bool fftprocessing = isfftprocessingActive(redshiftsTplFit.size());

  bool currentSampling = tplCatalog->m_logsampling;

  Log.LogInfo("COperatorLineModel::PrecomputeContinuumFit: fftprocessing = %d",
              fftprocessing);
  Log.LogDetail(
      "COperatorLineModel::PrecomputeContinuumFit: fitContinuum_dustfit = %d",
      m_opt_tplfit_dustFit);
  Log.LogDetail(
      "COperatorLineModel::PrecomputeContinuumFit: fitContinuum_igm = %d",
      m_opt_tplfit_extinction);

  if (fftprocessing) {
    if (m_templateFittingOperator == nullptr ||
        !m_templateFittingOperator->IsFFTProcessing()) // else reuse the shared
                                                       // pointer for secondpass
      m_templateFittingOperator = std::make_shared<COperatorTemplateFittingLog>(
          spectrum, logSampledSpectrum, *(Context.GetLambdaRange()),
          redshiftsTplFit);
    tplCatalog->m_logsampling = true;

  } else {
    if (m_templateFittingOperator == nullptr ||
        m_templateFittingOperator
            ->IsFFTProcessing()) { // else reuse the shared
                                   // pointer for secondpass
      if (m_opt_tplfit_use_photometry)
        m_templateFittingOperator =
            std::make_shared<COperatorTemplateFittingPhot>(
                spectrum, *(Context.GetLambdaRange()), photBandCat,
                ps->GetScoped<Float64>("linemodel.photometry.weight"),
                redshiftsTplFit);
      else
        m_templateFittingOperator = std::make_shared<COperatorTemplateFitting>(
            spectrum, *(Context.GetLambdaRange()), redshiftsTplFit);
    }

    tplCatalog->m_logsampling = false;
  }

  if (fftprocessing && ignoreLinesSupport == true) {
    ignoreLinesSupport = false;
    Flag.warning(Flag.IGNORELINESSUPPORT_DISABLED_FFT,
                 Formatter() << "  COperatorLineModel::" << __func__
                             << ": unable to ignoreLinesSupport if "
                                "fftprocessing. ignoreLinesSupport disabled");
  }
  std::vector<CMask> maskList;
  if (ignoreLinesSupport) {
    boost::chrono::thread_clock::time_point start_tplfitmaskprep =
        boost::chrono::thread_clock::now();

    maskList.resize(redshiftsTplFit.size());
    for (Int32 i = 0; i < redshiftsTplFit.size(); i++) {
      m_fittingManager->initModelAtZ(redshiftsTplFit[i],
                                     fftprocessing
                                         ? logSampledSpectrum.GetSpectralAxis()
                                         : spectrum.GetSpectralAxis());
      maskList[i] = m_fittingManager->getOutsideLinesMask();
    }

    boost::chrono::thread_clock::time_point stop_tplfitmaskprep =
        boost::chrono::thread_clock::now();
    Float64 duration_tplfitmaskprep =
        boost::chrono::duration_cast<boost::chrono::microseconds>(
            stop_tplfitmaskprep - start_tplfitmaskprep)
            .count();
    Float64 duration_tplfitmaskprep_seconds = duration_tplfitmaskprep / 1e6;
    Log.LogInfo("COperatorLineModel::PrecomputeContinuumFit: tplfit-precompute "
                "mask prep done in %.4e sec",
                duration_tplfitmaskprep_seconds);
  }

  std::vector<std::shared_ptr<CTemplateFittingResult>> chisquareResultsAllTpl;
  TStringList chisquareResultsTplName;
  fitContinuumTemplates(candidateIdx, redshiftsTplFit, maskList,
                        chisquareResultsAllTpl, chisquareResultsTplName);

  // fill the fit store with fitted values: only the best fitted values FOR
  // EACH TEMPLATE are used
  Float64 bestTplFitSNR = 0.0;
  Int32 nredshiftsTplFitResults = redshiftsTplFit.size();
  for (Int32 i = 0; i < nredshiftsTplFitResults; i++) {
    Float64 redshift = redshiftsTplFit[i];

    for (Int32 j = 0; j < chisquareResultsAllTpl.size(); j++) {
      const auto &chisquareResult = chisquareResultsAllTpl[j];

      bool retAdd = tplfitStore->Add(
          chisquareResultsTplName[j], chisquareResult->FitEbmvCoeff[i],
          chisquareResult->FitMeiksinIdx[i], redshift,
          chisquareResult->ChiSquare[i], chisquareResult->ChiSquarePhot[i],
          chisquareResult->FitAmplitude[i],
          chisquareResult->FitAmplitudeError[i],
          chisquareResult->FitAmplitudeSigma[i], chisquareResult->FitDtM[i],
          chisquareResult->FitMtM[i], chisquareResult->LogPrior[i]);

      if (!retAdd)
        THROWG(INTERNAL_ERROR, "Failed to add continuum fit to store");

      Float64 tplfitsnr = -1.;
      if (chisquareResult->FitMtM[i] > 0.)
        tplfitsnr =
            chisquareResult->FitDtM[i] / std::sqrt(chisquareResult->FitMtM[i]);

      if (tplfitsnr > bestTplFitSNR)
        bestTplFitSNR = tplfitsnr;
    }
  }
  tplfitStore->m_fitContinuum_tplFitSNRMax = bestTplFitSNR;
  Log.LogDetail("COperatorLineModel::PrecomputeContinuumFit: "
                "fitcontinuum_snrMAX set to %f",
                bestTplFitSNR);
  Log.LogDetail(
      "COperatorLineModel::PrecomputeContinuumFit: continuumcount set to %d",
      tplfitStore->GetContinuumCount());

  if (tplfitStore->GetContinuumCount() < m_opt_fitcontinuum_maxN) {
    THROWG(INTERNAL_ERROR, "Couldn't compute the required continuum count");
  }
  tplfitStore->m_opt_fitcontinuum_maxCount = m_opt_fitcontinuum_maxN;
  Log.LogInfo("COperatorLineModel::PrecomputeContinuumFit: fitStore with "
              "fitcontinuum_maxCount set to %f",
              tplfitStore->m_opt_fitcontinuum_maxCount);

  m_fittingManager->m_continuumManager->SetFitContinuum_FitStore(tplfitStore);
  tplCatalog->m_logsampling = currentSampling;

  boost::chrono::thread_clock::time_point stop_tplfitprecompute =
      boost::chrono::thread_clock::now();
  Float64 duration_tplfitprecompute =
      boost::chrono::duration_cast<boost::chrono::microseconds>(
          stop_tplfitprecompute - start_tplfitprecompute)
          .count();
  Float64 duration_tplfit_seconds = duration_tplfitprecompute / 1e6;
  Log.LogInfo("COperatorLineModel::PrecomputeContinuumFit: tplfit-precompute "
              "done in %.4e sec",
              duration_tplfit_seconds);
  Log.LogDetail("<proc-lm-tplfit><%d>", (Int32)duration_tplfit_seconds);

  evaluateContinuumAmplitude(tplfitStore);

  return tplfitStore;
}

void COperatorLineModel::evaluateContinuumAmplitude(
    const std::shared_ptr<CTemplatesFitStore> &tplfitStore) {
  // Check if best continuum amplitudes are negative fitted amplitudes at
  // all z
  CTemplatesFitStore::TemplateFitValues fitValues;
  Float64 max_fitamplitudeSigma_z = NAN;
  Float64 max_fitamplitudeSigma =
      tplfitStore->FindMaxAmplitudeSigma(max_fitamplitudeSigma_z, fitValues);
  if (max_fitamplitudeSigma < m_opt_continuum_neg_amp_threshold) {
    if (m_opt_continuumcomponent != "tplfitauto") {
      THROWG(INTERNAL_ERROR,
             Formatter() << "Negative "
                            "continuum amplitude found at z="
                         << max_fitamplitudeSigma_z << ": best continuum tpl "
                         << fitValues.tplName << ", amplitude/error = "
                         << fitValues.fitAmplitudeSigma
                         << " & error = " << fitValues.fitAmplitudeError);
    } else {
      Flag.warning(Flag.FORCE_FROMSPECTRUM_NEG_CONTINUUMAMP,
                   Formatter()
                       << ": Switching to spectrum continuum since Negative "
                          "continuum amplitude found at z="
                       << max_fitamplitudeSigma_z << ": best continuum tpl "
                       << fitValues.tplName
                       << ", amplitude/error = " << fitValues.fitAmplitudeSigma
                       << " & error = " << fitValues.fitAmplitudeError);
      m_opt_continuumcomponent = "fromspectrum";
      m_fittingManager->SetContinuumComponent("fromspectrum");
    }
  }
  // check if continuum is too weak comparing to the preset threshold
  if (std::abs(fitValues.fitAmplitudeSigma) <
      m_opt_continuum_null_amp_threshold) {
    Flag.warning(Flag.FORCE_NOCONTINUUM_WEAK_CONTINUUMAMP,
                 Formatter()
                     << ": Switching to nocontinuum since close to null "
                        "continuum amplitude found at z="
                     << max_fitamplitudeSigma_z << ": best continuum tpl "
                     << fitValues.tplName
                     << ", amplitude/error = " << fitValues.fitAmplitudeSigma
                     << " & error = " << fitValues.fitAmplitudeError);
    m_opt_continuumcomponent = "nocontinuum";
    m_fittingManager->SetContinuumComponent("nocontinuum");
  }
}
/**
 * @brief COperatorLineModel::SpanRedshiftWindow
 * @param z
 * @return extendedList
 */
TFloat64List COperatorLineModel::SpanRedshiftWindow(Float64 z) const {
  TFloat64List extendedList;

  const Float64 halfwindowsize_z = m_secondPass_halfwindowsize * (1. + z);
  TFloat64Range secondpass_window = {z - halfwindowsize_z,
                                     z + halfwindowsize_z};
  Int32 i_min, i_max;
  bool ret = secondpass_window.getClosedIntervalIndices(m_result->Redshifts,
                                                        i_min, i_max);
  if (!ret) {
    THROWG(INTERNAL_ERROR, "Second pass window outside z range");
  }
  for (Int32 i = i_min; i <= i_max; ++i) {
    extendedList.push_back(m_result->Redshifts[i]);
  }

  return extendedList;
}

Int32 COperatorLineModel::SetFirstPassCandidates(
    const TCandidateZbyRank &zCandidates) {
  m_firstpass_extremaResult =
      std::make_shared<CLineModelPassExtremaResult>(zCandidates.size());

  m_firstpass_extremaResult->m_ranked_candidates = zCandidates;

  // extend z around the extrema
  for (Int32 j = 0; j < zCandidates.size(); j++) {
    std::string Id = zCandidates[j].first;
    const std::shared_ptr<const TCandidateZ> &cand = zCandidates[j].second;

    Log.LogInfo("  Operator-Linemodel: Raw extr #%d, z_e.X=%f, m_e.Y=%e", j,
                cand->Redshift, cand->ValProba);

    const Float64 &x = cand->Redshift;
    const TFloat64List extendedList = SpanRedshiftWindow(x);
    m_firstpass_extremaResult->ExtendedRedshifts[j] = extendedList;
  }

  // now preparing the candidates extrema results
  for (Int32 i = 0; i < m_firstpass_extremaResult->size(); i++) {
    // find the index in the zaxis results
    Int32 idx = CIndexing<Float64>::getIndex(m_result->Redshifts,
                                             zCandidates[i].second->Redshift);

    // save basic fitting info from first pass
    m_firstpass_extremaResult->Elv[i] =
        m_result->LineModelSolutions[idx].EmissionVelocity;
    m_firstpass_extremaResult->Alv[i] =
        m_result->LineModelSolutions[idx].AbsorptionVelocity;

    // save the continuum fitting parameters from first pass
    m_firstpass_extremaResult->FittedTplName[i] =
        m_result->ContinuumModelSolutions[idx].tplName;
    m_firstpass_extremaResult->FittedTplAmplitude[i] =
        m_result->ContinuumModelSolutions[idx].tplAmplitude;
    m_firstpass_extremaResult->FittedTplAmplitudeError[i] =
        m_result->ContinuumModelSolutions[idx].tplAmplitudeError;
    m_firstpass_extremaResult->FittedTplMerit[i] =
        m_result->ContinuumModelSolutions[idx].tplMerit;
    m_firstpass_extremaResult->FittedTplMeritPhot[i] =
        m_result->ContinuumModelSolutions[idx].tplMeritPhot;
    m_firstpass_extremaResult->FittedTplEbmvCoeff[i] =
        m_result->ContinuumModelSolutions[idx].tplEbmvCoeff;
    m_firstpass_extremaResult->FittedTplMeiksinIdx[i] =
        m_result->ContinuumModelSolutions[idx].tplMeiksinIdx;
    m_firstpass_extremaResult->FittedTplRedshift[i] =
        m_result->ContinuumModelSolutions[idx].tplRedshift;
    m_firstpass_extremaResult->FittedTplDtm[i] =
        m_result->ContinuumModelSolutions[idx].tplDtm;
    m_firstpass_extremaResult->FittedTplMtm[i] =
        m_result->ContinuumModelSolutions[idx].tplMtm;
    m_firstpass_extremaResult->FittedTplLogPrior[i] =
        m_result->ContinuumModelSolutions[idx].tplLogPrior;
    m_firstpass_extremaResult->FittedTplpCoeffs[i] =
        m_result->ContinuumModelSolutions[idx].pCoeffs;

    //... TODO: more first pass results can be saved here if needed
  }

  return 0;
}

std::shared_ptr<const LineModelExtremaResult>
COperatorLineModel::buildFirstPassExtremaResults(
    const TCandidateZbyRank &zCandidates) {
  std::shared_ptr<LineModelExtremaResult> ExtremaResult =
      make_shared<LineModelExtremaResult>(zCandidates);

  for (Int32 i = 0; i < m_firstpass_extremaResult->size(); i++) {
    // find the index in the zaxis results
    Int32 idx = CIndexing<Float64>::getIndex(m_result->Redshifts,
                                             zCandidates[i].second->Redshift);

    CContinuumModelSolution csolution = m_result->ContinuumModelSolutions[idx];
    ExtremaResult->m_ranked_candidates[i]
        .second->updateFromContinuumModelSolution(csolution, true);
    // ExtremaResult->setCandidateFromContinuumSolution(i, csolution);

    // for saving velocities: use CLineModelSolution
    ExtremaResult->m_ranked_candidates[i].second->updateFromLineModelSolution(
        m_result->LineModelSolutions[idx]);

    ExtremaResult->m_savedModelFittingResults[i] =
        std::make_shared<CLineModelSolution>(m_result->LineModelSolutions[idx]);
  }

  return ExtremaResult;
}

/**
 * @brief Aggregate candidates from both firstpasses, while avoiding duplicate
 * candidates
 *
 * @param firstpass_results_b
 */
void COperatorLineModel::Combine_firstpass_candidates(
    std::shared_ptr<const CLineModelPassExtremaResult> firstpass_results_b) {
  Int32 startIdx = m_firstpass_extremaResult->size();
  TInt32List uniqueIdx_fpb =
      m_firstpass_extremaResult->getUniqueCandidates(firstpass_results_b);
  m_firstpass_extremaResult->Resize(m_firstpass_extremaResult->size() +
                                    uniqueIdx_fpb.size());

  for (Int32 keb = 0; keb < uniqueIdx_fpb.size(); keb++) {
    Int32 i = uniqueIdx_fpb[keb];
    // append the candidate to m_firstpass_extremaResult
    m_firstpass_extremaResult->m_ranked_candidates[startIdx + keb] =
        firstpass_results_b->m_ranked_candidates[i];

    // extend z around the extrema
    const Float64 &z_fpb = firstpass_results_b->Redshift(i);
    m_firstpass_extremaResult->ExtendedRedshifts[startIdx + keb] =
        (SpanRedshiftWindow(z_fpb));

    // save basic fitting info from first pass
    m_firstpass_extremaResult->Elv[startIdx + keb] =
        firstpass_results_b->Elv[i];
    m_firstpass_extremaResult->Alv[startIdx + keb] =
        firstpass_results_b->Alv[i];

    // save the continuum fitting parameters from first pass
    if (0) // cannot work since the fpb is linemodel without cont. tplfit...
    {
      m_firstpass_extremaResult->FittedTplName[startIdx + keb] =
          firstpass_results_b->FittedTplName[i];
      m_firstpass_extremaResult->FittedTplAmplitude[startIdx + keb] =
          firstpass_results_b->FittedTplAmplitude[i];
      m_firstpass_extremaResult->FittedTplAmplitudeError[startIdx + keb] =
          firstpass_results_b->FittedTplAmplitudeError[i];
      m_firstpass_extremaResult->FittedTplMerit[startIdx + keb] =
          firstpass_results_b->FittedTplMerit[i];
      m_firstpass_extremaResult->FittedTplMeritPhot[startIdx + keb] =
          firstpass_results_b->FittedTplMeritPhot[i];
      m_firstpass_extremaResult->FittedTplEbmvCoeff[startIdx + keb] =
          firstpass_results_b->FittedTplEbmvCoeff[i];
      m_firstpass_extremaResult->FittedTplMeiksinIdx[startIdx + keb] =
          firstpass_results_b->FittedTplMeiksinIdx[i];
      m_firstpass_extremaResult->FittedTplRedshift[startIdx + keb] =
          firstpass_results_b->FittedTplRedshift[i];
      m_firstpass_extremaResult->FittedTplDtm[startIdx + keb] =
          firstpass_results_b->FittedTplDtm[i];
      m_firstpass_extremaResult->FittedTplMtm[startIdx + keb] =
          firstpass_results_b->FittedTplMtm[i];
      m_firstpass_extremaResult->FittedTplLogPrior[startIdx + keb] =
          firstpass_results_b->FittedTplLogPrior[i];
      m_firstpass_extremaResult->FittedTplpCoeffs[startIdx + keb] =
          firstpass_results_b->FittedTplpCoeffs[i];
    } else {
      // find the index in the zaxis results

      Int32 idx = CIndexing<Float64>::getIndex(m_result->Redshifts, z_fpb);

      // save the continuum fitting parameters from first pass
      m_firstpass_extremaResult->FittedTplName[startIdx + keb] =
          m_result->ContinuumModelSolutions[idx].tplName;
      m_firstpass_extremaResult->FittedTplAmplitude[startIdx + keb] =
          m_result->ContinuumModelSolutions[idx].tplAmplitude;
      m_firstpass_extremaResult->FittedTplAmplitudeError[startIdx + keb] =
          m_result->ContinuumModelSolutions[idx].tplAmplitudeError;
      m_firstpass_extremaResult->FittedTplMerit[startIdx + keb] =
          m_result->ContinuumModelSolutions[idx].tplMerit;
      m_firstpass_extremaResult->FittedTplMeritPhot[startIdx + keb] =
          m_result->ContinuumModelSolutions[idx].tplMeritPhot;
      m_firstpass_extremaResult->FittedTplEbmvCoeff[startIdx + keb] =
          m_result->ContinuumModelSolutions[idx].tplEbmvCoeff;
      m_firstpass_extremaResult->FittedTplMeiksinIdx[startIdx + keb] =
          m_result->ContinuumModelSolutions[idx].tplMeiksinIdx;
      m_firstpass_extremaResult->FittedTplRedshift[startIdx + keb] =
          m_result->ContinuumModelSolutions[idx].tplRedshift;
      m_firstpass_extremaResult->FittedTplDtm[startIdx + keb] =
          m_result->ContinuumModelSolutions[idx].tplDtm;
      m_firstpass_extremaResult->FittedTplMtm[startIdx + keb] =
          m_result->ContinuumModelSolutions[idx].tplMtm;
      m_firstpass_extremaResult->FittedTplLogPrior[startIdx + keb] =
          m_result->ContinuumModelSolutions[idx].tplLogPrior;
      m_firstpass_extremaResult->FittedTplpCoeffs[startIdx + keb] =
          m_result->ContinuumModelSolutions[idx].pCoeffs;

      if (m_result->ContinuumModelSolutions[idx].tplName == "") {
        Flag.warning(Flag.TPL_NAME_EMPTY,
                     Formatter() << "COperatorLineModel::" << __func__
                                 << ": Saving first pass extremum w. "
                                    "ContinuumModelSolutions tplname="
                                 << m_result->ContinuumModelSolutions[idx]
                                        .tplName.c_str());
        Flag.warning(Flag.TPL_NAME_EMPTY,
                     Formatter()
                         << "COperatorLineModel::" << __func__
                         << ": Saving first pass extremum w. result idx=" << idx
                         << ", w. m_result->Redshifts[idx]="
                         << m_result->Redshifts[idx]);
      }
    }
  }
}

Int32 COperatorLineModel::ComputeSecondPass(
    const std::shared_ptr<const LineModelExtremaResult> &firstpassResults) {

  const CSpectrum &spectrum = *(Context.GetSpectrum());
  const CSpectrum &logSampledSpectrum =
      *(Context.GetRebinnedSpectrum()); // this is temporary
  std::shared_ptr<const CTemplateCatalog> tplCatalog =
      Context.GetTemplateCatalog();
  const CLineCatalogsTplRatio &tplRatioCatalog =
      *(Context.GetTplRatioCatalog());
  const std::shared_ptr<const CPhotBandCatalog> &photBandCat =
      Context.GetPhotBandCatalog();
  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();

  boost::chrono::thread_clock::time_point start_secondpass =
      boost::chrono::thread_clock::now();
  // Set model parameters to SECOND-PASS
  m_fittingManager->setPassMode(2);
  Int32 savedFitContinuumOption =
      m_fittingManager->m_continuumManager
          ->GetFitContinuum_Option(); // the first time was set in
                                      // precomputeContinuumFit
  Log.LogInfo(
      "  Operator-Linemodel: ---------- ---------- ---------- ----------");
  Log.LogInfo("  Operator-Linemodel: now computing second-pass");
  Log.LogInfo(
      "  Operator-Linemodel: ---------- ---------- ---------- ----------");

  // init lmfit variables
  mlmfit_modelInfoSave = false;
  mlmfit_savedModelSpectrumResults_lmfit.clear();
  mlmfit_savedModelFittingResults_lmfit.clear();
  mlmfit_savedModelRulesResults_lmfit.clear();
  mlmfit_savedBaselineResult_lmfit.clear();

  std::string opt_continuumfit_method =
      ps->GetScoped<std::string>("linemodel.secondpass.continuumfit");
  std::string opt_continuumreest =
      ps->GetScoped<std::string>("linemodel.continuumreestimation");
  std::string opt_fittingmethod =
      ps->GetScoped<std::string>("linemodel.fittingmethod");
  m_continnuum_fit_option = 0;
  if (opt_continuumfit_method == "fromfirstpass") {
    m_continnuum_fit_option = 2;
  } else if (opt_continuumfit_method == "retryall") {
    m_continnuum_fit_option = 0;
  } else if (opt_continuumfit_method == "refitfirstpass") {
    m_continnuum_fit_option = 3;
  } else {
    // TODO this should be a parameterException thrown at parameter setting
    // stage
    THROWG(INTERNAL_ERROR, Formatter() << "Invalid continnuum_fit_option: "
                                       << m_continnuum_fit_option);
  }
  if (m_opt_continuumcomponent == "tplfit" ||
      m_opt_continuumcomponent == "tplfitauto") {
    // precompute only whenever required and whenever the result can be a
    // tplfitStore
    if (m_continnuum_fit_option == 0 || m_continnuum_fit_option == 3) {
      m_tplfitStore_secondpass.resize(m_firstpass_extremaResult->size());
      tplCatalog->m_orthogonal = 1;
      for (Int32 i = 0; i < m_firstpass_extremaResult->size(); i++) {
        m_tplfitStore_secondpass[i] = PrecomputeContinuumFit(
            m_firstpass_extremaResult->ExtendedRedshifts[i], i);
        if (m_opt_continuumcomponent == "fromspectrum")
          break; // when set to "fromspectrum" by PrecomputeContinuumFit
                 // because negative continuum with tplfitauto
      }
      tplCatalog->m_orthogonal = 0; // finish using orthogTemplates
    } else {
      // since precompute is not called all the time, secondpass candidates do
      // not have systematically a tplfitstore_secondpass copy the firstpass
      // tplfitstore into the secondpass tplfitstore
      m_tplfitStore_secondpass.resize(1);
      if (m_continnuum_fit_option == 1 || m_continnuum_fit_option == 2)
        m_tplfitStore_secondpass[0] = m_tplfitStore_firstpass;
      m_fittingManager->m_continuumManager->SetFitContinuum_FitStore(
          m_tplfitStore_firstpass);
      // duplicated: double make sure that these info are present in the
      // modelElement
      m_fittingManager->m_continuumManager->SetFitContinuum_SNRMax(
          m_tplfitStore_firstpass->m_fitContinuum_tplFitSNRMax);
      m_fittingManager->m_continuumManager->m_opt_fitcontinuum_maxCount =
          m_tplfitStore_firstpass->m_opt_fitcontinuum_maxCount;
    }
  }

  m_secondpass_parameters_extremaResult.Resize(
      m_firstpass_extremaResult->size());
  if (m_continnuum_fit_option == 2) {
    for (Int32 i = 0; i < m_firstpass_extremaResult->size(); i++) {
      m_secondpass_parameters_extremaResult.FittedTplName[i] =
          m_firstpass_extremaResult->FittedTplName[i];
      m_secondpass_parameters_extremaResult.FittedTplAmplitude[i] =
          m_firstpass_extremaResult->FittedTplAmplitude[i];
      m_secondpass_parameters_extremaResult.FittedTplAmplitudeError[i] =
          m_firstpass_extremaResult->FittedTplAmplitudeError[i];
      m_secondpass_parameters_extremaResult.FittedTplMerit[i] =
          m_firstpass_extremaResult->FittedTplMerit[i];
      m_secondpass_parameters_extremaResult.FittedTplEbmvCoeff[i] =
          m_firstpass_extremaResult->FittedTplEbmvCoeff[i];
      m_secondpass_parameters_extremaResult.FittedTplMeiksinIdx[i] =
          m_firstpass_extremaResult->FittedTplMeiksinIdx[i];
      m_secondpass_parameters_extremaResult.FittedTplDtm[i] =
          m_firstpass_extremaResult->FittedTplDtm[i];
      m_secondpass_parameters_extremaResult.FittedTplMtm[i] =
          m_firstpass_extremaResult->FittedTplMtm[i];
      m_secondpass_parameters_extremaResult.FittedTplLogPrior[i] =
          m_firstpass_extremaResult->FittedTplLogPrior[i];
      m_secondpass_parameters_extremaResult.FittedTplRedshift[i] =
          m_firstpass_extremaResult->FittedTplRedshift[i];
      m_secondpass_parameters_extremaResult.FittedTplpCoeffs[i] =
          m_firstpass_extremaResult->FittedTplpCoeffs[i];
    }
  }

  // upcast LineModelExtremaResult to TCandidateZ
  m_secondpass_parameters_extremaResult.m_ranked_candidates.assign(
      firstpassResults->m_ranked_candidates.cbegin(),
      firstpassResults->m_ranked_candidates.cend());

  // now that we recomputed what should be recomputed, we define once for all
  // the secondpass
  //  estimate second pass parameters (mainly elv, alv...)
  EstimateSecondPassParameters(spectrum, *(Context.GetClampedLambdaRange()));

  // recompute the fine grid results around the extrema
  Int32 ret = RecomputeAroundCandidates(
      opt_continuumreest,
      m_continnuum_fit_option); // 0: retry all cont. templates at this stage

  // additional fitting with fittingmethod=svdlcp2
  if (m_opt_secondpasslcfittingmethod == "svdlc" ||
      m_opt_secondpasslcfittingmethod == "svdlcp2") {
    Log.LogInfo("  Operator-Linemodel: now computing second-pass %s on each "
                "secondpass candidate (n=%d)",
                m_opt_secondpasslcfittingmethod.c_str(),
                m_secondpass_parameters_extremaResult.size());
    bool useSecondPassRedshiftValue = true;
    if (useSecondPassRedshiftValue) {
      for (Int32 i = 0; i < m_secondpass_parameters_extremaResult.size(); i++) {
        m_secondpass_parameters_extremaResult.FittedTplRedshift[i] =
            m_secondpass_parameters_extremaResult.Redshift(i);
      }
    }
    m_fittingManager->SetFittingMethod(m_opt_secondpasslcfittingmethod);
    RecomputeAroundCandidates(opt_continuumreest, 2, true);
    m_fittingManager->SetFittingMethod(opt_fittingmethod);

    Log.LogInfo("  Operator-Linemodel: now re-computing the final chi2 for "
                "each candidate");
    RecomputeAroundCandidates(opt_continuumreest, 2);
  }

  boost::chrono::thread_clock::time_point stop_secondpass =
      boost::chrono::thread_clock::now();
  Float64 duration_secondpass =
      boost::chrono::duration_cast<boost::chrono::microseconds>(
          stop_secondpass - start_secondpass)
          .count();
  Float64 duration_secondpass_seconds = duration_secondpass / 1e6;
  Log.LogInfo("  Operator-Linemodel: second pass done in %.4e sec",
              duration_secondpass_seconds);
  Log.LogInfo("<proc-lm-secondpass><%d>", (Int32)duration_secondpass_seconds);

  m_fittingManager->m_continuumManager->SetFitContinuum_Option(
      savedFitContinuumOption);

  return 0;
}

std::shared_ptr<LineModelExtremaResult>
COperatorLineModel::buildExtremaResults(const CSpectrum &spectrum,
                                        const TFloat64Range &lambdaRange,
                                        const TCandidateZbyRank &zCandidates,
                                        const std::string &opt_continuumreest) {
  Int32 savedFitContinuumOption =
      m_fittingManager->m_continuumManager->GetFitContinuum_Option();
  Log.LogInfo("  Operator-Linemodel: Now storing extrema results");

  Int32 extremumCount = zCandidates.size();
  if (extremumCount > m_maxModelSaveCount) {
    THROWG(INTERNAL_ERROR, Formatter() << "ExtremumCount " << extremumCount
                                       << " is greater than maxModelSaveCount "
                                       << m_maxModelSaveCount);
  }

  std::shared_ptr<LineModelExtremaResult> ExtremaResult =
      make_shared<LineModelExtremaResult>(zCandidates);

  Int32 savedModels = 0;

  Log.LogDetail("  Operator-Linemodel: N extrema results will be saved : %d",
                extremumCount);
  for (Int32 i = 0; i < extremumCount; i++) {
    std::string Id = zCandidates[i].first;
    std::string parentId = zCandidates[i].second->ParentId; // retrieve parentID
    Float64 z = zCandidates[i].second->Redshift;

    // find the index in the zaxis results
    Int32 idx = CIndexing<Float64>::getIndex(m_result->Redshifts, z);
    Float64 m = m_result->ChiSquare[idx];

    Int32 i_2pass = -1;
    for (Int32 j = 0; j != m_secondpass_parameters_extremaResult.size(); ++j)
      if (m_secondpass_parameters_extremaResult.ID(j) == parentId)
        i_2pass = j;
    if (i_2pass == -1) {
      THROWG(INTERNAL_ERROR, Formatter()
                                 << "Impossible to find the first pass extrema "
                                    "id corresponding to 2nd pass extrema "
                                 << Id.c_str());
    }

    Log.LogInfo("");
    Log.LogInfo(
        "  Operator-Linemodel: Saving candidate #%d, idx=%d, z=%f, m=%f", i,
        idx, m_result->Redshifts[idx], m);

    m_fittingManager->m_continuumManager->SetFitContinuum_FitValues(
        m_result->ContinuumModelSolutions[idx].tplName,
        m_result->ContinuumModelSolutions[idx].tplAmplitude,
        m_result->ContinuumModelSolutions[idx].tplAmplitudeError,
        m_result->ContinuumModelSolutions[idx].tplMerit,
        m_result->ContinuumModelSolutions[idx].tplMeritPhot,
        m_result->ContinuumModelSolutions[idx].tplEbmvCoeff,
        m_result->ContinuumModelSolutions[idx].tplMeiksinIdx,
        m_result->ContinuumModelSolutions[idx].tplRedshift,
        m_result->ContinuumModelSolutions[idx].tplDtm,
        m_result->ContinuumModelSolutions[idx].tplMtm,
        m_result->ContinuumModelSolutions[idx].tplLogPrior,
        m_result->ContinuumModelSolutions[idx].pCoeffs);

    m_fittingManager->m_continuumManager->SetFitContinuum_Option(2);

    // reestimate the model (eventually with continuum reestimation) on the
    // extrema selected
    Int32 contreest_iterations = 0;
    if (opt_continuumreest == "always" || opt_continuumreest == "onlyextrema") {
      contreest_iterations = 8; // 4
    }

    if (m_enableWidthFitByGroups) {
      std::vector<TInt32List> idxVelfitGroups;
      // absorption
      idxVelfitGroups.clear();
      idxVelfitGroups = m_fittingManager->m_Elements.GetModelVelfitGroups(
          CLine::nType_Absorption);
      std::string alv_list_str = "";
      for (Int32 kgroup = 0; kgroup < idxVelfitGroups.size(); kgroup++) {
        for (Int32 ke = 0; ke < idxVelfitGroups[kgroup].size(); ke++) {
          m_fittingManager->SetVelocityAbsorptionOneElement(
              m_secondpass_parameters_extremaResult.GroupsALv[i_2pass][kgroup],
              idxVelfitGroups[kgroup][ke]);
        }
        alv_list_str.append(boost::str(
            boost::format("%.2f, ") %
            m_secondpass_parameters_extremaResult.GroupsALv[i_2pass][kgroup]));
      }
      Log.LogInfo("    Operator-Linemodel: saveResults with groups alv=%s",
                  alv_list_str.c_str());
      // emission
      idxVelfitGroups.clear();
      idxVelfitGroups = m_fittingManager->m_Elements.GetModelVelfitGroups(
          CLine::nType_Emission);
      std::string elv_list_str = "";
      for (Int32 kgroup = 0; kgroup < idxVelfitGroups.size(); kgroup++) {
        for (Int32 ke = 0; ke < idxVelfitGroups[kgroup].size(); ke++) {
          m_fittingManager->SetVelocityEmissionOneElement(
              m_secondpass_parameters_extremaResult.GroupsELv[i_2pass][kgroup],
              idxVelfitGroups[kgroup][ke]);
        }
        elv_list_str.append(boost::str(
            boost::format("%.2f") %
            m_secondpass_parameters_extremaResult.GroupsELv[i_2pass][kgroup]));
      }
      Log.LogInfo("    Operator-Linemodel: saveResults with groups elv=%s",
                  elv_list_str.c_str());

    } else {
      // m_fittingManager->SetVelocityEmission(m_secondpass_parameters_extremaResult.Elv[i_2pass]);
      // m_fittingManager->SetVelocityAbsorption(m_secondpass_parameters_extremaResult.Alv[i_2pass]);
      m_fittingManager->SetVelocityEmission(
          m_result->LineModelSolutions[idx].EmissionVelocity);
      m_fittingManager->SetVelocityAbsorption(
          m_result->LineModelSolutions[idx].AbsorptionVelocity);
    }

    if (!mlmfit_modelInfoSave) {
      m_result->ChiSquare[idx] = m_fittingManager->fit(
          m_result->Redshifts[idx], m_result->LineModelSolutions[idx],
          m_result->ContinuumModelSolutions[idx], contreest_iterations, true);
      m_result->ScaleMargCorrection[idx] =
          m_fittingManager->getScaleMargCorrection();
      if (m_fittingManager->getLineRatioType() == "tplratio")
        m_result->SetChisquareTplratioResult(
            idx, std::dynamic_pointer_cast<CTplratioManager>(
                     m_fittingManager->m_lineRatioManager));
      if (!m_estimateLeastSquareFast) {
        m_result->ChiSquareContinuum[idx] =
            m_fittingManager->getLeastSquareContinuumMerit();
      } else {
        m_result->ChiSquareContinuum[idx] =
            m_fittingManager->getLeastSquareContinuumMeritFast();
      }
      m_result->ScaleMargCorrectionContinuum[idx] =
          m_fittingManager->m_continuumManager
              ->getContinuumScaleMargCorrection();
    }
    if (m != m_result->ChiSquare[idx]) {
      Log.LogWarning("  Operator-Linemodel: m (%f for idx=%d) !=chi2 (%f) ", m,
                     idx, m_result->ChiSquare[idx]);
    }
    m = m_result->ChiSquare[idx]; // m_result->ChiSquare[idx];

    // save the model result
    // WARNING: saving results TODO: this is currently wrong !! the model
    // saved corresponds to the bestchi2 model. PDFs should be combined
    // prior to exporting the best model for each extrema...
    Int32 maxModelSave = std::min(m_maxModelSaveCount, extremumCount);
    Int32 maxSaveNLinemodelContinua = maxModelSave;
    if (savedModels < maxModelSave) {
      if (mlmfit_modelInfoSave) {

        Log.LogInfo("Save model store during lm_fit");
        ExtremaResult->m_savedModelSpectrumResults[i] =
            mlmfit_savedModelSpectrumResults_lmfit[i_2pass];
        ExtremaResult->m_savedModelFittingResults[i] =
            mlmfit_savedModelFittingResults_lmfit[i_2pass];
        ExtremaResult->m_savedModelRulesResults[i] =
            mlmfit_savedModelRulesResults_lmfit[i_2pass];
        if (savedModels < maxSaveNLinemodelContinua &&
            contreest_iterations > 0) {
          ExtremaResult->m_savedModelContinuumSpectrumResults[i] =
              mlmfit_savedBaselineResult_lmfit[i_2pass];
        }

      } else {
        // CModelSpectrumResult
        std::shared_ptr<CModelSpectrumResult> resultspcmodel;
        Int32 overrideModelSavedType = 0;
        // 0=save model, (DEFAULT)
        // 1=save model with lines removed,
        // 2=save model with only Em. lines removed.
        if (overrideModelSavedType == 0) {
          resultspcmodel = std::make_shared<CModelSpectrumResult>(
              m_fittingManager->getSpectrumModel()->GetModelSpectrum());
        } else if (overrideModelSavedType == 1 || overrideModelSavedType == 2) {
          Int32 lineTypeFilter = -1;
          if (overrideModelSavedType == 1) {
            lineTypeFilter = -1;
          } else if (overrideModelSavedType == 2) {
            lineTypeFilter = CLine::nType_Emission;
          }
          resultspcmodel = std::make_shared<CModelSpectrumResult>(
              m_fittingManager->getSpectrumModel()
                  ->GetObservedSpectrumWithLinesRemoved(lineTypeFilter));
        }
        // std::shared_ptr<CModelSpectrumResult>  resultspcmodel =
        // std::shared_ptr<CModelSpectrumResult>( new
        // CModelSpectrumResult(m_fittingManager->GetSpectrumModelContinuum())
        // );

        ExtremaResult->m_savedModelSpectrumResults[i] = resultspcmodel;

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

        // Save the reestimated continuum, only the first
        // n=maxSaveNLinemodelContinua extrema
        std::shared_ptr<CSpectraFluxResult> baselineResult =
            (std::shared_ptr<CSpectraFluxResult>)new CSpectraFluxResult();
        const CSpectrumFluxAxis &modelContinuumFluxAxis =
            m_fittingManager->getSpectrumModel()->GetModelContinuum();
        Int32 len = modelContinuumFluxAxis.GetSamplesCount();

        baselineResult->fluxes.resize(len);
        baselineResult->wavel.resize(len);
        for (Int32 k = 0; k < len; k++) {
          baselineResult->fluxes[k] = modelContinuumFluxAxis[k];
          baselineResult->wavel[k] = (spectrum.GetSpectralAxis())[k];
        }
        ExtremaResult->m_savedModelContinuumSpectrumResults[i] = baselineResult;
      }
      savedModels++;
    }

    // code here has been moved to TLineModelResult::updateFromModel
    ExtremaResult->m_ranked_candidates[i].second->updateFromModel(
        m_fittingManager, m_result, m_estimateLeastSquareFast, idx, i_2pass);

    // save the continuum tpl fitting results
    ExtremaResult->m_ranked_candidates[i].second->updateContinuumFromModel(
        m_fittingManager);

    CContinuumModelSolution csolution =
        m_fittingManager->m_continuumManager->GetContinuumModelSolution();
    ExtremaResult->m_ranked_candidates[i]
        .second->updateFromContinuumModelSolution(csolution, false);
    if (m_fittingManager->getLineRatioType() == "tplratio")
      ExtremaResult->m_ranked_candidates[i].second->updateTplRatioFromModel(
          std::dynamic_pointer_cast<CTplratioManager>(
              m_fittingManager->m_lineRatioManager));
    // save the tplcorr/tplratio results
  }

  // ComputeArea2(ExtremaResult);

  m_fittingManager->m_continuumManager->SetFitContinuum_Option(
      savedFitContinuumOption);

  return ExtremaResult;
}

/**
 * @brief COperatorLineModel::estimateSecondPassParameters
 * - Estimates best parameters: elv and alv
 * - Store parameters for further use into:
 *
 *
 * @return
 */
Int32 COperatorLineModel::EstimateSecondPassParameters(
    const CSpectrum &spectrum, const TFloat64Range &lambdaRange) {
  // setup velocity fitting

  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();

  const bool &enableVelocityFitting =
      ps->GetScoped<bool>("linemodel.velocityfit");
  const Float64 &velfitMinE =
      ps->GetScoped<Float64>("linemodel.emvelocityfitmin");
  const Float64 &velfitMaxE =
      ps->GetScoped<Float64>("linemodel.emvelocityfitmax");
  const Float64 &velfitStepE =
      ps->GetScoped<Float64>("linemodel.emvelocityfitstep");
  const Float64 &velfitMinA =
      ps->GetScoped<Float64>("linemodel.absvelocityfitmin");
  const Float64 &velfitMaxA =
      ps->GetScoped<Float64>("linemodel.absvelocityfitmax");
  const Float64 &velfitStepA =
      ps->GetScoped<Float64>("linemodel.absvelocityfitstep");
  const std::string &opt_continuumreest =
      ps->GetScoped<std::string>("linemodel.continuumreestimation");
  const std::string opt_fittingmethod =
      ps->GetScoped<std::string>("linemodel.fittingmethod");
  const std::string opt_lineRatioType =
      ps->GetScoped<std::string>("linemodel.lineRatioType");
  // HARDCODED - override: no-velocityfitting for abs
  // velfitMinA = opt_velocityAbsorption;
  // velfitMaxA = opt_velocityAbsorption;
  if (enableVelocityFitting) {
    Log.LogInfo(
        "  Operator-Linemodel: "
        "velocity fitting bounds for Emission: min=%.1f - max=%.1f - step=%.1f",
        velfitMinE, velfitMaxE, velfitStepE);
    Log.LogInfo("  Operator-Linemodel: "
                "velocity fitting bounds for Absorption: min=%.1f - max=%.1f - "
                "step=%.1f",
                velfitMinA, velfitMaxA, velfitStepA);
  }

  // enable/disable fit by groups. Once enabled, the velocity fitting groups
  // are defined in the line catalog from v4.0 on.
  m_enableWidthFitByGroups = true;
  m_fittingManager->logParameters();
  for (Int32 i = 0; i < m_firstpass_extremaResult->size(); i++) {
    Log.LogInfo("");
    Log.LogInfo("  Operator-Linemodel: Second pass - estimate parameters for "
                "candidate #%d",
                i);
    Float64 z = m_firstpass_extremaResult->Redshift(i);
    Float64 m = m_firstpass_extremaResult->ValProba(i);
    Log.LogInfo("  Operator-Linemodel: redshift=%.2f proba=%2f", z, m);

    m_secondpass_parameters_extremaResult.ExtendedRedshifts[i] =
        m_firstpass_extremaResult->ExtendedRedshifts[i];
    if (m_opt_continuumcomponent == "tplfit" ||
        m_opt_continuumcomponent == "tplfitauto") {
      // inject continuumFitValues of current candidate
      if (m_continnuum_fit_option == 0 || m_continnuum_fit_option == 3)
        m_fittingManager->m_continuumManager->SetFitContinuum_FitStore(
            m_tplfitStore_secondpass[i]);
      else {
        m_fittingManager->m_continuumManager->SetFitContinuum_FitStore(nullptr);
        m_fittingManager->m_continuumManager->SetFitContinuum_FitValues(
            m_firstpass_extremaResult->FittedTplName[i],
            m_firstpass_extremaResult->FittedTplAmplitude[i],
            m_firstpass_extremaResult->FittedTplAmplitudeError[i],
            m_firstpass_extremaResult->FittedTplMerit[i],
            m_firstpass_extremaResult->FittedTplMeritPhot[i],
            m_firstpass_extremaResult->FittedTplEbmvCoeff[i],
            m_firstpass_extremaResult->FittedTplMeiksinIdx[i],
            m_firstpass_extremaResult->FittedTplRedshift[i],
            m_firstpass_extremaResult->FittedTplDtm[i],
            m_firstpass_extremaResult->FittedTplMtm[i],
            m_firstpass_extremaResult->FittedTplLogPrior[i],
            m_firstpass_extremaResult->FittedTplpCoeffs[i]);
        m_fittingManager->m_continuumManager->SetFitContinuum_Option(2);
      }
    }
    // find the index in the zaxis results
    Int32 idx = -1;
    idx = CIndexing<Float64>::getIndex(m_result->Redshifts, z);

    // reestimate the model (eventually with continuum reestimation) on
    // the extrema selected
    Int32 contreest_iterations = 0;
    if (opt_continuumreest == "always") {
      contreest_iterations = 1;
    }

    // model.LoadModelSolution(m_result->LineModelSolutions[idx]);
    m_fittingManager->fit(
        m_result->Redshifts[idx], m_result->LineModelSolutions[idx],
        m_result->ContinuumModelSolutions[idx], contreest_iterations, false);
    // m = m_result->ChiSquare[idx];

    if (enableVelocityFitting) {
      // fit the emission and absorption width by minimizing the
      // linemodel merit with linemodel "hybrid" fitting method
      m_fittingManager->SetFittingMethod("hybrid");
      if (opt_lineRatioType == "tplratio") {
        m_fittingManager->SetFittingMethod("individual");
        std::dynamic_pointer_cast<CTplratioManager>(
            m_fittingManager->m_lineRatioManager)
            ->SetForcedisableTplratioISMfit(
                std::dynamic_pointer_cast<CTplratioManager>(
                    m_fittingManager->m_lineRatioManager)
                    ->m_opt_firstpass_forcedisableTplratioISMfit); // TODO: add
      }
      // new param
      // for this ?
      // m_fittingManager->m_enableAmplitudeOffsets = true;
      // contreest_iterations = 1;
      std::vector<TInt32List> idxVelfitGroups;
      for (Int32 iLineType = 0; iLineType < 2; iLineType++) {
        Float64 vInfLim;
        Float64 vSupLim;
        Float64 vStep;

        if (iLineType == 0) {
          Log.LogDetail("  Operator-Linemodel: manualStep velocity fit "
                        "ABSORPTION, for z = %.6f",
                        m_result->Redshifts[idx]);
          vInfLim = velfitMinA;
          vSupLim = velfitMaxA;
          vStep = velfitStepA;
          if (m_enableWidthFitByGroups) {
            idxVelfitGroups.clear();
            idxVelfitGroups = m_fittingManager->m_Elements.GetModelVelfitGroups(
                CLine::nType_Absorption);
            Log.LogDetail(
                "  Operator-Linemodel: VelfitGroups ABSORPTION - n = %d",
                idxVelfitGroups.size());
            if (m_firstpass_extremaResult->size() > 1 &&
                idxVelfitGroups.size() > 1) {
              Log.LogError(
                  "  Operator-Linemodel: not allowed to use more than 1 "
                  "group per E/A for "
                  "more than 1 extremum (see .json "
                  "linemodel.extremacount)");
            }
          }
        } else {
          Log.LogInfo("  Operator-Linemodel: manualStep velocity fit "
                      "EMISSION, for z = %.6f",
                      m_result->Redshifts[idx]);
          vInfLim = velfitMinE;
          vSupLim = velfitMaxE;
          vStep = velfitStepE;
          if (m_enableWidthFitByGroups) {
            idxVelfitGroups.clear();
            idxVelfitGroups = m_fittingManager->m_Elements.GetModelVelfitGroups(
                CLine::nType_Emission);
            Log.LogDetail(
                "  Operator-Linemodel: VelfitGroups EMISSION - n = %d",
                idxVelfitGroups.size());
            if (m_firstpass_extremaResult->size() > 1 &&
                idxVelfitGroups.size() > 1) {
              Log.LogError("  Operator-Linemodel: not allowed to use more "
                           "than 1 group per E/A for more than 1 extremum "
                           "(see .json linemodel.extremacount)");
            }
          }
        }

        // Prepare velocity grid to be checked
        TFloat64List velfitlist;
        Int32 optVelfit = 0; // lin
        // Int32 optVelfit = 1; //log todo ?
        if (optVelfit == 0) {
          Int32 nStepsLin = (int)((vSupLim - vInfLim) / vStep);
          for (Int32 kv = 0; kv < nStepsLin; kv++) {
            velfitlist.push_back(vInfLim + kv * vStep);
          }
        }
        Int32 nVelSteps = velfitlist.size();
        /*
          for (Int32 kv = 0; kv < nVelSteps; kv++)
          {
          Log.LogDetail("  Operator-Linemodel: velstep %d = %f", kv,
          velfitlist[kv]);
          }
        //*/

        for (Int32 kgroup = 0; kgroup < idxVelfitGroups.size(); kgroup++) {
          Log.LogDetail("  Operator-Linemodel: manualStep fitting group=%d",
                        kgroup);

          Float64 meritMin = DBL_MAX;
          Float64 vOptim = -1.0;
          Float64 z_vOptim = -1.0;
          const Int32 half_nb_zsteps = 6;
          Float64 z_front =
              m_secondpass_parameters_extremaResult.ExtendedRedshifts[i]
                  .front();
          Float64 z_back =
              m_secondpass_parameters_extremaResult.ExtendedRedshifts[i].back();
          Int32 idx_begin =
              CIndexing<Float64>::getIndex(m_result->Redshifts, z_front);
          Int32 idx_end =
              CIndexing<Float64>::getIndex(m_result->Redshifts, z_back);
          Int32 lowerzIdx = std::max(idx_begin, idx - half_nb_zsteps);
          Int32 higherzIdx = std::min(idx_end, idx + half_nb_zsteps);
          for (Int32 idzTest = lowerzIdx; idzTest <= higherzIdx; ++idzTest) {
            for (Int32 kv = 0; kv < nVelSteps; kv++) {
              Float64 vTest = velfitlist[kv];
              if (iLineType == 0) {
                if (m_enableWidthFitByGroups) {
                  for (Int32 ke = 0; ke < idxVelfitGroups[kgroup].size();
                       ke++) {
                    m_fittingManager->SetVelocityAbsorptionOneElement(
                        vTest, idxVelfitGroups[kgroup][ke]);
                  }
                } else {
                  m_fittingManager->SetVelocityAbsorption(vTest);
                }
              } else {
                if (m_enableWidthFitByGroups) {
                  for (Int32 ke = 0; ke < idxVelfitGroups[kgroup].size();
                       ke++) {
                    m_fittingManager->SetVelocityEmissionOneElement(
                        vTest, idxVelfitGroups[kgroup][ke]);
                  }
                } else {
                  m_fittingManager->SetVelocityEmission(vTest);
                }
              }

              // Log.LogInfo( "  Operator-Linemodel:
              // testing v=%f", vTest);
              Float64 meritv;
              Float64 zTest = m_result->Redshifts[idzTest];
              meritv = m_fittingManager->fit(
                  zTest,
                  m_result->LineModelSolutions[idx], // maybe this member result
                  // should be replaced by an
                  // unused variable
                  m_result->ContinuumModelSolutions
                      [idx], // maybe this member result should be replaced
                  // by an unused variable
                  contreest_iterations, false);

              //                                    if(m_enableWidthFitByGroups)
              //                                    {
              //                                        meritv
              //                                        =
              //                                        0.0;
              //                                        for(Int32
              //                                        ke=0;
              //                                        ke<idxVelfitGroups[kgroup].size();
              //                                        ke++)
              //                                        {
              //                                            meritv +=
              //                                            m_fittingManager->getModelErrorUnderElement(idxVelfitGroups[kgroup][ke]);
              //                                        }
              //                                    }

              Log.LogDebug("  Operator-Linemodel: testing velocity: "
                           "merit=%.3e for velocity = %.1f",
                           meritv, vTest);
              if (meritMin > meritv) {
                meritMin = meritv;
                if (iLineType == 0) {
                  vOptim = m_fittingManager->GetVelocityAbsorption();
                  z_vOptim = zTest;
                } else {
                  vOptim = m_fittingManager->GetVelocityEmission();
                  z_vOptim = zTest;
                }
              }
            }
          }

          if (vOptim != -1.0) {
            Log.LogDetail("  Operator-Linemodel: best Velocity found = %.1f",
                          vOptim);
            m_result->ChiSquare[idx] = meritMin;
            if (iLineType == 0) {
              if (m_enableWidthFitByGroups) {
                for (Int32 ke = 0; ke < idxVelfitGroups[kgroup].size(); ke++) {
                  m_fittingManager->SetVelocityAbsorptionOneElement(
                      vOptim, idxVelfitGroups[kgroup][ke]);
                }
                m_secondpass_parameters_extremaResult.GroupsALv[i][kgroup] =
                    vOptim;
              } else {
                m_fittingManager->SetVelocityAbsorption(vOptim);
              }

              m_secondpass_parameters_extremaResult.Alv[i] = vOptim;
              Log.LogDebug("    Operator-Linemodel: secondpass_parameters "
                           "extrema #%d set: alv=%.1f",
                           i, vOptim);
            } else {
              if (m_enableWidthFitByGroups) {
                for (Int32 ke = 0; ke < idxVelfitGroups[kgroup].size(); ke++) {
                  m_fittingManager->SetVelocityEmissionOneElement(
                      vOptim, idxVelfitGroups[kgroup][ke]);
                }
                m_secondpass_parameters_extremaResult.GroupsELv[i][kgroup] =
                    vOptim;
              } else {
                m_fittingManager->SetVelocityEmission(vOptim);
              }
              m_secondpass_parameters_extremaResult.Elv[i] = vOptim;
              Log.LogDebug("    Operator-Linemodel: secondpass_parameters "
                           "extrema #%d set: elv=%.1f (for z-optim=%.6f",
                           i, vOptim, z_vOptim);
            }
          }
        }
      }

      // restore some params
      m_fittingManager->SetFittingMethod(opt_fittingmethod);
      // m_fittingManager->m_enableAmplitudeOffsets = false;
      if (m_fittingManager->getLineRatioType() == "tplratio") {
        std::dynamic_pointer_cast<CTplratioManager>(
            m_fittingManager->m_lineRatioManager)
            ->SetForcedisableTplratioISMfit(
                false); // TODO: coordinate with SetPassMode() ?
      }
    } else {
      m_secondpass_parameters_extremaResult.Elv[i] =
          m_fittingManager->GetVelocityEmission();
      m_secondpass_parameters_extremaResult.Alv[i] =
          m_fittingManager->GetVelocityAbsorption();
      for (Int32 kg = 0;
           kg < m_secondpass_parameters_extremaResult.GroupsELv[i].size();
           kg++) {
        m_secondpass_parameters_extremaResult.GroupsELv[i][kg] =
            m_secondpass_parameters_extremaResult.Elv[i];
      }
      for (Int32 kg = 0;
           kg < m_secondpass_parameters_extremaResult.GroupsALv[i].size();
           kg++) {
        m_secondpass_parameters_extremaResult.GroupsALv[i][kg] =
            m_secondpass_parameters_extremaResult.Alv[i];
      }
    }
  }

  return 0;
}

Int32 COperatorLineModel::RecomputeAroundCandidates(
    const std::string &opt_continuumreest, const Int32 tplfit_option,
    const bool overrideRecomputeOnlyOnTheCandidate) {
  CLineModelPassExtremaResult &extremaResult =
      m_secondpass_parameters_extremaResult;
  if (extremaResult.size() < 1) {
    THROWG(INTERNAL_ERROR, "ExtremaResult is empty");
  }

  Log.LogInfo("");
  Log.LogInfo("  Operator-Linemodel: Second pass - recomputing around n=%d "
              "candidates",
              extremaResult.size());

  for (Int32 i = 0; i < extremaResult.size(); i++) {
    Log.LogInfo("");
    Log.LogInfo(
        "  Operator-Linemodel: Second pass - recompute around Candidate #%d",
        i);
    Log.LogInfo("  Operator-Linemodel: ---------- /\\ ---------- ---------- "
                "---------- Candidate #%d",
                i);
    Float64 Z = extremaResult.Redshift(i);

    if (m_enableWidthFitByGroups) {
      std::vector<TInt32List> idxVelfitGroups;
      // absorption
      idxVelfitGroups.clear();
      idxVelfitGroups = m_fittingManager->m_Elements.GetModelVelfitGroups(
          CLine::nType_Absorption);
      std::string alv_list_str = "";
      for (Int32 kgroup = 0; kgroup < idxVelfitGroups.size(); kgroup++) {
        for (Int32 ke = 0; ke < idxVelfitGroups[kgroup].size(); ke++) {
          m_fittingManager->SetVelocityAbsorptionOneElement(
              extremaResult.GroupsALv[i][kgroup], idxVelfitGroups[kgroup][ke]);
        }
        alv_list_str.append(boost::str(boost::format("%.2f, ") %
                                       extremaResult.GroupsALv[i][kgroup]));
      }
      Log.LogInfo("    Operator-Linemodel: recompute with groups alv=%s",
                  alv_list_str.c_str());
      // emission
      idxVelfitGroups.clear();
      idxVelfitGroups = m_fittingManager->m_Elements.GetModelVelfitGroups(
          CLine::nType_Emission);
      std::string elv_list_str = "";
      for (Int32 kgroup = 0; kgroup < idxVelfitGroups.size(); kgroup++) {
        for (Int32 ke = 0; ke < idxVelfitGroups[kgroup].size(); ke++) {
          m_fittingManager->SetVelocityEmissionOneElement(
              extremaResult.GroupsELv[i][kgroup], idxVelfitGroups[kgroup][ke]);
        }
        elv_list_str.append(boost::str(boost::format("%.2f") %
                                       extremaResult.GroupsELv[i][kgroup]));
      }
      Log.LogInfo("    Operator-Linemodel: recompute with groups elv=%s",
                  elv_list_str.c_str());

    } else {
      m_fittingManager->SetVelocityEmission(extremaResult.Elv[i]);
      m_fittingManager->SetVelocityAbsorption(extremaResult.Alv[i]);
      Log.LogInfo("    Operator-Linemodel: recompute with elv=%.1f, alv=%.1f",
                  m_fittingManager->GetVelocityEmission(),
                  m_fittingManager->GetVelocityAbsorption());
    }

    if (m_opt_continuumcomponent == "tplfit" ||
        m_opt_continuumcomponent == "tplfitauto") {
      // fix some fitcontinuum values for this extremum
      if (tplfit_option == 2) {
        m_fittingManager->m_continuumManager->SetFitContinuum_FitStore(nullptr);
        m_fittingManager->m_continuumManager->SetFitContinuum_FitValues(
            extremaResult.FittedTplName[i], extremaResult.FittedTplAmplitude[i],
            extremaResult.FittedTplAmplitudeError[i],
            extremaResult.FittedTplMerit[i],
            extremaResult.FittedTplMeritPhot[i],
            extremaResult.FittedTplEbmvCoeff[i],
            extremaResult.FittedTplMeiksinIdx[i],
            extremaResult.FittedTplRedshift[i], extremaResult.FittedTplDtm[i],
            extremaResult.FittedTplMtm[i], extremaResult.FittedTplLogPrior[i],
            extremaResult.FittedTplpCoeffs[i]);
        m_fittingManager->m_continuumManager->SetFitContinuum_Option(
            tplfit_option);
      } else if (tplfit_option == 0 ||
                 tplfit_option == 3) // for these cases we called precompute in
                                     // secondpass, so we have new fitstore
        m_fittingManager->m_continuumManager->SetFitContinuum_FitStore(
            m_tplfitStore_secondpass[i]);
      else if (tplfit_option == 1) {
        // nothing to do cause we already injected the fitStore for cases 1
        // and 2
        m_fittingManager->m_continuumManager->SetFitContinuum_FitStore(
            m_tplfitStore_firstpass); // 1
      }
    }

    // moved here to override the previously set option value
    // since all
    // m_fittingManager->SetFitContinuum_Option(tplfit_option);
    Log.LogInfo("    Operator-Linemodel: recompute with tplfit_option=%d",
                tplfit_option);

    // find the index in the zaxis results
    const Int32 idx = CIndexing<Float64>::getIndex(m_result->Redshifts, Z);

    // reestimate the model (eventually with continuum reestimation) on
    // the extrema selected
    Int32 contreest_iterations = 0;
    if (opt_continuumreest == "always") {
      contreest_iterations = 1;
    }

    // finally compute the redshifts on the z-range around the extremum
    // m_fittingManager->SetFittingMethod("nofit");
    Int32 n_progresssteps = extremaResult.ExtendedRedshifts[i].size();
    Log.LogInfo("    Operator-Linemodel: Fit n=%d values for z in [%.6f; %.6f]",
                n_progresssteps, extremaResult.ExtendedRedshifts[i].front(),
                extremaResult.ExtendedRedshifts[i].back());
    m_fittingManager->logParameters();
    for (const Float64 z : extremaResult.ExtendedRedshifts[i]) {
      const Int32 iz = CIndexing<Float64>::getIndex(m_result->Redshifts, z);
      Log.LogDetail("Fit for Extended redshift %d, z = %f", iz, z);

      m_result->ChiSquare[iz] = m_fittingManager->fit(
          m_result->Redshifts[iz], m_result->LineModelSolutions[iz],
          m_result->ContinuumModelSolutions[iz], contreest_iterations, false);
      m_result->ScaleMargCorrection[iz] =
          m_fittingManager->getScaleMargCorrection();
      if (m_opt_continuumcomponent == "tplfit" ||
          m_opt_continuumcomponent == "tplfitauto") {
        if (tplfit_option == 0 ||
            tplfit_option == 3) // retryall & refitfirstpass
          m_result->SetChisquareTplContinuumResult(
              iz,
              m_fittingManager->m_continuumManager->GetFitContinuum_FitStore());
        // nothing to do when fromfirstpass: keep
        // m_result->ChiSquareTplContinuum from first pass
      }
      if (m_fittingManager->getLineRatioType() == "tplratio")
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
          m_fittingManager->m_continuumManager
              ->getContinuumScaleMargCorrection();
    }
  }

  return 0;
}

Int32 COperatorLineModel::Init(const TFloat64List &redshifts) {

  m_tplCategoryList = {Context.GetCurrentCategory()};
  // initialize empty results so that it can be returned anyway in case of an
  // error
  m_result = std::make_shared<CLineModelResult>();

  if (Context.GetSpectrum()->GetSpectralAxis().IsInLinearScale() == false) {
    THROWG(INTERNAL_ERROR, "input spectrum is not in linear scale.");
  }
  std::shared_ptr<const CParameterStore> ps = Context.GetParameterStore();

  CAutoScope autoscope(Context.m_ScopeStack, "linemodel");

  m_opt_continuumcomponent = ps->GetScoped<std::string>("continuumcomponent");

  // sort the redshifts
  m_sortedRedshifts = redshifts;
  std::sort(m_sortedRedshifts.begin(), m_sortedRedshifts.end());

  if (Context.GetCurrentMethod() == "LineModelSolve") {

    m_secondPass_halfwindowsize =
        ps->GetScoped<Float64>("secondpass.halfwindowsize");

    m_opt_firstpass_fittingmethod =
        ps->GetScoped<std::string>("firstpass.fittingmethod");
  }
  //
  if (m_opt_continuumcomponent == "tplfit" ||
      m_opt_continuumcomponent == "tplfitauto") {
    m_opt_fitcontinuum_maxN = ps->GetScoped<Int32>("continuumfit.count");
    Log.LogDetail(
        "  method Linemodel wit tplfit: fitcontinuum_maxN set to %.0f",
        m_opt_fitcontinuum_maxN);

    m_opt_tplfit_fftprocessing =
        ps->GetScoped<bool>("continuumfit.fftprocessing");
    m_opt_tplfit_fftprocessing_secondpass =
        m_opt_tplfit_fftprocessing; // TODO add a real parameter or remove this
                                    // member
    if (ps->HasScoped<bool>("enablephotometry"))
      m_opt_tplfit_use_photometry = ps->GetScoped<bool>("enablephotometry");
    m_opt_tplfit_dustFit = ps->GetScoped<bool>("continuumfit.ismfit");
    m_opt_tplfit_extinction = ps->GetScoped<bool>("continuumfit.igmfit");

    m_opt_tplfit_ignoreLinesSupport =
        ps->GetScoped<bool>("continuumfit.ignorelinesupport");
    if (Context.GetCurrentMethod() == "LineModelSolve") {

      m_opt_secondpasslcfittingmethod =
          ps->GetScoped<std::string>("secondpasslcfittingmethod");
    }
    m_opt_continuum_neg_amp_threshold =
        ps->GetScoped<Float64>("continuumfit.negativethreshold");

    m_opt_continuum_null_amp_threshold =
        ps->GetScoped<Float64>("continuumfit.nullthreshold");
  }

  return 0;
}

std::shared_ptr<COperatorResult> COperatorLineModel::getResult() {
  return m_result;
}

/**
 * @brief COperatorLineModel::initContaminant
 * prepare the contaminant to be used when instantiating a multimodel in
 * compute()
 * @return
 */
Int32 COperatorLineModel::initContaminant(
    std::shared_ptr<CModelSpectrumResult> contModelSpectrum,
    Int32 iRollContaminated, Float64 contLambdaOffset) {
  /*
    Log.LogInfo("  Operator-Linemodel: Initializing contaminant for roll #%d,
    " "with offset=%.2f", iRollContaminated, contLambdaOffset);

    m_iRollContaminated = iRollContaminated;
    m_contLambdaOffset = contLambdaOffset;
    const std::string &category = "emission";
    m_tplContaminant =
        std::shared_ptr<CTemplate>(new CTemplate("contaminant", category));
    Int32 length = contModelSpectrum->ModelFlux.size();

    CSpectrumFluxAxis  spcFluxAxis =
    contModelSpectrum->GetSpectrum().GetFluxAxis(); CSpectrumSpectralAxis
    spcSpectralAxis = contModelSpectrum->GetSpectrum().GetSpectralAxis();

    // applying offset for ra/dec distance between main source and contaminant
    spcSpectralAxis.ApplyOffset(m_contLambdaOffset);

    m_tplContaminant->SetSpectralAndFluxAxes(std::move(spcSpectralAxis),
    std::move(spcFluxAxis)); m_enableLoadContTemplate = true;
  */
  return 0;
}

std::shared_ptr<CModelSpectrumResult>
COperatorLineModel::GetContaminantSpectrumResult() {
  return m_savedContaminantSpectrumResult;
}

Int32 COperatorLineModel::interpolateLargeGridOnFineGrid(
    const TFloat64List &redshiftsLargeGrid,
    const TFloat64List &redshiftsFineGrid, const TFloat64List &meritLargeGrid,
    TFloat64List &meritFineGrid) const {
  //* // GSL method LIN
  Log.LogDetail("  Operator-Linemodel: First-Pass - interp FROM large grid "
                "z0=%f to zEnd=%f (n=%d)",
                redshiftsLargeGrid[0],
                redshiftsLargeGrid[redshiftsLargeGrid.size() - 1],
                redshiftsLargeGrid.size());
  Log.LogDetail("  Operator-Linemodel: First-Pass - interp TO fine grid "
                "z0=%f to zEnd=%f (n=%d)",
                redshiftsFineGrid[0],
                redshiftsFineGrid[redshiftsFineGrid.size() - 1],
                redshiftsFineGrid.size());

  // initialise and allocate the gsl objects
  // lin
  gsl_interp *interpolation =
      gsl_interp_alloc(gsl_interp_linear, meritLargeGrid.size());
  gsl_interp_init(interpolation, &(redshiftsLargeGrid.front()),
                  &(meritLargeGrid.front()), meritLargeGrid.size());
  gsl_interp_accel *accelerator = gsl_interp_accel_alloc();

  for (Int32 j = 0; j < redshiftsFineGrid.size(); j++) {
    Float64 Xrebin = redshiftsFineGrid[j];
    if (Xrebin < redshiftsLargeGrid[0] ||
        Xrebin > redshiftsLargeGrid[redshiftsLargeGrid.size() - 1]) {
      continue;
    }
    meritFineGrid[j] =
        gsl_interp_eval(interpolation, &redshiftsLargeGrid.front(),
                        &meritLargeGrid.front(), Xrebin, accelerator); // lin
  }

  gsl_interp_free(interpolation);
  gsl_interp_accel_free(accelerator);
  //*/

  return 0;
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
  Float64 A = boost::numeric::bounds<float>::lowest();
  const Float64 *flux = fluxAxis.GetSamples();
  const Float64 *spectral = spectralAxis.GetSamples();
  // const Float64* error = fluxAxis.GetError();

  // A = max, good value ?
  for (Int32 i = start; i < end; i++) {
    Float64 y = flux[i];
    if (y > A) {
      A = y;
    }
  }

  if (A <= 0) {
    return 0.0;
  }
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
      Float64 x = spectral[i];
      Float64 Yi = A * exp(-1. * (x - mu) * (x - mu) / (2 * c * c));
      sum2 += pow(Yi - flux[i], 2.0);
      // sum2 += pow( Yi - flux[i] , 2.0 ) / pow( error[i], 2.0 );
    }
    if (sum2 < minsum2) {
      minc = c;
      minsum2 = sum2;
    }
    icmpt++;
    c = c + cstepup;
  }

  if (minc < 0) {
    minc = 0;
  }
  return minc;
}

CLineModelSolution COperatorLineModel::fitWidthByGroups(
    std::shared_ptr<const CInputContext> context, Float64 redshift) {
  /*  CDataStore &datastore = context.GetDataStore();
  const TFloat64Range &lambdaRange = context.GetLambdaRange();
  Float64 redshift_min =
  datastore.GetFloat64Range("redshiftrange").GetBegin(); Float64 redshift_max
  = datastore.GetFloat64Range("redshiftrange").GetEnd(); CLineModelSolution
  modelSolution; CLineModelSolution bestModelSolution;

  CContinuumModelSolution continuumModelSolution;

  m_fittingManager->SetFittingMethod("hybrid");

    if (opt_lineRatioType == "tplratio")
    {
    m_fittingManager->SetFittingMethod("individual");
    }

  m_fittingManager->SetForcedisableTplratioISMfit(true);//m_fittingManager->m_opt_firstpass_forcedisableTplratioISMfit);
  //todo, add new param for this ?

  //TODO these params should be in a dedicated scope "velocityfit"
  datastore.PushScope("linemodel");
  const Float64&velfitMinE =
  datastore.GetScopedFloat64Param("emvelocityfitmin"); const
  Float64&velfitMaxE = datastore.GetScopedFloat64Param("emvelocityfitmax");
  const Float64&velfitStepE =
  datastore.GetScopedFloat64Param("emvelocityfitstep"); const
  Float64&velfitMinA = datastore.GetScopedFloat64Param("absvelocityfitmin");
  const Float64&velfitMaxA =
  datastore.GetScopedFloat64Param("absvelocityfitmax"); const
  Float64&velfitStepA = datastore.GetScopedFloat64Param("absvelocityfitstep");
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
  TFloat64Range(velfitMinE,velfitMaxE).SpreadOver(velfitStepE); TFloat64List
  velFitAList = TFloat64Range(velfitMinA,velfitMaxA).SpreadOver(velfitStepA);

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
  TFloat64Range(redshift-opt_manvelfit_dzmin,redshift+opt_manvelfit_dzmax).SpreadOver(opt_manvelfit_dzstep);

  fitVelocityByGroups(velFitEList,zList,CLine::nType_Emission);
  fitVelocityByGroups(velFitAList,zList,CLine::nType_Absorption);

*/
  CLineModelSolution clms;
  return clms;
}

void COperatorLineModel::fitVelocityByGroups(TFloat64List velfitlist,
                                             TFloat64List zfitlist,
                                             Int32 lineType) {
  Int32 nVelSteps = velfitlist.size();

  std::vector<TInt32List> idxVelfitGroups =
      m_fittingManager->m_Elements.GetModelVelfitGroups(lineType);
  TFloat64List GroupsV;

  for (Int32 kgroup = 0; kgroup < idxVelfitGroups.size(); kgroup++) {
    Log.LogDetail("  Operator-Linemodel: manualStep fitting group=%d", kgroup);

    Float64 meritMin = DBL_MAX;
    Float64 vOptim = -1.0;
    Float64 z_vOptim = -1.0;
    for (Int32 kdz = 0; kdz < zfitlist.size(); kdz++) {
      Float64 dzTest = zfitlist[kdz];
      for (Int32 kv = 0; kv < nVelSteps; kv++) {
        Float64 vTest = velfitlist[kv];
        for (Int32 ke = 0; ke < idxVelfitGroups[kgroup].size(); ke++) {
          m_fittingManager->setVelocity(vTest, idxVelfitGroups[kgroup][ke],
                                        lineType);
        }
      }
    }
  }
}

CLineModelSolution COperatorLineModel::computeForLineMeas(
    std::shared_ptr<const CInputContext> inputContext,
    const TFloat64List &redshiftsGrid, Float64 &bestz) {
  std::shared_ptr<const CParameterStore> params =
      inputContext->GetParameterStore();
  if (params->GetScoped<bool>("linemodel.velocityfit"))
    THROWG(INTERNAL_ERROR, "velocityfit not implemented yet");

  Int32 amplitudeOffsetsDegree =
      params->GetScoped<Int32>("linemodel.polynomialdegree");
  if (amplitudeOffsetsDegree < 0 || amplitudeOffsetsDegree > 2)
    THROWG(INTERNAL_ERROR, "the polynomial degree "
                           "parameter should be between 0 and 2");

  m_fittingManager = std::make_shared<CLineModelFitting>();

  m_fittingManager->setPassMode(
      3); // does m_fittingManager->m_enableAmplitudeOffsets = true;

  // init catalog offsets

  m_estimateLeastSquareFast = 0;
  m_fittingManager->m_lineRatioManager->SetLeastSquareFastEstimationEnabled(
      m_estimateLeastSquareFast);

  m_fittingManager->initDtd();

  CLineModelSolution modelSolution(m_fittingManager->m_RestLineList);
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
      bestModelSolution = modelSolution;
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
  m_fittingManager->getSpectrumModel()->refreshModel();
  return m_fittingManager->getSpectrumModel()->GetModelSpectrum();
}
