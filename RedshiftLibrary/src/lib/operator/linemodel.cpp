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
#include "RedshiftLibrary/linemodel/templatesfitstore.h"
#include "RedshiftLibrary/linemodel/templatesortho.h"
#include "RedshiftLibrary/operator/spectraFluxResult.h"
#include "RedshiftLibrary/operator/templatefitting.h"
#include "RedshiftLibrary/operator/templatefittinglog.h"
#include "RedshiftLibrary/operator/templatefittingresult.h"
#include "RedshiftLibrary/operator/templatefittingwithphot.h"
#include "RedshiftLibrary/spectrum/axis.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "RedshiftLibrary/spectrum/tools.h"
#include "RedshiftLibrary/statistics/deltaz.h"
#include "RedshiftLibrary/statistics/priorhelper.h"

#include "RedshiftLibrary/common/flag.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/common/indexing.h"
#include "RedshiftLibrary/log/log.h"

#include <boost/chrono/thread_clock.hpp>
#include <boost/format.hpp>

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
Int32 COperatorLineModel::ComputeFirstPass(
    const CSpectrum &spectrum,
    const CSpectrum &logSampledSpectrum, // this is temporary
    const CTemplateCatalog &tplCatalog,
    const CLineCatalogsTplShape &tplRatioCatalog,
    const TFloat64Range &lambdaRange,
    const std::shared_ptr<const CPhotBandCatalog> &photBandCat,
    const Float64 photo_weight, const std::string &opt_fittingmethod,
    const std::string &opt_lineWidthType, const Float64 opt_velocityEmission,
    const Float64 opt_velocityAbsorption, const std::string &opt_continuumreest,
    const std::string &opt_rules, const bool &opt_velocityFitting,
    const Int32 &opt_twosteplargegridstep_ratio,
    const string &opt_twosteplargegridsampling, const std::string &opt_rigidity,
    const Float64 opt_haprior) {
  TFloat64List largeGridRedshifts;
  // redefine redshift grid
  m_enableFastFitLargeGrid = 0;
  if (opt_twosteplargegridstep_ratio > 1) {
    m_enableFastFitLargeGrid = 1;
    CreateRedshiftLargeGrid(opt_twosteplargegridstep_ratio, largeGridRedshifts);
  } else {
    largeGridRedshifts = m_sortedRedshifts;
  }

  TFloat64Range clampedlambdaRange;
  spectrum.GetSpectralAxis().ClampLambdaRange(lambdaRange, clampedlambdaRange);
  m_model = std::make_shared<CLineModelFitting>(
      spectrum, clampedlambdaRange, tplCatalog, m_tplCategoryList,
      m_RestLineList, opt_fittingmethod, m_opt_continuumcomponent,
      m_opt_continuum_neg_amp_threshold, opt_lineWidthType,
      m_linesmodel_nsigmasupport, opt_velocityEmission, opt_velocityAbsorption,
      opt_rules, opt_rigidity);
  m_model->setHaPriorOption(opt_haprior);

  /*
  CMultiRollModel model( spectrum,
                         tplCatalog,//orthoTplCatalog,
                         m_tplCategoryList,
                         restLineList,
                         opt_fittingmethod,
                         opt_continuumcomponent,
                         opt_lineWidthType,
                         opt_velocityEmission,
                         opt_velocityAbsorption,
                         opt_rules,
                         opt_rigidity);

  bool enableLoadContTemplateOverride = false; //manual switch for hardcoded
  contaminant bypass load if(enableLoadContTemplateOverride)
  {
      m_enableLoadContTemplate = false; //disable stock contaminant in case of
  overrided contaminant Int32 iContTemplate = 0; //idx of the roll being
  contaminated std::shared_ptr<CTemplate> tplContaminant;

      Log.LogInfo( "  Operator-Linemodel: OVERRIDDEN loading contaminant for
  roll #%d", iContTemplate );

      //hardcoded load from file on disk: of the contaminant for first model
      std::string templatePath =
  "/home/aschmitt/data/euclid/simulation2017-SC3_test_zweiroll/amazed/output_rolls_source3/euc_testsc3_zweiroll_source3_roll0_F_i0/linemodelsolve.linemodel_spc_extrema_0.txt";
      const std::string& category = "emission";
      tplContaminant = std::shared_ptr<CTemplate>( new CTemplate(
  "contaminant", category ) ); CSpectrumIOGenericReader asciiReader; if(
  !asciiReader.Read( templatePath.c_str(), *tplContaminant ) ) { Log.LogError(
  "Fail to read contaminant template: %s", templatePath.c_str() ); return -1;
      }else{
          Log.LogInfo( "Successfully loaded contaminant template: %s",
  templatePath.c_str() );
      }
      //debug:
      //FILE* f = fopen( "contaminantLoaded.txt", "w+" );
      //for(Int32 k=0; k<tplContaminant->GetSampleCount(); k++)
      //{
      //    fprintf( f, "%f\t%e\n", tplContaminant->GetSpectralAxis()[k],
  tplContaminant->GetFluxAxis()[k]);
      //}
      //fclose( f );
      //

      Float64 lambdaOffset_forSourceSpatialOffsetInDispersionDirection=2000;
  //hardcoded
      tplContaminant->GetSpectralAxis().ApplyOffset(lambdaOffset_forSourceSpatialOffsetInDispersionDirection);
      //debug:
      FILE* f2 = fopen( "contaminantShifted.txt", "w+" );
      for(Int32 k=0; k<tplContaminant->GetSampleCount(); k++)
      {
          fprintf( f2, "%f\t%e\n", tplContaminant->GetSpectralAxis()[k],
  tplContaminant->GetFluxAxis()[k]);
      }
      fclose( f2 );

      //tplContaminant
      model.LoadFitContaminantTemplate(iContTemplate, *tplContaminant,
  lambdaRange);
  }
  if(m_enableLoadContTemplate)
  {
      Log.LogInfo( "  Operator-Linemodel: loading contaminant for roll #%d",
  m_iRollContaminated ); if(m_tplContaminant==0)
      {
          Log.LogError( "  Operator-Linemodel: Contaminant data is invalid...
  aborting." ); }else{
          //
          if(0)
          {
          //debug:
          FILE* f2 = fopen( "contaminantShifted.txt", "w+" );
          for(Int32 k=0; k<m_tplContaminant->GetSampleCount(); k++)
          {
              fprintf( f2, "%f\t%e\n", m_tplContaminant->GetSpectralAxis()[k],
  m_tplContaminant->GetFluxAxis()[k]);
          }
          fclose( f2 );
          //
          }

          //apply contamination to the multiroll model
          model.LoadFitContaminantTemplate(m_iRollContaminated,
  *m_tplContaminant, lambdaRange); m_savedContaminantSpectrumResult =
  model.GetContaminantSpectrumResult(m_iRollContaminated);
      }
  }
  //*/

  // set some model parameters
  m_model->m_opt_firstpass_fittingmethod = m_opt_firstpass_fittingmethod;
  m_model->m_opt_secondpass_fittingmethod = opt_fittingmethod;

  Int32 opt_tplfit_integer_chi2_ebmv = -1;
  if (m_opt_tplfit_dustFit) {
    opt_tplfit_integer_chi2_ebmv = -10;
  }
  // should be replaced with option passed in param.json?
  Int32 observedFrame = 0;
  // passing the ignorelinesSupport option to the secondpass; //before it was
  // hardcoded to 0
  m_model->SetSecondpassContinuumFitPrms(
      opt_tplfit_integer_chi2_ebmv, m_opt_tplfit_extinction,
      m_opt_tplfit_ignoreLinesSupport, observedFrame);

  m_model->m_opt_lya_forcefit = m_opt_lya_forcefit;
  m_model->m_opt_lya_forcedisablefit = m_opt_lya_forcedisablefit;
  m_model->m_opt_lya_fit_asym_min = m_opt_lya_fit_asym_min;
  m_model->m_opt_lya_fit_asym_max = m_opt_lya_fit_asym_max;
  m_model->m_opt_lya_fit_asym_step = m_opt_lya_fit_asym_step;
  m_model->m_opt_lya_fit_width_min = m_opt_lya_fit_width_min;
  m_model->m_opt_lya_fit_width_max = m_opt_lya_fit_width_max;
  m_model->m_opt_lya_fit_width_step = m_opt_lya_fit_width_step;
  m_model->m_opt_lya_fit_delta_min = m_opt_lya_fit_delta_min;
  m_model->m_opt_lya_fit_delta_max = m_opt_lya_fit_delta_max;
  m_model->m_opt_lya_fit_delta_step = m_opt_lya_fit_delta_step;

  m_model->m_opt_enable_improveBalmerFit = m_opt_enableImproveBalmerFit;

  if (opt_rigidity == "tplshape") {
    // init catalog tplratios
    Log.LogInfo("  Operator-Linemodel: Tpl-ratios init");
    m_model->m_CatalogTplShape = tplRatioCatalog; // pass tplRatioCatalog
    bool tplratioInitRet = m_model->initTplratioCatalogs(m_opt_tplratio_ismFit);
    if (!tplratioInitRet) {
      THROWG(INTERNAL_ERROR,
             "  Operator-Linemodel: Failed to init tpl-ratios. aborting...");
    }

    m_model->m_opt_firstpass_forcedisableTplratioISMfit =
        !m_opt_firstpass_tplratio_ismFit;

    InitTplratioPriors();
  }

  // TODO: check option tplfit
  Int32 nfitcontinuum = 0;
  if (m_opt_continuumcomponent == "tplfit" ||
      m_opt_continuumcomponent == "tplfitauto")
    for (const auto &category : m_tplCategoryList)
      nfitcontinuum += tplCatalog.GetTemplateCount(category);
  m_result->Init(m_sortedRedshifts, m_RestLineList, nfitcontinuum,
                 m_model->getTplshape_count(), m_model->getTplshape_priors());

  Log.LogInfo("  Operator-Linemodel: initialized");

  // commom between firstpass and secondpass processes
  m_phelperContinuum = std::make_shared<CPriorHelper>();
  m_phelperContinuum->Init(m_opt_tplfit_continuumprior_dirpath.c_str(), 0);
  m_phelperContinuum->SetBetaA(m_opt_tplfit_continuumprior_betaA);
  m_phelperContinuum->SetBetaTE(m_opt_tplfit_continuumprior_betaTE);
  m_phelperContinuum->SetBetaZ(m_opt_tplfit_continuumprior_betaZ);
  m_model->SetFitContinuum_PriorHelper(m_phelperContinuum);

  Log.LogInfo("  Operator-Linemodel: start processing");

  // Set model parameters to FIRST-PASS
  m_model->setPassMode(1);
  Log.LogInfo(
      "  Operator-Linemodel: ---------- ---------- ---------- ----------");
  Log.LogInfo("  Operator-Linemodel: now computing first-pass");
  Log.LogInfo(
      "  Operator-Linemodel: ---------- ---------- ---------- ----------");

  // fit continuum
  ////////////////////
  if (m_opt_continuumcomponent == "tplfit" ||
      m_opt_continuumcomponent == "tplfitauto") {
    tplCatalog.m_orthogonal = 1;
    m_tplfitStore_firstpass =
        PrecomputeContinuumFit(spectrum, logSampledSpectrum, tplCatalog,
                               lambdaRange, largeGridRedshifts, photBandCat,
                               photo_weight, m_opt_tplfit_ignoreLinesSupport);
    tplCatalog.m_orthogonal = 0;
  } else {
    if (m_opt_continuumcomponent == "tplfit" ||
        m_opt_continuumcomponent == "tplfitauto") {
      m_model->m_opt_fitcontinuum_maxCount = m_opt_fitcontinuum_maxN;
    }
  }
  if (m_opt_continuumcomponent == "tplfit" ||
      m_opt_continuumcomponent == "tplfitauto") {
    m_model->m_opt_firstpass_forcedisableMultipleContinuumfit =
        m_opt_firstpass_multiplecontinuumfit_disable;
  }

  //    //hack, zero outside of the support
  //    ///////////////////////////////////////////////////////////////////////////////////////////
  //    model.setModelSpcObservedOnSupportZeroOutside(lambdaRange);

  //    std::shared_ptr<CSpectraFluxResult> baselineResult =
  //    (std::shared_ptr<CSpectraFluxResult>) new CSpectraFluxResult();
  //    baselineResult->m_optio = 0;
  //    const CSpectrum& modelSpc = model.GetModelSpectrum();
  //    Int32 len = modelSpc.GetSampleCount();

  //    baselineResult->fluxes.resize(len);
  //    baselineResult->wavel.resize(len);
  //    for( Int32 k=0; k<len; k++ )
  //    {
  //        baselineResult->fluxes[k] = modelSpc.GetFluxAxis()[k];
  //        baselineResult->wavel[k]  = (spectrum.GetSpectralAxis())[k];
  //    }

  //    std::string nameBaselineStr =
  //    (boost::format("linemodel_template_zeroedoutsidelines")).str();
  //    dataStore.StoreScopedGlobalResult(nameBaselineStr.c_str(),
  //    baselineResult); return NULL;
  //    // end of hack
  //    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  m_result->nSpcSamples = m_model->getSpcNSamples();
  m_result->dTransposeD = m_model->getDTransposeD();
  m_result->cstLog = m_model->getLikelihood_cstLog();

  Int32 contreest_iterations = 0;
  if (opt_continuumreest == "always") {
    contreest_iterations = 1;
  }

  // Set model parameter: abs lines limit
  Float64 absLinesLimit = 1.0; //-1 to disable, 1.0 is typical
  m_model->SetAbsLinesLimit(absLinesLimit);
  Log.LogInfo("  Operator-Linemodel: set abs lines limit to %f (ex: -1 means "
              "disabled)",
              absLinesLimit);

  // Set model parameter: continuum least-square estimation fast
  // note: this fast method requires continuum templates and linemodels to be
  // orthogonal. The velfit option turns this trickier...
  m_estimateLeastSquareFast = 0;
  m_model->SetLeastSquareFastEstimationEnabled(m_estimateLeastSquareFast);
  Log.LogInfo("  Operator-Linemodel: set estimateLeastSquareFast to %d (ex: "
              "0 means disabled)",
              m_estimateLeastSquareFast);

  //
  TBoolList allAmplitudesZero;
  Int32 indexLargeGrid = 0;
  TFloat64List calculatedLargeGridRedshifts;
  TFloat64List calculatedLargeGridMerits;
  std::vector<TFloat64List> calculatedChiSquareTplshapes(
      m_result->ChiSquareTplshapes.size());
  std::vector<TFloat64List> calculatedChisquareTplContinuum(
      m_result->ChiSquareTplContinuum.size());

  boost::chrono::thread_clock::time_point start_mainloop =
      boost::chrono::thread_clock::now();

  //#pragma omp parallel for
  for (Int32 i = 0; i < m_result->Redshifts.size(); i++) {
    if (m_enableFastFitLargeGrid == 0 ||
        m_result->Redshifts[i] == largeGridRedshifts[indexLargeGrid]) {
      m_result->ChiSquare[i] = m_model->fit(
          m_result->Redshifts[i], m_result->LineModelSolutions[i],
          m_result->ContinuumModelSolutions[i], contreest_iterations, false);
      calculatedLargeGridRedshifts.push_back(m_result->Redshifts[i]);
      calculatedLargeGridMerits.push_back(m_result->ChiSquare[i]);
      m_result->ScaleMargCorrection[i] = m_model->getScaleMargCorrection();
      if (m_opt_continuumcomponent == "tplfit" ||
          m_opt_continuumcomponent == "tplfitauto") {
        m_result->SetChisquareTplContinuumResult(i, m_tplfitStore_firstpass);
        for (Int32 k = 0, s = m_result->ChiSquareTplContinuum.size(); k < s;
             k++) {
          calculatedChisquareTplContinuum[k].push_back(
              m_result->ChiSquareTplContinuum[k][i]);
        }
      }
      m_result->SetChisquareTplshapeResult(
          i, m_model->GetChisquareTplshape(), m_model->GetScaleMargTplshape(),
          m_model->GetStrongELPresentTplshape(),
          m_model->getHaELPresentTplshape(),
          m_model->GetNLinesAboveSNRTplshape(),
          m_model->GetPriorLinesTplshape());

      for (Int32 k = 0, s = m_result->ChiSquareTplshapes.size(); k < s; k++) {
        calculatedChiSquareTplshapes[k].push_back(
            m_result->ChiSquareTplshapes[k][i]);
      }
      if (!m_estimateLeastSquareFast) {
        m_result->ChiSquareContinuum[i] =
            m_model->getLeastSquareContinuumMerit();
      } else {
        m_result->ChiSquareContinuum[i] =
            m_model->getLeastSquareContinuumMeritFast();
      }
      m_result->ScaleMargCorrectionContinuum[i] =
          m_model->getContinuumScaleMargCorrection();
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
      m_result->SetChisquareTplshapeResult(
          i, m_result->getChisquareTplshapeResult(i - 1),
          m_result->getScaleMargCorrTplshapeResult(i - 1),
          m_result->getStrongELPresentTplshapeResult(i - 1),
          m_result->getHaELPresentTplshapeResult(i - 1),
          m_result->getNLinesAboveSNRTplshapeResult(i - 1),
          m_result->getPriorLinesTplshapeResult(i - 1));
      if (!m_estimateLeastSquareFast) {
        m_result->ChiSquareContinuum[i] =
            m_model->getLeastSquareContinuumMerit();
      } else {
        m_result->ChiSquareContinuum[i] =
            m_model->getLeastSquareContinuumMeritFast();
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
    THROWG(INTERNAL_ERROR, "  Operator-Linemodel: All amplitudes (continuum & "
                           "model) are zero for all z. Aborting...");
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

  for (Int32 kts = 0, s = m_result->ChiSquareTplshapes.size(); kts < s; kts++) {
    if (m_result->Redshifts.size() > calculatedChiSquareTplshapes[kts].size() &&
        calculatedChiSquareTplshapes[kts].size() > 1) {
      interpolateLargeGridOnFineGrid(
          calculatedLargeGridRedshifts, m_result->Redshifts,
          calculatedChiSquareTplshapes[kts], m_result->ChiSquareTplshapes[kts]);
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

/**
 * Estimate once for all the continuum amplitude which is only dependent from
 * the tplName, ism and igm indexes. This is useful when the model fitting
 * option corresponds to fitting separately the continuum and the lines. In such
 * case, playing with (fit) Lines parameters (velocity, line offsets, lines
 * amps, etc.) do not affect continuum amplitudes.. thus we can save time  by
 * fitting once-for-all the continuum amplitudes, prior to fitting the lines.
 * @candidateIdx@ is also an indicator of pass mode
 * */
std::shared_ptr<CTemplatesFitStore> COperatorLineModel::PrecomputeContinuumFit(
    const CSpectrum &spectrum, const CSpectrum &logSampledSpectrum,
    const CTemplateCatalog &tplCatalog, const TFloat64Range &lambdaRange,
    const TFloat64List &redshifts,
    const std::shared_ptr<const CPhotBandCatalog> &photBandCat,
    const Float64 photometry_weight, bool ignoreLinesSupport,
    Int32 candidateIdx) {
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

  for (Int32 kztplfit = 0;
       kztplfit < std::min(Int32(redshiftsTplFit.size()), Int32(10));
       kztplfit++) {
    Log.LogDebug("COperatorLineModel::PrecomputeContinuumFit: continuum tpl "
                 "redshift list[%d] = %f",
                 kztplfit, redshiftsTplFit[kztplfit]);
  }

  std::vector<std::shared_ptr<CTemplateFittingResult>> chisquareResultsAllTpl;
  TStringList chisquareResultsTplName;

  bool fftprocessing = m_model->GetPassNumber() == 1
                           ? m_opt_tplfit_fftprocessing
                           : m_opt_tplfit_fftprocessing_secondpass;
  Log.LogDebug(Formatter()
               << "COperatorLineModel::PrecomputeContinuumFit: redshtplfitsize "
               << redshiftsTplFit.size());
  const Int32 fftprocessing_min_sample_nb =
      100; // warning arbitrary number of redshift samples threshold below which
           // fftprocessing is considered slower
           // TB revised
  if ((redshiftsTplFit.size() < fftprocessing_min_sample_nb) && fftprocessing) {
    fftprocessing = false;
    Log.LogInfo("COperatorLineModel::PrecomputeContinuumFit: auto deselect fft "
                "processing"
                " (faster when only few redshifts calc. points)");
  }

  bool currentSampling = tplCatalog.m_logsampling;
  std::string opt_interp = "precomputedfinegrid"; //"lin";
  Log.LogInfo("COperatorLineModel::PrecomputeContinuumFit: fftprocessing = %d",
              fftprocessing);
  Log.LogDetail(
      "COperatorLineModel::PrecomputeContinuumFit: fitContinuum_dustfit = %d",
      m_opt_tplfit_dustFit);
  Log.LogDetail(
      "COperatorLineModel::PrecomputeContinuumFit: fitContinuum_igm = %d",
      m_opt_tplfit_extinction);
  Log.LogDetail("COperatorLineModel::PrecomputeContinuumFit: fitContinuum "
                "opt_interp = %s",
                opt_interp.c_str());

  TFloat64Range clampedlambdaRange;
  if (fftprocessing) {
    if (m_templateFittingOperator == nullptr ||
        !m_templateFittingOperator->IsFFTProcessing()) // else reuse the shared
                                                       // pointer for secondpass
      m_templateFittingOperator = std::make_shared<COperatorTemplateFittingLog>(
          spectrum, logSampledSpectrum, lambdaRange, redshiftsTplFit);
    logSampledSpectrum.GetSpectralAxis().ClampLambdaRange(lambdaRange,
                                                          clampedlambdaRange);
    tplCatalog.m_logsampling = true;
  } else {
    if (m_templateFittingOperator == nullptr ||
        m_templateFittingOperator
            ->IsFFTProcessing()) // else reuse the shared pointer for secondpass
      if (m_opt_tplfit_use_photometry) {
        m_templateFittingOperator =
            std::make_shared<COperatorTemplateFittingPhot>(
                spectrum, lambdaRange, photBandCat, photometry_weight,
                redshiftsTplFit);
      } else {
        m_templateFittingOperator = std::make_shared<COperatorTemplateFitting>(
            spectrum, lambdaRange, redshiftsTplFit);
      }
    spectrum.GetSpectralAxis().ClampLambdaRange(lambdaRange,
                                                clampedlambdaRange);
    tplCatalog.m_logsampling = false;
  }

  Float64 overlapThreshold = 1.0;
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
    for (Int32 kztplfit = 0; kztplfit < redshiftsTplFit.size(); kztplfit++) {
      m_model->initModelAtZ(redshiftsTplFit[kztplfit],
                            fftprocessing ? logSampledSpectrum.GetSpectralAxis()
                                          : spectrum.GetSpectralAxis());
      maskList[kztplfit] = m_model->getOutsideLinesMask();
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

  Int32 opt_tplfit_integer_chi2_ebmv = -1;
  if (m_opt_tplfit_dustFit) {
    opt_tplfit_integer_chi2_ebmv = -10;
  }

  Float64 EbmvCoeff = -1.;
  Int32 meiksinIdx = -1;
  bool keepismigm = false;
  if (m_model->GetPassNumber() == 2) { // if we are in secondpass
    if (m_continnuum_fit_option == 3 &&
        (m_opt_tplfit_dustFit || m_opt_tplfit_extinction)) { // refitfirstpass
      if (candidateIdx < 0 ||
          candidateIdx > m_firstpass_extremaResult->size() - 1) {
        THROWG(INTERNAL_ERROR, "Candidate index is out of range");
      }
      // using one template per Z with fixed values for ism/igm (if Z changes,
      // template change;)
      keepismigm = true; // telling that we want to keep the ism and igm indexes
      meiksinIdx = m_firstpass_extremaResult->FittedTplMeiksinIdx[candidateIdx];
      EbmvCoeff = m_firstpass_extremaResult->FittedTplEbmvCoeff[candidateIdx];
      // access any template and retrieve the ismcorrection object
      opt_tplfit_integer_chi2_ebmv =
          tplCatalog.GetTemplate(m_tplCategoryList[0], 0)
              ->m_ismCorrectionCalzetti->GetEbmvIndex(EbmvCoeff);
    }
  }

  bool found = false;
  for (Int32 i = 0; i < m_tplCategoryList.size(); i++) {
    std::string category = m_tplCategoryList[i];
    Log.LogDebug(Formatter()
                 << "Processing " << tplCatalog.GetTemplateCount(category)
                 << " templates");

    for (Int32 j = 0; j < tplCatalog.GetTemplateCount(category); j++) {
      std::shared_ptr<const CTemplate> tpl =
          tplCatalog.GetTemplate(category, j);
      Log.LogDebug(Formatter() << "Processing tpl " << tpl->GetName());
      // case where we only want to refit using one template:
      if (m_continnuum_fit_option == 3) {
        if (tpl->GetName() !=
            m_firstpass_extremaResult->FittedTplName[candidateIdx])
          continue;
        else
          found = true;
      }

      CPriorHelper::TPriorZEList zePriorData;

      bool retGetPrior = m_phelperContinuum->GetTplPriorData(
          tpl->GetName(), redshiftsTplFit, zePriorData);
      if (retGetPrior == false) {
        THROWG(INTERNAL_ERROR, "Failed to get prior "
                               "for chi2 continuum precomp fit. aborting...");
      }

      m_templateFittingOperator->SetRedshifts(redshiftsTplFit);
      auto templatefittingResult =
          std::dynamic_pointer_cast<CTemplateFittingResult>(
              m_templateFittingOperator->Compute(
                  tpl, overlapThreshold, maskList, opt_interp,
                  m_opt_tplfit_extinction, opt_tplfit_integer_chi2_ebmv,
                  zePriorData, keepismigm, EbmvCoeff, meiksinIdx));

      if (!templatefittingResult) {
        THROWG(INTERNAL_ERROR, Formatter()
                                   << "failed "
                                      "to compute chisquare value for tpl=%"
                                   << tpl->GetName());
      } else {
        chisquareResultsAllTpl.push_back(templatefittingResult);
        chisquareResultsTplName.push_back(tpl->GetName());
      }
      if (found)
        break;
    }
    if (found)
      break;
  }

  // fill the fit store with fitted values: only the best fitted values FOR EACH
  // TEMPLATE are used
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
      // Log.LogInfo("  Operator-Linemodel: check prior data, tplfitStore->Add
      // logprior = %e", chisquareResult->LogPrior[i]);

      if (!retAdd) {
        THROWG(INTERNAL_ERROR, "Failed to add "
                               "continuum fit to store. aborting...");
      }

      Float64 tplfitsnr = -1.;
      if (chisquareResult->FitMtM[i] > 0.) {
        tplfitsnr =
            chisquareResult->FitDtM[i] / std::sqrt(chisquareResult->FitMtM[i]);
      }
      if (tplfitsnr > bestTplFitSNR) {
        bestTplFitSNR = tplfitsnr;
      }
    }
  }
  tplfitStore->m_fitContinuum_tplFitSNRMax = bestTplFitSNR;
  Log.LogDetail("COperatorLineModel::PrecomputeContinuumFit: "
                "fitcontinuum_snrMAX set to %f",
                bestTplFitSNR);
  Log.LogDetail(
      "COperatorLineModel::PrecomputeContinuumFit: continuumcount set to %d",
      tplfitStore->GetContinuumCount());

  Int32 v = std::min(m_opt_fitcontinuum_maxN, tplfitStore->GetContinuumCount());
  // TODO if v < a_value_to_compute log warning or throw error (if v = 1, we'll
  // have 1st pass redshift=nan)
  tplfitStore->m_opt_fitcontinuum_maxCount =
      (m_opt_fitcontinuum_maxN == -1 ? tplfitStore->GetContinuumCount() : v);
  Log.LogInfo("COperatorLineModel::PrecomputeContinuumFit: fitStore with "
              "fitcontinuum_maxCount set to %f",
              tplfitStore->m_opt_fitcontinuum_maxCount);

  m_model->SetFitContinuum_FitStore(tplfitStore);
  tplCatalog.m_logsampling = currentSampling;

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

  // Check if best continuum amplitudes are negative fitted amplitudes at all z
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
      m_model->SetContinuumComponent("fromspectrum");
    }
  }

  return tplfitStore;
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
    THROWG(INTERNAL_ERROR, "Second pass "
                           "window outside z range");
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

    m_result->LineModelSolutions[idx].fillLineIds();
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
  TInt32List uniqueIdx_fpb =
      m_firstpass_extremaResult->getUniqueCandidates(firstpass_results_b);
  m_firstpass_extremaResult->Resize(m_firstpass_extremaResult->size() +
                                    uniqueIdx_fpb.size());

  Int32 startIdx = m_firstpass_extremaResult->size();
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
    const CSpectrum &spectrum, const CSpectrum &logSampledSpectrum,
    const CTemplateCatalog &tplCatalog, const TFloat64Range &lambdaRange,
    const std::shared_ptr<const CPhotBandCatalog> &photBandCat,
    const std::shared_ptr<const LineModelExtremaResult> &firstpassResults,
    const Float64 photo_weight, const std::string &opt_fittingmethod,
    const std::string &opt_lineWidthType, const Float64 opt_velocityEmission,
    const Float64 opt_velocityAbsorption, const std::string &opt_continuumreest,
    const std::string &opt_rules, const bool &opt_velocityFitting,
    const std::string &opt_rigidity, const Float64 &opt_emvelocityfitmin,
    const Float64 &opt_emvelocityfitmax, const Float64 &opt_emvelocityfitstep,
    const Float64 &opt_absvelocityfitmin, const Float64 &opt_absvelocityfitmax,
    const Float64 &opt_absvelocityfitstep,
    const std::string &opt_continuumfit_method) {
  boost::chrono::thread_clock::time_point start_secondpass =
      boost::chrono::thread_clock::now();
  // Set model parameters to SECOND-PASS
  m_model->setPassMode(2);
  Int32 savedFitContinuumOption =
      m_model->GetFitContinuum_Option(); // the first time was set in
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
    THROWG(
        INTERNAL_ERROR,
        Formatter() << "Operator-Linemodel: continnuum_fit_option not found: "
                    << m_continnuum_fit_option);
  }
  if (m_opt_continuumcomponent == "tplfit" ||
      m_opt_continuumcomponent == "tplfitauto") {
    // precompute only whenever required and whenever the result can be a
    // tplfitStore
    if (m_continnuum_fit_option == 0 || m_continnuum_fit_option == 3) {
      m_tplfitStore_secondpass.resize(m_firstpass_extremaResult->size());
      tplCatalog.m_orthogonal = 1;
      for (Int32 i = 0; i < m_firstpass_extremaResult->size(); i++) {
        m_tplfitStore_secondpass[i] = PrecomputeContinuumFit(
            spectrum, logSampledSpectrum, tplCatalog, lambdaRange,
            m_firstpass_extremaResult->ExtendedRedshifts[i], photBandCat,
            photo_weight, m_opt_tplfit_ignoreLinesSupport, i);
        if (m_opt_continuumcomponent == "fromspectrum")
          break; // when set to "fromspectrum" by PrecomputeContinuumFit because
                 // negative continuum with tplfitauto
      }
      tplCatalog.m_orthogonal = 0; // finish using orthogTemplates
    } else {
      // since precompute is not called all the time, secondpass candidates do
      // not have systematically a tplfitstore_secondpass copy the firstpass
      // tplfitstore into the secondpass tplfitstore
      m_tplfitStore_secondpass.resize(1);
      if (m_continnuum_fit_option == 1 || m_continnuum_fit_option == 2)
        m_tplfitStore_secondpass[0] = m_tplfitStore_firstpass;
      // belwo line is redundant since in the first call to precompute we
      // already set the fitStore
      m_model->SetFitContinuum_FitStore(m_tplfitStore_firstpass);
      // duplicated: double make sure that these info are present in the
      // modelElement
      m_model->SetFitContinuum_SNRMax(
          m_tplfitStore_firstpass->m_fitContinuum_tplFitSNRMax);
      m_model->m_opt_fitcontinuum_maxCount =
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
  TFloat64Range clampedlambdaRange;
  spectrum.GetSpectralAxis().ClampLambdaRange(lambdaRange, clampedlambdaRange);

  // upcast LineModelExtremaResult to TCandidateZ
  m_secondpass_parameters_extremaResult.m_ranked_candidates.assign(
      firstpassResults->m_ranked_candidates.cbegin(),
      firstpassResults->m_ranked_candidates.cend());

  // now that we recomputed what should be recomputed, we define once for all
  // the secondpass
  //  estimate second pass parameters (mainly elv, alv...)
  EstimateSecondPassParameters(
      spectrum, clampedlambdaRange, opt_continuumreest, opt_fittingmethod,
      opt_rigidity, opt_velocityFitting, opt_emvelocityfitmin,
      opt_emvelocityfitmax, opt_emvelocityfitstep, opt_absvelocityfitmin,
      opt_absvelocityfitmax, opt_absvelocityfitstep);

  // recompute the fine grid results around the extrema
  Int32 ret = RecomputeAroundCandidates(
      clampedlambdaRange, opt_continuumreest,
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
    m_model->SetFittingMethod(m_opt_secondpasslcfittingmethod);
    RecomputeAroundCandidates(clampedlambdaRange, opt_continuumreest, 2, true);
    m_model->SetFittingMethod(opt_fittingmethod);

    Log.LogInfo("  Operator-Linemodel: now re-computing the final chi2 for "
                "each candidate");
    RecomputeAroundCandidates(clampedlambdaRange, opt_continuumreest, 2);
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

  m_model->SetFitContinuum_Option(savedFitContinuumOption);

  return 0;
}

std::shared_ptr<LineModelExtremaResult>
COperatorLineModel::buildExtremaResults(const CSpectrum &spectrum,
                                        const TFloat64Range &lambdaRange,
                                        const TCandidateZbyRank &zCandidates,
                                        const std::string &opt_continuumreest) {
  Int32 savedFitContinuumOption = m_model->GetFitContinuum_Option();
  Log.LogInfo("  Operator-Linemodel: Now storing extrema results");

  Int32 extremumCount = zCandidates.size();
  if (extremumCount > m_maxModelSaveCount) {
    THROWG(INTERNAL_ERROR, Formatter() << "ExtremumCount " << extremumCount
                                       << " is greater the maxModelSaveCount "
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

    m_model->SetFitContinuum_FitValues(
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

    m_model->SetFitContinuum_Option(2);

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
      idxVelfitGroups =
          m_model->m_Elements.GetModelVelfitGroups(CLine::nType_Absorption);
      std::string alv_list_str = "";
      for (Int32 kgroup = 0; kgroup < idxVelfitGroups.size(); kgroup++) {
        for (Int32 ke = 0; ke < idxVelfitGroups[kgroup].size(); ke++) {
          m_model->SetVelocityAbsorptionOneElement(
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
      idxVelfitGroups =
          m_model->m_Elements.GetModelVelfitGroups(CLine::nType_Emission);
      std::string elv_list_str = "";
      for (Int32 kgroup = 0; kgroup < idxVelfitGroups.size(); kgroup++) {
        for (Int32 ke = 0; ke < idxVelfitGroups[kgroup].size(); ke++) {
          m_model->SetVelocityEmissionOneElement(
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
      // m_model->SetVelocityEmission(m_secondpass_parameters_extremaResult.Elv[i_2pass]);
      // m_model->SetVelocityAbsorption(m_secondpass_parameters_extremaResult.Alv[i_2pass]);
      m_model->SetVelocityEmission(
          m_result->LineModelSolutions[idx].EmissionVelocity);
      m_model->SetVelocityAbsorption(
          m_result->LineModelSolutions[idx].AbsorptionVelocity);
    }

    if (!mlmfit_modelInfoSave) {
      m_result->ChiSquare[idx] = m_model->fit(
          m_result->Redshifts[idx], m_result->LineModelSolutions[idx],
          m_result->ContinuumModelSolutions[idx], contreest_iterations, true);
      m_result->ScaleMargCorrection[idx] = m_model->getScaleMargCorrection();
      m_result->SetChisquareTplshapeResult(
          idx, m_model->GetChisquareTplshape(), m_model->GetScaleMargTplshape(),
          m_model->GetStrongELPresentTplshape(),
          m_model->getHaELPresentTplshape(),
          m_model->GetNLinesAboveSNRTplshape(),
          m_model->GetPriorLinesTplshape());
      if (!m_estimateLeastSquareFast) {
        m_result->ChiSquareContinuum[idx] =
            m_model->getLeastSquareContinuumMerit();
      } else {
        m_result->ChiSquareContinuum[idx] =
            m_model->getLeastSquareContinuumMeritFast();
      }
      m_result->ScaleMargCorrectionContinuum[idx] =
          m_model->getContinuumScaleMargCorrection();
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
              m_model->GetModelSpectrum());
        } else if (overrideModelSavedType == 1 || overrideModelSavedType == 2) {
          Int32 lineTypeFilter = -1;
          if (overrideModelSavedType == 1) {
            lineTypeFilter = -1;
          } else if (overrideModelSavedType == 2) {
            lineTypeFilter = CLine::nType_Emission;
          }
          resultspcmodel = std::make_shared<CModelSpectrumResult>(
              m_model->GetObservedSpectrumWithLinesRemoved(lineTypeFilter));
        }
        // std::shared_ptr<CModelSpectrumResult>  resultspcmodel =
        // std::shared_ptr<CModelSpectrumResult>( new
        // CModelSpectrumResult(m_model->GetSpectrumModelContinuum()) );

        ExtremaResult->m_savedModelSpectrumResults[i] = resultspcmodel;

        m_result->LineModelSolutions[idx].fillLineIds();
        ExtremaResult->m_savedModelFittingResults[i] =
            std::make_shared<CLineModelSolution>(
                m_result->LineModelSolutions[idx]);

        // CModelRulesResult
        std::shared_ptr<CModelRulesResult> resultrulesmodel =
            std::make_shared<CModelRulesResult>(m_model->GetModelRulesLog());
        ExtremaResult->m_savedModelRulesResults[i] = resultrulesmodel;

        if (savedModels < maxSaveNLinemodelContinua) {
          // Save the reestimated continuum, only the first
          // n=maxSaveNLinemodelContinua extrema
          std::shared_ptr<CSpectraFluxResult> baselineResult =
              (std::shared_ptr<CSpectraFluxResult>)new CSpectraFluxResult();
          const CSpectrumFluxAxis &modelContinuumFluxAxis =
              m_model->GetModelContinuum();
          Int32 len = modelContinuumFluxAxis.GetSamplesCount();

          baselineResult->fluxes.resize(len);
          baselineResult->wavel.resize(len);
          for (Int32 k = 0; k < len; k++) {
            baselineResult->fluxes[k] = modelContinuumFluxAxis[k];
            baselineResult->wavel[k] = (spectrum.GetSpectralAxis())[k];
          }
          ExtremaResult->m_savedModelContinuumSpectrumResults[i] =
              baselineResult;
        }
      }
      savedModels++;
    }

    // code here has been moved to TLineModelResult::updateFromModel
    ExtremaResult->m_ranked_candidates[i].second->updateFromModel(
        m_model, m_result, m_estimateLeastSquareFast, idx, i_2pass);

    // save the continuum tpl fitting results
    ExtremaResult->m_ranked_candidates[i].second->updateContinuumFromModel(
        m_model);

    CContinuumModelSolution csolution = m_model->GetContinuumModelSolution();
    ExtremaResult->m_ranked_candidates[i]
        .second->updateFromContinuumModelSolution(csolution, false);
    ExtremaResult->m_ranked_candidates[i].second->updateTplRatioFromModel(
        m_model);
    // save the tplcorr/tplratio results
  }

  // ComputeArea2(ExtremaResult);

  m_model->SetFitContinuum_Option(savedFitContinuumOption);

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
    const CSpectrum &spectrum, const TFloat64Range &lambdaRange,
    const std::string &opt_continuumreest, const std::string &opt_fittingmethod,
    const string &opt_rigidity, const bool &opt_velocityFitting,
    const Float64 &opt_emvelocityfitmin, const Float64 &opt_emvelocityfitmax,
    const Float64 &opt_emvelocityfitstep, const Float64 &opt_absvelocityfitmin,
    const Float64 &opt_absvelocityfitmax,
    const Float64 &opt_absvelocityfitstep) {
  // setup velocity fitting
  bool enableVelocityFitting = true;
  Float64 velfitMinE = opt_emvelocityfitmin;
  Float64 velfitMaxE = opt_emvelocityfitmax;
  Float64 velfitStepE = opt_emvelocityfitstep;
  Float64 velfitMinA = opt_absvelocityfitmin;
  Float64 velfitMaxA = opt_absvelocityfitmax;
  Float64 velfitStepA = opt_absvelocityfitstep;
  // HARDCODED - override: no-velocityfitting for abs
  // velfitMinA = opt_velocityAbsorption;
  // velfitMaxA = opt_velocityAbsorption;
  if (!opt_velocityFitting) {
    enableVelocityFitting = false;
  } else {
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
        m_model->SetFitContinuum_FitStore(m_tplfitStore_secondpass[i]);
      else {
        m_model->SetFitContinuum_FitValues(
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
        m_model->SetFitContinuum_Option(2);
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
    m_model->fit(m_result->Redshifts[idx], m_result->LineModelSolutions[idx],
                 m_result->ContinuumModelSolutions[idx], contreest_iterations,
                 false);
    // m = m_result->ChiSquare[idx];
    if (enableVelocityFitting) {
      bool enableManualStepVelocityFit = true;
      bool enableLMVelocityFit = false;
      if (enableLMVelocityFit) {
        // fit the emission and absorption width using the linemodel
        // lmfit strategy
        m_model->SetFittingMethod("lmfit");
        // m_model->SetElementIndexesDisabledAuto();

        Log.LogInfo("  Operator-Linemodel: Lm fit for extrema %d", i);
        m_model->fit(
            m_result->Redshifts[idx], m_result->LineModelSolutions[idx],
            m_result->ContinuumModelSolutions[idx], contreest_iterations, true);
        mlmfit_modelInfoSave = true;
        // CModelSpectrumResult
        std::shared_ptr<CModelSpectrumResult> resultspcmodel =
            std::make_shared<CModelSpectrumResult>(m_model->GetModelSpectrum());

        mlmfit_savedModelSpectrumResults_lmfit.push_back(resultspcmodel);

        std::shared_ptr<CLineModelSolution> resultfitmodel =
            std::make_shared<CLineModelSolution>(
                m_result->LineModelSolutions[idx]);
        mlmfit_savedModelFittingResults_lmfit.push_back(resultfitmodel);

        // CModelRulesResult
        std::shared_ptr<CModelRulesResult> resultrulesmodel =
            std::shared_ptr<CModelRulesResult>(
                new CModelRulesResult(m_model->GetModelRulesLog()));
        mlmfit_savedModelRulesResults_lmfit.push_back(resultrulesmodel);

        std::shared_ptr<CSpectraFluxResult> baselineResult_lmfit =
            (std::shared_ptr<CSpectraFluxResult>)new CSpectraFluxResult();
        const CSpectrumFluxAxis &modelContinuumFluxAxis =
            m_model->GetModelContinuum();
        Int32 len = modelContinuumFluxAxis.GetSamplesCount();

        baselineResult_lmfit->fluxes.resize(len);
        baselineResult_lmfit->wavel.resize(len);
        for (Int32 k = 0; k < len; k++) {
          baselineResult_lmfit->fluxes[k] = modelContinuumFluxAxis[k];
          baselineResult_lmfit->wavel[k] = (spectrum.GetSpectralAxis())[k];
        }
        mlmfit_savedBaselineResult_lmfit.push_back(baselineResult_lmfit);

        z = m_result->LineModelSolutions[idx].Redshift;
        m_secondpass_parameters_extremaResult.Redshift_lmfit.push_back(z);

        m_model->SetFittingMethod(opt_fittingmethod);
        m_model->m_Elements.ResetElementIndexesDisabled();
        Int32 velocityHasBeenReset =
            m_model->ApplyVelocityBound(velfitMinE, velfitMaxE);
        enableManualStepVelocityFit = velocityHasBeenReset;
      }

      if (enableManualStepVelocityFit) {
        // fit the emission and absorption width by minimizing the
        // linemodel merit with linemodel "hybrid" fitting method
        m_model->SetFittingMethod("hybrid");
        if (opt_rigidity == "tplshape") {
          m_model->SetFittingMethod("individual");
        }
        m_model->SetForcedisableTplratioISMfit(
            m_model->m_opt_firstpass_forcedisableTplratioISMfit); // TODO: add
                                                                  // new param
                                                                  // for this ?
        // m_model->m_enableAmplitudeOffsets = true;
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
              idxVelfitGroups = m_model->m_Elements.GetModelVelfitGroups(
                  CLine::nType_Absorption);
              Log.LogDetail(
                  "  Operator-Linemodel: VelfitGroups ABSORPTION - n = %d",
                  idxVelfitGroups.size());
              if (m_firstpass_extremaResult->size() > 1 &&
                  idxVelfitGroups.size() > 1) {
                Log.LogError(
                    "  Operator-Linemodel: not allowed to use more than 1 "
                    "group per E/A for "
                    "more than 1 extremum (see .json linemodel.extremacount)");
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
              idxVelfitGroups = m_model->m_Elements.GetModelVelfitGroups(
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
                m_secondpass_parameters_extremaResult.ExtendedRedshifts[i]
                    .back();
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
                      m_model->SetVelocityAbsorptionOneElement(
                          vTest, idxVelfitGroups[kgroup][ke]);
                    }
                  } else {
                    m_model->SetVelocityAbsorption(vTest);
                  }
                } else {
                  if (m_enableWidthFitByGroups) {
                    for (Int32 ke = 0; ke < idxVelfitGroups[kgroup].size();
                         ke++) {
                      m_model->SetVelocityEmissionOneElement(
                          vTest, idxVelfitGroups[kgroup][ke]);
                    }
                  } else {
                    m_model->SetVelocityEmission(vTest);
                  }
                }

                // Log.LogInfo( "  Operator-Linemodel:
                // testing v=%f", vTest);
                Float64 meritv;
                Float64 zTest = m_result->Redshifts[idzTest];
                meritv = m_model->fit(
                    zTest,
                    m_result
                        ->LineModelSolutions[idx], // maybe this member result
                                                   // should be replaced by an
                                                   // unused variable
                    m_result->ContinuumModelSolutions
                        [idx], // maybe this member result should be replaced by
                               // an unused variable
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
                //                                            m_model->getModelErrorUnderElement(idxVelfitGroups[kgroup][ke]);
                //                                        }
                //                                    }

                Log.LogDebug("  Operator-Linemodel: testing velocity: "
                             "merit=%.3e for velocity = %.1f",
                             meritv, vTest);
                if (meritMin > meritv) {
                  meritMin = meritv;
                  if (iLineType == 0) {
                    vOptim = m_model->GetVelocityAbsorption();
                    z_vOptim = zTest;
                  } else {
                    vOptim = m_model->GetVelocityEmission();
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
                  for (Int32 ke = 0; ke < idxVelfitGroups[kgroup].size();
                       ke++) {
                    m_model->SetVelocityAbsorptionOneElement(
                        vOptim, idxVelfitGroups[kgroup][ke]);
                  }
                  m_secondpass_parameters_extremaResult.GroupsALv[i][kgroup] =
                      vOptim;
                } else {
                  m_model->SetVelocityAbsorption(vOptim);
                }

                m_secondpass_parameters_extremaResult.Alv[i] = vOptim;
                Log.LogDebug("    Operator-Linemodel: secondpass_parameters "
                             "extrema #%d set: alv=%.1f",
                             i, vOptim);
              } else {
                if (m_enableWidthFitByGroups) {
                  for (Int32 ke = 0; ke < idxVelfitGroups[kgroup].size();
                       ke++) {
                    m_model->SetVelocityEmissionOneElement(
                        vOptim, idxVelfitGroups[kgroup][ke]);
                  }
                  m_secondpass_parameters_extremaResult.GroupsELv[i][kgroup] =
                      vOptim;
                } else {
                  m_model->SetVelocityEmission(vOptim);
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
        m_model->SetFittingMethod(opt_fittingmethod);
        // m_model->m_enableAmplitudeOffsets = false;
        m_model->SetForcedisableTplratioISMfit(
            false); // TODO: coordinate with SetPassMode() ?
      }
    } else {
      m_secondpass_parameters_extremaResult.Elv[i] =
          m_model->GetVelocityEmission();
      m_secondpass_parameters_extremaResult.Alv[i] =
          m_model->GetVelocityAbsorption();
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
    const TFloat64Range &lambdaRange, const std::string &opt_continuumreest,
    const Int32 tplfit_option, const bool overrideRecomputeOnlyOnTheCandidate) {
  CLineModelPassExtremaResult &extremaResult =
      m_secondpass_parameters_extremaResult;
  if (extremaResult.size() < 1) {
    THROWG(INTERNAL_ERROR,
           "  Operator-Linemodel: RecomputeAroundCandidates n<1...");
  }

  Log.LogInfo("");
  Log.LogInfo(
      "  Operator-Linemodel: Second pass - recomputing around n=%d candidates",
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
      idxVelfitGroups =
          m_model->m_Elements.GetModelVelfitGroups(CLine::nType_Absorption);
      std::string alv_list_str = "";
      for (Int32 kgroup = 0; kgroup < idxVelfitGroups.size(); kgroup++) {
        for (Int32 ke = 0; ke < idxVelfitGroups[kgroup].size(); ke++) {
          m_model->SetVelocityAbsorptionOneElement(
              extremaResult.GroupsALv[i][kgroup], idxVelfitGroups[kgroup][ke]);
        }
        alv_list_str.append(boost::str(boost::format("%.2f, ") %
                                       extremaResult.GroupsALv[i][kgroup]));
      }
      Log.LogInfo("    Operator-Linemodel: recompute with groups alv=%s",
                  alv_list_str.c_str());
      // emission
      idxVelfitGroups.clear();
      idxVelfitGroups =
          m_model->m_Elements.GetModelVelfitGroups(CLine::nType_Emission);
      std::string elv_list_str = "";
      for (Int32 kgroup = 0; kgroup < idxVelfitGroups.size(); kgroup++) {
        for (Int32 ke = 0; ke < idxVelfitGroups[kgroup].size(); ke++) {
          m_model->SetVelocityEmissionOneElement(
              extremaResult.GroupsELv[i][kgroup], idxVelfitGroups[kgroup][ke]);
        }
        elv_list_str.append(boost::str(boost::format("%.2f") %
                                       extremaResult.GroupsELv[i][kgroup]));
      }
      Log.LogInfo("    Operator-Linemodel: recompute with groups elv=%s",
                  elv_list_str.c_str());

    } else {
      m_model->SetVelocityEmission(extremaResult.Elv[i]);
      m_model->SetVelocityAbsorption(extremaResult.Alv[i]);
      Log.LogInfo("    Operator-Linemodel: recompute with elv=%.1f, alv=%.1f",
                  m_model->GetVelocityEmission(),
                  m_model->GetVelocityAbsorption());
    }

    if (m_opt_continuumcomponent == "tplfit" ||
        m_opt_continuumcomponent == "tplfitauto") {
      // fix some fitcontinuum values for this extremum
      if (tplfit_option == 2) {
        m_model->SetFitContinuum_FitValues(
            extremaResult.FittedTplName[i], extremaResult.FittedTplAmplitude[i],
            extremaResult.FittedTplAmplitudeError[i],
            extremaResult.FittedTplMerit[i],
            extremaResult.FittedTplMeritPhot[i],
            extremaResult.FittedTplEbmvCoeff[i],
            extremaResult.FittedTplMeiksinIdx[i],
            extremaResult.FittedTplRedshift[i], extremaResult.FittedTplDtm[i],
            extremaResult.FittedTplMtm[i], extremaResult.FittedTplLogPrior[i],
            extremaResult.FittedTplpCoeffs[i]);
        m_model->SetFitContinuum_Option(tplfit_option);
      }
      if (tplfit_option == 0 ||
          tplfit_option == 3) // for these cases we called precompute in
                              // secondpass, so we have new fitstore
        m_model->SetFitContinuum_FitStore(m_tplfitStore_secondpass[i]);
      else {
        if (tplfit_option == 1) {
          // nothing to do cause we already injected the fitStore for cases 1
          // and 2
          m_model->SetFitContinuum_FitStore(m_tplfitStore_firstpass); // 1
        }
      }
    }

    // moved here to override the previously set option value
    // since all
    // m_model->SetFitContinuum_Option(tplfit_option);
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
    // m_model->SetFittingMethod("nofit");
    Int32 n_progresssteps = extremaResult.ExtendedRedshifts[i].size();
    Log.LogInfo("    Operator-Linemodel: Fit n=%d values for z in [%.6f; %.6f]",
                n_progresssteps, extremaResult.ExtendedRedshifts[i].front(),
                extremaResult.ExtendedRedshifts[i].back());
    for (const Float64 z : extremaResult.ExtendedRedshifts[i]) {
      const Int32 iz = CIndexing<Float64>::getIndex(m_result->Redshifts, z);
      Log.LogDetail("Fit for Extended redshift %d, z = %f", iz, z);

      m_result->ChiSquare[iz] = m_model->fit(
          m_result->Redshifts[iz], m_result->LineModelSolutions[iz],
          m_result->ContinuumModelSolutions[iz], contreest_iterations, false);
      m_result->ScaleMargCorrection[iz] = m_model->getScaleMargCorrection();
      if (m_opt_continuumcomponent == "tplfit" ||
          m_opt_continuumcomponent == "tplfitauto")
        m_result->SetChisquareTplContinuumResult(
            iz, m_model->GetFitContinuum_FitStore());
      m_result->SetChisquareTplshapeResult(
          iz, m_model->GetChisquareTplshape(), m_model->GetScaleMargTplshape(),
          m_model->GetStrongELPresentTplshape(),
          m_model->getHaELPresentTplshape(),
          m_model->GetNLinesAboveSNRTplshape(),
          m_model->GetPriorLinesTplshape());
      if (!m_estimateLeastSquareFast) {
        m_result->ChiSquareContinuum[iz] =
            m_model->getLeastSquareContinuumMerit();
      } else {
        m_result->ChiSquareContinuum[iz] =
            m_model->getLeastSquareContinuumMeritFast();
      }
      m_result->ScaleMargCorrectionContinuum[iz] =
          m_model->getContinuumScaleMargCorrection();
    }
  }

  return 0;
}

Int32 COperatorLineModel::Init(const CSpectrum &spectrum,
                               const TFloat64List &redshifts,
                               const CLineCatalog::TLineVector restLineList,
                               const TStringList &tplCategoryList,
                               const std::string &opt_continuumcomponent,
                               const Float64 nsigmasupport,
                               const Float64 halfwdwsize,
                               const Float64 radius) {
  m_RestLineList = std::move(restLineList);
  m_tplCategoryList = tplCategoryList;
  // initialize empty results so that it can be returned anyway in case of an
  // error
  m_result = std::make_shared<CLineModelResult>();

  if (spectrum.GetSpectralAxis().IsInLinearScale() == false) {
    THROWG(INTERNAL_ERROR, "  Operator-Linemodel: input spectrum is not in "
                           "linear scale (ignored).");
  }
  m_opt_continuumcomponent = opt_continuumcomponent;

  // sort the redshifts
  m_sortedRedshifts = redshifts;
  std::sort(m_sortedRedshifts.begin(), m_sortedRedshifts.end());

  // set the nsigmasupport
  m_linesmodel_nsigmasupport = nsigmasupport;
  m_secondPass_halfwindowsize = halfwdwsize;
  m_extremaRedshiftSeparation = radius;

  return 0;
}

void COperatorLineModel::InitTplratioPriors() {
  std::shared_ptr<CPriorHelper> phelperLines = make_shared<CPriorHelper>();
  phelperLines->Init(m_opt_tplratio_prior_dirpath.c_str(), 1);
  phelperLines->SetBetaA(m_opt_tplratio_prior_betaA);
  phelperLines->SetBetaTE(m_opt_tplratio_prior_betaTE);
  phelperLines->SetBetaZ(m_opt_tplratio_prior_betaZ);

  m_model->SetTplshape_PriorHelper(phelperLines);
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
    Log.LogInfo("  Operator-Linemodel: Initializing contaminant for roll #%d, "
                "with offset=%.2f",
                iRollContaminated, contLambdaOffset);

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
 *return zero. For each value of c within the range: Sum the squared difference
 *between the flux and the maximum amplitude with a exponential decay
 *parameterized by c. Save the minimal result. If the result is not greater than
 *zero, return zero. Return the result.
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
  Float64 redshift_min = datastore.GetFloat64Range("redshiftrange").GetBegin();
  Float64 redshift_max = datastore.GetFloat64Range("redshiftrange").GetEnd();
  CLineModelSolution modelSolution;
  CLineModelSolution bestModelSolution;

  CContinuumModelSolution continuumModelSolution;

  m_model->SetFittingMethod("hybrid");

    if (opt_rigidity == "tplshape")
    {
    m_model->SetFittingMethod("individual");
    }

  m_model->SetForcedisableTplratioISMfit(true);//m_model->m_opt_firstpass_forcedisableTplratioISMfit);
  //todo, add new param for this ?

  //TODO these params should be in a dedicated scope "velocityfit"
  datastore.PushScope("linemodel");
  const Float64&velfitMinE =
  datastore.GetScopedFloat64Param("emvelocityfitmin"); const Float64&velfitMaxE
  = datastore.GetScopedFloat64Param("emvelocityfitmax"); const
  Float64&velfitStepE = datastore.GetScopedFloat64Param("emvelocityfitstep");
  const Float64&velfitMinA =
  datastore.GetScopedFloat64Param("absvelocityfitmin"); const Float64&velfitMaxA
  = datastore.GetScopedFloat64Param("absvelocityfitmax"); const
  Float64&velfitStepA = datastore.GetScopedFloat64Param("absvelocityfitstep");
  const Float64&opt_manvelfit_dzmin =
  datastore.GetScopedFloat64Param("manvelocityfitdzmin"); const
  Float64&opt_manvelfit_dzmax =
  datastore.GetScopedFloat64Param("manvelocityfitdzmax"); const
  Float64&opt_manvelfit_dzstep =
  datastore.GetScopedFloat64Param("manvelocityfitdzstep"); datastore.PopScope();

  //TODO there should be a mapping between this variable and opt_continuumreest
  Int32 contreest_iterations = 0;
  // contreest_iterations = 1;


  //output
  TFloat64List GroupsElv;
  TFloat64List GroupsAlv;
  Float64 Elv,Alv;

  Float64 dzInfLim = roundf(opt_manvelfit_dzmin*10000)/10000;//set precision to
  10^4 Float64 dzStep = opt_manvelfit_dzstep; Float64 dzSupLim =
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
      m_model->m_Elements.GetModelVelfitGroups(lineType);
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
          m_model->setVelocity(vTest, idxVelfitGroups[kgroup][ke], lineType);
        }
      }
    }
  }
}

CLineModelSolution COperatorLineModel::computeForLineMeas(
    std::shared_ptr<const CInputContext> inputContext,
    const TFloat64List &redshiftsGrid, Float64 &bestz) {
  const CSpectrum &spc = *(inputContext->GetSpectrum());
  const CSpectrum &rebinnedSpc = *(inputContext->GetRebinnedSpectrum());
  const CTemplateCatalog &tplCatalog = *(inputContext->GetTemplateCatalog());
  std::shared_ptr<const CParameterStore> params =
      inputContext->GetParameterStore();

  std::shared_ptr<const CLSF> lsf = inputContext->GetSpectrum()->GetLSF();

  //  std::string opt_fittingmethod_ortho =
  //  params->GetScoped<std::string>("continuumfit.fittingmethod");
  std::string opt_lineWidthType =
      params->GetScoped<std::string>("linewidthtype");

  Float64 opt_velocityEmission = params->GetScoped<Float64>(
      "velocityemission"); // set by client, not in parameters.json
  Float64 opt_velocityAbsorption = params->GetScoped<Float64>(
      "velocityabsorption"); // set by client, not in parameters.json
  std::string opt_rules = params->GetScoped<std::string>("rules");
  std::string opt_rigidity = params->GetScoped<std::string>("rigidity");
  bool velocityfit = params->GetScoped<bool>("velocityfit");
  if (velocityfit)
    THROWG(INTERNAL_ERROR,"velocityfit not implemented yet");

  Int32 amplitudeOffsetsDegree = params->GetScoped<Int32>("polynomialdegree");
  if (amplitudeOffsetsDegree < 0 || amplitudeOffsetsDegree > 2)
    THROWG(INTERNAL_ERROR,"the polynomial degree "
        "parameter should be between 0 and 2");

  const TFloat64Range &lambdaRange = inputContext->m_lambdaRange;
  // bool opt_tplfit_ignoreLinesSupport =
  // params->GetScoped<std::string>("continuumfit.ignorelinesupport");
  const std::string opt_fittingmethod =
      "hybrid"; // params->GetScoped<std::string>("fittingmethod");
  const std::string &opt_continuumcomponent =
      "nocontinuum"; // params->GetScoped<std::string>("continuumcomponent");
  spc.SetContinuumEstimationMethod("zero");
  m_opt_continuum_neg_amp_threshold = -1; // unused params->GetScoped<Float64>(
                                          // "continuumfit.negativethreshold");

  TFloat64Range clampedlambdaRange;
  spc.GetSpectralAxis().ClampLambdaRange(lambdaRange, clampedlambdaRange);

  m_model = std::make_shared<CLineModelFitting>(
      spc, clampedlambdaRange, tplCatalog, m_tplCategoryList, m_RestLineList,
      opt_fittingmethod, opt_continuumcomponent,
      m_opt_continuum_neg_amp_threshold, opt_lineWidthType,
      m_linesmodel_nsigmasupport, opt_velocityEmission, opt_velocityAbsorption,
      opt_rules, opt_rigidity, amplitudeOffsetsDegree);

  m_model->setPassMode(3); // does m_model->m_enableAmplitudeOffsets = true;

  // init catalog offsets

  // commom between firstpass and secondpass processes
  m_phelperContinuum = std::make_shared<CPriorHelper>();
  m_phelperContinuum->Init(m_opt_tplfit_continuumprior_dirpath.c_str(), 0);
  m_phelperContinuum->SetBetaA(m_opt_tplfit_continuumprior_betaA);
  m_phelperContinuum->SetBetaTE(m_opt_tplfit_continuumprior_betaTE);
  m_phelperContinuum->SetBetaZ(m_opt_tplfit_continuumprior_betaZ);
  m_model->SetFitContinuum_PriorHelper(m_phelperContinuum);

  m_estimateLeastSquareFast = 0;
  m_model->SetLeastSquareFastEstimationEnabled(m_estimateLeastSquareFast);

  m_model->initDtd();

  CLineModelSolution modelSolution;
  CContinuumModelSolution continuumModelSolution;
  CLineModelSolution bestModelSolution;

  Float64 bestScore = DBL_MAX;
  bestz = NAN;
  for (const Float64 &z : redshiftsGrid) {
    Log.LogDebug(Formatter() << "test with z=" << z);

    Float64 score =
        m_model->fit(z, modelSolution, continuumModelSolution, 0, true);

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
  m_model->LoadModelSolution(bestModelSolution);
  m_model->refreshModel();
  return m_model->GetModelSpectrum();
}
