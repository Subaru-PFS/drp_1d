#include <RedshiftLibrary/common/mask.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/linemodel/templatesfitstore.h>
#include <RedshiftLibrary/linemodel/templatesortho.h>
#include <RedshiftLibrary/linemodel/templatesorthostore.h>
#include <RedshiftLibrary/operator/chisquare2.h>
#include <RedshiftLibrary/operator/chisquareloglambda.h>
#include <RedshiftLibrary/operator/chisquareresult.h>
#include <RedshiftLibrary/operator/linemodel.h>
#include <RedshiftLibrary/operator/spectraFluxResult.h>
#include <RedshiftLibrary/spectrum/axis.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/tools.h>
#include <RedshiftLibrary/statistics/deltaz.h>
#include <RedshiftLibrary/statistics/priorhelpercontinuum.h>

#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/spectrum/io/fitswriter.h>

#include "boost/format.hpp"
#include <boost/chrono/thread_clock.hpp>
#include <boost/numeric/conversion/bounds.hpp>
#include <boost/progress.hpp>

#include <RedshiftLibrary/processflow/datastore.h>
#include <RedshiftLibrary/spectrum/io/genericreader.h>

#include <algorithm> // std::sort
#include <assert.h>
#include <float.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_spline.h>
#include <math.h>

#define NOT_OVERLAP_VALUE NAN
#include <stdio.h>

using namespace NSEpic;
using namespace std;

/**
 **/
COperatorLineModel::COperatorLineModel() { m_maxModelSaveCount = 10; }

/**
 * \brief Empty destructor.
 **/
COperatorLineModel::~COperatorLineModel() {}

/**
 * @brief COperatorLineModel::ComputeFirstPass
 * @return 0=no errors, -1=error
 */
Int32 COperatorLineModel::ComputeFirstPass(CDataStore &dataStore,
                                           const CSpectrum &spectrum,
                                           const CSpectrum &spectrumContinuum,
                                           const CTemplateCatalog &tplCatalog,
                                           const TStringList &tplCategoryList,
                                           const std::string opt_calibrationPath,
                                           const CRayCatalog &restraycatalog,
                                           const std::string &opt_lineTypeFilter,
                                           const std::string &opt_lineForceFilter,
                                           const TFloat64Range &lambdaRange,
                                           const std::string &opt_fittingmethod,
                                           const std::string &opt_continuumcomponent,
                                           const std::string &opt_lineWidthType,
                                           const Float64 opt_resolution,
                                           const Float64 opt_velocityEmission,
                                           const Float64 opt_velocityAbsorption,
                                           const std::string &opt_continuumreest,
                                           const std::string &opt_rules,
                                           const std::string &opt_velocityFitting,
                                           const Float64 &opt_twosteplargegridstep,
                                           const string &opt_twosteplargegridsampling, const std::string &opt_rigidity,
                                           const std::string &opt_tplratioCatRelPath,
                                           const std::string &opt_offsetCatRelPath)
{

    // redefine redshift grid
    m_enableFastFitLargeGrid = 0;
    if (opt_twosteplargegridstep > 0.0)
    {
        m_enableFastFitLargeGrid = 1;
    }

    TFloat64List largeGridRedshifts;
    //*
    if (m_enableFastFitLargeGrid == 1)
    {
        // Log.LogInfo("Line Model, Fast Fit Large Grid enabled");
        // calculate on a wider grid, defined by a minimum step
        Float64 dz_thres = opt_twosteplargegridstep;
        std::vector<Int32> removed_inds;
        Int32 lastKeptInd = 0;
        for (Int32 i = 1; i < m_sortedRedshifts.size() - 1; i++)
        {
            Bool conditionKeepSample = true;
            if (opt_twosteplargegridsampling == "log")
            {
                Float64 onepz = abs(1. + m_sortedRedshifts[i]);
                conditionKeepSample =
                    (abs(m_sortedRedshifts[i] -
                         m_sortedRedshifts[lastKeptInd]) < onepz * dz_thres) &&
                    (abs(m_sortedRedshifts[i + 1] - m_sortedRedshifts[i]) <
                     onepz * dz_thres);
            } else
            {
                conditionKeepSample =
                    (abs(m_sortedRedshifts[i] -
                         m_sortedRedshifts[lastKeptInd]) < dz_thres) &&
                    (abs(m_sortedRedshifts[i + 1] - m_sortedRedshifts[i]) <
                     dz_thres);
            }
            if (conditionKeepSample)
            {
                removed_inds.push_back(i);
            } else
            {
                lastKeptInd = i;
            }
        }
        Int32 rmInd = 0;
        for (Int32 i = 1; i < m_sortedRedshifts.size() - 1; i++)
        {
            bool addToLargeGrid = true;
            if (removed_inds.size() > 0)
            {
                if (removed_inds[rmInd] == i)
                {
                    rmInd++;
                    addToLargeGrid = false;
                }
            }
            if (addToLargeGrid)
            {
                largeGridRedshifts.push_back(m_sortedRedshifts[i]);
            }
        }
        if (largeGridRedshifts.size() < 1)
        {
            m_enableFastFitLargeGrid = 0;
            Log.LogInfo("  Operator-Linemodel: FastFitLargeGrid auto disabled: "
                        "raw %d redshifts will be calculated",
                        m_sortedRedshifts.size());
        } else
        {
            Log.LogInfo(
                "  Operator-Linemodel: FastFitLargeGrid enabled: %d redshifts "
                "will be calculated on the large grid (%d initially)",
                largeGridRedshifts.size(), m_sortedRedshifts.size());
        }
    }

    Int32 typeFilter = -1;
    if (opt_lineTypeFilter == "A")
    {
        typeFilter = CRay::nType_Absorption;
    } else if (opt_lineTypeFilter == "E")
    {
        typeFilter = CRay::nType_Emission;
    }

    Int32 forceFilter = -1; // CRay::nForce_Strong;
    if (opt_lineForceFilter == "S")
    {
        forceFilter = CRay::nForce_Strong;
    }
    Log.LogDebug("restRayList force filter = %d", forceFilter);
    CRayCatalog::TRayVector restRayList =
        restraycatalog.GetFilteredList(typeFilter, forceFilter);
    Log.LogDebug("restRayList.size() = %d", restRayList.size());

    //*
    //tpl orthogonalization
    bool enableOrtho = (opt_continuumcomponent == "tplfit" || opt_continuumcomponent == "tplfitauto");
    Log.LogInfo("  Operator-Linemodel: TemplatesOrthogonalization enabled = %d", enableOrtho);

    // prepare continuum templates catalog
    CTemplatesOrthogonalization tplOrtho(
                tplCatalog,
                tplCategoryList,
                opt_calibrationPath,
                restRayList,
                opt_fittingmethod,
                opt_continuumcomponent,
                opt_lineWidthType,
                opt_resolution,
                opt_velocityEmission,
                opt_velocityAbsorption,
                opt_rules,
                opt_rigidity,
                enableOrtho);

    // CTemplateCatalog orthoTplCatalog = tplOrtho.getOrthogonalTplCatalog();
    CTemplatesOrthoStore orthoTplStore = tplOrtho.getOrthogonalTplStore();
    Int32 ctlgIdx = 0; // only one ortho config for now
    std::shared_ptr<CTemplateCatalog> orthoTplCatalog = orthoTplStore.getTplCatalog(ctlgIdx);
    Log.LogInfo("  Operator-Linemodel: Templates store prepared.");
    //*/

    m_model = std::shared_ptr<CLineModelElementList>(new CLineModelElementList(
                                                         spectrum,
                                                         spectrumContinuum,
                                                         tplCatalog,
                                                         tplCategoryList,
                                                         opt_calibrationPath,
                                                         restRayList,
                                                         opt_fittingmethod,
                                                         opt_continuumcomponent,
                                                         opt_lineWidthType,
                                                         opt_resolution,
                                                         opt_velocityEmission,
                                                         opt_velocityAbsorption,
                                                         opt_rules,
                                                         opt_rigidity));
    Float64 setssSizeInit = 0.1;
    m_model->SetSourcesizeDispersion(setssSizeInit);
    Log.LogInfo("  Operator-Linemodel: sourcesize init to: ss=%.2f",
                setssSizeInit);

    /*
    CMultiRollModel model( spectrum,
                                 spectrumContinuum,
                                 tplCatalog,//orthoTplCatalog,
                                 tplCategoryList,
                                 opt_calibrationPath,
                                 restRayList,
                                 opt_fittingmethod,
                                 opt_continuumcomponent,
                                 opt_lineWidthType,
                                 opt_resolution,
                                 opt_velocityEmission,
                                 opt_velocityAbsorption,
                                 opt_rules,
                                 opt_rigidity);


    Bool enableLoadContTemplateOverride = false; //manual switch for hardcoded
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

    //set some model parameters
    m_model->m_opt_firstpass_fittingmethod = m_opt_firstpass_fittingmethod;
    m_model->m_opt_secondpass_fittingmethod = opt_fittingmethod;


    m_model->m_opt_lya_forcefit=m_opt_lya_forcefit=="yes";
    m_model->m_opt_lya_forcedisablefit=m_opt_lya_forcedisablefit=="yes";
    m_model->m_opt_lya_fit_asym_min=m_opt_lya_fit_asym_min;
    m_model->m_opt_lya_fit_asym_max=m_opt_lya_fit_asym_max;
    m_model->m_opt_lya_fit_asym_step=m_opt_lya_fit_asym_step;
    m_model->m_opt_lya_fit_width_min=m_opt_lya_fit_width_min;
    m_model->m_opt_lya_fit_width_max=m_opt_lya_fit_width_max;
    m_model->m_opt_lya_fit_width_step=m_opt_lya_fit_width_step;
    m_model->m_opt_lya_fit_delta_min=m_opt_lya_fit_delta_min;
    m_model->m_opt_lya_fit_delta_max=m_opt_lya_fit_delta_max;
    m_model->m_opt_lya_fit_delta_step=m_opt_lya_fit_delta_step;


    if (opt_rigidity == "tplshape")
    {
        // init catalog tplratios
        Log.LogInfo("  Operator-Linemodel: Tpl-ratios init");
        bool tplratioInitRet =
            m_model->initTplratioCatalogs(opt_tplratioCatRelPath, m_opt_tplratio_ismFit);
        if (!tplratioInitRet)
        {
            Log.LogError(
                "  Operator-Linemodel: Failed to init tpl-ratios. aborting...");
            throw runtime_error("  Operator-Linemodel: Failed to init tpl-ratios. aborting...");
            return -1;
        }

        m_model->m_opt_firstpass_forcedisableTplratioISMfit = !m_opt_firstpass_tplratio_ismFit;
    }

    // init catalog offsets
    Log.LogInfo("  Operator-Linemodel: Lambda offsets init");
    try
    {
        m_model->initLambdaOffsets(opt_offsetCatRelPath);
    } catch (std::exception const &e)
    {
        Log.LogError("  Operator-Linemodel: Failed to init lambda offsets. "
                     "Continuing without offsets...");
    }

    Int32 resultInitRet = m_result->Init(m_sortedRedshifts, restRayList,
                                         m_model->getTplshape_count(),
                                         m_model->getTplshape_priors(),
                                         m_model->getTplshape_priorsPz());
    if (resultInitRet != 0)
    {
        Log.LogError("  Operator-Linemodel: ERROR while initializing linemodel "
                     "result (ret=%d)",
                     resultInitRet);
        return -1;
    }

    Log.LogInfo("  Operator-Linemodel: initialized");

    // fit continuum
    bool enableFitContinuumPrecomputed = true;
    if (enableFitContinuumPrecomputed && (opt_continuumcomponent == "tplfit" || opt_continuumcomponent == "tplfitauto") )
    {
        Float64 redshiftStepForContinuumFit = opt_twosteplargegridstep;
        if(redshiftStepForContinuumFit<=0.)
        {
            Log.LogError("  Operator-Linemodel: opt_firstpass_zgridstep must be defined for this fit continuum procedure.");
            throw runtime_error("  Operator-Linemodel: opt_firstpass_zgridstep must be defined (!= -1) for this fit continuum procedure.");
        }
        PrecomputeContinuumFit(spectrum,
                               spectrumContinuum,
                               *orthoTplCatalog,
                               tplCategoryList,
                               opt_calibrationPath,
                               lambdaRange,
                               redshiftStepForContinuumFit,
                               opt_twosteplargegridsampling,
                               m_opt_tplfit_ignoreLinesSupport);
    }else{
        if(opt_continuumcomponent == "tplfit" || opt_continuumcomponent == "tplfitauto")
        {
            m_model->m_opt_fitcontinuum_maxCount = m_opt_fitcontinuum_maxN;
        }
    }
    if(opt_continuumcomponent == "tplfit" || opt_continuumcomponent == "tplfitauto")
    {
        m_model->m_opt_firstpass_forcedisableMultipleContinuumfit = m_opt_firstpass_multiplecontinuumfit_disable;
    }

    //    //hack, zero outside of the support
    //    ///////////////////////////////////////////////////////////////////////////////////////////
    //    model.setModelSpcObservedOnSupportZeroOutside(lambdaRange);

    //    std::shared_ptr<CSpectraFluxResult> baselineResult =
    //    (std::shared_ptr<CSpectraFluxResult>) new CSpectraFluxResult();
    //    baselineResult->m_optio = 0;
    //    const CSpectrum& modelSpc = model.GetModelSpectrum();
    //    UInt32 len = modelSpc.GetSampleCount();

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

    m_result->nSpcSamples = m_model->getSpcNSamples(lambdaRange);
    m_result->dTransposeDNocontinuum =
        m_model->getDTransposeD(lambdaRange, "nocontinuum");
    m_result->dTransposeD = m_model->getDTransposeD(lambdaRange, "raw");
    m_result->cstLog = m_model->getLikelihood_cstLog(lambdaRange);

    Int32 contreest_iterations = 0;
    if (opt_continuumreest == "always")
    {
        contreest_iterations = 1;
    } else
    {
        contreest_iterations = 0;
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

    Log.LogInfo("  Operator-Linemodel: start processing");

    // Set model parameters to FIRST-PASS
    m_model->setPassMode(1);
    Log.LogInfo(
        "  Operator-Linemodel: ---------- ---------- ---------- ----------");
    Log.LogInfo("  Operator-Linemodel: now computing first-pass");
    Log.LogInfo(
        "  Operator-Linemodel: ---------- ---------- ---------- ----------");

    // WARNING: HACK, first pass with continuum from spectrum.
    // model.SetContinuumComponent("fromspectrum");
    //
    Int32 indexLargeGrid = 0;
    std::vector<Float64> calculatedLargeGridRedshifts;
    std::vector<Float64> calculatedLargeGridMerits;
    std::vector<TFloat64List> calculatedChiSquareTplshapes;
    for (Int32 k = 0; k < m_result->ChiSquareTplshapes.size(); k++)
    {
        TFloat64List _chi2tpl;
        calculatedChiSquareTplshapes.push_back(_chi2tpl);
    }
    boost::chrono::thread_clock::time_point start_mainloop =
        boost::chrono::thread_clock::now();

    boost::progress_display show_progress(m_result->Redshifts.size());
    //#pragma omp parallel for
    for (Int32 i = 0; i < m_result->Redshifts.size(); i++)
    {
        if (m_enableFastFitLargeGrid == 0 || i == 0 ||
                m_result->Redshifts[i] == largeGridRedshifts[indexLargeGrid])
        {
            m_result->ChiSquare[i] = m_model->fit(m_result->Redshifts[i],
                                                  lambdaRange,
                                                  m_result->LineModelSolutions[i],
                                                  m_result->ContinuumModelSolutions[i],
                                                  contreest_iterations,
                                                  false);
            calculatedLargeGridRedshifts.push_back(m_result->Redshifts[i]);
            calculatedLargeGridMerits.push_back(m_result->ChiSquare[i]);
            m_result->ScaleMargCorrection[i] = m_model->getScaleMargCorrection();
            m_result->SetChisquareTplshapeResult(i,
                                                 m_model->GetChisquareTplshape(),
                                                 m_model->GetScaleMargTplshape(),
                                                 m_model->GetStrongELPresentTplshape(),
                                                 m_model->GetNLinesAboveSNRTplshape());
            for (Int32 k = 0; k < m_result->ChiSquareTplshapes.size(); k++)
            {
                calculatedChiSquareTplshapes[k].push_back(
                            m_result->ChiSquareTplshapes[k][i]);
            }
            if (!m_estimateLeastSquareFast)
            {
                m_result->ChiSquareContinuum[i] =
                        m_model->getLeastSquareContinuumMerit(lambdaRange);
            } else
            {
                m_result->ChiSquareContinuum[i] =
                        m_model->getLeastSquareContinuumMeritFast();
            }
            m_result->ScaleMargCorrectionContinuum[i] =
                    m_model->getContinuumScaleMargCorrection();
            Log.LogDebug("  Operator-Linemodel: Z interval %d: Chi2 = %f", i,
                         m_result->ChiSquare[i]);
            indexLargeGrid++;
            // Log.LogInfo( "\nLineModel Infos: large grid step %d", i);
        } else
        {
            m_result->ChiSquare[i] =
                m_result->ChiSquare[i - 1] +
                1e-2; // these values will be replaced by the fine grid
                      // interpolation below...
            m_result->ScaleMargCorrection[i] =
                m_model->getScaleMargCorrection();
            m_result->LineModelSolutions[i] =
                m_result->LineModelSolutions[i - 1];
            m_result->ContinuumModelSolutions[i]=
                m_result->ContinuumModelSolutions[i - 1];
            m_result->SetChisquareTplshapeResult(
                i, m_result->GetChisquareTplshapeResult(i - 1),
                m_result->GetScaleMargCorrTplshapeResult(i - 1),
                m_result->GetStrongELPresentTplshapeResult(i - 1),
                m_result->GetNLinesAboveSNRTplshapeResult(i - 1));
            if (!m_estimateLeastSquareFast)
            {
                m_result->ChiSquareContinuum[i] =
                    m_model->getLeastSquareContinuumMerit(lambdaRange);
            } else
            {
                m_result->ChiSquareContinuum[i] =
                    m_model->getLeastSquareContinuumMeritFast();
            }
            m_result->ScaleMargCorrectionContinuum[i] =
                m_model->getContinuumScaleMargCorrection();
        }
        ++show_progress;
    }

    // now interpolate large grid merit results onto the fine grid
    if (m_result->Redshifts.size() > calculatedLargeGridMerits.size() &&
        calculatedLargeGridMerits.size() > 1)
    {
        interpolateLargeGridOnFineGrid(
            calculatedLargeGridRedshifts, m_result->Redshifts,
            calculatedLargeGridMerits, m_result->ChiSquare);
    }
    for (Int32 kts = 0; kts < m_result->ChiSquareTplshapes.size(); kts++)
    {
        if (m_result->Redshifts.size() >
                calculatedChiSquareTplshapes[kts].size() &&
            calculatedChiSquareTplshapes[kts].size() > 1)
        {
            interpolateLargeGridOnFineGrid(calculatedLargeGridRedshifts,
                                           m_result->Redshifts,
                                           calculatedChiSquareTplshapes[kts],
                                           m_result->ChiSquareTplshapes[kts]);
        }
    }

    // WARNING: HACK, first pass with continuum from spectrum.
    // model.SetContinuumComponent(opt_continuumcomponent);
    // model.InitFitContinuum();
    //
    boost::chrono::thread_clock::time_point stop_mainloop =
        boost::chrono::thread_clock::now();
    Float64 duration_mainloop =
        boost::chrono::duration_cast<boost::chrono::microseconds>(
            stop_mainloop - start_mainloop)
            .count();
    Float64 duration_firstpass_seconds = duration_mainloop / 1e6;
    Log.LogInfo("  Operator-Linemodel: first-pass done in %.4e sec",
                duration_firstpass_seconds);
    Log.LogInfo("  Operator-Linemodel: <proc-lm-firstpass><%d>",
                (Int32)duration_firstpass_seconds);

    return 0;
}


void COperatorLineModel::PrecomputeContinuumFit(const CSpectrum &spectrum,
                                                 const CSpectrum &spectrumContinuum,
                                                 const CTemplateCatalog &tplCatalog,
                                                 const TStringList &tplCategoryList,
                                                 const std::string opt_calibrationPath,
                                                 const TFloat64Range &lambdaRange,
                                                 const Float64 redshiftStep,
                                                 const string zsampling,
                                                 bool ignoreLinesSupport)
{
    boost::chrono::thread_clock::time_point start_tplfitprecompute =
        boost::chrono::thread_clock::now();
    Float64 minRedshift = m_sortedRedshifts[0];
    Float64 maxRedshift = m_sortedRedshifts[m_sortedRedshifts.size()-1];

    Log.LogInfo("  Operator-Linemodel: continuum tpl fitting: sampling:%s, "
                "step=%.5e, min=%.5e, max=%.5e",
                zsampling.c_str(),
                redshiftStep,
                minRedshift,
                maxRedshift);

    CTemplatesFitStore *tplfitStore = new CTemplatesFitStore(minRedshift,
                                                             maxRedshift,
                                                             redshiftStep,
                                                             zsampling);
    std::vector<Float64> redshiftsTplFit = tplfitStore->GetRedshiftList();
    Log.LogInfo("  Operator-Linemodel: continuum tpl redshift list n=%d",redshiftsTplFit.size());
    for(UInt32 kztplfit=0; kztplfit<std::min(Int32(redshiftsTplFit.size()), Int32(10)); kztplfit++)
    {
        Log.LogDebug("  Operator-Linemodel: continuum tpl redshift list[%d] = %f",
                    kztplfit,
                    redshiftsTplFit[kztplfit]);
    }
    std::vector<std::shared_ptr<CChisquareResult>> chisquareResultsAllTpl;
    std::vector<std::string> chisquareResultsTplName;

    std::string opt_chi2operator = "chisquarelog"; //"chisquare2"; //
    if (redshiftsTplFit.size() < 100 && opt_chi2operator != "chisquare2")
        // warning arbitrary number of redshifts threshold
        // to consider chisquare2 faster than chisquarelog
    {
        opt_chi2operator = "chisquare2";
        Log.LogInfo("  Operator-Linemodel: precomputing- auto select chisquare2 operator"
                    " (faster when only few redshifts calc. points)");
    }
    std::string opt_interp = "lin"; //"precomputedfinegrid"; //
    Log.LogInfo("  Operator-Linemodel: precomputing- with operator = %s",
                opt_chi2operator.c_str());
    Log.LogInfo("  Operator-Linemodel: precomputing-fitContinuum_dustfit = %d",
                m_opt_tplfit_dustFit);
    Log.LogInfo("  Operator-Linemodel: precomputing-fitContinuum_igm = %d",
                m_opt_tplfit_extinction);
    Log.LogInfo("  Operator-Linemodel: precomputing-fitContinuum opt_interp = %s",
                opt_interp.c_str());

    std::shared_ptr<COperator> chiSquareOperator;
    if (opt_chi2operator == "chisquarelog")
    {
        // COperatorChiSquareLogLambda* chiSquareOperator;
        bool enableLogRebin = true;
        chiSquareOperator = std::shared_ptr<COperatorChiSquareLogLambda>(
            new COperatorChiSquareLogLambda(opt_calibrationPath));
        std::shared_ptr<COperatorChiSquareLogLambda> chiSquareLogOperator =
            std::dynamic_pointer_cast<COperatorChiSquareLogLambda>(chiSquareOperator);
        chiSquareLogOperator->enableSpcLogRebin(enableLogRebin);
    } else if (opt_chi2operator == "chisquare2")
    {
        chiSquareOperator = std::shared_ptr<COperatorChiSquare2>(
            new COperatorChiSquare2(opt_calibrationPath));
    } else
    {
        Log.LogError("  Operator-Linemodel: unable to parse chisquare continuum fit operator");
    }

    Float64 overlapThreshold = 1.0;
    if (opt_chi2operator != "chisquare2" && ignoreLinesSupport==true)
    {
        ignoreLinesSupport=false;
        Log.LogWarning("  Operator-Linemodel: unable to ignoreLinesSupport if NOT chisquare2-operator is used. Disabled");
    }
    std::vector<CMask> maskList;
    if(ignoreLinesSupport)
    {
        boost::chrono::thread_clock::time_point start_tplfitmaskprep =
                boost::chrono::thread_clock::now();

        maskList.resize(redshiftsTplFit.size());
        for(UInt32 kztplfit=0; kztplfit<redshiftsTplFit.size(); kztplfit++)
        {
            m_model->initModelAtZ(redshiftsTplFit[kztplfit], lambdaRange, spectrum.GetSpectralAxis());
            maskList[kztplfit]=m_model->getOutsideLinesMask();
        }

        boost::chrono::thread_clock::time_point stop_tplfitmaskprep =
                boost::chrono::thread_clock::now();
        Float64 duration_tplfitmaskprep =
            boost::chrono::duration_cast<boost::chrono::microseconds>(
                stop_tplfitmaskprep - start_tplfitmaskprep).count();
        Float64 duration_tplfitmaskprep_seconds = duration_tplfitmaskprep / 1e6;
        Log.LogInfo("  Operator-Linemodel: tplfit-precompute mask prep done in %.4e sec",
                    duration_tplfitmaskprep_seconds);
    }

    CPriorHelperContinuum *phelperContinuum = new CPriorHelperContinuum();
    phelperContinuum->Init(m_opt_tplfit_continuumprior_reldirpath.c_str());
    phelperContinuum->SetBeta(m_opt_tplfit_continuumprior_beta);

    for (UInt32 i = 0; i < tplCategoryList.size(); i++)
    {
        std::string category = tplCategoryList[i];
        //for (UInt32 j = 0; j < orthoTplCatalog->GetTemplateCount(category); j++)
        for (UInt32 j = 0; j < tplCatalog.GetTemplateCount(category); j++)
        {
            //const CTemplate &tpl = orthoTplCatalog->GetTemplate(category, j);
            const CTemplate &tpl = tplCatalog.GetTemplate(category, j);

            Int32 opt_tplfit_integer_chi2_dustfit = -1;
            if(m_opt_tplfit_dustFit)
            {
                opt_tplfit_integer_chi2_dustfit=-10;
            }

            CPriorHelperContinuum::TPriorZEList zePriorData;
            //*
            bool retGetPrior = phelperContinuum->GetTplPriorData(tpl.GetName(), redshiftsTplFit, zePriorData);
            if(retGetPrior==false)
            {
                Log.LogError("  Operator-Linemodel: Failed to get prior for chi2 continuum precomp fit. aborting...");
                throw runtime_error("  Operator-Linemodel: Failed to get prior for chi2 continuum precomp fit. aborting...");
            }
            //Float64 priorDataLogCheck = zePriorData[0][0].logpriorTZE;
            //Log.LogInfo("  Operator-Linemodel: check prior data, zePriorData[0][0].logpriorTZE = %e", priorDataLogCheck);
            //*/


            auto chisquareResult =
                std::dynamic_pointer_cast<CChisquareResult>(
                    chiSquareOperator->Compute(
                            spectrum,
                            tpl,
                            lambdaRange,
                            redshiftsTplFit,
                            overlapThreshold,
                            maskList,
                            opt_interp,
                            m_opt_tplfit_extinction,
                            opt_tplfit_integer_chi2_dustfit,
                            zePriorData));
            if (!chisquareResult)
            {
                Log.LogInfo("  Operator-Linemodel failed to compute chi "
                            "square value for tpl=%s",
                            tpl.GetName().c_str());
            } else
            {
                chisquareResultsAllTpl.push_back(chisquareResult);
                chisquareResultsTplName.push_back(tpl.GetName());
            }
        }
    }
    chiSquareOperator.reset();

    // fill the fit store with fitted values: only the best fitted values FOR EACH TEMPLATE are used
    Float64 bestTplFitSNR = 0.0;
    Int32 nredshiftsTplFitResults = redshiftsTplFit.size();
    for (Int32 i = 0; i < nredshiftsTplFitResults; i++)
    {
        Float64 redshift = redshiftsTplFit[i];

        for (UInt32 j = 0; j < chisquareResultsAllTpl.size(); j++)
        {
            auto chisquareResult =
                std::dynamic_pointer_cast<CChisquareResult>(
                    chisquareResultsAllTpl[j]);

            Int32 retAdd = tplfitStore->Add(chisquareResultsTplName[j],
                             chisquareResult->FitDustCoeff[i],
                             chisquareResult->FitMeiksinIdx[i],
                             redshift,
                             chisquareResult->ChiSquare[i],
                             chisquareResult->FitAmplitude[i],
                             chisquareResult->FitDtM[i],
                             chisquareResult->FitMtM[i],
                             chisquareResult->LogPrior[i]);
            //Log.LogInfo("  Operator-Linemodel: check prior data, tplfitStore->Add logprior = %e", chisquareResult->LogPrior[i]);


           if(!retAdd)
           {
               Log.LogError("  Operator-Linemodel: Failed to add continuum fit to store. aborting...");
               throw runtime_error("  Operator-Linemodel: Failed to add continuum fit to store. aborting...");
           }

           Float64 tplfitsnr = -1;
           if(chisquareResult->FitMtM[i]>0.)
           {
               tplfitsnr = chisquareResult->FitDtM[i]/std::sqrt(chisquareResult->FitMtM[i]);
           }
           if(tplfitsnr>bestTplFitSNR)
           {
               bestTplFitSNR = tplfitsnr;
           }
        }
    }
    m_model->SetFitContinuum_SNRMax(bestTplFitSNR);
    Log.LogInfo("  Operator-Linemodel: fitcontinuum_snrMAX set to %f", bestTplFitSNR);

    // Set tplFitStore if needed
    m_model->SetFitContinuum_FitStore(tplfitStore);
    m_model->SetFitContinuum_PriorHelper(phelperContinuum);

    boost::chrono::thread_clock::time_point stop_tplfitprecompute =
        boost::chrono::thread_clock::now();
    Float64 duration_tplfitprecompute =
        boost::chrono::duration_cast<boost::chrono::microseconds>(
            stop_tplfitprecompute - start_tplfitprecompute).count();
    Float64 duration_tplfit_seconds = duration_tplfitprecompute / 1e6;
    Log.LogInfo("  Operator-Linemodel: tplfit-precompute done in %.4e sec",
                duration_tplfit_seconds);
    Log.LogInfo("<proc-lm-tplfit><%d>", (Int32)duration_tplfit_seconds);

    if(m_opt_fitcontinuum_maxN==-1)
    {
        m_model->m_opt_fitcontinuum_maxCount = tplfitStore->GetContinuumCount();
    }else{
        Int32 effectiveContinuumCount = std::min(m_opt_fitcontinuum_maxN, tplfitStore->GetContinuumCount());
        m_model->m_opt_fitcontinuum_maxCount = effectiveContinuumCount;
    }
    Log.LogInfo("  Operator-Linemodel: fitcontinuum_maxCount set to %d",
                m_model->m_opt_fitcontinuum_maxCount);

}

/**
 * @brief COperatorLineModel::ComputeCandidates
 * @param opt_extremacount
 * @param meritCut: optionally cut the number of candidates by merit (-1=disabled)
 * @return
 */
Int32 COperatorLineModel::ComputeCandidates(const Int32 opt_extremacount,
                                            const Int32 opt_sign,
                                            const std::vector<Float64> floatValues,
                                            const Float64 meritCut)
{
    Log.LogDebug("  Operator-Linemodel: opt_extremacount = %d", opt_extremacount);
    TFloat64Range redshiftsRange(
        m_result->Redshifts[0],
        m_result->Redshifts[m_result->Redshifts.size() - 1]);
    Log.LogDebug("  Operator-Linemodel: redshiftsRange.GetBegin() = %f, "
                 "redshiftsRange.GetEnd() = %f",
                 redshiftsRange.GetBegin(), redshiftsRange.GetEnd());

    if (m_result->Redshifts.size() == 1)
    {
        m_firstpass_extremumList.push_back(
            SPoint(m_result->Redshifts[0], m_result->ChiSquare[0]));
        Log.LogInfo("  Operator-Linemodel: found only 1 redshift calculated, "
                    "thus using only 1 extremum");
    } else if (opt_extremacount == -1)
    {
        for (Int32 ke = 0; ke < m_result->Redshifts.size(); ke++)
        {
            m_firstpass_extremumList.push_back(
                SPoint(m_result->Redshifts[ke], m_result->ChiSquare[ke]));
        }
        Log.LogInfo("  Operator-Linemodel: all initial redshifts considered as "
                    "extrema");
    } else
    {
        Log.LogInfo("  Operator-Linemodel: ChiSquare min val = %e",
                    m_result->GetMinChiSquare());
        Log.LogInfo("  Operator-Linemodel: ChiSquare max val = %e",
                    m_result->GetMaxChiSquare());
        Bool invertForMinSearch = true;
        if (opt_sign == 1)
        {
            invertForMinSearch = false;
        }
        CExtremum extremum(redshiftsRange, opt_extremacount, invertForMinSearch,
                           2);
        extremum.Find(m_result->Redshifts, floatValues, m_firstpass_extremumList);
        if (m_firstpass_extremumList.size() == 0)
        {
            Log.LogError("  Operator-Linemodel: Extremum find method failed");
            return -1;
        }

        Log.LogInfo("  Operator-Linemodel: found %d extrema",
                    m_firstpass_extremumList.size());
    }

    // remove extrema with merit threshold (input floatValues MUST be log-proba !)
    if(meritCut>0.0){
        Float64 meritThres = meritCut; //30=default logProba value for 1% missed values on PFS-cosmo noOiiDoublet
        Int32 keepMinN = 2;
        Int32 nExtrema = m_firstpass_extremumList.size();
        Int32 iExtremumFinalList = 0;
        for (Int32 i = 0; i < nExtrema; i++)
        {
            Float64 meritDiff = m_firstpass_extremumList[0].Y-m_firstpass_extremumList[iExtremumFinalList].Y;
            if(meritDiff>meritThres && i>=keepMinN)
            {
                Log.LogInfo("  Operator-Linemodel: Candidates selection by proba cut: removing i=%d, final_i=%d, e.X=%f, e.Y=%e",
                            i,
                            iExtremumFinalList,
                            m_firstpass_extremumList[iExtremumFinalList].X,
                            m_firstpass_extremumList[iExtremumFinalList].Y);
                m_firstpass_extremumList.erase(m_firstpass_extremumList.begin() + iExtremumFinalList);
            }else{
                iExtremumFinalList++;
            }
        }
    }

    /*
    // Refine Extremum with a second maximum search around the z candidates:
    // This corresponds to the finer xcorrelation in EZ Pandora (in standard_DP
    fctn in SolveKernel.py) Float64 radius = 0.001; for( Int32 i=0;
    i<m_firstpass_extremumList.size(); i++ )
    {
        Float64 x = m_firstpass_extremumList[i].X;
        Float64 left_border = max(redshiftsRange.GetBegin(), x-radius);
        Float64 right_border=min(redshiftsRange.GetEnd(), x+radius);

        TPointList m_extremumListFine;
        TFloat64Range rangeFine = TFloat64Range( left_border, right_border );
        CExtremum extremumFine( rangeFine , 1, true);
        extremumFine.Find( m_result->Redshifts, m_result->ChiSquare,
    m_extremumListFine ); if(m_extremumListFine.size()>0){ m_firstpass_extremumList[i] =
    m_extremumListFine[0];
        }
    }
    //*/



    //*
    // extend z around the extrema
    for (Int32 i = 0; i < m_firstpass_extremumList.size(); i++)
    {
        Log.LogInfo("  Operator-Linemodel: Raw extr #%d, z_e.X=%f, m_e.Y=%e", i,
                    m_firstpass_extremumList[i].X, m_firstpass_extremumList[i].Y);
        Float64 x = m_firstpass_extremumList[i].X;
        Float64 left_border =
            max(redshiftsRange.GetBegin(), x - m_secondPass_extensionradius*(1.+x));
        Float64 right_border =
            min(redshiftsRange.GetEnd(), x + m_secondPass_extensionradius*(1.+x));

        for (Int32 i = 0; i < m_result->Redshifts.size(); i++)
        {
            if (m_result->Redshifts[i] >= left_border &&
                m_result->Redshifts[i] <= right_border)
            {
                m_result->ExtremaResult.ExtremaExtendedRedshifts.push_back(
                    m_result->Redshifts[i]);
            }
        }
    }
    //*/
    // todo: remove duplicate redshifts from the extended extrema list


    //now preparing the candidates extrema results
    m_firstpass_extremaResult.Resize(m_firstpass_extremumList.size());
    for (Int32 i = 0; i < m_firstpass_extremumList.size(); i++)
    {
        Float64 z = m_firstpass_extremumList[i].X;
        Float64 m = m_firstpass_extremumList[i].Y;
        // find the index in the zaxis results
        Int32 idx = -1;
        for (UInt32 i2 = 0; i2 < m_result->Redshifts.size(); i2++)
        {
            if (m_result->Redshifts[i2] == z)
            {
                idx = i2;
                break;
            }
        }
        if (idx == -1)
        {
            Log.LogInfo("Problem. could not find extrema solution index...");
            continue;
        }

        //save basic fitting info from first pass
        m_firstpass_extremaResult.Extrema[i] = z;
        m_firstpass_extremaResult.ExtremaMerit[i] = m;
        m_firstpass_extremaResult.Elv[i] = m_result->LineModelSolutions[idx].EmissionVelocity;
        m_firstpass_extremaResult.Alv[i] = m_result->LineModelSolutions[idx].AbsorptionVelocity;

        //save the continuum fitting parameters from first pass
        m_firstpass_extremaResult.FittedTplName[i] = m_result->ContinuumModelSolutions[idx].tplName;
        m_firstpass_extremaResult.FittedTplAmplitude[i] = m_result->ContinuumModelSolutions[idx].tplAmplitude;
        m_firstpass_extremaResult.FittedTplMerit[i] = m_result->ContinuumModelSolutions[idx].tplMerit;
        m_firstpass_extremaResult.FittedTplDustCoeff[i] = m_result->ContinuumModelSolutions[idx].tplDustCoeff;
        m_firstpass_extremaResult.FittedTplMeiksinIdx[i] = m_result->ContinuumModelSolutions[idx].tplMeiksinIdx;
        m_firstpass_extremaResult.FittedTplRedshift[i] = m_result->ContinuumModelSolutions[idx].tplRedshift;
        m_firstpass_extremaResult.FittedTplDtm[i] = m_result->ContinuumModelSolutions[idx].tplDtm;
        m_firstpass_extremaResult.FittedTplMtm[i] = m_result->ContinuumModelSolutions[idx].tplMtm;
        m_firstpass_extremaResult.FittedTplLogPrior[i] = m_result->ContinuumModelSolutions[idx].tplLogPrior;
        m_firstpass_extremaResult.FittedTplpCoeffs[i] = m_result->ContinuumModelSolutions[idx].pCoeffs;

        if(m_result->ContinuumModelSolutions[idx].tplName=="")
        {
            Log.LogError(" Saving first pass extremum w. ContinuumModelSolutions tplname=%s", m_result->ContinuumModelSolutions[idx].tplName.c_str());
            Log.LogError(" Saving first pass extremum w. result idx=%d, w. m_result->Redshifts[idx]=%f", idx, m_result->Redshifts[idx]);
        }

        //... todo: more first pass results can be saved here if needed
    }


    return 0;
}

Int32 COperatorLineModel::Combine_firstpass_candidates(std::shared_ptr<CLineModelExtremaResult> firstpass_results_b)
{
    Int32 retval = 0;
    Float64 skip_thres_absdiffz = 5e-4; //threshold to remove duplicate extrema/candidates


    for (Int32 keb = 0; keb < firstpass_results_b->Extrema.size(); keb++)
    {
        Float64 z_fpb = firstpass_results_b->Extrema[keb];
        Float64 m_fpb = firstpass_results_b->ExtremaMerit[keb];
        //skip if z_fpb is nearly the same as any z_fp
        Float64 minAbsDiffz = DBL_MAX;
        for (Int32 ke = 0; ke < m_firstpass_extremaResult.Extrema.size(); ke++)
        {
            Float64 z_diff = z_fpb - m_firstpass_extremaResult.Extrema[ke];
            if(std::abs(z_diff) < minAbsDiffz)
            {
                minAbsDiffz = std::abs(z_diff);
            }
        }
        if(minAbsDiffz<skip_thres_absdiffz)
        {
            Log.LogInfo(" Combine firstpass results: dropping cand B, idx=%d, z_fpb=%f", keb, firstpass_results_b->Extrema[keb]);
            continue;
        }

        //append the candidate to m_firstpass_extremumList and m_firstpass_extremaResult
        m_firstpass_extremumList.push_back(SPoint(z_fpb, m_fpb));

        //*
        // extend z around the extrema
        TFloat64Range redshiftsRange(
            m_result->Redshifts[0],
            m_result->Redshifts[m_result->Redshifts.size() - 1]);
        Float64 left_border =
                max(redshiftsRange.GetBegin(), z_fpb - m_secondPass_extensionradius*(1.+z_fpb));
        Float64 right_border =
                min(redshiftsRange.GetEnd(), z_fpb + m_secondPass_extensionradius*(1.+z_fpb));
        for (Int32 i = 0; i < m_result->Redshifts.size(); i++)
        {
            if (m_result->Redshifts[i] >= left_border &&
                    m_result->Redshifts[i] <= right_border)
            {
                m_result->ExtremaResult.ExtremaExtendedRedshifts.push_back(
                            m_result->Redshifts[i]);
            }
        }
        //*/
        // todo: remove duplicate redshifts from the extended extrema list

        //save basic fitting info from first pass
        m_firstpass_extremaResult.Extrema.push_back(z_fpb);
        m_firstpass_extremaResult.ExtremaMerit.push_back(m_fpb);
        m_firstpass_extremaResult.Elv.push_back(firstpass_results_b->Elv[keb]);
        m_firstpass_extremaResult.Alv.push_back(firstpass_results_b->Alv[keb]);

        //save the continuum fitting parameters from first pass
        if(0) //cannot work since the fpb is linemodel without cont. tplfit...
        {
            m_firstpass_extremaResult.FittedTplName.push_back(firstpass_results_b->FittedTplName[keb]);
            m_firstpass_extremaResult.FittedTplAmplitude.push_back(firstpass_results_b->FittedTplAmplitude[keb]);
            m_firstpass_extremaResult.FittedTplMerit.push_back(firstpass_results_b->FittedTplMerit[keb]);
            m_firstpass_extremaResult.FittedTplDustCoeff.push_back(firstpass_results_b->FittedTplDustCoeff[keb]);
            m_firstpass_extremaResult.FittedTplMeiksinIdx.push_back(firstpass_results_b->FittedTplMeiksinIdx[keb]);
            m_firstpass_extremaResult.FittedTplRedshift.push_back(firstpass_results_b->FittedTplRedshift[keb]);
            m_firstpass_extremaResult.FittedTplDtm.push_back(firstpass_results_b->FittedTplDtm[keb]);
            m_firstpass_extremaResult.FittedTplMtm.push_back(firstpass_results_b->FittedTplMtm[keb]);
            m_firstpass_extremaResult.FittedTplLogPrior.push_back(firstpass_results_b->FittedTplLogPrior[keb]);
            m_firstpass_extremaResult.FittedTplpCoeffs.push_back(firstpass_results_b->FittedTplpCoeffs[keb]);
        }else{
            // find the index in the zaxis results
            Int32 idx = -1;
            for (UInt32 i2 = 0; i2 < m_result->Redshifts.size(); i2++)
            {
                if (m_result->Redshifts[i2] == z_fpb)
                {
                    idx = i2;
                    break;
                }
            }
            if (idx == -1)
            {
                Log.LogInfo("Problem. could not find fpb extrema solution index...");
                continue;
            }
            //save the continuum fitting parameters from first pass
            m_firstpass_extremaResult.FittedTplName.push_back(m_result->ContinuumModelSolutions[idx].tplName);
            m_firstpass_extremaResult.FittedTplAmplitude.push_back(m_result->ContinuumModelSolutions[idx].tplAmplitude);
            m_firstpass_extremaResult.FittedTplMerit.push_back(m_result->ContinuumModelSolutions[idx].tplMerit);
            m_firstpass_extremaResult.FittedTplDustCoeff.push_back(m_result->ContinuumModelSolutions[idx].tplDustCoeff);
            m_firstpass_extremaResult.FittedTplMeiksinIdx.push_back(m_result->ContinuumModelSolutions[idx].tplMeiksinIdx);
            m_firstpass_extremaResult.FittedTplRedshift.push_back(m_result->ContinuumModelSolutions[idx].tplRedshift);
            m_firstpass_extremaResult.FittedTplDtm.push_back(m_result->ContinuumModelSolutions[idx].tplDtm);
            m_firstpass_extremaResult.FittedTplMtm.push_back(m_result->ContinuumModelSolutions[idx].tplMtm);
            m_firstpass_extremaResult.FittedTplLogPrior.push_back(m_result->ContinuumModelSolutions[idx].tplLogPrior);
            m_firstpass_extremaResult.FittedTplpCoeffs.push_back(m_result->ContinuumModelSolutions[idx].pCoeffs);

            if(m_result->ContinuumModelSolutions[idx].tplName=="")
            {
                Log.LogError(" Saving first pass extremum w. ContinuumModelSolutions tplname=%s", m_result->ContinuumModelSolutions[idx].tplName.c_str());
                Log.LogError(" Saving first pass extremum w. result idx=%d, w. m_result->Redshifts[idx]=%f", idx, m_result->Redshifts[idx]);
            }

        }
    }


    return retval;
}

Int32 COperatorLineModel::ComputeSecondPass(CDataStore &dataStore,
                                            const CSpectrum &spectrum,
                                            const CSpectrum &spectrumContinuum,
                                            const CTemplateCatalog &tplCatalog,
                                            const TStringList &tplCategoryList,
                                            const std::string opt_calibrationPath,
                                            const CRayCatalog &restraycatalog,
                                            const std::string &opt_lineTypeFilter,
                                            const std::string &opt_lineForceFilter,
                                            const TFloat64Range &lambdaRange,
                                            const std::string &opt_fittingmethod,
                                            const std::string &opt_continuumcomponent,
                                            const std::string &opt_lineWidthType,
                                            const Float64 opt_resolution,
                                            const Float64 opt_velocityEmission,
                                            const Float64 opt_velocityAbsorption,
                                            const std::string &opt_continuumreest,
                                            const std::string &opt_rules,
                                            const std::string &opt_velocityFitting,
                                            const std::string &opt_rigidity,
                                            const Float64 &opt_emvelocityfitmin,
                                            const Float64 &opt_emvelocityfitmax,
                                            const Float64 &opt_emvelocityfitstep,
                                            const Float64 &opt_absvelocityfitmin,
                                            const Float64 &opt_absvelocityfitmax,
                                            const Float64 &opt_absvelocityfitstep,
                                            const Float64 &opt_manvelfit_dzmin,
                                            const Float64 &opt_manvelfit_dzmax,
                                            const Float64 &opt_manvelfit_dzstep,
                                            const std::string &opt_continuumfit_method)
{
    boost::chrono::thread_clock::time_point start_secondpass =
            boost::chrono::thread_clock::now();

    // Set model parameters to SECOND-PASS
    m_model->setPassMode(2);
    Int32 savedFitContinuumOption = m_model->GetFitContinuum_Option();

    Log.LogInfo("  Operator-Linemodel: ---------- ---------- ---------- ----------");
    Log.LogInfo("  Operator-Linemodel: now computing second-pass");
    Log.LogInfo("  Operator-Linemodel: ---------- ---------- ---------- ----------");

    //init lmfit variables
    mlmfit_modelInfoSave = false;
    mlmfit_savedModelSpectrumResults_lmfit.clear();
    mlmfit_savedModelFittingResults_lmfit.clear();
    mlmfit_savedModelRulesResults_lmfit.clear();
    mlmfit_savedBaselineResult_lmfit.clear();

    // estimate second pass parameters (mainly elv, alv...)
    EstimateSecondPassParameters(spectrum,
                                 lambdaRange,
                                 opt_continuumreest,
                                 opt_fittingmethod,
                                 opt_rigidity,
                                 opt_velocityFitting,
                                 opt_emvelocityfitmin,
                                 opt_emvelocityfitmax,
                                 opt_emvelocityfitstep,
                                 opt_absvelocityfitmin,
                                 opt_absvelocityfitmax,
                                 opt_absvelocityfitstep,
                                 opt_manvelfit_dzmin,
                                 opt_manvelfit_dzmax,
                                 opt_manvelfit_dzstep);

    // recompute the fine grid results around the extrema
    Int32 continnuum_fit_option = 0;
    if(opt_continuumfit_method=="fromfirstpass")
    {
        continnuum_fit_option=2;
    }else if(opt_continuumfit_method=="retryall")
    {
        continnuum_fit_option=0;
    }else{
        Log.LogError("  Operator-Linemodel: continnuum_fit_option not found: %d", continnuum_fit_option);
        throw std::runtime_error("  Operator-Linemodel: continnuum_fit_option not found");
    }
    RecomputeAroundCandidates(m_firstpass_extremumList,
                              lambdaRange,
                              opt_continuumreest,
                              continnuum_fit_option); //0: retry all cont. templates at this stage

    // additional fitting with fittingmethod=svdlcp2
    if(m_opt_secondpasslcfittingmethod=="svdlc" || m_opt_secondpasslcfittingmethod=="svdlcp2")
    {
        TPointList _secondpass_extremumList;
        for (Int32 i = 0; i < m_secondpass_parameters_extremaResult.Extrema.size(); i++)
        {
            SPoint _point(m_secondpass_parameters_extremaResult.Extrema[i],
                          m_secondpass_parameters_extremaResult.ExtremaMerit[i]);
            _secondpass_extremumList.push_back(_point);
        }
        Log.LogInfo("  Operator-Linemodel: now computing second-pass %s on each secondpass candidate (n=%d)",
                    m_opt_secondpasslcfittingmethod.c_str(),
                    _secondpass_extremumList.size());
        bool useSecondPassRedshiftValue = true;
        if(useSecondPassRedshiftValue)
        {
            for (Int32 i = 0; i < m_secondpass_parameters_extremaResult.Extrema.size(); i++)
            {
                m_secondpass_parameters_extremaResult.FittedTplRedshift[i] = m_secondpass_parameters_extremaResult.Extrema[i];
            }
        }
        m_model->SetFittingMethod(m_opt_secondpasslcfittingmethod);
        RecomputeAroundCandidates(_secondpass_extremumList,
                                  lambdaRange,
                                  opt_continuumreest,
                                  2,
                                  true);
        m_model->SetFittingMethod(opt_fittingmethod);

        Log.LogInfo("  Operator-Linemodel: now re-computing the final chi2 for each candidate");
        RecomputeAroundCandidates(_secondpass_extremumList,
                                  lambdaRange,
                                  opt_continuumreest,
                                  2);

    }

    boost::chrono::thread_clock::time_point stop_secondpass =
        boost::chrono::thread_clock::now();
    Float64 duration_secondpass =
        boost::chrono::duration_cast<boost::chrono::microseconds>(
            stop_secondpass - start_secondpass).count();
    Float64 duration_secondpass_seconds = duration_secondpass / 1e6;
    Log.LogInfo("  Operator-Linemodel: second pass done in %.4e sec",
                duration_secondpass_seconds);
    Log.LogInfo("<proc-lm-secondpass><%d>", (Int32)duration_secondpass_seconds);



    Log.LogInfo("  Operator-Linemodel: Now storing extrema results");
    Int32 extremumCount = m_secondpass_parameters_extremaResult.Extrema.size();
    m_result->ExtremaResult.Resize(extremumCount);

    // Int32 start =
    // spectrum.GetSpectralAxis().GetIndexAtWaveLength(lambdaRange.GetBegin());
    // Int32 end =
    // spectrum.GetSpectralAxis().GetIndexAtWaveLength(lambdaRange.GetEnd());
    // Int32 nsamples = end - start + 1;
    Int32 savedModels = 0;
    m_savedModelSpectrumResults.clear();
    m_savedModelFittingResults.clear();
    m_savedModelRulesResults.clear();
    m_savedModelContinuumSpectrumResults.clear();

    for (Int32 i = 0; i < extremumCount; i++)
    {
        Int32 index_extremum = m_secondpass_indiceSortedCandidatesList[i];
        Float64 z = m_secondpass_parameters_extremaResult.Extrema[index_extremum];
        Float64 m = m_secondpass_parameters_extremaResult.ExtremaMerit[index_extremum];

        m_model->SetFitContinuum_FitValues(m_secondpass_parameters_extremaResult.FittedTplName[index_extremum],
                                           m_secondpass_parameters_extremaResult.FittedTplAmplitude[index_extremum],
                                           m_secondpass_parameters_extremaResult.FittedTplMerit[index_extremum],
                                           m_secondpass_parameters_extremaResult.FittedTplDustCoeff[index_extremum],
                                           m_secondpass_parameters_extremaResult.FittedTplMeiksinIdx[index_extremum],
                                           m_secondpass_parameters_extremaResult.FittedTplRedshift[index_extremum],
                                           m_secondpass_parameters_extremaResult.FittedTplDtm[index_extremum],
                                           m_secondpass_parameters_extremaResult.FittedTplMtm[index_extremum],
                                           m_secondpass_parameters_extremaResult.FittedTplLogPrior[index_extremum],
                                           m_secondpass_parameters_extremaResult.FittedTplpCoeffs[index_extremum]);
        m_model->SetFitContinuum_Option(2);


        // find the index in the zaxis results
        Int32 idx = -1;
        for (UInt32 i2 = 0; i2 < m_result->Redshifts.size(); i2++)
        {
            if (m_result->Redshifts[i2] == z)
            {
                idx = i2;
                break;
            }
        }
        if (idx == -1)
        {
            Log.LogInfo("Problem. could not find extrema solution index...");
            continue;
        }
        Log.LogInfo("  Operator-Linemodel: Saving candidate #%d, idx=%d, z=%f, m=%f",
                    index_extremum, idx, m_result->Redshifts[idx], m_result->ChiSquare[idx]);

        // reestimate the model (eventually with continuum reestimation) on the
        // extrema selected
        Int32 contreest_iterations = 0;
        if (opt_continuumreest == "always" ||
            opt_continuumreest == "onlyextrema")
        {
            contreest_iterations = 8; // 4
        } else
        {
            contreest_iterations = 0;
        }

        m_model->SetVelocityEmission(m_secondpass_parameters_extremaResult.Elv[index_extremum]);
        m_model->SetVelocityAbsorption(m_secondpass_parameters_extremaResult.Alv[index_extremum]);


        if (!mlmfit_modelInfoSave)
        {
            m_result->ChiSquare[idx] = m_model->fit(m_result->Redshifts[idx],
                                                    lambdaRange,
                                                    m_result->LineModelSolutions[idx],
                                                    m_result->ContinuumModelSolutions[idx],
                                                    contreest_iterations,
                                                    true);
            m_result->ScaleMargCorrection[idx] = m_model->getScaleMargCorrection();
            m_result->SetChisquareTplshapeResult(idx,
                                                 m_model->GetChisquareTplshape(),
                                                 m_model->GetScaleMargTplshape(),
                                                 m_model->GetStrongELPresentTplshape(),
                                                 m_model->GetNLinesAboveSNRTplshape());
            if (!m_estimateLeastSquareFast)
            {
                m_result->ChiSquareContinuum[idx] = m_model->getLeastSquareContinuumMerit(lambdaRange);
            }else
            {
                m_result->ChiSquareContinuum[idx] = m_model->getLeastSquareContinuumMeritFast();
            }
            m_result->ScaleMargCorrectionContinuum[idx] = m_model->getContinuumScaleMargCorrection();
        }
        if (m != m_result->ChiSquare[idx])
        {
            Log.LogInfo("  Operator-Linemodel: m (%f for idx=%d) !=chi2 (%f) ",
                        m, idx, m_result->ChiSquare[idx]);
        }
        m = m_result->ChiSquare[idx]; // m_result->ChiSquare[idx];

        // save the model result
        // WARNING: saving results TODO: this is currently wrong !! the model
        // saved corresponds to the bestchi2 model. PDFs should be combined
        // prior to exporting the best model for each extrema...
        static Int32 maxModelSave = std::min(m_maxModelSaveCount, extremumCount);
        Int32 maxSaveNLinemodelContinua = maxModelSave;
        if (savedModels < maxModelSave)
        {
            if (mlmfit_modelInfoSave)
            {
                Log.LogInfo("Save model store during lm_fit");
                m_savedModelSpectrumResults.push_back(mlmfit_savedModelSpectrumResults_lmfit[index_extremum]);
                m_savedModelFittingResults.push_back(mlmfit_savedModelFittingResults_lmfit[index_extremum]);
                m_savedModelRulesResults.push_back(mlmfit_savedModelRulesResults_lmfit[index_extremum]);
                if (savedModels < maxSaveNLinemodelContinua && contreest_iterations > 0)
                {
                    m_savedModelContinuumSpectrumResults.push_back(mlmfit_savedBaselineResult_lmfit[index_extremum]);
                }
            } else
            {
                // CModelSpectrumResult
                std::shared_ptr<CModelSpectrumResult> resultspcmodel;
                Int32 overrideModelSavedType = 0;
                // 0=save model, (DEFAULT)
                // 1=save model with lines removed,
                // 2=save model with only Em. lines removed.
                if (overrideModelSavedType == 0)
                {
                    resultspcmodel = std::shared_ptr<CModelSpectrumResult>(
                        new CModelSpectrumResult(m_model->GetModelSpectrum()));
                } else if (overrideModelSavedType == 1 || overrideModelSavedType == 2)
                {
                    Int32 lineTypeFilter = -1;
                    if (overrideModelSavedType == 1)
                    {
                        lineTypeFilter = -1;
                    } else if (overrideModelSavedType == 2)
                    {
                        lineTypeFilter = CRay::nType_Emission;
                    }
                    resultspcmodel = std::shared_ptr<CModelSpectrumResult>(
                        new CModelSpectrumResult(m_model->GetObservedSpectrumWithLinesRemoved(lineTypeFilter)));
                }
                // std::shared_ptr<CModelSpectrumResult>  resultspcmodel =
                // std::shared_ptr<CModelSpectrumResult>( new
                // CModelSpectrumResult(m_model->GetSpectrumModelContinuum()) );

                m_savedModelSpectrumResults.push_back(resultspcmodel);

                // CModelFittingResult
                std::shared_ptr<CModelFittingResult> resultfitmodel =
                    std::shared_ptr<CModelFittingResult>(
                        new CModelFittingResult(
                            m_result->LineModelSolutions[idx],
                            m_result->Redshifts[idx], m_result->ChiSquare[idx],
                            m_result->restRayList,
                            m_model->GetVelocityEmission(),
                            m_model->GetVelocityAbsorption()));
                m_savedModelFittingResults.push_back(resultfitmodel);

                // CModelContinuumFittingResult
                std::shared_ptr<CModelContinuumFittingResult>
                    resultfitcontinuummodel =
                        std::shared_ptr<CModelContinuumFittingResult>(
                            new CModelContinuumFittingResult(
                                m_result->Redshifts[idx],
                                m_model->getFitContinuum_tplName(),
                                m_model->getFitContinuum_tplMerit(),
                                m_model->getFitContinuum_tplAmplitude(),
                                m_model->getFitContinuum_tplIsmDustCoeff(),
                                m_model->getFitContinuum_tplIgmMeiksinIdx(),
                                m_model->getFitContinuum_snr()));
                m_savedModelContinuumFittingResults.push_back(
                    resultfitcontinuummodel);

                // CModelRulesResult
                std::shared_ptr<CModelRulesResult> resultrulesmodel =
                    std::shared_ptr<CModelRulesResult>(
                        new CModelRulesResult(m_model->GetModelRulesLog()));
                m_savedModelRulesResults.push_back(resultrulesmodel);

                if (savedModels < maxSaveNLinemodelContinua)
                {
                    // Save the reestimated continuum, only the first
                    // n=maxSaveNLinemodelContinua extrema
                    std::shared_ptr<CSpectraFluxResult> baselineResult =
                        (std::shared_ptr<
                            CSpectraFluxResult>)new CSpectraFluxResult();
                    baselineResult->m_optio = 0;
                    const CSpectrumFluxAxis &modelContinuumFluxAxis =
                        m_model->GetModelContinuum();
                    UInt32 len = modelContinuumFluxAxis.GetSamplesCount();

                    baselineResult->fluxes.resize(len);
                    baselineResult->wavel.resize(len);
                    for (Int32 k = 0; k < len; k++)
                    {
                        baselineResult->fluxes[k] = modelContinuumFluxAxis[k];
                        baselineResult->wavel[k] =
                            (spectrum.GetSpectralAxis())[k];
                    }
                    m_savedModelContinuumSpectrumResults.push_back(
                        baselineResult);
                }
            }
            savedModels++;
        }

        m_result->ExtremaResult.Extrema[i] = z;
        m_result->ExtremaResult.ExtremaMerit[i] = m;
        m_result->ExtremaResult.Elv[i] = m_model->GetVelocityEmission();
        m_result->ExtremaResult.Alv[i] = m_model->GetVelocityAbsorption();

        if (!m_estimateLeastSquareFast)
        {
            m_result->ExtremaResult.ExtremaMeritContinuum[i] =
                m_model->getLeastSquareContinuumMerit(lambdaRange);
        } else
        {
            m_result->ExtremaResult.ExtremaMeritContinuum[i] =
                m_model->getLeastSquareContinuumMeritFast();
        }

        m_result->ExtremaResult.ExtremaLastPass[i] =
            z; // refined extremum is initialized here.

        // computing errz (or deltaz, dz...): should probably be computed in
        // linemodelresult.cpp instead ?
        Float64 dz = -1.;
        if (m_result->Redshifts.size() > 1)
        {
            CDeltaz deltaz;
            Float64 zRangeHalf = 0.005;
            TFloat64Range range = TFloat64Range(z - zRangeHalf, z + zRangeHalf);
            // Int32 ret = deltaz.Compute(m_result->ChiSquare,
            // m_result->Redshifts, z, range, dz);
            Int32 ret = deltaz.Compute3ddl(m_result->ChiSquare,
                                           m_result->Redshifts, z, range, dz);
            if (ret != 0)
            {
                Log.LogError("  Operator-Linemodel: Deltaz computation failed");
            }
        }
        m_result->ExtremaResult.DeltaZ[i] = dz;

        // store model Ha SNR & Flux
        m_result->ExtremaResult.snrHa[i] =
            m_result->LineModelSolutions[idx].snrHa;
        m_result->ExtremaResult.lfHa[i] =
            m_result->LineModelSolutions[idx].lfHa;

        // store model OII SNR & Flux
        m_result->ExtremaResult.snrOII[i] =
            m_result->LineModelSolutions[idx].snrOII;
        m_result->ExtremaResult.lfOII[i] =
            m_result->LineModelSolutions[idx].lfOII;

        // store the model norm
        m_result->ExtremaResult.mTransposeM[i] =
            m_model->EstimateMTransposeM(lambdaRange);

        // scale marginalization correction
        Float64 corrScaleMarg = m_model->getScaleMargCorrection(); //
        m_result->ExtremaResult.CorrScaleMarg[i] = corrScaleMarg;

        static Float64 cutThres = 3.0;
        Int32 nValidLines =
            m_result->GetNLinesOverCutThreshold(i, cutThres, cutThres);
        m_result->ExtremaResult.Posterior[i] =
            nValidLines; // m/Float64(1+nValidLines);
        Float64 cumulStrongELSNR = m_model->getCumulSNRStrongEL(); // getStrongerMultipleELAmpCoeff(); //
        m_result->ExtremaResult.StrongELSNR[i] = cumulStrongELSNR;

        std::vector<std::string> strongELSNRAboveCut = m_model->getLinesAboveSNR(3.5);
        m_result->ExtremaResult.StrongELSNRAboveCut[i] = strongELSNRAboveCut;


        m_result->ExtremaResult.LogArea[i] = -DBL_MAX;
        m_result->ExtremaResult.LogAreaCorrectedExtrema[i] = -1.0;

        Int32 nddl = m_model->GetNElements(); // get the total number of
                                              // elements in the model
        nddl = m_result->LineModelSolutions[idx]
                   .nDDL; // override nddl by the actual number of elements in
                          // the fitted model
        m_result->ExtremaResult.NDof[i] =
            m_model->GetModelNonZeroElementsNDdl();

        Float64 bic = m + nddl * log(m_result->nSpcSamples); // BIC
        // Float64 aic = m + 2*nddl; //AIC
        m_result->ExtremaResult.bic[i] = bic;
        // m_result->bic[i] = aic + (2*nddl*(nddl+1) )/(nsamples-nddl-1);
        // //AICc, better when nsamples small

        // compute continuum indexes
        CContinuumIndexes continuumIndexes;
        m_result->ExtremaResult.ContinuumIndexes[i] =
            continuumIndexes.getIndexes(spectrum, z);

        // save the outsideLinesMask
        m_result->ExtremaResult.OutsideLinesMask[i] =
            m_model->getOutsideLinesMask();

        m_result->ExtremaResult.OutsideLinesSTDFlux[i] = m_model->getOutsideLinesSTD(1, lambdaRange);
        m_result->ExtremaResult.OutsideLinesSTDError[i] = m_model->getOutsideLinesSTD(2, lambdaRange);
        Float64 ratioSTD = -1;
        if(m_result->ExtremaResult.OutsideLinesSTDError[i]>0.0)
        {
            ratioSTD = m_result->ExtremaResult.OutsideLinesSTDFlux[i]/m_result->ExtremaResult.OutsideLinesSTDError[i];
            Float64 ratio_thres = 1.5;
            if(abs(ratioSTD)>ratio_thres || abs(ratioSTD)<1./ratio_thres)
            {
                Log.LogWarning( "  Operator-Linemodel: STD estimations outside lines do not match: ratio=%e, flux-STD=%e, error-std=%e", ratioSTD, m_result->ExtremaResult.OutsideLinesSTDFlux[i], m_result->ExtremaResult.OutsideLinesSTDError[i]);
            }else{
                Log.LogInfo( "  Operator-Linemodel: STD estimations outside lines found matching: ratio=%e, flux-STD=%e, error-std=%e", ratioSTD, m_result->ExtremaResult.OutsideLinesSTDFlux[i], m_result->ExtremaResult.OutsideLinesSTDError[i]);
            }
        }else{
            Log.LogWarning( "  Operator-Linemodel: unable to get STD estimations..." );
        }


        // save the continuum tpl fitting results
        m_result->ExtremaResult.FittedTplName[i] =
            m_model->getFitContinuum_tplName();
        m_result->ExtremaResult.FittedTplAmplitude[i] =
            m_model->getFitContinuum_tplAmplitude();
        m_result->ExtremaResult.FittedTplMerit[i] =
            m_model->getFitContinuum_tplMerit();
        m_result->ExtremaResult.FittedTplDustCoeff[i] =
            m_model->getFitContinuum_tplIsmDustCoeff();
        m_result->ExtremaResult.FittedTplMeiksinIdx[i] =
            m_model->getFitContinuum_tplIgmMeiksinIdx();

        CContinuumModelSolution csolution = m_model->GetContinuumModelSolution();
        m_result->ExtremaResult.FittedTplRedshift[i] = csolution.tplRedshift;
        m_result->ExtremaResult.FittedTplpCoeffs[i] = csolution.pCoeffs;

        m_result->ExtremaResult.FittedTplDtm[i] = csolution.tplDtm;
        m_result->ExtremaResult.FittedTplMtm[i] = csolution.tplMtm;
        m_result->ExtremaResult.FittedTplLogPrior[i] = csolution.tplLogPrior;

        // save the tplcorr/tplratio results
        m_result->ExtremaResult.FittedTplshapeName[i] =
            m_model->getTplshape_bestTplName();
        m_result->ExtremaResult.FittedTplshapeIsmCoeff[i] =
            m_model->getTplshape_bestTplIsmCoeff();
        m_result->ExtremaResult.FittedTplshapeAmplitude[i] =
            m_model->getTplshape_bestAmplitude();
    }

    // ComputeArea2(*m_result);

    m_model->SetFitContinuum_Option(savedFitContinuumOption);

    return 0;
}

/**
 * @brief COperatorLineModel::estimateSecondPassParameters
 * - Estimates best parameters: elv and alv
 * - Store parameters for further use into:
 *
 *
 * @return
 */
Int32 COperatorLineModel::EstimateSecondPassParameters(const CSpectrum &spectrum,
                                                       const TFloat64Range &lambdaRange,
                                                       const std::string &opt_continuumreest,
                                                       const std::string &opt_fittingmethod,
                                                       const string &opt_rigidity,
                                                       const std::string &opt_velocityFitting,
                                                       const Float64 &opt_emvelocityfitmin,
                                                       const Float64 &opt_emvelocityfitmax,
                                                       const Float64 &opt_emvelocityfitstep,
                                                       const Float64 &opt_absvelocityfitmin,
                                                       const Float64 &opt_absvelocityfitmax,
                                                       const Float64 &opt_absvelocityfitstep,
                                                       const Float64 &opt_manvelfit_dzmin,
                                                       const Float64 &opt_manvelfit_dzmax,
                                                       const Float64 &opt_manvelfit_dzstep)
{
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
    if (opt_velocityFitting != "yes")
    {
        enableVelocityFitting = false;
    } else
    {
        Log.LogInfo("  Operator-Linemodel: "
                    "velocity fitting bounds for Emission: min=%.1f - max=%.1f - step=%.1f",
                    velfitMinE, velfitMaxE, velfitStepE);
        Log.LogInfo("  Operator-Linemodel: "
                    "velocity fitting bounds for Absorption: min=%.1f - max=%.1f - step=%.1f",
                    velfitMinA, velfitMaxA, velfitStepA);
    }

    // enable/disable fit by groups. Once enabled, the velocity fitting groups
    // are defined in the line catalog from v4.0 on.
    m_enableWidthFitByGroups = true;




    bool enable_secondpass_parameters_estimation = true;
    if (!enable_secondpass_parameters_estimation) //?
    {
        m_secondpass_parameters_extremaResult = m_firstpass_extremaResult;
        return 0;
    }

    m_secondpass_parameters_extremaResult.Resize(m_firstpass_extremumList.size());
    for (Int32 i = 0; i < m_firstpass_extremumList.size(); i++)
    {
        Log.LogInfo("");
        Log.LogInfo("  Operator-Linemodel: Second pass - estimate parameters for candidate #%d", i);
        Log.LogInfo("  Operator-Linemodel: ---------- /\\ ---------- ---------- ---------- Candidate #%d", i);
        Float64 z = m_firstpass_extremumList[i].X;
        Float64 m = m_firstpass_extremumList[i].Y;

        m_secondpass_parameters_extremaResult.Extrema[i] = z;
        m_secondpass_parameters_extremaResult.ExtremaMerit[i] = m;

        // fix the fitcontinuum values for this extremum : keep the 1st pass values
        if(m_opt_secondpass_estimateParms_tplfit_fixfromfirstpass)
        {
            m_secondpass_parameters_extremaResult.FittedTplName[i] = m_firstpass_extremaResult.FittedTplName[i];
            m_secondpass_parameters_extremaResult.FittedTplAmplitude[i] = m_firstpass_extremaResult.FittedTplAmplitude[i];
            m_secondpass_parameters_extremaResult.FittedTplMerit[i] = m_firstpass_extremaResult.FittedTplMerit[i];
            m_secondpass_parameters_extremaResult.FittedTplDustCoeff[i] = m_firstpass_extremaResult.FittedTplDustCoeff[i];
            m_secondpass_parameters_extremaResult.FittedTplMeiksinIdx[i] = m_firstpass_extremaResult.FittedTplMeiksinIdx[i];
            m_secondpass_parameters_extremaResult.FittedTplDtm[i] = m_firstpass_extremaResult.FittedTplDtm[i];
            m_secondpass_parameters_extremaResult.FittedTplMtm[i] = m_firstpass_extremaResult.FittedTplMtm[i];
            m_secondpass_parameters_extremaResult.FittedTplLogPrior[i] = m_firstpass_extremaResult.FittedTplLogPrior[i];
            m_secondpass_parameters_extremaResult.FittedTplRedshift[i] = m_firstpass_extremaResult.FittedTplRedshift[i];
            m_secondpass_parameters_extremaResult.FittedTplpCoeffs[i] = m_firstpass_extremaResult.FittedTplpCoeffs[i] ;

            m_model->SetFitContinuum_FitValues(m_firstpass_extremaResult.FittedTplName[i],
                                               m_firstpass_extremaResult.FittedTplAmplitude[i],
                                               m_firstpass_extremaResult.FittedTplMerit[i],
                                               m_firstpass_extremaResult.FittedTplDustCoeff[i],
                                               m_firstpass_extremaResult.FittedTplMeiksinIdx[i],
                                               m_firstpass_extremaResult.FittedTplRedshift[i],
                                               m_firstpass_extremaResult.FittedTplDtm[i],
                                               m_firstpass_extremaResult.FittedTplMtm[i],
                                               m_firstpass_extremaResult.FittedTplLogPrior[i],
                                               m_firstpass_extremaResult.FittedTplpCoeffs[i]);
            m_model->SetFitContinuum_Option(2);
        }


        // find the index in the zaxis results
        Int32 idx = -1;
        for (UInt32 i2 = 0; i2 < m_result->Redshifts.size(); i2++)
        {
            if (m_result->Redshifts[i2] == z)
            {
                idx = i2;
                break;
            }
        }
        if (idx == -1)
        {
            Log.LogInfo("  Operator-Linemodel: Problem. could not find "
                        "extrema solution index...");
            continue;
        }

        // reestimate the model (eventually with continuum reestimation) on
        // the extrema selected
        Int32 contreest_iterations = 0;
        if (opt_continuumreest == "always")
        {
            contreest_iterations = 1;
        } else
        {
            contreest_iterations = 0;
        }

        // model.LoadModelSolution(m_result->LineModelSolutions[idx]);
        m_model->fit(m_result->Redshifts[idx],
                     lambdaRange,
                     m_result->LineModelSolutions[idx],
                     m_result->ContinuumModelSolutions[idx],
                     contreest_iterations, false);
        // m = m_result->ChiSquare[idx];
        if (enableVelocityFitting)
        {
            Bool enableManualStepVelocityFit = true;
            Bool enableLMVelocityFit = false;
            if (enableLMVelocityFit)
            {
                // fit the emission and absorption width using the linemodel
                // lmfit strategy
                m_model->SetFittingMethod("lmfit");
                // m_model->SetElementIndexesDisabledAuto();

                Log.LogInfo("  Operator-Linemodel: Lm fit for extrema %d",
                            i);
                m_model->fit(m_result->Redshifts[idx], lambdaRange,
                             m_result->LineModelSolutions[idx],
                             m_result->ContinuumModelSolutions[idx],
                             contreest_iterations, true);
                mlmfit_modelInfoSave = true;
                // CModelSpectrumResult
                std::shared_ptr<CModelSpectrumResult> resultspcmodel =
                        std::shared_ptr<CModelSpectrumResult>(
                            new CModelSpectrumResult(
                                m_model->GetModelSpectrum()));
                // std::shared_ptr<CModelSpectrumResult>  resultspcmodel =
                // std::shared_ptr<CModelSpectrumResult>( new
                // CModelSpectrumResult(m_model->GetObservedSpectrumWithLinesRemoved())
                // );
                mlmfit_savedModelSpectrumResults_lmfit.push_back(resultspcmodel);
                // CModelFittingResult
                std::shared_ptr<CModelFittingResult> resultfitmodel =
                        std::shared_ptr<CModelFittingResult>(
                            new CModelFittingResult(
                                m_result->LineModelSolutions[idx],
                                m_result->Redshifts[idx],
                                m_result->ChiSquare[idx], m_result->restRayList,
                                m_model->GetVelocityEmission(),
                                m_model->GetVelocityAbsorption()));
                mlmfit_savedModelFittingResults_lmfit.push_back(resultfitmodel);
                // CModelRulesResult
                std::shared_ptr<CModelRulesResult> resultrulesmodel =
                        std::shared_ptr<CModelRulesResult>(
                            new CModelRulesResult(m_model->GetModelRulesLog()));
                mlmfit_savedModelRulesResults_lmfit.push_back(resultrulesmodel);

                std::shared_ptr<CSpectraFluxResult> baselineResult_lmfit =
                        (std::shared_ptr<
                         CSpectraFluxResult>)new CSpectraFluxResult();
                baselineResult_lmfit->m_optio = 0;
                const CSpectrumFluxAxis &modelContinuumFluxAxis =
                        m_model->GetModelContinuum();
                UInt32 len = modelContinuumFluxAxis.GetSamplesCount();

                baselineResult_lmfit->fluxes.resize(len);
                baselineResult_lmfit->wavel.resize(len);
                for (Int32 k = 0; k < len; k++)
                {
                    baselineResult_lmfit->fluxes[k] =
                            modelContinuumFluxAxis[k];
                    baselineResult_lmfit->wavel[k] =
                            (spectrum.GetSpectralAxis())[k];
                }
                mlmfit_savedBaselineResult_lmfit.push_back(baselineResult_lmfit);

                z = m_result->LineModelSolutions[idx].Redshift;
                m_result->ExtremaResult.lmfitPass.push_back(z);
                // m_result->Redshifts[idx] = z;

                m_model->SetFittingMethod(opt_fittingmethod);
                m_model->ResetElementIndexesDisabled();
                Int32 velocityHasBeenReset =
                        m_model->ApplyVelocityBound(velfitMinE, velfitMaxE);
                enableManualStepVelocityFit = velocityHasBeenReset;
            }

            if (enableManualStepVelocityFit)
            {
                // fit the emission and absorption width by minimizing the
                // linemodel merit with linemodel "hybrid" fitting method
                m_model->SetFittingMethod("hybrid");
                if (opt_rigidity == "tplshape")
                {
                    m_model->SetFittingMethod("individual");
                }
                // m_model->m_enableAmplitudeOffsets = true;
                // contreest_iterations = 1;
                std::vector<std::vector<Int32>> idxVelfitGroups;

                for (Int32 iLineType = 0; iLineType < 2; iLineType++)
                {
                    Float64 vInfLim;
                    Float64 vSupLim;
                    Float64 vStep;

                    Float64 dzInfLim = opt_manvelfit_dzmin;
                    Float64 dzStep = opt_manvelfit_dzstep;
                    Float64 dzSupLim = opt_manvelfit_dzmax;

                    if (iLineType == 0)
                    {
                        Log.LogInfo("  Operator-Linemodel: manualStep velocity fit ABSORPTION, for z = %.6f",
                                    m_result->Redshifts[idx]);
                        vInfLim = velfitMinA;
                        vSupLim = velfitMaxA;
                        vStep = velfitStepA;
                        if (m_enableWidthFitByGroups)
                        {
                            idxVelfitGroups.clear();
                            idxVelfitGroups = m_model->GetModelVelfitGroups(CRay::nType_Absorption);
                            Log.LogInfo("  Operator-Linemodel: VelfitGroups ABSORPTION - n = %d",
                                        idxVelfitGroups.size());
                            if (m_firstpass_extremumList.size() > 1 && idxVelfitGroups.size() > 1)
                            {
                                Log.LogError(
                                            "  Operator-Linemodel: not allowed to "
                                            "use more than 1 group per E/A for "
                                            "more than 1 extremum (see .json "
                                            "linemodel.extremacount)");
                            }
                        }
                    } else
                    {
                        Log.LogInfo("  Operator-Linemodel: manualStep velocity fit EMISSION, for z = %.6f",
                                    m_result->Redshifts[idx]);
                        vInfLim = velfitMinE;
                        vSupLim = velfitMaxE;
                        vStep = velfitStepE;
                        if (m_enableWidthFitByGroups)
                        {
                            idxVelfitGroups.clear();
                            idxVelfitGroups = m_model->GetModelVelfitGroups(
                                        CRay::nType_Emission);
                            Log.LogInfo("  Operator-Linemodel: VelfitGroups EMISSION - n = %d",
                                        idxVelfitGroups.size());
                            if (m_firstpass_extremumList.size() > 1 && idxVelfitGroups.size() > 1)
                            {
                                Log.LogError(
                                            "  Operator-Linemodel: not allowed to "
                                            "use more than 1 group per E/A for "
                                            "more than 1 extremum (see .json "
                                            "linemodel.extremacount)");
                            }
                        }
                    }

                    // Prepare velocity grid to be checked
                    Int32 nSteps = (int)((vSupLim - vInfLim) / vStep);

                    if (m_result->Redshifts[idx] + dzInfLim <
                            m_result->Redshifts[0])
                    {
                        dzInfLim = m_result->Redshifts[0] -
                                m_result->Redshifts[idx];
                    }
                    if (m_result->Redshifts[idx] + dzSupLim >
                            m_result->Redshifts[m_result->Redshifts.size() - 1])
                    {
                        dzSupLim =
                                m_result->Redshifts[m_result->Redshifts.size() -
                                1] -
                                m_result->Redshifts[idx];
                    }

                    Int32 nDzSteps = (int)((dzSupLim - dzInfLim) / dzStep);
                    if (nDzSteps == 0)
                    {
                        nDzSteps = 1;
                        dzInfLim = 0.;
                        dzSupLim = 0.;
                    }else
                    {
                        Log.LogInfo("  Operator-Linemodel: dzInfLim n=%e", dzInfLim);
                        Log.LogInfo("  Operator-Linemodel: dzSupLim n=%e", dzSupLim);
                        Log.LogInfo("  Operator-Linemodel: manualStep n=%d", nDzSteps);
                    }

                    Int32 n_progresssteps = idxVelfitGroups.size() * nDzSteps * nSteps;
                    boost::progress_display show_progress(n_progresssteps);
                    for (Int32 kgroup = 0; kgroup < idxVelfitGroups.size(); kgroup++)
                    {
                        Log.LogInfo("  Operator-Linemodel: manualStep fitting group=%d", kgroup);

                        Float64 meritMin = DBL_MAX;
                        Float64 vOptim = -1.0;
                        Float64 z_vOptim = -1.0;
                        for (Int32 kdz = 0; kdz < nDzSteps; kdz++)
                        {
                            Float64 dzTest = dzInfLim + kdz * dzStep;
                            for (Int32 kv = 0; kv < nSteps; kv++)
                            {
                                Float64 vTest = vInfLim + kv * vStep;
                                if (iLineType == 0)
                                {
                                    if (m_enableWidthFitByGroups)
                                    {
                                        for (Int32 ke = 0; ke < idxVelfitGroups[kgroup].size(); ke++)
                                        {
                                            m_model->SetVelocityAbsorptionOneElement(vTest,
                                                                                     idxVelfitGroups[kgroup][ke]);
                                        }
                                    } else
                                    {
                                        m_model->SetVelocityAbsorption(vTest);
                                    }
                                } else
                                {
                                    if (m_enableWidthFitByGroups)
                                    {
                                        for (Int32 ke = 0; ke < idxVelfitGroups[kgroup].size(); ke++)
                                        {
                                            m_model->SetVelocityEmissionOneElement(vTest,
                                                                                   idxVelfitGroups[kgroup][ke]);
                                        }
                                    } else
                                    {
                                        m_model->SetVelocityEmission(vTest);
                                    }
                                }

                                // Log.LogInfo( "  Operator-Linemodel:
                                // testing v=%f", vTest);
                                Float64 meritv;
                                Float64 zTest = m_result->Redshifts[idx] + dzTest*(1.+m_result->Redshifts[idx]);
                                meritv = m_model->fit(zTest,
                                                      lambdaRange,
                                                      m_result->LineModelSolutions[idx], //maybe this member result should be replaced by an unused variable
                                                      m_result->ContinuumModelSolutions[idx], //maybe this member result should be replaced by an unused variable
                                                      contreest_iterations,
                                                      false);

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
                                //                                            meritv += m_model->getModelErrorUnderElement(idxVelfitGroups[kgroup][ke]);
                                //                                        }
                                //                                    }


                                Log.LogDebug("  Operator-Linemodel: testing velocity: merit=%.3e for velocity = %.1f", meritv, vTest);
                                if (meritMin > meritv)
                                {
                                    meritMin = meritv;
                                    if (iLineType == 0)
                                    {
                                        vOptim = m_model->GetVelocityAbsorption();
                                        z_vOptim = zTest;
                                    } else
                                    {
                                        vOptim = m_model->GetVelocityEmission();
                                        z_vOptim = zTest;
                                    }
                                }
                                ++show_progress;
                            }
                        }
                        if (vOptim != -1.0)
                        {
                            Log.LogInfo("  Operator-Linemodel: best Velocity found = %.1f", vOptim);
                            m_result->ChiSquare[idx] = meritMin;
                            if (iLineType == 0)
                            {
                                if (m_enableWidthFitByGroups)
                                {
                                    for (Int32 ke = 0; ke < idxVelfitGroups[kgroup].size(); ke++)
                                    {
                                        m_model ->SetVelocityAbsorptionOneElement( vOptim,
                                                                                   idxVelfitGroups[kgroup][ke]);
                                    }
                                    //todo: elv/alv fitting per fitting groups: should be saved in m_secondpass_parameters_extremaResult[i].GroupsLv
                                } else
                                {
                                    m_model->SetVelocityAbsorption(vOptim);
                                }

                                m_secondpass_parameters_extremaResult.Alv[i] = vOptim;
                                Log.LogDebug("    Operator-Linemodel: secondpass_parameters extrema #%d set: alv=%.1f", i, vOptim);
                            } else
                            {
                                if (m_enableWidthFitByGroups)
                                {
                                    for (Int32 ke = 0; ke < idxVelfitGroups[kgroup].size(); ke++)
                                    {
                                        m_model->SetVelocityEmissionOneElement( vOptim,
                                                                                idxVelfitGroups[kgroup][ke]);
                                    }
                                    //todo: elv/alv fitting per fitting groups: should be saved in m_secondpass_parameters_extremaResult[i].GroupsLv
                                } else
                                {
                                    m_model->SetVelocityEmission(vOptim);
                                }
                                m_secondpass_parameters_extremaResult.Elv[i] = vOptim;
                                Log.LogInfo("    Operator-Linemodel: secondpass_parameters extrema #%d set: elv=%.1f (for z-optim=%.6f", i, vOptim, z_vOptim);
                            }
                        }
                    }
                }
                m_model->SetFittingMethod(opt_fittingmethod);
                // m_model->m_enableAmplitudeOffsets = false;
            }
        }else{
            m_secondpass_parameters_extremaResult.Elv[i] = m_model->GetVelocityEmission();
            m_secondpass_parameters_extremaResult.Alv[i] = m_model->GetVelocityAbsorption();

        }
    }

    return 0;
}

Int32 COperatorLineModel::RecomputeAroundCandidates(TPointList input_extremumList,
                                                    const TFloat64Range &lambdaRange,
                                                    const string &opt_continuumreest,
                                                    const Int32 tplfit_option,
                                                    const bool overrideRecomputeOnlyOnTheCandidate)
{
    if(input_extremumList.size()<1)
    {
        Log.LogError("  Operator-Linemodel: RecomputeAroundCandidates n<1...");
        throw runtime_error("  Operator-Linemodel: RecomputeAroundCandidates n<1...");
        return -1;
    }

    TPointList _secondpass_recomputed_extremumList;
    _secondpass_recomputed_extremumList.resize(input_extremumList.size());
    Log.LogInfo("  Operator-Linemodel: Second pass - recomputing around n=%d candidates", input_extremumList.size());

    bool enable_recompute_around_candidate = true;
    if (enable_recompute_around_candidate)
    {
        m_secondpass_parameters_extremaResult.ExtremaExtendedRedshifts.clear();
        for (Int32 i = 0; i < input_extremumList.size(); i++)
        {
            Log.LogInfo("");
            Log.LogInfo("  Operator-Linemodel: Second pass - recompute around Candidate #%d", i);
            Log.LogInfo("  Operator-Linemodel: ---------- /\\ ---------- ---------- ---------- Candidate #%d", i);
            Float64 z = input_extremumList[i].X;
            m_model->SetVelocityEmission(m_secondpass_parameters_extremaResult.Elv[i]);
            m_model->SetVelocityAbsorption(m_secondpass_parameters_extremaResult.Alv[i]);
            Log.LogInfo("    Operator-Linemodel: recompute with elv=%.1f, alv=%.1f",
                        m_model->GetVelocityEmission(),
                        m_model->GetVelocityAbsorption());

            // fix some fitcontinuum values for this extremum
            if(tplfit_option==2 || tplfit_option==3)
            {
                m_model->SetFitContinuum_FitValues(m_secondpass_parameters_extremaResult.FittedTplName[i],
                                                   m_secondpass_parameters_extremaResult.FittedTplAmplitude[i],
                                                   m_secondpass_parameters_extremaResult.FittedTplMerit[i],
                                                   m_secondpass_parameters_extremaResult.FittedTplDustCoeff[i],
                                                   m_secondpass_parameters_extremaResult.FittedTplMeiksinIdx[i],
                                                   m_secondpass_parameters_extremaResult.FittedTplRedshift[i],
                                                   m_secondpass_parameters_extremaResult.FittedTplDtm[i],
                                                   m_secondpass_parameters_extremaResult.FittedTplMtm[i],
                                                   m_secondpass_parameters_extremaResult.FittedTplLogPrior[i],
                                                   m_secondpass_parameters_extremaResult.FittedTplpCoeffs[i]);
            }
            m_model->SetFitContinuum_Option(tplfit_option);
            Log.LogInfo("    Operator-Linemodel: recompute with tplfit_option=%d", tplfit_option);

            // find the index in the zaxis results
            Int32 idx = -1;
            for (UInt32 i2 = 0; i2 < m_result->Redshifts.size(); i2++)
            {
                if (m_result->Redshifts[i2] == z)
                {
                    idx = i2;
                    break;
                }
            }
            if (idx == -1)
            {
                Log.LogInfo("  Operator-Linemodel: Problem. could not find "
                            "extrema solution index...");
                continue;
            }

            // reestimate the model (eventually with continuum reestimation) on
            // the extrema selected
            Int32 contreest_iterations = 0;
            if (opt_continuumreest == "always")
            {
                contreest_iterations = 1;
            } else
            {
                contreest_iterations = 0;
            }

            // finally compute the redshifts on the z-range around the extremum
            TFloat64Range redshiftsRange(
                m_result->Redshifts[0],
                m_result->Redshifts[m_result->Redshifts.size() - 1]);
            Float64 left_border = max(redshiftsRange.GetBegin(),
                                      z - m_secondPass_extensionradius*(1.+z));
            Float64 right_border = min(redshiftsRange.GetEnd(),
                                       z + m_secondPass_extensionradius*(1.+z));
            // m_model->SetFittingMethod("nofit");
            _secondpass_recomputed_extremumList[i].Y = DBL_MAX;
            _secondpass_recomputed_extremumList[i].X = m_result->Redshifts[idx];
            Int32 idx2 = idx;

            //find the candidate z-range min/max indexes
            Int32 izmin_cand = -1;
            Int32 izmax_cand = -1;
            if(!overrideRecomputeOnlyOnTheCandidate)
            {
                izmin_cand = m_result->Redshifts.size();
                izmax_cand = -1;
                for (Int32 iz = 0; iz < m_result->Redshifts.size(); iz++)
                {
                    if (m_result->Redshifts[iz] >= left_border &&
                            m_result->Redshifts[iz] <= right_border)
                    {
                        if(izmin_cand>iz)
                        {
                            izmin_cand = iz;
                        }
                        if(izmax_cand<iz)
                        {
                            izmax_cand = iz;
                        }
                    }
                }
            }else{
                izmin_cand = idx;
                izmax_cand = idx;
            }

            Int32 n_progresssteps = izmax_cand-izmin_cand+1;
            Log.LogInfo("    Operator-Linemodel: Fit n=%d values for z in [%.6f; %.6f]",
                        n_progresssteps,
                        m_result->Redshifts[izmin_cand],
                        m_result->Redshifts[izmax_cand]);
            boost::progress_display show_progress(n_progresssteps);
            for (Int32 iz = izmin_cand; iz <= izmax_cand; iz++)
            {
                Log.LogDetail("Fit for Extended redshift %d, z = %f", iz, m_result->Redshifts[iz]);
                m_secondpass_parameters_extremaResult.ExtremaExtendedRedshifts.push_back(m_result->Redshifts[iz]);

                m_result->ChiSquare[iz] =
                        m_model->fit(m_result->Redshifts[iz],
                                     lambdaRange,
                                     m_result->LineModelSolutions[iz],
                                     m_result->ContinuumModelSolutions[iz],
                                     contreest_iterations, false);
                m_result->ScaleMargCorrection[iz] =
                        m_model->getScaleMargCorrection();
                m_result->SetChisquareTplshapeResult(iz,
                                                     m_model->GetChisquareTplshape(),
                                                     m_model->GetScaleMargTplshape(),
                                                     m_model->GetStrongELPresentTplshape(),
                                                     m_model->GetNLinesAboveSNRTplshape());
                if (!m_estimateLeastSquareFast)
                {
                    m_result->ChiSquareContinuum[iz] =
                            m_model->getLeastSquareContinuumMerit(lambdaRange);
                } else
                {
                    m_result->ChiSquareContinuum[iz] =
                            m_model->getLeastSquareContinuumMeritFast();
                }
                m_result->ScaleMargCorrectionContinuum[iz] =
                        m_model->getContinuumScaleMargCorrection();
                if (m_result->ChiSquare[iz] < _secondpass_recomputed_extremumList[i].Y)
                {
                    _secondpass_recomputed_extremumList[i].X = m_result->Redshifts[iz];
                    _secondpass_recomputed_extremumList[i].Y = m_result->ChiSquare[iz];

                    // set the second pass parameters used in the model export procedure in computeSecondPass()
                    m_secondpass_parameters_extremaResult.Extrema[i] = m_result->Redshifts[iz];
                    m_secondpass_parameters_extremaResult.ExtremaMerit[i] = m_result->ChiSquare[iz];
                    CContinuumModelSolution csolution = m_model->GetContinuumModelSolution();
                    m_secondpass_parameters_extremaResult.FittedTplName[i] = csolution.tplName;
                    m_secondpass_parameters_extremaResult.FittedTplAmplitude[i] = csolution.tplAmplitude;
                    m_secondpass_parameters_extremaResult.FittedTplMerit[i] = csolution.tplMerit;
                    m_secondpass_parameters_extremaResult.FittedTplDustCoeff[i] = csolution.tplDustCoeff;
                    m_secondpass_parameters_extremaResult.FittedTplMeiksinIdx[i] = csolution.tplMeiksinIdx;
                    m_secondpass_parameters_extremaResult.FittedTplRedshift[i] = csolution.tplRedshift;
                    m_secondpass_parameters_extremaResult.FittedTplDtm[i] = csolution.tplDtm;
                    m_secondpass_parameters_extremaResult.FittedTplMtm[i] = csolution.tplMtm;
                    m_secondpass_parameters_extremaResult.FittedTplLogPrior[i] = csolution.tplLogPrior;
                    m_secondpass_parameters_extremaResult.FittedTplpCoeffs[i] = csolution.pCoeffs;


                    idx2 = iz;
                }
                ++show_progress;
            }
            // m_model->SetFittingMethod(opt_fittingmethod);


            Log.LogInfo("  Operator-Linemodel: Recomputed extr #%d, idx=%d, "
                        "z_e.X=%f, m_e.Y=%f",
                        i,
                        idx2,
                        _secondpass_recomputed_extremumList[i].X,
                        _secondpass_recomputed_extremumList[i].Y);
            Log.LogInfo("  Operator-Linemodel: Recomputed extr #%d, FittedTplName=%s",
                        i,
                        m_secondpass_parameters_extremaResult.FittedTplName[i].c_str());
            Log.LogInfo("  Operator-Linemodel: Recomputed extr #%d, FittedTplAmplitude=%.4e",
                        i,
                        m_secondpass_parameters_extremaResult.FittedTplAmplitude[i]);
            Log.LogInfo("  Operator-Linemodel: Recomputed extr #%d, FittedTplDustCoeff=%f, FittedTplMeiksinIdx=%d",
                        i,
                        m_secondpass_parameters_extremaResult.FittedTplDustCoeff[i],
                        m_secondpass_parameters_extremaResult.FittedTplMeiksinIdx[i]);
            Log.LogInfo("  Operator-Linemodel: Recomputed extr #%d, FittedTplLogPrior=%e",
                        i,
                        m_secondpass_parameters_extremaResult.FittedTplLogPrior[i]);
            Log.LogInfo("  Operator-Linemodel: Recomputed extr #%d, FittedTplRedshiftf=%.6f",
                        i,
                        m_secondpass_parameters_extremaResult.FittedTplRedshift[i]);

            Float64 pCoeff0 = -1;
            Float64 pCoeff1 = -1;
            Float64 pCoeff2 = -1;
            if(m_secondpass_parameters_extremaResult.FittedTplpCoeffs[i].size()>0)
            {
                pCoeff0=m_secondpass_parameters_extremaResult.FittedTplpCoeffs[i][0];
            }
            if(m_secondpass_parameters_extremaResult.FittedTplpCoeffs[i].size()>1)
            {
                pCoeff1=m_secondpass_parameters_extremaResult.FittedTplpCoeffs[i][1];
            }
            if(m_secondpass_parameters_extremaResult.FittedTplpCoeffs[i].size()>2)
            {
                pCoeff2=m_secondpass_parameters_extremaResult.FittedTplpCoeffs[i][2];
            }
            Log.LogInfo("  Operator-Linemodel: Recomputed extr #%d, FittedTplpCoeffs_0=%.4e, "
                        "FittedTplpCoeffs_1=%.4e, FittedTplpCoeffs_2=%.4e",
                        i,
                        pCoeff0,
                        pCoeff1,
                        pCoeff2);
        }
    } else
    {
        for (Int32 i = 0; i < input_extremumList.size(); i++)
        {
            _secondpass_recomputed_extremumList[i].X = input_extremumList[i].X;
            _secondpass_recomputed_extremumList[i].Y = input_extremumList[i].Y;
        }
    }

    // sort extremumList using merit values : smallest to highest
    // m_secondpass_indiceSortedCandidatesList will contain the indexes order
    // todo: recode using map and sort from std lib
    Int32 extremumCount = m_secondpass_parameters_extremaResult.Extrema.size();
    m_secondpass_indiceSortedCandidatesList.clear();
    for (Int32 ie = 0; ie < extremumCount; ie++)
    {
        Int32 iYmin = 0;
        Float64 YMin = DBL_MAX;
        for (Int32 ie2 = 0; ie2 < _secondpass_recomputed_extremumList.size(); ie2++)
        {
            if (YMin > _secondpass_recomputed_extremumList[ie2].Y)
            {
                YMin = _secondpass_recomputed_extremumList[ie2].Y;
                iYmin = ie2;
            }
        }
        _secondpass_recomputed_extremumList.erase(_secondpass_recomputed_extremumList.begin() + iYmin);

        // find the initial index in _secondpass_recomputed_extremumList
        for (Int32 ie2 = 0; ie2 < m_secondpass_indiceSortedCandidatesList.size(); ie2++)
        {
            if (iYmin >= m_secondpass_indiceSortedCandidatesList[ie2])
            {
                iYmin++;
            }
        }
        m_secondpass_indiceSortedCandidatesList.push_back(iYmin);
    }
    if (mlmfit_modelInfoSave)
    {
        TFloat64List OrderedLMZ;
        for (Int32 ie2 = 0; ie2 < m_secondpass_indiceSortedCandidatesList.size(); ie2++)
        {
            OrderedLMZ.push_back(
                m_result->ExtremaResult.lmfitPass[m_secondpass_indiceSortedCandidatesList[ie2]]);
        }
        m_result->ExtremaResult.lmfitPass = OrderedLMZ;
    }
    if (m_secondpass_indiceSortedCandidatesList.size() == 0)
    {
        Log.LogError("  Operator-Linemodel: Extremum Ordering failed");
    }

    return 0;
}

Int32 COperatorLineModel::Init(const CSpectrum &spectrum,
                               const TFloat64List &redshifts)
{
    // initialize empty results so that it can be returned anyway in case of an
    // error
    m_result = std::shared_ptr<CLineModelResult>(new CLineModelResult());

    if (spectrum.GetSpectralAxis().IsInLinearScale() == false)
    {
        Log.LogError(
            "Line Model, input spectrum is not in linear scale (ignored).");
        return -1;
    }

    // sort the redshifts
    m_sortedRedshifts = redshifts;
    std::sort(m_sortedRedshifts.begin(), m_sortedRedshifts.end());

    return 0;
}

std::shared_ptr<COperatorResult> COperatorLineModel::getResult()
{
    return m_result;
}


std::shared_ptr<CLineModelExtremaResult> COperatorLineModel::GetFirstpassExtremaResult() const
{
    std::shared_ptr<CLineModelExtremaResult> extremaresult = std::shared_ptr<CLineModelExtremaResult>(new CLineModelExtremaResult(m_firstpass_extremaResult));
    return extremaresult;
}

std::shared_ptr<COperatorResult> COperatorLineModel::Compute(
        CDataStore &dataStore,
        const CSpectrum &spectrum,
        const CSpectrum &spectrumContinuum,
        const CTemplateCatalog &tplCatalog,
        const TStringList &tplCategoryList,
        const std::string opt_calibrationPath,
        const CRayCatalog &restraycatalog,
        const std::string &opt_lineTypeFilter,
        const std::string &opt_lineForceFilter,
        const TFloat64Range &lambdaRange,
        const TFloat64List &redshifts,
        const Int32 opt_extremacount,
        const std::string &opt_fittingmethod,
        const std::string &opt_continuumcomponent,
        const std::string &opt_lineWidthType,
        const Float64 opt_resolution,
        const Float64 opt_velocityEmission,
        const Float64 opt_velocityAbsorption,
        const std::string &opt_continuumreest,
        const std::string &opt_rules,
        const std::string &opt_velocityFitting,
        const Float64 &opt_twosteplargegridstep,
        const string &opt_twosteplargegridsampling,
        const std::string &opt_rigidity,
        const string &opt_tplratioCatRelPath,
        const string &opt_offsetCatRelPath,
        const Float64 &opt_emvelocityfitmin,
        const Float64 &opt_emvelocityfitmax,
        const Float64 &opt_emvelocityfitstep,
        const Float64 &opt_absvelocityfitmin,
        const Float64 &opt_absvelocityfitmax,
        const Float64 &opt_absvelocityfitstep,
        const Float64 &opt_manvelfit_dzmin,
        const Float64 &opt_manvelfit_dzmax,
        const Float64 &opt_manvelfit_dzstep)
{
    // initialize empty results so that it can be returned anyway in case of an
    // error
    m_result = std::shared_ptr<CLineModelResult>(new CLineModelResult());

    if (spectrum.GetSpectralAxis().IsInLinearScale() == false)
    {
        Log.LogError(
            "Line Model, input spectrum is not in linear scale (ignored).");
    }

    // sort the redshifts
    m_sortedRedshifts = redshifts;
    std::sort(m_sortedRedshifts.begin(), m_sortedRedshifts.end());

    //**************************************************
    // FIRST PASS
    //**************************************************
    Int32 retFirstPass = ComputeFirstPass(dataStore,
                                          spectrum,
                                          spectrumContinuum,
                                          tplCatalog,
                                          tplCategoryList,
                                          opt_calibrationPath,
                                          restraycatalog,
                                          opt_lineTypeFilter,
                                          opt_lineForceFilter,
                                          lambdaRange,
                                          opt_fittingmethod,
                                          opt_continuumcomponent,
                                          opt_lineWidthType,
                                          opt_resolution,
                                          opt_velocityEmission,
                                          opt_velocityAbsorption,
                                          opt_continuumreest,
                                          opt_rules, opt_velocityFitting,
                                          opt_twosteplargegridstep,
                                          opt_twosteplargegridsampling,
                                          opt_rigidity,
                                          opt_tplratioCatRelPath,
                                          opt_offsetCatRelPath);
    if (retFirstPass != 0)
    {
        Log.LogError("Line Model, first pass failed. Aborting");
        return m_result;
    }

    //**************************************************
    // Compute z-candidates
    //**************************************************

    Log.LogInfo("Line Model, compute z-candidates from best-chisquare values "
                "(nb. no priors)");
    Int32 retCandidates =
        ComputeCandidates(opt_extremacount, -1, m_result->ChiSquare, -1);
    if (retCandidates != 0)
    {
        Log.LogError("Line Model, compute z-candidates failed. Aborting");
        return m_result;
    }

    //**************************************************
    // SECOND PASS
    //**************************************************
    Int32 retSecondPass = ComputeSecondPass(
                dataStore,
                spectrum,
                spectrumContinuum,
                tplCatalog,
                tplCategoryList,
                opt_calibrationPath,
                restraycatalog,
                opt_lineTypeFilter,
                opt_lineForceFilter,
                lambdaRange,
                opt_fittingmethod,
                opt_continuumcomponent,
                opt_lineWidthType,
                opt_resolution,
                opt_velocityEmission,
                opt_velocityAbsorption,
                opt_continuumreest,
                opt_rules,
                opt_velocityFitting,
                opt_rigidity,
                opt_emvelocityfitmin,
                opt_emvelocityfitmax,
                opt_emvelocityfitstep,
                opt_absvelocityfitmin,
                opt_absvelocityfitmax,
                opt_absvelocityfitstep,
                opt_manvelfit_dzmin,
                opt_manvelfit_dzmax,
                opt_manvelfit_dzstep);
    if (retSecondPass != 0)
    {
        Log.LogError("Line Model, second pass failed. Aborting");
        return m_result;
    }

    return m_result;
}

/**
 * @brief COperatorLineModel::processPass
 * @return
 */
std::shared_ptr<COperatorResult> COperatorLineModel::computeWithUltimPass(
    CDataStore &dataStore, const CSpectrum &spectrum,
    const CSpectrum &spectrumContinuum, const CTemplateCatalog &tplCatalog,
    const TStringList &tplCategoryList, const std::string opt_calibrationPath,
    const CRayCatalog &restraycatalog, const std::string &opt_lineTypeFilter,
    const std::string &opt_lineForceFilter, const TFloat64Range &lambdaRange,
    const TFloat64List &redshifts, const Int32 opt_extremacount,
    const std::string &opt_fittingmethod,
    const std::string &opt_continuumcomponent,
    const std::string &opt_lineWidthType, const Float64 opt_resolution,
    const Float64 opt_velocityEmission, const Float64 opt_velocityAbsorption,
    const std::string &opt_continuumreest, const std::string &opt_rules,
    const std::string &opt_velocityFitting,
    const Float64 &opt_twosteplargegridstep,
    const string &opt_twosteplargegridsampling, const std::string &opt_rigidity,
    const string &opt_tplratioCatRelPath, const string &opt_offsetCatRelPath,
    const Float64 &opt_emvelocityfitmin, const Float64 &opt_emvelocityfitmax,
    const Float64 &opt_emvelocityfitstep, const Float64 &opt_absvelocityfitmin,
    const Float64 &opt_absvelocityfitmax, const Float64 &opt_absvelocityfitstep)
{
    auto result = std::dynamic_pointer_cast<CLineModelResult>(Compute(
        dataStore, spectrum, spectrumContinuum, tplCatalog, tplCategoryList,
        opt_calibrationPath, restraycatalog, opt_lineTypeFilter,
        opt_lineForceFilter, lambdaRange, redshifts, opt_extremacount,
        opt_fittingmethod, opt_continuumcomponent, opt_lineWidthType,
        opt_resolution, opt_velocityEmission, opt_velocityAbsorption,
        opt_continuumreest, opt_rules, opt_velocityFitting,
        opt_twosteplargegridstep, opt_twosteplargegridsampling, opt_rigidity,
        opt_tplratioCatRelPath, opt_offsetCatRelPath, opt_emvelocityfitmin,
        opt_emvelocityfitmax, opt_emvelocityfitstep, opt_absvelocityfitmin,
        opt_absvelocityfitmax, opt_absvelocityfitstep));

    if (result && opt_rigidity == "tplshape")
    {
        Log.LogInfo("  Operator-Linemodel - Last Pass: begin");
        //
        // do the last pass on the 1st extremum range
        //
        Float64 halfRange = 1e-3;
        Float64 lastPassStep = 1e-4;
        // find the best extremum redhift
        Float64 bestRedshift = -1;
        Float64 bestMerit = DBL_MAX;
        Int32 bestIndex = -1;
        for (Int32 k = 0; k < result->ExtremaResult.Extrema.size(); k++)
        {
            Log.LogInfo("  Operator-Linemodel - Last Pass: k = %d", k);
            Log.LogInfo(
                "  Operator-Linemodel - Last Pass: result->Extrema[k] = %.5f",
                result->ExtremaResult.Extrema[k]);
            Log.LogInfo("  Operator-Linemodel - Last Pass: "
                        "result->ExtremaMerit[k] = %.5f",
                        result->ExtremaResult.ExtremaMerit[k]);

            if (bestMerit > result->ExtremaResult.ExtremaMerit[k])
            {
                bestMerit = result->ExtremaResult.ExtremaMerit[k];
                bestRedshift = result->ExtremaResult.Extrema[k];
                bestIndex = k;
            }
        }
        Log.LogInfo("  Operator-Linemodel - Last Pass: around extrema z = %.5f",
                    bestRedshift);
        Float64 z = bestRedshift - halfRange;
        std::vector<Float64> lastPassRedshifts;
        if (z < redshifts[0])
        {
            z = redshifts[0];
        }
        while (z < bestRedshift + halfRange)
        {
            lastPassRedshifts.push_back(z);
            z += lastPassStep;
        }
        Log.LogInfo(
            "  Operator-Linemodel - Last Pass: range zmin=%.5f, zmax=%.5f",
            lastPassRedshifts[0],
            lastPassRedshifts[lastPassRedshifts.size() - 1]);

        std::string opt_rigidity_lastPass = "rules";
        Int32 opt_extremacount_lastPass = 1;

        Int32 maxSaveBackup = m_maxModelSaveCount;
        m_maxModelSaveCount = 0;
        auto lastPassResult = std::dynamic_pointer_cast<CLineModelResult>(
            Compute(dataStore, spectrum, spectrumContinuum, tplCatalog,
                    tplCategoryList, opt_calibrationPath, restraycatalog,
                    opt_lineTypeFilter, opt_lineForceFilter, lambdaRange,
                    lastPassRedshifts, opt_extremacount_lastPass,
                    opt_fittingmethod, opt_continuumcomponent,
                    opt_lineWidthType, opt_resolution, opt_velocityEmission,
                    opt_velocityAbsorption, opt_continuumreest, opt_rules,
                    opt_velocityFitting, opt_twosteplargegridstep,
                    opt_twosteplargegridsampling, opt_rigidity_lastPass,
                    opt_tplratioCatRelPath, opt_offsetCatRelPath,
                    opt_emvelocityfitmin, opt_emvelocityfitmax,
                    opt_emvelocityfitstep, opt_absvelocityfitmin,
                    opt_absvelocityfitmax, opt_absvelocityfitstep));

        m_maxModelSaveCount = maxSaveBackup;
        Float64 refinedExtremum = lastPassResult->ExtremaResult.Extrema[0];
        Log.LogInfo("  Operator-Linemodel - Last Pass: found refined z=%.5f",
                    refinedExtremum);

        result->ExtremaResult.ExtremaLastPass[bestIndex] = refinedExtremum;
    } else
    {
        Log.LogInfo("  Operator-Linemodel - Last Pass: failed to do linemodel, "
                    "rigidity=%s",
                    opt_rigidity.c_str());
    }

    return result;
}

/**
 * @brief COperatorLineModel::initContaminant
 * prepare the contaminant to be used when instantiating a multimodel in
 * compute()
 * @return
 */
Int32 COperatorLineModel::initContaminant(
    std::shared_ptr<CModelSpectrumResult> contModelSpectrum,
    Int32 iRollContaminated, Float64 contLambdaOffset)
{
    //
    Log.LogInfo("  Operator-Linemodel: Initializing contaminant for roll #%d, "
                "with offset=%.2f",
                iRollContaminated, contLambdaOffset);

    m_iRollContaminated = iRollContaminated;
    m_contLambdaOffset = contLambdaOffset;
    const std::string &category = "emission";
    m_tplContaminant =
        std::shared_ptr<CTemplate>(new CTemplate("contaminant", category));
    Int32 length = contModelSpectrum->GetSpectrum().GetSampleCount();
    CSpectrumAxis &contModelFluxAxis =
        contModelSpectrum->GetSpectrum().GetFluxAxis();
    CSpectrumAxis &contModelSpectralAxis =
        contModelSpectrum->GetSpectrum().GetSpectralAxis();

    CSpectrumAxis &spcFluxAxis = m_tplContaminant->GetFluxAxis();
    spcFluxAxis.SetSize(length);
    CSpectrumAxis &spcSpectralAxis = m_tplContaminant->GetSpectralAxis();
    spcSpectralAxis.SetSize(length);

    for (Int32 k = 0; k < length; k++)
    {
        spcFluxAxis[k] = contModelFluxAxis[k];
        spcSpectralAxis[k] = contModelSpectralAxis[k];
    }

    // applying offset for ra/dec distance between main source and contaminant
    m_tplContaminant->GetSpectralAxis().ApplyOffset(m_contLambdaOffset);

    m_enableLoadContTemplate = true;
    return 0;
}

///
/// \brief COperatorLineModel::storeGlobalModelResults
/// stores the linemodel results as global results in the datastore
///
void COperatorLineModel::storeGlobalModelResults(CDataStore &dataStore)
{
    Int32 nResults = m_savedModelSpectrumResults.size();
    if (nResults > m_savedModelFittingResults.size())
    {
        Log.LogError("Line Model, not as many model fitting results as model "
                     "spectrum results, (nspc = %d, nfit = %d)",
                     m_savedModelSpectrumResults.size(),
                     m_savedModelFittingResults.size());
        nResults = m_savedModelFittingResults.size();
    }

    for (Int32 k = 0; k < nResults; k++)
    {
        std::string fname_spc =
            (boost::format("linemodel_spc_extrema_%1%") % k).str();
        dataStore.StoreScopedGlobalResult(fname_spc.c_str(),
                                          m_savedModelSpectrumResults[k]);

        std::string fname_fit =
            (boost::format("linemodel_fit_extrema_%1%") % k).str();
        dataStore.StoreScopedGlobalResult(fname_fit.c_str(),
                                          m_savedModelFittingResults[k]);

        std::string fname_fitcontinuum =
            (boost::format("linemodel_fitcontinuum_extrema_%1%") % k).str();
        dataStore.StoreScopedGlobalResult(
            fname_fitcontinuum.c_str(), m_savedModelContinuumFittingResults[k]);

        std::string fname_rules =
            (boost::format("linemodel_rules_extrema_%1%") % k).str();
        dataStore.StoreScopedGlobalResult(fname_rules.c_str(),
                                          m_savedModelRulesResults[k]);
    }

    for (Int32 k = 0; k < m_savedModelContinuumSpectrumResults.size(); k++)
    {
        std::string nameBaselineStr =
            (boost::format("linemodel_continuum_extrema_%1%") % k).str();
        dataStore.StoreScopedGlobalResult(
            nameBaselineStr.c_str(), m_savedModelContinuumSpectrumResults[k]);
    }
}

///
/// \brief COperatorLineModel::storePerTemplateModelResults
/// stores the linemodel results as per template results in the datastore
///
void COperatorLineModel::storePerTemplateModelResults(CDataStore &dataStore,
                                                      const CTemplate &tpl)
{
    Int32 nResults = m_savedModelSpectrumResults.size();
    if (nResults > m_savedModelFittingResults.size())
    {
        Log.LogError("Line Model, not as many model fitting results as model "
                     "spectrum results, (nspc = %d, nfit = %d)",
                     m_savedModelSpectrumResults.size(),
                     m_savedModelFittingResults.size());
        nResults = m_savedModelFittingResults.size();
    }

    for (Int32 k = 0; k < nResults; k++)
    {
        std::string fname_spc =
            (boost::format("linemodel_spc_extrema_%1%") % k).str();
        dataStore.StoreScopedPerTemplateResult(tpl, fname_spc.c_str(),
                                               m_savedModelSpectrumResults[k]);

        std::string fname_fit =
            (boost::format("linemodel_fit_extrema_%1%") % k).str();
        dataStore.StoreScopedPerTemplateResult(tpl, fname_fit.c_str(),
                                               m_savedModelFittingResults[k]);

        std::string fname_fitcontinuum =
            (boost::format("linemodel_fitcontinuum_extrema_%1%") % k).str();
        dataStore.StoreScopedPerTemplateResult(
            tpl, fname_fitcontinuum.c_str(),
            m_savedModelContinuumFittingResults[k]);

        std::string fname_rules =
            (boost::format("linemodel_rules_extrema_%1%") % k).str();
        dataStore.StoreScopedPerTemplateResult(tpl, fname_rules.c_str(),
                                               m_savedModelRulesResults[k]);
    }
}

std::shared_ptr<CModelSpectrumResult>
COperatorLineModel::GetModelSpectrumResult(Int32 idx)
{
    Int32 nResults = m_savedModelSpectrumResults.size();
    if (idx >= nResults)
    {
        return NULL;
    } else
    {
        return m_savedModelSpectrumResults[idx];
    }
}

std::shared_ptr<CSpectraFluxResult>
COperatorLineModel::GetModelSpectrumContinuumResult(Int32 idx)
{
    Int32 nResults = m_savedModelContinuumSpectrumResults.size();
    if (idx >= nResults)
    {
        return NULL;
    } else
    {
        return m_savedModelContinuumSpectrumResults[idx];
    }
}

std::shared_ptr<CModelSpectrumResult>
COperatorLineModel::GetContaminantSpectrumResult()
{
    return m_savedContaminantSpectrumResult;
}

Int32 COperatorLineModel::interpolateLargeGridOnFineGrid(
        TFloat64List redshiftsLargeGrid,
        TFloat64List redshiftsFineGrid,
        TFloat64List meritLargeGrid,
        TFloat64List &meritFineGrid)
{
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

    for (Int32 j = 0; j < redshiftsFineGrid.size(); j++)
    {
        Float64 Xrebin = redshiftsFineGrid[j];
        if (Xrebin < redshiftsLargeGrid[0] ||
            Xrebin > redshiftsLargeGrid[redshiftsLargeGrid.size() - 1])
        {
            continue;
        }
        meritFineGrid[j] = gsl_interp_eval(
            interpolation, &redshiftsLargeGrid.front(), &meritLargeGrid.front(),
            Xrebin, accelerator); // lin
    }

    gsl_interp_free(interpolation);
    gsl_interp_accel_free(accelerator);
    //*/

    return 0;
}

void COperatorLineModel::ComputeArea1(CLineModelResult &results)
{
    // prepare p
    Float64 maxp = DBL_MIN;
    CSpectrum pspc;
    CSpectrumFluxAxis &spcFluxAxis = pspc.GetFluxAxis();
    spcFluxAxis.SetSize(results.Redshifts.size());
    CSpectrumSpectralAxis &spcSpectralAxis = pspc.GetSpectralAxis();
    spcSpectralAxis.SetSize(results.Redshifts.size());
    for (Int32 i2 = 0; i2 < results.Redshifts.size(); i2++)
    {
        if (maxp < results.ChiSquare[i2])
        {
            maxp = results.ChiSquare[i2];
        }
    }
    for (Int32 i2 = 0; i2 < results.Redshifts.size(); i2++)
    {
        spcFluxAxis[i2] = expm1(-(results.ChiSquare[i2] - maxp) / 2.0);
        spcSpectralAxis[i2] = results.Redshifts[i2];
    }

    /*//debug:
    FILE* f = fopen( "getbestredshiftbayes_dbg.txt", "w+" );
    for( Int32 i=0; i<spcFluxAxis.GetSamplesCount(); i++ )
    {
        if( spcFluxAxis[i] < 0.0001 ){
            fprintf( f, "%e %e\n", spcSpectralAxis[i], spcFluxAxis[i]);
        }else{
            fprintf( f, "%f %f\n", spcSpectralAxis[i], spcFluxAxis[i]);
        }
    }
    fclose( f );
    //*/

    Float64 winsize = 0.0025;
    Float64 inclusionThresRatio = 0.25;
    Int32 iz0 = 0;
    for (Int32 i = 0; i < results.ExtremaResult.Extrema.size(); i++)
    {
        // find iz, izmin and izmax
        Int32 izmin = -1;
        Int32 iz = -1;
        Int32 izmax = -1;
        for (Int32 i2 = 0; i2 < results.Redshifts.size(); i2++)
        {
            if (iz == -1 &&
                (results.ExtremaResult.Extrema[i]) <= results.Redshifts[i2])
            {
                iz = i2;
                if (i == 0)
                {
                    iz0 = iz;
                }
            }
            if (izmin == -1 && (results.ExtremaResult.Extrema[i] -
                                winsize / 2.0) <= results.Redshifts[i2])
            {
                izmin = i2;
            }
            if (izmax == -1 && (results.ExtremaResult.Extrema[i] +
                                winsize / 2.0) <= results.Redshifts[i2])
            {
                izmax = i2;
                break;
            }
        }
        Float64 di = abs(results.ChiSquare[iz] - maxp);
        Float64 d0 = abs(results.ChiSquare[iz0] - maxp);
        if (di < inclusionThresRatio * d0)
        {
            continue;
        }

        /*
        CGaussianFitSimple fitter;
        CGaussianFitSimple::EStatus status = fitter.Compute( pspc, TInt32Range(
        izmin, izmax ) );
        if(status!=NSEpic::CGaussianFitSimple::nStatus_Success){
            continue;
        }

        Float64 gaussAmp;
        Float64 gaussPos;
        Float64 gaussWidth;
        fitter.GetResults( gaussAmp, gaussPos, gaussWidth );
        Float64 gaussAmpErr;
        Float64 gaussPosErr;
        Float64 gaussWidthErr;
        fitter.GetResultsError( gaussAmpErr, gaussPosErr, gaussWidthErr );
        */
        Float64 gaussWidth =
            FitBayesWidth(spcSpectralAxis, spcFluxAxis,
                          results.ExtremaResult.Extrema[i], izmin, izmax);
        Float64 gaussAmp = spcFluxAxis[iz];

        Float64 area = 0.0;
        for (Int32 i2 = izmin; i2 < izmax; i2++)
        {
            Float64 x = spcSpectralAxis[i2];
            Float64 Yi =
                gaussAmp * exp(-1. * (x - results.ExtremaResult.Extrema[i]) *
                               (x - results.ExtremaResult.Extrema[i]) /
                               (2 * gaussWidth * gaussWidth));
            area += Yi;
        }

        // Float64 area = gaussAmp*gaussWidth*sqrt(2.0*3.141592654);
        results.ExtremaResult.LogArea[i] = area;
    }
}

///
/// \brief COperatorLineModel::ComputeArea2
/// computes the Laplace approx for a given Chi2 result around the N best
/// extrema
///
void COperatorLineModel::ComputeArea2(CLineModelResult &results)
{
    Float64 maxp = DBL_MIN;
    for (Int32 i2 = 0; i2 < results.Redshifts.size(); i2++)
    {
        if (maxp < results.ChiSquare[i2])
        {
            maxp = results.ChiSquare[i2];
        }
    }
    Float64 winsize = 0.001;
    Float64 inclusionThresRatio = 0.01;
    Int32 iz0 = 0;
    for (Int32 indz = 0; indz < results.ExtremaResult.Extrema.size(); indz++)
    {
        // find iz, izmin and izmax
        Int32 izmin = -1;
        Int32 iz = -1;
        Int32 izmax = -1;
        for (Int32 i2 = 0; i2 < results.Redshifts.size(); i2++)
        {
            if (iz == -1 &&
                (results.ExtremaResult.Extrema[indz]) <= results.Redshifts[i2])
            {
                iz = i2;
                if (indz == 0)
                {
                    iz0 = iz;
                }
            }
            if (izmin == -1 && (results.ExtremaResult.Extrema[indz] -
                                winsize / 2.0) <= results.Redshifts[i2])
            {
                izmin = i2;
            }
            if (izmax == -1 && (results.ExtremaResult.Extrema[indz] +
                                winsize / 2.0) <= results.Redshifts[i2])
            {
                izmax = i2;
                break;
            }
        }
        Float64 di = abs(results.ChiSquare[iz] - maxp);
        Float64 d0 = abs(results.ChiSquare[iz0] - maxp);
        if (di < inclusionThresRatio * d0)
        {
            continue;
        }
        if (izmin == -1 || izmax == -1)
        {
            continue;
        }

        // quadratic fit
        int i, n;
        double chisq;
        gsl_matrix *X, *cov;
        gsl_vector *y, *w, *c;

        n = izmax - izmin + 1;

        X = gsl_matrix_alloc(n, 3);
        y = gsl_vector_alloc(n);
        w = gsl_vector_alloc(n);

        c = gsl_vector_alloc(3);
        cov = gsl_matrix_alloc(3, 3);

        double x0 = results.ExtremaResult.Extrema[indz];
        for (i = 0; i < n; i++)
        {
            double xi, yi, ei;
            xi = results.Redshifts[i + izmin];
            yi = results.ChiSquare[i + izmin];
            ei = 1.0; // todo, estimate weighting ?
            gsl_matrix_set(X, i, 0, 1.0);
            gsl_matrix_set(X, i, 1, xi - x0);
            gsl_matrix_set(X, i, 2, (xi - x0) * (xi - x0));

            gsl_vector_set(y, i, yi);
            gsl_vector_set(w, i, 1.0 / (ei * ei));
        }

        {
            gsl_multifit_linear_workspace *work =
                gsl_multifit_linear_alloc(n, 3);
            gsl_multifit_wlinear(X, w, y, c, cov, &chisq, work);
            gsl_multifit_linear_free(work);
        }

#define C(i) (gsl_vector_get(c, (i)))
#define COV(i, j) (gsl_matrix_get(cov, (i), (j)))

        double zcorr = x0 - C(1) / (2.0 * C(2));
        double sigma = sqrt(1.0 / C(2));
        Float64 a = (Float64)(C(0));
        Float64 b2sur4c = (Float64)(C(1) * C(1) / ((Float64)(4.0 * C(2))));
        Float64 logK = (-(a - b2sur4c) / 2.0);
        Float64 logarea = log(sigma) + logK + log(2.0 * M_PI);
        if (0)
        {
            Log.LogInfo("Extrema: %g", results.ExtremaResult.Extrema[indz]);
            Log.LogInfo("# best fit: Y = %g + %g X + %g X^2", C(0), C(1), C(2));
            if (false) // debug
            {
                Log.LogInfo("# covariance matrix:\n");
                Log.LogInfo("[ %+.5e, %+.5e, %+.5e  \n", COV(0, 0), COV(0, 1),
                            COV(0, 2));
                Log.LogInfo("  %+.5e, %+.5e, %+.5e  \n", COV(1, 0), COV(1, 1),
                            COV(1, 2));
                Log.LogInfo("  %+.5e, %+.5e, %+.5e ]\n", COV(2, 0), COV(2, 1),
                            COV(2, 2));
            }
            Log.LogInfo("# chisq/n = %g", chisq / n);
            Log.LogInfo("# zcorr = %g", zcorr);
            Log.LogInfo("# sigma = %g", sigma);
            Log.LogInfo("# logarea = %g", logarea);
            Log.LogInfo("\n");
        }

        gsl_matrix_free(X);
        gsl_vector_free(y);
        gsl_vector_free(w);
        gsl_vector_free(c);
        gsl_matrix_free(cov);

        results.ExtremaResult.LogArea[indz] = logarea;
        results.ExtremaResult.SigmaZ[indz] = sigma;
        results.ExtremaResult.LogAreaCorrectedExtrema[indz] = zcorr;
    }
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
Float64 COperatorLineModel::FitBayesWidth(CSpectrumSpectralAxis &spectralAxis,
                                          CSpectrumFluxAxis &fluxAxis,
                                          Float64 z, Int32 start, Int32 end)
{
    Float64 A = boost::numeric::bounds<float>::lowest();
    const Float64 *flux = fluxAxis.GetSamples();
    const Float64 *spectral = spectralAxis.GetSamples();
    // const Float64* error = fluxAxis.GetError();

    // A = max, good value ?
    for (Int32 i = start; i < end; i++)
    {
        Float64 y = flux[i];
        if (y > A)
        {
            A = y;
        }
    }

    if (A <= 0)
    {
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
    while (icmpt < maxIteration)
    {
        sum2 = 0.0;
        for (Int32 i = start; i < end; i++)
        {
            Float64 x = spectral[i];
            Float64 Yi = A * exp(-1. * (x - mu) * (x - mu) / (2 * c * c));
            sum2 += pow(Yi - flux[i], 2.0);
            // sum2 += pow( Yi - flux[i] , 2.0 ) / pow( error[i], 2.0 );
        }
        if (sum2 < minsum2)
        {
            minc = c;
            minsum2 = sum2;
        }
        icmpt++;
        c = c + cstepup;
    }

    if (minc < 0)
    {
        minc = 0;
    }
    return minc;
}
