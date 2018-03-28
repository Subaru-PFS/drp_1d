#include <RedshiftLibrary/operator/linemodel.h>
#include <RedshiftLibrary/linemodel/templatesortho.h>
#include <RedshiftLibrary/linemodel/templatesorthostore.h>
#include <RedshiftLibrary/spectrum/axis.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/tools.h>
#include <RedshiftLibrary/common/mask.h>
#include <RedshiftLibrary/operator/chisquare2.h>
#include <RedshiftLibrary/operator/chisquareloglambda.h>
#include <RedshiftLibrary/operator/chisquareresult.h>
#include <RedshiftLibrary/operator/spectraFluxResult.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/statistics/deltaz.h>
#include <RedshiftLibrary/linemodel/templatesfitstore.h>

#include <RedshiftLibrary/spectrum/io/fitswriter.h>
#include <RedshiftLibrary/log/log.h>

#include <boost/numeric/conversion/bounds.hpp>
#include "boost/format.hpp"
#include <boost/chrono/thread_clock.hpp>

#include <RedshiftLibrary/processflow/datastore.h>
#include <RedshiftLibrary/spectrum/io/genericreader.h>

#include <math.h>
#include <float.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <algorithm>    // std::sort
#include <gsl/gsl_multifit.h>
#include <assert.h>

#define NOT_OVERLAP_VALUE NAN
#include <stdio.h>

using namespace NSEpic;
using namespace std;

/**
 **/
COperatorLineModel::COperatorLineModel()
{
    m_maxModelSaveCount = 10;
}

/**
 * \brief Empty destructor.
 **/
COperatorLineModel::~COperatorLineModel()
{

}

/**
 * \brief Call the Linemodel to fit models and select the parts of the results that are relevant.
 * Print an error if the spectral axis is not in linear scale.
 * Sort the list of redshifts.
 * Get the list of rest lines filtered by type and force.
 * Create a CLineModelResult object, resize its members and populate it with the filtered list.
 * Create an ElementList object.
 * For each sorted redshift value, call ModelFit.
 * Find the extrema in the set of per-redshift results.
 * Refine Extremum with a second maximum search around the z candidates.
 * Filter out peaks outside intervals centered on the refined extrema.
 * Save the results.
 **/
std::shared_ptr<COperatorResult> COperatorLineModel::Compute(CDataStore &dataStore,
                                  const CSpectrum& spectrum,
                                  const CSpectrum& spectrumContinuum,
                                  const CTemplateCatalog& tplCatalog,
                                  const TStringList& tplCategoryList,
                                  const std::string opt_calibrationPath,
                                  const CRayCatalog& restraycatalog,
                                  const std::string& opt_lineTypeFilter,
                                  const std::string& opt_lineForceFilter,
                                  const TFloat64Range& lambdaRange,
                                  const TFloat64List& redshifts,
                                  const Int32 opt_extremacount,
                                  const std::string& opt_fittingmethod,
                                  const std::string& opt_continuumcomponent,
                                  const std::string& opt_lineWidthType,
                                  const Float64 opt_resolution,
                                  const Float64 opt_velocityEmission,
                                  const Float64 opt_velocityAbsorption,
                                  const std::string& opt_continuumreest,
                                  const std::string& opt_rules,
                                  const std::string& opt_velocityFitting,
                                  const Float64 &opt_twosteplargegridstep,
                                  const std::string& opt_rigidity,
                                 const Float64 &opt_emvelocityfitmin,
                                 const Float64 &opt_emvelocityfitmax,
                                 const Float64 &opt_absvelocityfitmin,
                                 const Float64 &opt_absvelocityfitmax)
{
    // initialize empty results so that it can be returned anyway in case of an error
    std::shared_ptr<CLineModelResult> result = std::shared_ptr<CLineModelResult>( new CLineModelResult() );

    if( spectrum.GetSpectralAxis().IsInLinearScale()==false )
    {
        Log.LogError( "Line Model, input spectrum is not in linear scale (ignored)." );
    }

    TFloat64List sortedRedshifts = redshifts;
    std::sort(sortedRedshifts.begin(), sortedRedshifts.end());

    // redefine redshift grid
    Int32 enableFastFitLargeGrid = 0;
    if(opt_twosteplargegridstep > 0.0)
    {
        enableFastFitLargeGrid = 1;
    }

    TFloat64List largeGridRedshifts;
    //*
    if(enableFastFitLargeGrid==1){
        //Log.LogInfo("Line Model, Fast Fit Large Grid enabled");
        //calculate on a wider grid, defined by a minimum step
        Float64 dz_thres = opt_twosteplargegridstep;
        std::vector<Int32> removed_inds;
        Int32 lastKeptInd = 0;
        for(Int32 i=1; i<sortedRedshifts.size()-1; i++){
            if( abs(sortedRedshifts[i]-sortedRedshifts[lastKeptInd])<dz_thres && abs(sortedRedshifts[i+1]-sortedRedshifts[i])<dz_thres ){
                removed_inds.push_back(i);
            }else{
                lastKeptInd=i;
            }
        }
        Int32 rmInd = 0;
        for(Int32 i=1; i<sortedRedshifts.size()-1; i++){
            bool addToLargeGrid = true;
            if(removed_inds.size()>0)
            {
                if(removed_inds[rmInd]==i){
                    rmInd++;
                    addToLargeGrid = false;
                }
            }
            if(addToLargeGrid){
                largeGridRedshifts.push_back(sortedRedshifts[i]);
            }
        }
        if(largeGridRedshifts.size()<1)
        {
            enableFastFitLargeGrid = 0;
            Log.LogInfo( "  Operator-Linemodel: FastFitLargeGrid auto disabled: raw %d redshifts will be calculated", sortedRedshifts.size());
        }else{
            Log.LogInfo( "  Operator-Linemodel: FastFitLargeGrid enabled: %d redshifts will be calculated on the large grid (%d initially)", largeGridRedshifts.size(), sortedRedshifts.size());
        }
    }

    Int32 typeFilter = -1;
    if(opt_lineTypeFilter == "A"){
        typeFilter = CRay::nType_Absorption;
    }else if(opt_lineTypeFilter == "E"){
        typeFilter = CRay::nType_Emission;
    }

    Int32 forceFilter = -1;//CRay::nForce_Strong;
    if(opt_lineForceFilter == "S"){
        forceFilter = CRay::nForce_Strong;
    }
    CRayCatalog::TRayVector restRayList = restraycatalog.GetFilteredList( typeFilter, forceFilter);
    Log.LogDebug( "restRayList.size() = %d", restRayList.size() );

    bool enableOrtho = true;
    Log.LogInfo( "  Operator-Linemodel: TemplatesOrthogonalization enabled = %d", enableOrtho );

    //prepare continuum templates catalog
    CTemplatesOrthogonalization tplOrtho(tplCatalog,
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
    //CTemplateCatalog orthoTplCatalog = tplOrtho.getOrthogonalTplCatalog();
    CTemplatesOrthoStore orthoTplStore = tplOrtho.getOrthogonalTplStore();
    Int32 ctlgIdx = 0; //only one ortho config for now
    std::shared_ptr<CTemplateCatalog> orthoTplCatalog = orthoTplStore.getTplCatalog(ctlgIdx);
    Log.LogInfo( "  Operator-Linemodel: Templates store prepared.");

    //*
    CLineModelElementList model( spectrum,
                                 spectrumContinuum,
                                 tplCatalog,//*orthoTplCatalog,//
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
    Float64 setssSizeInit = 0.1;
    model.SetSourcesizeDispersion(setssSizeInit);
    Log.LogInfo( "  Operator-Linemodel: sourcesize init to: ss=%.1f", setssSizeInit);
    //*/


    /*
    CMultiRollModel model( spectrum,
                                 spectrumContinuum,
                                 tplCatalog,//*orthoTplCatalog,//
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


    Bool enableLoadContTemplateOverride = false; //manual switch for hardcoded contaminant bypass load
    if(enableLoadContTemplateOverride)
    {
        m_enableLoadContTemplate = false; //disable stock contaminant in case of overrided contaminant
        Int32 iContTemplate = 0; //idx of the roll being contaminated
        std::shared_ptr<CTemplate> tplContaminant;

        Log.LogInfo( "  Operator-Linemodel: OVERRIDDEN loading contaminant for roll #%d", iContTemplate );

        //hardcoded load from file on disk: of the contaminant for first model
        std::string templatePath = "/home/aschmitt/data/euclid/simulation2017-SC3_test_zweiroll/amazed/output_rolls_source3/euc_testsc3_zweiroll_source3_roll0_F_i0/linemodelsolve.linemodel_spc_extrema_0.txt";
        const std::string& category = "emission";
        tplContaminant = std::shared_ptr<CTemplate>( new CTemplate( "contaminant", category ) );
        CSpectrumIOGenericReader asciiReader;
        if( !asciiReader.Read( templatePath.c_str(), *tplContaminant ) ) {
            Log.LogError( "Fail to read contaminant template: %s", templatePath.c_str() );
            return result;
        }else{
            Log.LogInfo( "Successfully loaded contaminant template: %s", templatePath.c_str() );
        }
        //debug:
        //FILE* f = fopen( "contaminantLoaded.txt", "w+" );
        //for(Int32 k=0; k<tplContaminant->GetSampleCount(); k++)
        //{
        //    fprintf( f, "%f\t%e\n", tplContaminant->GetSpectralAxis()[k], tplContaminant->GetFluxAxis()[k]);
        //}
        //fclose( f );
        //

        Float64 lambdaOffset_forSourceSpatialOffsetInDispersionDirection=2000; //hardcoded
        tplContaminant->GetSpectralAxis().ApplyOffset(lambdaOffset_forSourceSpatialOffsetInDispersionDirection);
        //debug:
        FILE* f2 = fopen( "contaminantShifted.txt", "w+" );
        for(Int32 k=0; k<tplContaminant->GetSampleCount(); k++)
        {
            fprintf( f2, "%f\t%e\n", tplContaminant->GetSpectralAxis()[k], tplContaminant->GetFluxAxis()[k]);
        }
        fclose( f2 );

        //tplContaminant
        model.LoadFitContaminantTemplate(iContTemplate, *tplContaminant, lambdaRange);
    }
    if(m_enableLoadContTemplate)
    {
        Log.LogInfo( "  Operator-Linemodel: loading contaminant for roll #%d", m_iRollContaminated );
        if(m_tplContaminant==0)
        {
            Log.LogError( "  Operator-Linemodel: Contaminant data is invalid... aborting." );
        }else{
            //
            if(0)
            {
            //debug:
            FILE* f2 = fopen( "contaminantShifted.txt", "w+" );
            for(Int32 k=0; k<m_tplContaminant->GetSampleCount(); k++)
            {
                fprintf( f2, "%f\t%e\n", m_tplContaminant->GetSpectralAxis()[k], m_tplContaminant->GetFluxAxis()[k]);
            }
            fclose( f2 );
            //
            }

            //apply contamination to the multiroll model
            model.LoadFitContaminantTemplate(m_iRollContaminated, *m_tplContaminant, lambdaRange);
            m_savedContaminantSpectrumResult = model.GetContaminantSpectrumResult(m_iRollContaminated);
        }
    }
    //*/


    //init catalog tplratios
    Log.LogInfo( "  Operator-Linemodel: Tpl-ratios init");
    bool tplratioInitRet = model.initTplratioCatalogs();
    if(!tplratioInitRet)
    {
        Log.LogError( "  Operator-Linemodel: Failed to init tpl-ratios. aborting...");
        return result;
    }

    //init catalog offsets
    Log.LogInfo( "  Operator-Linemodel: Lambda offsets init");
    bool offsetsInitRet = model.initLambdaOffsets();
    if(!offsetsInitRet)
    {
        Log.LogError( "  Operator-Linemodel: Failed to init lambda offsets. Continuing without offsets...");
    }

    Int32 resultInitRet = result->Init( sortedRedshifts, restRayList, model.getTplshape_count(), model.getTplshape_priors());
    if(resultInitRet!=0)
    {
        Log.LogError( "  Operator-Linemodel: ERROR while initializing linemodel result (ret=%d)", resultInitRet);
        return result;
    }


    Log.LogInfo( "  Operator-Linemodel: initialized");

    //setup velocity fitting
    bool enableVelocityFitting = true;
    Float64 velfitMinE = opt_emvelocityfitmin;
    Float64 velfitMaxE = opt_emvelocityfitmax;
    Float64 velfitMinA = opt_absvelocityfitmin;
    Float64 velfitMaxA = opt_absvelocityfitmax;
    //HARDCODED - override: no-velocityfitting for abs
    //velfitMinA = opt_velocityAbsorption;
    //velfitMaxA = opt_velocityAbsorption;
    if(opt_velocityFitting != "yes"){
        enableVelocityFitting = false;
    }else{        
        Log.LogInfo( "  Operator-Linemodel: velocity fitting bounds for Emission: min=%.1f - max=%.1f", velfitMinE, velfitMaxE);
        Log.LogInfo( "  Operator-Linemodel: velocity fitting bounds for Absorption: min=%.1f - max=%.1f", velfitMinA, velfitMaxA);
    }

    //enable/disable fit by groups. Once enabled, the velocity fitting groups are defined in the line catalog from v4.0 on.
    m_enableWidthFitByGroups = true;

    //fit continuum
    bool enableFitContinuumPrecomputed = true;
    if(enableFitContinuumPrecomputed && opt_continuumcomponent == "tplfit")
    {
        boost::chrono::thread_clock::time_point start_tplfitprecompute = boost::chrono::thread_clock::now();
        Float64 redshiftStep = 1.5e-4;
        Float64 minRedshift = sortedRedshifts[0];
        Float64 maxRedshift = sortedRedshifts[sortedRedshifts.size()-1];
        Log.LogInfo( "  Operator-Linemodel: continuum tpl fitting: step=%.5f, min=%.1f, max=%.1f", redshiftStep, minRedshift, maxRedshift);

        CTemplatesFitStore* tplfitStore = new CTemplatesFitStore(minRedshift, maxRedshift, redshiftStep);
        std::vector<Float64> redshiftsTplFit = tplfitStore->GetRedshiftList();
        Log.LogInfo( "  Operator-Linemodel: continuum tpl redshift list n=%d", redshiftsTplFit.size());
        std::vector<std::shared_ptr<CChisquareResult>> chisquareResultsAllTpl;
        std::vector<std::string> chisquareResultsTplName;

        std::string opt_chi2operator = "chisquarelog"; //"chisquare2"; //
        if(redshiftsTplFit.size()<100 && opt_chi2operator != "chisquare2") // warning arbitrary number of redshifts threshold to consider chisquare2 faster than chisquarelog
        {
            opt_chi2operator = "chisquare2";
            Log.LogInfo( "  Operator-Linemodel: precomputing- auto select chisquare2 operator (faster for few redshifts calc. points)" );
        }
        std::string opt_interp = "precomputedfinegrid"; // "lin"; //
        Int32 opt_dustFit = 1;
        Int32 opt_extinction = 1;
        Log.LogInfo( "  Operator-Linemodel: precomputing- with operator = %s", opt_chi2operator.c_str() );
        Log.LogInfo( "  Operator-Linemodel: precomputing-fitContinuum_dustfit = %d", opt_dustFit );
        Log.LogInfo( "  Operator-Linemodel: precomputing-fitContinuum_igm = %d", opt_extinction );

        std::shared_ptr<COperator> chiSquareOperator;
        if(opt_chi2operator=="chisquarelog")
        {
            //COperatorChiSquareLogLambda* chiSquareOperator;
            bool enableLogRebin = true;
            chiSquareOperator = std::shared_ptr<COperatorChiSquareLogLambda>( new COperatorChiSquareLogLambda(opt_calibrationPath));
            std::shared_ptr<COperatorChiSquareLogLambda> chiSquareLogOperator = std::dynamic_pointer_cast<COperatorChiSquareLogLambda>( chiSquareOperator );
            chiSquareLogOperator->enableSpcLogRebin(enableLogRebin);
        }
        else if(opt_chi2operator=="chisquare2"){
            chiSquareOperator = std::shared_ptr<COperatorChiSquare2>( new COperatorChiSquare2(opt_calibrationPath));
        }else{
            Log.LogError( "  Operator-Linemodel: unable to parce chisquare continuum fit operator");
        }


        Float64 overlapThreshold = 1.0;
        bool ignoreLinesSupport=0;
        std::vector<CMask> maskList; //can't get the lines support BEFORE initializing the model

        for( UInt32 i=0; i<tplCategoryList.size(); i++ )
        {
            std::string category = tplCategoryList[i];
            for( UInt32 j=0; j<orthoTplCatalog->GetTemplateCount( category ); j++ )
            {
                const CTemplate& tpl = orthoTplCatalog->GetTemplate( category, j );

                auto  chisquareResult = std::dynamic_pointer_cast<CChisquareResult>( chiSquareOperator->Compute( spectrum, tpl, lambdaRange, redshiftsTplFit, overlapThreshold, maskList, opt_interp, opt_extinction, opt_dustFit ) );
                if( !chisquareResult )
                {
                    Log.LogInfo( "  Operator-Linemodel failed to compute chi square value for tpl=%s", tpl.GetName().c_str());
                }else{
                    chisquareResultsAllTpl.push_back(chisquareResult);
                    chisquareResultsTplName.push_back(tpl.GetName());
                }
            }
        }

        chiSquareOperator.reset();

        //fill the results with Best Values
        Int32 nTplFitResults = redshiftsTplFit.size();
        for (Int32 i=0;i<nTplFitResults;i++)
        {
            Float64 redshift = redshiftsTplFit[i];
            Float64 bestMerit = DBL_MAX;
            Float64 bestFitAmplitude;
            Float64 bestFitDustCoeff;
            Int32   bestFitMeiksinIdx;
            Float64 bestFitDtM;
            Float64 bestFitMtM;
            std::string bestTplName;

            for( UInt32 j=0; j<chisquareResultsAllTpl.size(); j++ )
            {
                auto chisquareResult = std::dynamic_pointer_cast<CChisquareResult>( chisquareResultsAllTpl[j] );
                Float64 merit = chisquareResult->ChiSquare[i];


                if(merit<bestMerit)
                {
                    bestMerit = merit;
                    bestFitAmplitude = chisquareResult->FitAmplitude[i];
                    bestFitDustCoeff = chisquareResult->FitDustCoeff[i];
                    bestFitMeiksinIdx = chisquareResult->FitMeiksinIdx[i];
                    bestFitDtM = chisquareResult->FitDtM[i];
                    bestFitMtM = chisquareResult->FitMtM[i];
                    bestTplName = chisquareResultsTplName[j];
                }
            }
            tplfitStore->Add(redshift, bestMerit, bestFitAmplitude, bestFitDustCoeff, bestFitMeiksinIdx, bestFitDtM, bestFitMtM, bestTplName);
        }

        //Set tplFitStore if needed
        model.SetFitContinuum_FitStore(tplfitStore);

        boost::chrono::thread_clock::time_point stop_tplfitprecompute = boost::chrono::thread_clock::now();
        Float64 duration_tplfitprecompute = boost::chrono::duration_cast<boost::chrono::microseconds>(stop_tplfitprecompute - start_tplfitprecompute).count();
        Float64 duration_tplfit_seconds = duration_tplfitprecompute/1e6;
        Log.LogInfo( "  Operator-Linemodel: tplfit-precompute done in %.4e sec", duration_tplfit_seconds);
        Log.LogInfo("<proc-lm-tplfit><%d>", (Int32)duration_tplfit_seconds);

    }

//    //hack, zero outside of the support  ///////////////////////////////////////////////////////////////////////////////////////////
//    model.setModelSpcObservedOnSupportZeroOutside(lambdaRange);


//    std::shared_ptr<CSpectraFluxResult> baselineResult = (std::shared_ptr<CSpectraFluxResult>) new CSpectraFluxResult();
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

//    std::string nameBaselineStr = (boost::format("linemodel_template_zeroedoutsidelines")).str();
//    dataStore.StoreScopedGlobalResult(nameBaselineStr.c_str(), baselineResult);
//    return NULL;
//    // end of hack //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    result->nSpcSamples = model.getSpcNSamples(lambdaRange);
    result->dTransposeDNocontinuum = model.getDTransposeD(lambdaRange, "nocontinuum");
    result->dTransposeD = model.getDTransposeD(lambdaRange, "raw");
    result->cstLog = model.getLikelihood_cstLog(lambdaRange);

    Int32 contreest_iterations = 0;
    if( opt_continuumreest == "always" )
    {
        contreest_iterations = 1;
    }
    else
    {
        contreest_iterations = 0;
    }

    //Set model parameter: abs lines limit
    Float64 absLinesLimit = 1.0; //-1 to disable, 1.0 is typical
    model.SetAbsLinesLimit(absLinesLimit);
    Log.LogInfo( "  Operator-Linemodel: set abs lines limit to %f (ex: -1 means disabled)", absLinesLimit);

    //Set model parameter: continuum least-square estimation fast
    //note: this fast method requires continuum templates and linemodels to be orthogonal. The velfit option turns this trickier...
    Int32 estimateLeastSquareFast = 0;
    model.SetLeastSquareFastEstimationEnabled(estimateLeastSquareFast);
    Log.LogInfo( "  Operator-Linemodel: set estimateLeastSquareFast to %d (ex: 0 means disabled)", estimateLeastSquareFast);

    Log.LogInfo( "  Operator-Linemodel: start processing");

    //Set model parameters to FIRST-PASS
    model.setPassMode(1);
    Log.LogInfo( "  Operator-Linemodel: ---------- ---------- ---------- ----------" );
    Log.LogInfo( "  Operator-Linemodel: now computing first-pass");
    Log.LogInfo( "  Operator-Linemodel: ---------- ---------- ---------- ----------" );

    //WARNING: HACK, first pass with continuum from spectrum.
    //model.SetContinuumComponent("fromspectrum");
    //
    Int32 indexLargeGrid = 0;
    std::vector<Float64> calculatedLargeGridRedshifts;
    std::vector<Float64> calculatedLargeGridMerits;
    std::vector<TFloat64List> calculatedChiSquareTplshapes;
    for(Int32 k=0; k<result->ChiSquareTplshapes.size(); k++)
    {
        TFloat64List _chi2tpl;
        calculatedChiSquareTplshapes.push_back(_chi2tpl);
    }
    boost::chrono::thread_clock::time_point start_mainloop = boost::chrono::thread_clock::now();

    //#pragma omp parallel for
    for (Int32 i=0;i<result->Redshifts.size();i++)
    {
        if(enableFastFitLargeGrid==0 || i==0 || result->Redshifts[i] == largeGridRedshifts[indexLargeGrid])
        {
            result->ChiSquare[i] = model.fit( result->Redshifts[i], lambdaRange, result->LineModelSolutions[i], contreest_iterations, false );
            calculatedLargeGridRedshifts.push_back(result->Redshifts[i]);
            calculatedLargeGridMerits.push_back(result->ChiSquare[i]);
            result->ScaleMargCorrection[i] = model.getScaleMargCorrection();
            result->SetChisquareTplshapeResult(i , model.GetChisquareTplshape(), model.GetScaleMargTplshape(), model.GetStrongELPresentTplshape());
            for(Int32 k=0; k<result->ChiSquareTplshapes.size(); k++)
            {
                calculatedChiSquareTplshapes[k].push_back(result->ChiSquareTplshapes[k][i]);
            }
            if(estimateLeastSquareFast)
            {
                result->ChiSquareContinuum[i] = model.getLeastSquareContinuumMerit(lambdaRange);
            }else{
                result->ChiSquareContinuum[i] = model.getLeastSquareContinuumMeritFast();
            }
            result->ScaleMargCorrectionContinuum[i] = model.getContinuumScaleMargCorrection();
            Log.LogDebug( "Z interval %d: Chi2 = %f", i, result->ChiSquare[i] );
            indexLargeGrid++;
            //Log.LogInfo( "\nLineModel Infos: large grid step %d", i);
        }else{
            result->ChiSquare[i] = result->ChiSquare[i-1] + 1e-2; //these values will be replaced by the fine grid interpolation below...
            result->ScaleMargCorrection[i] = model.getScaleMargCorrection();
            result->LineModelSolutions[i] = result->LineModelSolutions[i-1];
            result->SetChisquareTplshapeResult( i, result->GetChisquareTplshapeResult(i-1), result->GetScaleMargCorrTplshapeResult(i-1), result->GetStrongELPresentTplshapeResult(i-1));
            if(estimateLeastSquareFast)
            {
                result->ChiSquareContinuum[i] = model.getLeastSquareContinuumMerit(lambdaRange);
            }else{
                result->ChiSquareContinuum[i] = model.getLeastSquareContinuumMeritFast();
            }
            result->ScaleMargCorrectionContinuum[i] = model.getContinuumScaleMargCorrection();
        }
    }


    //now interpolate large grid merit results onto the fine grid
    if(result->Redshifts.size()>calculatedLargeGridMerits.size() && calculatedLargeGridMerits.size()>1)
    {
        interpolateLargeGridOnFineGrid( calculatedLargeGridRedshifts, result->Redshifts, calculatedLargeGridMerits, result->ChiSquare);
    }
    for(Int32 kts=0; kts<result->ChiSquareTplshapes.size(); kts++)
    {
        if(result->Redshifts.size()>calculatedChiSquareTplshapes[kts].size() && calculatedChiSquareTplshapes[kts].size()>1)
        {
            interpolateLargeGridOnFineGrid( calculatedLargeGridRedshifts, result->Redshifts, calculatedChiSquareTplshapes[kts], result->ChiSquareTplshapes[kts]);
        }
    }


    //WARNING: HACK, first pass with continuum from spectrum.
    //model.SetContinuumComponent(opt_continuumcomponent);
    //model.InitFitContinuum();
    //
    boost::chrono::thread_clock::time_point stop_mainloop = boost::chrono::thread_clock::now();
    Float64 duration_mainloop = boost::chrono::duration_cast<boost::chrono::microseconds>(stop_mainloop - start_mainloop).count();
    Float64 duration_firstpass_seconds = duration_mainloop/1e6;
    Log.LogInfo( "  Operator-Linemodel: first-pass done in %.4e sec", duration_firstpass_seconds);
    Log.LogInfo("  Operator-Linemodel: <proc-lm-fistpass><%d>", (Int32)duration_firstpass_seconds);

    // extrema
    Int32 extremumCount = opt_extremacount;
    Log.LogDebug( "  Operator-Linemodel: opt_extremacount = %d", opt_extremacount );
    TPointList extremumList;
    TFloat64Range redshiftsRange(result->Redshifts[0], result->Redshifts[result->Redshifts.size()-1]);
    Log.LogDebug( "  Operator-Linemodel: redshiftsRange.GetBegin() = %f, redshiftsRange.GetEnd() = %f",
		redshiftsRange.GetBegin(), redshiftsRange.GetEnd() );

    if(result->Redshifts.size() == 1)
    {
        extremumList.push_back(SPoint( result->Redshifts[0], result->ChiSquare[0] ));
        Log.LogInfo("  Operator-Linemodel: only 1 redshift calculated, only 1 extremum");
    }else if(opt_extremacount==-1)
    {
        for(Int32 ke=0; ke<result->Redshifts.size(); ke++)
        {
            extremumList.push_back(SPoint( result->Redshifts[ke], result->ChiSquare[ke] ));
        }
        Log.LogInfo("  Operator-Linemodel: all initial redshifts considered as extrema");
    }else
    {
        Log.LogInfo("  Operator-Linemodel: ChiSquare min val = %e", result->GetMinChiSquare() );
        Log.LogInfo("  Operator-Linemodel: ChiSquare max val = %e", result->GetMaxChiSquare() );
        CExtremum extremum( redshiftsRange, extremumCount, true, 2);
        extremum.Find( result->Redshifts, result->ChiSquare, extremumList );
        if( extremumList.size() == 0 )
        {
            Log.LogError("  Operator-Linemodel: Extremum find method failed");
            return result;
        }

        Log.LogInfo( "  Operator-Linemodel: found %d extrema", extremumList.size());
    }
    /*
    // Refine Extremum with a second maximum search around the z candidates:
    // This corresponds to the finer xcorrelation in EZ Pandora (in standard_DP fctn in SolveKernel.py)
    Float64 radius = 0.001;
    for( Int32 i=0; i<extremumList.size(); i++ )
    {
        Float64 x = extremumList[i].X;
        Float64 left_border = max(redshiftsRange.GetBegin(), x-radius);
        Float64 right_border=min(redshiftsRange.GetEnd(), x+radius);

        TPointList extremumListFine;
        TFloat64Range rangeFine = TFloat64Range( left_border, right_border );
        CExtremum extremumFine( rangeFine , 1, true);
        extremumFine.Find( result->Redshifts, result->ChiSquare, extremumListFine );
        if(extremumListFine.size()>0){
            extremumList[i] = extremumListFine[0];
        }
    }
    //*/

    //*
    // extend z around the extrema
    Float64 extensionradius = m_secondPass_extensionradius;
    for( Int32 i=0; i<extremumList.size(); i++ )
    {
        Log.LogInfo("  Operator-Linemodel: Raw extr #%d, z_e.X=%f, m_e.Y=%e", i, extremumList[i].X, extremumList[i].Y);
        Float64 x = extremumList[i].X;
        Float64 left_border = max(redshiftsRange.GetBegin(), x-extensionradius);
        Float64 right_border=min(redshiftsRange.GetEnd(), x+extensionradius);

        for (Int32 i=0;i<result->Redshifts.size();i++)
        {
            if(result->Redshifts[i] >= left_border && result->Redshifts[i] <= right_border){
                result->ExtremaResult.ExtremaExtendedRedshifts.push_back(result->Redshifts[i]);
            }
        }
    }
    //*/
   //todo: remove duplicate redshifts from the extended extrema list


    //Set model parameters to SECOND-PASS
    model.setPassMode(2);
    Log.LogInfo( "  Operator-Linemodel: ---------- ---------- ---------- ----------" );
    Log.LogInfo( "  Operator-Linemodel: now computing second-pass");
    Log.LogInfo( "  Operator-Linemodel: ---------- ---------- ---------- ----------" );

    std::vector<Float64> extrema_velocityEL;
    std::vector<Float64> extrema_velocityAL;
    extrema_velocityEL.resize(extremumList.size());
    extrema_velocityAL.resize(extremumList.size());

    TPointList extremumList2;
    extremumList2.resize(extremumList.size());

    bool modelInfoSave = false;
    std::vector<std::shared_ptr<CModelSpectrumResult>  > savedModelSpectrumResults_lmfit;
    std::vector<std::shared_ptr<CModelFittingResult>  > savedModelFittingResults_lmfit;
    std::vector<std::shared_ptr<CModelRulesResult>  > savedModelRulesResults_lmfit;
    std::vector<  std::shared_ptr<CSpectraFluxResult>> savedBaselineResult_lmfit;

    //recompute the fine grid results around the extrema
    if(enableFastFitLargeGrid==1 || enableVelocityFitting)
    {
        for( Int32 i=0; i<extremumList.size(); i++ )
        {
            Log.LogInfo( "  Operator-Linemodel: Computing second pass for extremum #%d", i);
            Log.LogInfo( "  Operator-Linemodel: ---------- /\\ ---------- ---------- ---------- %d", i );
            Float64 z = extremumList[i].X;
            Float64 m = extremumList[i].Y;

            //find the index in the zaxis results
            Int32 idx=-1;
            for ( UInt32 i2=0; i2<result->Redshifts.size(); i2++)
            {
                if(result->Redshifts[i2] == z){
                    idx = i2;
                    break;
                }
            }
            if(idx==-1){
                Log.LogInfo( "  Operator-Linemodel: Problem. could not find extrema solution index...");
                continue;
            }


            // reestimate the model (eventually with continuum reestimation) on the extrema selected
            if(opt_continuumreest == "always"){
                contreest_iterations = 1;
            }else{
                contreest_iterations  = 0;
            }

            //model.LoadModelSolution(result->LineModelSolutions[idx]);
            model.fit( result->Redshifts[idx], lambdaRange, result->LineModelSolutions[idx], contreest_iterations, false );
            m = result->ChiSquare[idx];
            if(enableVelocityFitting)
            {
                Bool enableManualStepVelocityFit = true;
                Bool enableLMVelocityFit = false;
                Bool enableLBFGSVelocityFit = false;
                if(enableLMVelocityFit)
                {
                    //fit the emission and absorption width using the linemodel lmfit strategy
                    model.SetFittingMethod("lmfit");
                    //model.SetElementIndexesDisabledAuto();
                    Float64 meritTmp;
                    Log.LogInfo("  Operator-Linemodel: Lm fit for extrema %d", i);
                    meritTmp = model.fit( result->Redshifts[idx], lambdaRange, result->LineModelSolutions[idx], contreest_iterations, true );
                    modelInfoSave = true;
                    // CModelSpectrumResult
                    std::shared_ptr<CModelSpectrumResult>  resultspcmodel = std::shared_ptr<CModelSpectrumResult>( new CModelSpectrumResult(model.GetModelSpectrum()) );
                    //std::shared_ptr<CModelSpectrumResult>  resultspcmodel = std::shared_ptr<CModelSpectrumResult>( new CModelSpectrumResult(model.GetObservedSpectrumWithLinesRemoved()) );
                    savedModelSpectrumResults_lmfit.push_back(resultspcmodel);
                    // CModelFittingResult
                    std::shared_ptr<CModelFittingResult>  resultfitmodel = std::shared_ptr<CModelFittingResult>( new CModelFittingResult(result->LineModelSolutions[idx], result->Redshifts[idx], result->ChiSquare[idx], result->restRayList, model.GetVelocityEmission(), model.GetVelocityAbsorption()) );
                    savedModelFittingResults_lmfit.push_back(resultfitmodel);
                    // CModelRulesResult
                    std::shared_ptr<CModelRulesResult>  resultrulesmodel = std::shared_ptr<CModelRulesResult>( new CModelRulesResult( model.GetModelRulesLog() ));
                    savedModelRulesResults_lmfit.push_back(resultrulesmodel);


                    std::shared_ptr<CSpectraFluxResult> baselineResult_lmfit = (std::shared_ptr<CSpectraFluxResult>) new CSpectraFluxResult();
                    baselineResult_lmfit->m_optio = 0;
                    const CSpectrumFluxAxis& modelContinuumFluxAxis = model.GetModelContinuum();
                    UInt32 len = modelContinuumFluxAxis.GetSamplesCount();

                    baselineResult_lmfit->fluxes.resize(len);
                    baselineResult_lmfit->wavel.resize(len);
                    for( Int32 k=0; k<len; k++ )
                    {
                        baselineResult_lmfit->fluxes[k] = modelContinuumFluxAxis[k];
                        baselineResult_lmfit->wavel[k]  = (spectrum.GetSpectralAxis())[k];
                    }
                    savedBaselineResult_lmfit.push_back(baselineResult_lmfit);

                    z = result->LineModelSolutions[idx].Redshift;
                    result->ExtremaResult.lmfitPass.push_back(z);
                    //result->Redshifts[idx] = z;

                    model.SetFittingMethod(opt_fittingmethod);
                    model.ResetElementIndexesDisabled();
                    Int32 velocityHasBeenReset = model.ApplyVelocityBound(velfitMinE, velfitMaxE);
                    enableManualStepVelocityFit = velocityHasBeenReset;
                }

                if(enableLBFGSVelocityFit){
                    //fit the emission and absorption width using the linemodel lbfgsfit strategy
                    model.SetFittingMethod("lbfgsfit");
                    //model.SetElementIndexesDisabledAuto();
                    Float64 meritTmp;
                    meritTmp = model.fit( result->Redshifts[idx], lambdaRange, result->LineModelSolutions[idx], contreest_iterations, false );

                    model.SetFittingMethod(opt_fittingmethod);
                    model.ResetElementIndexesDisabled();
                    Int32 velocityHasBeenReset = model.ApplyVelocityBound(velfitMinE, velfitMaxE);
                    enableManualStepVelocityFit = velocityHasBeenReset;
                }

                if(enableManualStepVelocityFit){
                    //fit the emission and absorption width by minimizing the linemodel merit with linemodel "hybrid" fitting method
                    model.SetFittingMethod("hybrid");
                    //model.m_enableAmplitudeOffsets = true;
                    //contreest_iterations = 1;
                    std::vector<std::vector<Int32>> idxVelfitGroups;

                    for(Int32 iLineType = 0; iLineType<2; iLineType++)
                    {
                        Float64 vInfLim;
                        Float64 vSupLim;

                        if(iLineType==0)
                        {
                            Log.LogInfo( "  Operator-Linemodel: manualStep velocity fit ABSORPTION, for z = %.6f", result->Redshifts[idx]);
                            vInfLim = velfitMinA;
                            vSupLim = velfitMaxA;
                            if(m_enableWidthFitByGroups)
                            {
                                idxVelfitGroups.clear();
                                idxVelfitGroups = model.GetModelVelfitGroups(CRay::nType_Absorption);
                                Log.LogInfo( "  Operator-Linemodel: VelfitGroups ABSORPTION - n = %d", idxVelfitGroups.size());
                                if(extremumList.size()>1 && idxVelfitGroups.size()>1)
                                {
                                    Log.LogError("  Operator-Linemodel: not allowed to use more than 1 group per E/A for more than 1 extremum (see .json linemodel.extremacount)");
                                }
                            }
                        }else{
                            Log.LogInfo( "  Operator-Linemodel: manualStep velocity fit EMISSION, for z = %.6f", result->Redshifts[idx]);
                            vInfLim = velfitMinE;
                            vSupLim = velfitMaxE;
                            if(m_enableWidthFitByGroups)
                            {
                                idxVelfitGroups.clear();
                                idxVelfitGroups = model.GetModelVelfitGroups(CRay::nType_Emission);
                                Log.LogInfo( "  Operator-Linemodel: VelfitGroups EMISSION - n = %d", idxVelfitGroups.size());
                                if(extremumList.size()>1 && idxVelfitGroups.size()>1)
                                {
                                    Log.LogError("  Operator-Linemodel: not allowed to use more than 1 group per E/A for more than 1 extremum (see .json linemodel.extremacount)");
                                }
                            }
                        }

                        //Prepare velocity grid to be checked
                        Float64 vStep = m_secondPass_velfit_vStep;;
                        Int32 nSteps = (int)((vSupLim-vInfLim)/vStep);


                        Float64 dzInfLim = m_secondPass_velfit_dzInfLim;
                        if(result->Redshifts[idx]+dzInfLim<result->Redshifts[0])
                        {
                            dzInfLim = result->Redshifts[0]-result->Redshifts[idx];
                        }
                        Float64 dzSupLim = m_secondPass_velfit_dzSupLim;
                        if(result->Redshifts[idx]+dzSupLim>result->Redshifts[result->Redshifts.size()-1])
                        {
                            dzSupLim = result->Redshifts[result->Redshifts.size()-1]-result->Redshifts[idx];
                        }

                        Float64 dzStep = m_secondPass_velfit_dzStep;
                        Int32 nDzSteps = (int)((dzSupLim-dzInfLim)/dzStep);
                        if(nDzSteps==0)
                        {
                            nDzSteps=1;
                            dzInfLim = 0.;
                            dzSupLim = 0.;
                        }
                        Log.LogInfo( "  Operator-Linemodel: dzInfLim n=%e", dzInfLim);
                        Log.LogInfo( "  Operator-Linemodel: dzSupLim n=%e", dzSupLim);
                        Log.LogInfo( "  Operator-Linemodel: manualStep n=%d", nDzSteps);


                        for(Int32 kgroup=0; kgroup<idxVelfitGroups.size(); kgroup++)
                        {
                            Log.LogInfo( "  Operator-Linemodel: manualStep fitting group=%d", kgroup);

                            Float64 meritMin = DBL_MAX;
                            Float64 vOptim = -1.0;
                            Float64 dzOptim = -1.0;
                            for(Int32 kdz=0; kdz<nDzSteps; kdz++)
                            {
                                Float64 dzTest = dzInfLim+kdz*dzStep;
                                for(Int32 kv=0; kv<nSteps; kv++)
                                {
                                    Float64 vTest = vInfLim+kv*vStep;
                                    if(iLineType==0)
                                    {
                                        if(m_enableWidthFitByGroups)
                                        {
                                            for(Int32 ke=0; ke<idxVelfitGroups[kgroup].size(); ke++)
                                            {
                                                model.SetVelocityAbsorptionOneElement(vTest, idxVelfitGroups[kgroup][ke]);
                                            }
                                        }else{
                                            model.SetVelocityAbsorption(vTest);
                                        }
                                    }else
                                    {
                                        if(m_enableWidthFitByGroups)
                                        {
                                            for(Int32 ke=0; ke<idxVelfitGroups[kgroup].size(); ke++)
                                            {
                                                model.SetVelocityEmissionOneElement(vTest, idxVelfitGroups[kgroup][ke]);
                                            }
                                        }else{
                                            model.SetVelocityEmission(vTest);
                                        }
                                    }

                                    //Log.LogInfo( "  Operator-Linemodel: testing v=%f", vTest);
                                    Float64 meritv;
                                    meritv = model.fit( result->Redshifts[idx]+dzTest, lambdaRange, result->LineModelSolutions[idx], contreest_iterations, false );

//                                    if(m_enableWidthFitByGroups)
//                                    {
//                                        meritv = 0.0;
//                                        for(Int32 ke=0; ke<idxVelfitGroups[kgroup].size(); ke++)
//                                        {
//                                            meritv += model.getModelErrorUnderElement(idxVelfitGroups[kgroup][ke]);
//                                        }
//                                    }

                                    //

                                    //Log.LogInfo("  Operator-Linemodel: testing velocity: merit=%.3e for velocity = %.1f", meritv, vTest);
                                    if(meritMin>meritv)
                                    {
                                        meritMin = meritv;
                                        if(iLineType==0)
                                        {
                                            vOptim = model.GetVelocityAbsorption();
                                        }else
                                        {
                                            vOptim = model.GetVelocityEmission();
                                        }
                                        dzOptim = dzTest;
                                    }
                                }
                            }
                            if(vOptim != -1.0)
                            {
                                Log.LogInfo( "  Operator-Linemodel: best Velocity found = %.1f", vOptim);
                                result->ChiSquare[idx] =  meritMin;
                                if(iLineType==0)
                                {
                                    if(m_enableWidthFitByGroups)
                                    {
                                        for(Int32 ke=0; ke<idxVelfitGroups[kgroup].size(); ke++)
                                        {
                                            model.SetVelocityAbsorptionOneElement(vOptim, idxVelfitGroups[kgroup][ke]);
                                        }
                                    }else{
                                        model.SetVelocityAbsorption(vOptim);
                                    }
                                }else
                                {
                                    if(m_enableWidthFitByGroups)
                                    {
                                        for(Int32 ke=0; ke<idxVelfitGroups[kgroup].size(); ke++)
                                        {
                                            model.SetVelocityEmissionOneElement(vOptim, idxVelfitGroups[kgroup][ke]);
                                        }
                                    }else{
                                        model.SetVelocityEmission(vOptim);
                                    }
                                }

                            }
                        }
                    }
                    model.SetFittingMethod(opt_fittingmethod);
                    //model.m_enableAmplitudeOffsets = false;
                }
            }
            extrema_velocityEL[i]=model.GetVelocityEmission();
            extrema_velocityAL[i]=model.GetVelocityAbsorption();

            //finally compute the redshifts on the z-range around the extremum
            Float64 left_border = max(redshiftsRange.GetBegin(), z-extensionradius);
            Float64 right_border=min(redshiftsRange.GetEnd(), z+extensionradius);
            //model.SetFittingMethod("nofit");
            extremumList2[i].Y = DBL_MAX;
            extremumList2[i].X = result->Redshifts[idx];
            Int32 idx2 = idx;
            for (Int32 iz=0;iz<result->Redshifts.size();iz++)
            {
                if(result->Redshifts[iz] >= left_border && result->Redshifts[iz] <= right_border){
                    //Log.LogInfo("Fit for Extended redshift %d, z = %f", iz, result->Redshifts[iz]);
                    result->ChiSquare[iz] = model.fit( result->Redshifts[iz], lambdaRange, result->LineModelSolutions[iz], contreest_iterations, false );
                    result->ScaleMargCorrection[iz] = model.getScaleMargCorrection();
                    result->SetChisquareTplshapeResult(iz , model.GetChisquareTplshape(), model.GetScaleMargTplshape(), model.GetStrongELPresentTplshape());
                    if(estimateLeastSquareFast)
                    {
                        result->ChiSquareContinuum[iz] = model.getLeastSquareContinuumMerit(lambdaRange);
                    }else{
                        result->ChiSquareContinuum[iz] = model.getLeastSquareContinuumMeritFast();
                    }
                    result->ScaleMargCorrectionContinuum[iz] = model.getContinuumScaleMargCorrection();
                    if(result->ChiSquare[iz]< extremumList2[i].Y)
                    {
                        extremumList2[i].X = result->Redshifts[iz];
                        extremumList2[i].Y = result->ChiSquare[iz];
                        idx2 = iz;
                    }
                }
            }
            //model.SetFittingMethod(opt_fittingmethod);

            Log.LogInfo("  Operator-Linemodel: Recomputed extr #%d, idx=%d, z_e.X=%f, m_e.Y=%f", i, idx2, extremumList2[i].X, extremumList2[i].Y);

        }
    }else
    {
        for( Int32 i=0; i<extremumList.size(); i++ )
        {
            extremumList2[i].X = extremumList[i].X;
            extremumList2[i].Y = extremumList[i].Y;
        }
    }

    //return result;
    //reorder extremumList using .Y values : smallest to highest
    TPointList extremumListOrdered;
    std::vector<Float64> extrema_velocityELOrdered;
    std::vector<Float64> extrema_velocityALOrdered;
    extremumCount = extremumList2.size();
    std::vector<Int32> indiceList2;

    for(Int32 ie=0; ie<extremumCount; ie++)
    {
        Int32 iYmin=0;
        Float64 YMin=DBL_MAX;
        for(Int32 ie2=0; ie2<extremumList2.size(); ie2++)
        {
            if(YMin>extremumList2[ie2].Y){
                YMin = extremumList2[ie2].Y;
                iYmin = ie2;
            }
        }

        extremumListOrdered.push_back(extremumList2[iYmin]);
        extremumList2.erase(extremumList2.begin() + iYmin);
        extrema_velocityELOrdered.push_back(extrema_velocityEL[iYmin]);
        extrema_velocityEL.erase(extrema_velocityEL.begin() + iYmin);
        extrema_velocityALOrdered.push_back(extrema_velocityAL[iYmin]);
        extrema_velocityAL.erase(extrema_velocityAL.begin() + iYmin);
        //founding the initial index in extrenumList2
        for(Int32 ie2=0; ie2<indiceList2.size(); ie2++)
        {
            if(iYmin>=indiceList2[ie2]){
                iYmin++;
            }
        }
        indiceList2.push_back(iYmin);
    }
    if(modelInfoSave){
        TFloat64List OrderedLMZ;
        for(Int32 ie2=0; ie2<indiceList2.size(); ie2++)
        {
            OrderedLMZ.push_back(result->ExtremaResult.lmfitPass[indiceList2[ie2]]);
        }
        result->ExtremaResult.lmfitPass = OrderedLMZ;
    }
    if( extremumListOrdered.size() == 0 )
    {
        Log.LogError("Line Model, Extremum Ordering failed");
    }

    Log.LogInfo("  Operator-Linemodel: Now storing extrema results");
    // store extrema results
    result->ExtremaResult.Resize(extremumCount);

    //Int32 start = spectrum.GetSpectralAxis().GetIndexAtWaveLength(lambdaRange.GetBegin());
    //Int32 end = spectrum.GetSpectralAxis().GetIndexAtWaveLength(lambdaRange.GetEnd());
    //Int32 nsamples = end - start + 1;
    Int32 savedModels = 0;
    m_savedModelSpectrumResults.clear();
    m_savedModelFittingResults.clear();
    m_savedModelRulesResults.clear();
    m_savedModelContinuumSpectrumResults.clear();

    for( Int32 i=0; i<extremumCount; i++ )
    {
        Float64 z = extremumListOrdered[i].X;
        Float64 m = extremumListOrdered[i].Y;

        //find the index in the zaxis results
        Int32 idx=-1;
        for ( UInt32 i2=0; i2<result->Redshifts.size(); i2++)
        {
            if(result->Redshifts[i2] == z){
                idx = i2;
                break;
            }
        }
        if(idx==-1){
            Log.LogInfo( "Problem. could not find extrema solution index...");
            continue;
        }
        Log.LogInfo("  Operator-Linemodel: Saving extr #%d, idx=%d, z=%f, m=%f", i, idx, result->Redshifts[idx], result->ChiSquare[idx]);


        // reestimate the model (eventually with continuum reestimation) on the extrema selected
        if(opt_continuumreest == "always" || opt_continuumreest == "onlyextrema"){
            contreest_iterations = 8; //4
        }else{
            contreest_iterations  = 0;
        }


        if(enableVelocityFitting && extremumCount>1){
//            ModelFit( model, lambdaRange, result->Redshifts[idx], result->ChiSquare[idx], result->LineModelSolutions[idx], contreest_iterations);
//            m = result->ChiSquare[idx];
//            //fit the emission and absorption width using the lindemodel lmfit strategy
//            model.SetFittingMethod("lmfit");
//            model.SetElementIndexesDisabledAuto();
//            ModelFit( model, lambdaRange, result->Redshifts[idx], result->ChiSquare[idx], result->LineModelSolutions[idx], contreest_iterations);


//            model.SetFittingMethod(opt_fittingmethod);
//            model.ResetElementIndexesDisabled();
//            model.ApplyVelocityBound();
//        }else{
            model.SetVelocityEmission(extrema_velocityELOrdered[i]);
            model.SetVelocityAbsorption(extrema_velocityALOrdered[i]);
        }


        if(!modelInfoSave){
            result->ChiSquare[idx] = model.fit( result->Redshifts[idx], lambdaRange, result->LineModelSolutions[idx], contreest_iterations, true );
            result->ScaleMargCorrection[idx] = model.getScaleMargCorrection();
            result->SetChisquareTplshapeResult(idx , model.GetChisquareTplshape(), model.GetScaleMargTplshape(), model.GetStrongELPresentTplshape());
            if(estimateLeastSquareFast)
            {
                result->ChiSquareContinuum[idx] = model.getLeastSquareContinuumMerit(lambdaRange);
            }else{
                result->ChiSquareContinuum[idx] = model.getLeastSquareContinuumMeritFast();
            }
            result->ScaleMargCorrectionContinuum[idx] = model.getContinuumScaleMargCorrection();
        }
        if(m!=result->ChiSquare[idx])
        {
            Log.LogInfo("  Operator-Linemodel: m (%f for idx=%d) !=chi2 (%f) ", m, idx, result->ChiSquare[idx]);
        }
        m = result->ChiSquare[idx];//result->ChiSquare[idx];

        //save the model result
        //WARNING: saving results TODO: this is currently wrong !! the model saved corresponds to the bestchi2 model. PDFs should be combined prior to exporting the best model for each extrema...
        static Int32 maxModelSave = std::min(m_maxModelSaveCount, extremumCount);
        Int32 maxSaveNLinemodelContinua = maxModelSave;
        if( savedModels<maxModelSave )
        {
          if(modelInfoSave){
            Log.LogInfo("Save model store during lm_fit");
            m_savedModelSpectrumResults.push_back(savedModelSpectrumResults_lmfit[indiceList2[i]]);
            m_savedModelFittingResults.push_back(savedModelFittingResults_lmfit[indiceList2[i]]);
            m_savedModelRulesResults.push_back(savedModelRulesResults_lmfit[indiceList2[i]]);
            if( savedModels < maxSaveNLinemodelContinua && contreest_iterations>0)
            {
                m_savedModelContinuumSpectrumResults.push_back(savedBaselineResult_lmfit[indiceList2[i]]);
            }
          }else{
            // CModelSpectrumResult
            std::shared_ptr<CModelSpectrumResult>  resultspcmodel;
            Int32 overrideModelSavedType = 0;
            //0=save model, (DEFAULT)
            //1=save model with lines removed,
            //2=save model with only Em. lines removed.
            if(overrideModelSavedType==0)
            {
                resultspcmodel = std::shared_ptr<CModelSpectrumResult>( new CModelSpectrumResult(model.GetModelSpectrum()) );
            }else if(overrideModelSavedType==1 || overrideModelSavedType==2)
            {
                Int32 lineTypeFilter=-1;
                if(overrideModelSavedType==1)
                {
                    lineTypeFilter=-1;
                }else if(overrideModelSavedType==2)
                {
                    lineTypeFilter = CRay::nType_Emission;
                }
                resultspcmodel = std::shared_ptr<CModelSpectrumResult>( new CModelSpectrumResult(model.GetObservedSpectrumWithLinesRemoved(lineTypeFilter)) );
            }
            //std::shared_ptr<CModelSpectrumResult>  resultspcmodel = std::shared_ptr<CModelSpectrumResult>( new CModelSpectrumResult(model.GetSpectrumModelContinuum()) );

            m_savedModelSpectrumResults.push_back(resultspcmodel);

            // CModelFittingResult
            std::shared_ptr<CModelFittingResult>  resultfitmodel = std::shared_ptr<CModelFittingResult>( new CModelFittingResult(result->LineModelSolutions[idx], result->Redshifts[idx], result->ChiSquare[idx], result->restRayList, model.GetVelocityEmission(), model.GetVelocityAbsorption()) );
            m_savedModelFittingResults.push_back(resultfitmodel);

            // CModelRulesResult
            std::shared_ptr<CModelRulesResult>  resultrulesmodel = std::shared_ptr<CModelRulesResult>( new CModelRulesResult( model.GetModelRulesLog() ));
            m_savedModelRulesResults.push_back(resultrulesmodel);



            if( savedModels < maxSaveNLinemodelContinua )
            {
                // Save the reestimated continuum, only the first n=maxSaveNLinemodelContinua extrema
                std::shared_ptr<CSpectraFluxResult> baselineResult = (std::shared_ptr<CSpectraFluxResult>) new CSpectraFluxResult();
                baselineResult->m_optio = 0;
                const CSpectrumFluxAxis& modelContinuumFluxAxis = model.GetModelContinuum();
                UInt32 len = modelContinuumFluxAxis.GetSamplesCount();

                baselineResult->fluxes.resize(len);
                baselineResult->wavel.resize(len);
                for( Int32 k=0; k<len; k++ )
                {
                    baselineResult->fluxes[k] = modelContinuumFluxAxis[k];
                    baselineResult->wavel[k]  = (spectrum.GetSpectralAxis())[k];
                }
                m_savedModelContinuumSpectrumResults.push_back(baselineResult);
            }
          }
          savedModels++;
        }


        result->ExtremaResult.Extrema[i] = z;
        result->ExtremaResult.ExtremaMerit[i] = m;
        if(estimateLeastSquareFast)
        {
            result->ExtremaResult.ExtremaMeritContinuum[i] = model.getLeastSquareContinuumMerit(lambdaRange);
        }else{
            result->ExtremaResult.ExtremaMeritContinuum[i] = model.getLeastSquareContinuumMeritFast();
        }

        result->ExtremaResult.ExtremaLastPass[i] = z; //refined extremum is initialized here.


        //computing errz (or deltaz, dz...): should probably be computed in linemodelresult.cpp instead ?
        Float64 dz=-1.;
        if(result->Redshifts.size()>1)
        {
            CDeltaz deltaz;
            Float64 zRangeHalf = 0.005;
            TFloat64Range range = TFloat64Range(z-zRangeHalf, z+zRangeHalf);
            //Int32 ret = deltaz.Compute(result->ChiSquare, result->Redshifts, z, range, dz);
            Int32 ret = deltaz.Compute3ddl(result->ChiSquare, result->Redshifts, z, range, dz);
            if(ret!=0)
            {
                Log.LogError("  Operator-Linemodel: Deltaz computation failed");
            }
        }
        result->ExtremaResult.DeltaZ[i] = dz;

        //store model Ha SNR
        result->ExtremaResult.snrHa[i] = result->LineModelSolutions[idx].snrHa;

        //store the model norm
        result->ExtremaResult.mTransposeM[i] = model.EstimateMTransposeM(lambdaRange);

        //scale marginalization correction
        Float64 corrScaleMarg = model.getScaleMargCorrection();//
        result->ExtremaResult.CorrScaleMarg[i] = corrScaleMarg;

        static Float64 cutThres = 3.0;
        Int32 nValidLines = result->GetNLinesOverCutThreshold(i, cutThres, cutThres);
        result->ExtremaResult.Posterior[i] = nValidLines;//m/Float64(1+nValidLines);
        Float64 cumulStrongELSNR = model.getCumulSNRStrongEL(); //getStrongerMultipleELAmpCoeff();
        result->ExtremaResult.StrongELSNR[i] = cumulStrongELSNR;

        result->ExtremaResult.LogArea[i] = -DBL_MAX;
        result->ExtremaResult.LogAreaCorrectedExtrema[i] = -1.0;


        Int32 nddl = model.GetNElements(); //get the total number of elements in the model
        nddl = result->LineModelSolutions[idx].nDDL; //override nddl by the actual number of elements in the fitted model
        result->ExtremaResult.NDof[i] = model.GetModelNonZeroElementsNDdl();

        Float64 bic = m + nddl*log(result->nSpcSamples); //BIC
        //Float64 aic = m + 2*nddl; //AIC
        result->ExtremaResult.bic[i] = bic;
        //result->bic[i] = aic + (2*nddl*(nddl+1) )/(nsamples-nddl-1);  //AICc, better when nsamples small

        //compute continuum indexes
        CContinuumIndexes continuumIndexes;
        result->ExtremaResult.ContinuumIndexes[i] = continuumIndexes.getIndexes( spectrum, z );

        //save the outsideLinesMask
        result->ExtremaResult.OutsideLinesMask[i] = model.getOutsideLinesMask();

        //save the continuum tpl fitting results
        result->ExtremaResult.FittedTplName[i] = model.getFitContinuum_tplName();
        result->ExtremaResult.FittedTplAmplitude[i] = model.getFitContinuum_tplAmplitude();
        result->ExtremaResult.FittedTplMerit[i] = model.getFitContinuum_tplMerit();
        result->ExtremaResult.FittedTplDustCoeff[i] = model.getFitContinuum_tplIsmDustCoeff();
        result->ExtremaResult.FittedTplMeiksinIdx[i] = model.getFitContinuum_tplIgmMeiksinIdx();

        //save the tplcorr results
        result->ExtremaResult.FittedTplshapeName[i] = model.getTplshape_bestTplName();
    }

    //ComputeArea2(*result);

    return result;

}


/**
 * @brief COperatorLineModel::processPass
 * @return
 */
std::shared_ptr<COperatorResult> COperatorLineModel::computeWithUltimPass(CDataStore &dataStore,
                                  const CSpectrum& spectrum,
                                  const CSpectrum& spectrumContinuum,
                                  const CTemplateCatalog& tplCatalog,
                                  const TStringList& tplCategoryList,
                                  const std::string opt_calibrationPath,
                                  const CRayCatalog& restraycatalog,
                                  const std::string& opt_lineTypeFilter,
                                  const std::string& opt_lineForceFilter,
                                  const TFloat64Range& lambdaRange,
                                  const TFloat64List& redshifts,
                                  const Int32 opt_extremacount,
                                  const std::string& opt_fittingmethod,
                                  const std::string& opt_continuumcomponent,
                                  const std::string& opt_lineWidthType,
                                  const Float64 opt_resolution,
                                  const Float64 opt_velocityEmission,
                                  const Float64 opt_velocityAbsorption,
                                  const std::string& opt_continuumreest,
                                  const std::string& opt_rules,
                                  const std::string& opt_velocityFitting,
                                  const Float64 &opt_twosteplargegridstep,
                                  const std::string& opt_rigidity,
                                  const Float64 &opt_velocityfitmin,
                                  const Float64 &opt_velocityfitmax)
{
    auto result = std::dynamic_pointer_cast<CLineModelResult>(Compute( dataStore,
                            spectrum,
                            spectrumContinuum,
                            tplCatalog,
                            tplCategoryList,
                            opt_calibrationPath,
                            restraycatalog,
                            opt_lineTypeFilter,
                            opt_lineForceFilter,
                            lambdaRange,
                            redshifts,
                            opt_extremacount,
                            opt_fittingmethod,
                            opt_continuumcomponent,
                            opt_lineWidthType,
                            opt_resolution,
                            opt_velocityEmission,
                            opt_velocityAbsorption,
                            opt_continuumreest,
                            opt_rules,
                            opt_velocityFitting,
                            opt_twosteplargegridstep,
                            opt_rigidity,
                            opt_velocityfitmin,
                            opt_velocityfitmax));

    if(result && opt_rigidity=="tplshape")
    {
       Log.LogInfo("  Operator-Linemodel - Last Pass: begin");
       //
       //do the last pass on the 1st extremum range
       //
       Float64 halfRange = 1e-3;
       Float64 lastPassStep = 1e-4;
       //find the best extremum redhift
       Float64 bestRedshift=-1;
       Float64 bestMerit=DBL_MAX;
       Int32 bestIndex=-1;
       for(Int32 k=0; k<result->ExtremaResult.Extrema.size(); k++)
       {
           Log.LogInfo("  Operator-Linemodel - Last Pass: k = %d", k);
           Log.LogInfo("  Operator-Linemodel - Last Pass: result->Extrema[k] = %.5f", result->ExtremaResult.Extrema[k]);
           Log.LogInfo("  Operator-Linemodel - Last Pass: result->ExtremaMerit[k] = %.5f", result->ExtremaResult.ExtremaMerit[k]);

            if(bestMerit>result->ExtremaResult.ExtremaMerit[k])
            {
                bestMerit = result->ExtremaResult.ExtremaMerit[k];
                bestRedshift = result->ExtremaResult.Extrema[k];
                bestIndex = k;
            }
       }
       Log.LogInfo("  Operator-Linemodel - Last Pass: around extrema z = %.5f", bestRedshift);
       Float64 z=bestRedshift-halfRange;
       std::vector<Float64> lastPassRedshifts;
       if(z<redshifts[0])
       {
           z=redshifts[0];
       }
       while(z<bestRedshift+halfRange)
       {
           lastPassRedshifts.push_back(z);
           z+=lastPassStep;
       }
       Log.LogInfo("  Operator-Linemodel - Last Pass: range zmin=%.5f, zmax=%.5f", lastPassRedshifts[0], lastPassRedshifts[lastPassRedshifts.size()-1]);


       std::string opt_rigidity_lastPass = "rules";
       Int32 opt_extremacount_lastPass = 1;

       Int32 maxSaveBackup = m_maxModelSaveCount;
       m_maxModelSaveCount = 0;
       auto lastPassResult = std::dynamic_pointer_cast<CLineModelResult>(Compute( dataStore,
                               spectrum,
                               spectrumContinuum,
                               tplCatalog,
                               tplCategoryList,
                               opt_calibrationPath,
                               restraycatalog,
                               opt_lineTypeFilter,
                               opt_lineForceFilter,
                               lambdaRange,
                               lastPassRedshifts,
                               opt_extremacount_lastPass,
                               opt_fittingmethod,
                               opt_continuumcomponent,
                               opt_lineWidthType,
                               opt_resolution,
                               opt_velocityEmission,
                               opt_velocityAbsorption,
                               opt_continuumreest,
                               opt_rules,
                               opt_velocityFitting,
                               opt_twosteplargegridstep,
                               opt_rigidity_lastPass,
                               opt_velocityfitmin,
                               opt_velocityfitmax));

        m_maxModelSaveCount = maxSaveBackup;
        Float64 refinedExtremum = lastPassResult->ExtremaResult.Extrema[0];
        Log.LogInfo("  Operator-Linemodel - Last Pass: found refined z=%.5f", refinedExtremum);

        result->ExtremaResult.ExtremaLastPass[bestIndex] = refinedExtremum;
    }else{
        Log.LogInfo("  Operator-Linemodel - Last Pass: failed to do linemodel, rigidity=%s", opt_rigidity.c_str());
    }

    return result;
}

/**
 * @brief COperatorLineModel::initContaminant
 * prepare the contaminant to be used when instantiating a multimodel in compute()
 * @return
 */
Int32 COperatorLineModel::initContaminant(std::shared_ptr<CModelSpectrumResult> contModelSpectrum, Int32 iRollContaminated, Float64 contLambdaOffset)
{
    //
    Log.LogInfo("  Operator-Linemodel: Initializing contaminant for roll #%d, with offset=%.2f", iRollContaminated, contLambdaOffset);

    m_iRollContaminated = iRollContaminated;
    m_contLambdaOffset = contLambdaOffset;
    const std::string& category = "emission";
    m_tplContaminant = std::shared_ptr<CTemplate>( new CTemplate( "contaminant", category ) );
    Int32 length = contModelSpectrum->GetSpectrum().GetSampleCount();
    CSpectrumAxis& contModelFluxAxis = contModelSpectrum->GetSpectrum().GetFluxAxis();
    CSpectrumAxis& contModelSpectralAxis = contModelSpectrum->GetSpectrum().GetSpectralAxis();

    CSpectrumAxis& spcFluxAxis = m_tplContaminant->GetFluxAxis();
    spcFluxAxis.SetSize( length );
    CSpectrumAxis& spcSpectralAxis = m_tplContaminant->GetSpectralAxis();
    spcSpectralAxis.SetSize( length );

    for(Int32 k=0; k<length; k++)
    {
        spcFluxAxis[k] = contModelFluxAxis[k];
        spcSpectralAxis[k] =  contModelSpectralAxis[k];
    }

    //applying offset for ra/dec distance between main source and contaminant
    m_tplContaminant->GetSpectralAxis().ApplyOffset(m_contLambdaOffset);

    m_enableLoadContTemplate = true;
    return 0;
}

///
/// \brief COperatorLineModel::storeGlobalModelResults
/// stores the linemodel results as global results in the datastore
///
void COperatorLineModel::storeGlobalModelResults( CDataStore &dataStore )
{
    Int32 nResults = m_savedModelSpectrumResults.size();
    if( nResults > m_savedModelFittingResults.size())
    {
        Log.LogError( "Line Model, not as many model fitting results as model spectrum results, (nspc = %d, nfit = %d)",  m_savedModelSpectrumResults.size(), m_savedModelFittingResults.size());
        nResults = m_savedModelFittingResults.size();
    }

    for(Int32 k=0; k<nResults; k++)
    {
        std::string fname_spc = (boost::format("linemodel_spc_extrema_%1%") % k).str();
        dataStore.StoreScopedGlobalResult( fname_spc.c_str(), m_savedModelSpectrumResults[k] );

        std::string fname_fit = (boost::format("linemodel_fit_extrema_%1%") % k).str();
        dataStore.StoreScopedGlobalResult( fname_fit.c_str(), m_savedModelFittingResults[k] );

        std::string fname_rules = (boost::format("linemodel_rules_extrema_%1%") % k).str();
        dataStore.StoreScopedGlobalResult( fname_rules.c_str(), m_savedModelRulesResults[k] );
    }


    for(Int32 k=0; k<m_savedModelContinuumSpectrumResults.size(); k++)
    {
        std::string nameBaselineStr = (boost::format("linemodel_continuum_extrema_%1%") % k).str();
        dataStore.StoreScopedGlobalResult(nameBaselineStr.c_str(), m_savedModelContinuumSpectrumResults[k]);
    }

}

///
/// \brief COperatorLineModel::storePerTemplateModelResults
/// stores the linemodel results as per template results in the datastore
///
void COperatorLineModel::storePerTemplateModelResults( CDataStore &dataStore, const CTemplate& tpl )
{
    Int32 nResults = m_savedModelSpectrumResults.size();
    if( nResults > m_savedModelFittingResults.size())
    {
        Log.LogError( "Line Model, not as many model fitting results as model spectrum results, (nspc = %d, nfit = %d)",  m_savedModelSpectrumResults.size(), m_savedModelFittingResults.size());
        nResults = m_savedModelFittingResults.size();
    }

    for(Int32 k=0; k<nResults; k++)
    {
        std::string fname_spc = (boost::format("linemodel_spc_extrema_%1%") % k).str();
        dataStore.StoreScopedPerTemplateResult(  tpl, fname_spc.c_str(), m_savedModelSpectrumResults[k] );

        std::string fname_fit = (boost::format("linemodel_fit_extrema_%1%") % k).str();
        dataStore.StoreScopedPerTemplateResult(  tpl, fname_fit.c_str(), m_savedModelFittingResults[k] );

        std::string fname_rules = (boost::format("linemodel_rules_extrema_%1%") % k).str();
        dataStore.StoreScopedPerTemplateResult(  tpl, fname_rules.c_str(), m_savedModelRulesResults[k] );
    }

}

std::shared_ptr<CModelSpectrumResult> COperatorLineModel::GetModelSpectrumResult(Int32 idx)
{
    Int32 nResults = m_savedModelSpectrumResults.size();
    if(idx>=nResults)
    {
        return NULL;
    }else{
        return m_savedModelSpectrumResults[idx];
    }
}


std::shared_ptr<CSpectraFluxResult> COperatorLineModel::GetModelSpectrumContinuumResult(Int32 idx)
{
    Int32 nResults = m_savedModelContinuumSpectrumResults.size();
    if(idx>=nResults)
    {
        return NULL;
    }else{
        return m_savedModelContinuumSpectrumResults[idx];
    }
}



std::shared_ptr<CModelSpectrumResult> COperatorLineModel::GetContaminantSpectrumResult()
{
    return m_savedContaminantSpectrumResult;
}


Int32 COperatorLineModel::interpolateLargeGridOnFineGrid(TFloat64List redshiftsLargeGrid, TFloat64List redshiftsFineGrid, TFloat64List meritLargeGrid, TFloat64List& meritFineGrid)
{
    //* // GSL method LIN
    Log.LogDetail( "  Operator-Linemodel: First-Pass - interp FROM large grid z0=%f to zEnd=%f (n=%d)",  redshiftsLargeGrid[0], redshiftsLargeGrid[redshiftsLargeGrid.size()-1], redshiftsLargeGrid.size());
    Log.LogDetail( "  Operator-Linemodel: First-Pass - interp TO fine grid z0=%f to zEnd=%f (n=%d)",  redshiftsFineGrid[0], redshiftsFineGrid[redshiftsFineGrid.size()-1], redshiftsFineGrid.size());

    //initialise and allocate the gsl objects
    //lin
    gsl_interp *interpolation = gsl_interp_alloc (gsl_interp_linear,meritLargeGrid.size());
    gsl_interp_init(interpolation, &(redshiftsLargeGrid.front()), &(meritLargeGrid.front()), meritLargeGrid.size());
    gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();

    for( Int32 j=0; j<redshiftsFineGrid.size(); j++ )
    {
        Float64 Xrebin = redshiftsFineGrid[j];
        if(Xrebin<redshiftsLargeGrid[0] || Xrebin>redshiftsLargeGrid[redshiftsLargeGrid.size()-1])
        {
            continue;
        }
        meritFineGrid[j] = gsl_interp_eval(interpolation, &redshiftsLargeGrid.front(), &meritLargeGrid.front(), Xrebin, accelerator); //lin
    }

    gsl_interp_free(interpolation);
    gsl_interp_accel_free (accelerator);
    //*/

    return 0;

}

void COperatorLineModel::ComputeArea1(CLineModelResult& results)
{
    //prepare p
    Float64 maxp = DBL_MIN;
    CSpectrum pspc;
    CSpectrumFluxAxis& spcFluxAxis = pspc.GetFluxAxis();
    spcFluxAxis.SetSize( results.Redshifts.size() );
    CSpectrumSpectralAxis& spcSpectralAxis = pspc.GetSpectralAxis();
    spcSpectralAxis.SetSize( results.Redshifts.size()  );
    for( Int32 i2=0; i2<results.Redshifts.size(); i2++ )
    {
        if( maxp < results.ChiSquare[i2]){
            maxp = results.ChiSquare[i2];
        }
    }
    for( Int32 i2=0; i2<results.Redshifts.size(); i2++ )
    {
        spcFluxAxis[i2] = exp(-(results.ChiSquare[i2]-maxp)/2.0)-1.0;
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
    Int32 iz0=0;
    for( Int32 i=0; i<results.ExtremaResult.Extrema.size(); i++ )
    {
        //find iz, izmin and izmax
        Int32 izmin= -1;
        Int32 iz= -1;
        Int32 izmax= -1;
        for( Int32 i2=0; i2<results.Redshifts.size(); i2++ )
        {
            if(iz == -1 && (results.ExtremaResult.Extrema[i]) <= results.Redshifts[i2]){
                iz = i2;
                if(i==0){
                    iz0=iz;
                }
            }
            if(izmin == -1 && (results.ExtremaResult.Extrema[i] - winsize/2.0) <= results.Redshifts[i2]){
                izmin = i2;
            }
            if(izmax == -1 && (results.ExtremaResult.Extrema[i] + winsize/2.0) <= results.Redshifts[i2]){
                izmax = i2;
                break;
            }
        }
        Float64 di = abs(results.ChiSquare[iz]-maxp) ;
        Float64 d0 = abs(results.ChiSquare[iz0]-maxp) ;
        if( di < inclusionThresRatio*d0){
            continue;
        }

        /*
        CGaussianFitSimple fitter;
        CGaussianFitSimple::EStatus status = fitter.Compute( pspc, TInt32Range( izmin, izmax ) );
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
        Float64 gaussWidth = FitBayesWidth( spcSpectralAxis, spcFluxAxis, results.ExtremaResult.Extrema[i], izmin, izmax);
        Float64 gaussAmp = spcFluxAxis[iz];

        Float64 intsize = 0.001;
        Float64 area=0.0;
        for( Int32 i2=izmin; i2<izmax; i2++ )
        {
            Float64 x = spcSpectralAxis[i2];
            Float64 Yi = gaussAmp * exp (-1.*(x-results.ExtremaResult.Extrema[i])*(x-results.ExtremaResult.Extrema[i])/(2*gaussWidth*gaussWidth));
            area += Yi;
        }

        //Float64 area = gaussAmp*gaussWidth*sqrt(2.0*3.141592654);
        results.ExtremaResult.LogArea[i] = area;
    }
}

///
/// \brief COperatorLineModel::ComputeArea2
/// computes the Laplace approx for a given Chi2 result around the N best extrema
///
void COperatorLineModel::ComputeArea2(CLineModelResult& results)
{
    Float64 maxp = DBL_MIN;
    for( Int32 i2=0; i2<results.Redshifts.size(); i2++ )
    {
        if( maxp < results.ChiSquare[i2]){
            maxp = results.ChiSquare[i2];
        }
    }
    Float64 winsize = 0.001;
    Float64 inclusionThresRatio = 0.01;
    Int32 iz0=0;
    for( Int32 indz=0; indz<results.ExtremaResult.Extrema.size(); indz++ )
    {
        //find iz, izmin and izmax
        Int32 izmin= -1;
        Int32 iz= -1;
        Int32 izmax= -1;
        for( Int32 i2=0; i2<results.Redshifts.size(); i2++ )
        {
            if(iz == -1 && (results.ExtremaResult.Extrema[indz]) <= results.Redshifts[i2]){
                iz = i2;
                if(indz==0){
                    iz0=iz;
                }
            }
            if(izmin == -1 && (results.ExtremaResult.Extrema[indz] - winsize/2.0) <= results.Redshifts[i2]){
                izmin = i2;
            }
            if(izmax == -1 && (results.ExtremaResult.Extrema[indz] + winsize/2.0) <= results.Redshifts[i2]){
                izmax = i2;
                break;
            }
        }
        Float64 di = abs(results.ChiSquare[iz]-maxp) ;
        Float64 d0 = abs(results.ChiSquare[iz0]-maxp) ;
        if( di < inclusionThresRatio*d0){
            continue;
        }
        if(izmin == -1 || izmax == -1){
            continue;
        }

        //quadratic fit
        int i, n;
        double xi, yi, ei, chisq;
        gsl_matrix *X, *cov;
        gsl_vector *y, *w, *c;

        n = izmax - izmin +1;

        X = gsl_matrix_alloc (n, 3);
        y = gsl_vector_alloc (n);
        w = gsl_vector_alloc (n);

        c = gsl_vector_alloc (3);
        cov = gsl_matrix_alloc (3, 3);

        double x0 = results.ExtremaResult.Extrema[indz];
        for (i = 0; i < n; i++)
        {
            xi = results.Redshifts[i+izmin];
            yi = results.ChiSquare[i+izmin];
            ei = 1.0; //todo, estimate weighting ?
            gsl_matrix_set (X, i, 0, 1.0);
            gsl_matrix_set (X, i, 1, xi-x0);
            gsl_matrix_set (X, i, 2, (xi-x0)*(xi-x0));

            gsl_vector_set (y, i, yi);
            gsl_vector_set (w, i, 1.0/(ei*ei));
        }

        {
          gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, 3);
          gsl_multifit_wlinear (X, w, y, c, cov, &chisq, work);
          gsl_multifit_linear_free (work);
        }

#define C(i) (gsl_vector_get(c,(i)))
#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))

        double zcorr = x0-C(1)/(2.0*C(2));
        double sigma = sqrt(1.0/C(2));
        Float64 a = (Float64)(C(0));
        Float64 b2sur4c = (Float64)(C(1)*C(1)/((Float64)(4.0*C(2))));
        Float64 logK = ( -(a - b2sur4c)/2.0 );
        Float64 logarea = log(sigma) + logK + log(2.0*M_PI);
        if(0){
            Log.LogInfo("Extrema: %g", results.ExtremaResult.Extrema[indz]);
            Log.LogInfo("# best fit: Y = %g + %g X + %g X^2", C(0), C(1), C(2));
	    if( false ) //debug
	      {
		Log.LogInfo("# covariance matrix:\n");
		Log.LogInfo("[ %+.5e, %+.5e, %+.5e  \n", COV(0,0), COV(0,1), COV(0,2));
		Log.LogInfo("  %+.5e, %+.5e, %+.5e  \n", COV(1,0), COV(1,1), COV(1,2));
		Log.LogInfo("  %+.5e, %+.5e, %+.5e ]\n", COV(2,0), COV(2,1), COV(2,2));
	      }
            Log.LogInfo("# chisq/n = %g", chisq/n);
            Log.LogInfo("# zcorr = %g", zcorr);
            Log.LogInfo("# sigma = %g", sigma);
            Log.LogInfo("# logarea = %g", logarea);
            Log.LogInfo("\n");
        }

        gsl_matrix_free (X);
        gsl_vector_free (y);
        gsl_vector_free (w);
        gsl_vector_free (c);
        gsl_matrix_free (cov);

        results.ExtremaResult.LogArea[indz] = logarea;
        results.ExtremaResult.SigmaZ[indz] = sigma;
        results.ExtremaResult.LogAreaCorrectedExtrema[indz] = zcorr;
    }
}

/**
 * \brief Returns a non-negative value for the width that yields the least squared difference between the flux and a exponentially decayed maximum amplitude.
 * Find the maximum flux amplitude. If this not greater than zero, return zero.
 * For each value of c within the range:
 *   Sum the squared difference between the flux and the maximum amplitude with a exponential decay parameterized by c.
 *   Save the minimal result.
 * If the result is not greater than zero, return zero.
 * Return the result.
 **/
Float64 COperatorLineModel::FitBayesWidth( CSpectrumSpectralAxis& spectralAxis, CSpectrumFluxAxis& fluxAxis, Float64 z, Int32 start, Int32 end)
{
    Float64 A = boost::numeric::bounds<float>::lowest();
    const Float64* flux = fluxAxis.GetSamples();
    const Float64* spectral = spectralAxis.GetSamples();
    //const Float64* error = fluxAxis.GetError();

    //A = max, good value ?
    for ( Int32 i = start; i < end; i++)
    {
        Float64 y = flux[i];
        if( y>A )
	  {
            A = y;
	  }
    }

    if( A<=0 )
      {
        return 0.0;
      }
    //c fitting iteration loop
    Float64 mu = z;
    Float64 c = 0.0001;
    Float64 cmax = 0.05;
    Int32 maxIteration = 500;
    Float64 cstepup = (cmax-c)/((Float64)(maxIteration+1));
    Float64 sum2 = boost::numeric::bounds<float>::highest();
    Float64 minsum2 = boost::numeric::bounds<float>::highest();
    Float64 minc = c;
    Int32 icmpt = 0;
    while( icmpt<maxIteration )
      {
        sum2 = 0.0;
        for ( Int32 i = start; i < end; i++)
        {
            Float64 x = spectral[i];
            Float64 Yi = A * exp (-1.*(x-mu)*(x-mu)/(2*c*c));
            sum2 += pow( Yi - flux[i] , 2.0 );
            //sum2 += pow( Yi - flux[i] , 2.0 ) / pow( error[i], 2.0 );
        }
        if(sum2<minsum2){
            minc = c;
            minsum2 = sum2;
        }
        icmpt++;
        c = c+cstepup;
    }

    if(minc<0){
        minc=0;
    }
    return minc;
}
