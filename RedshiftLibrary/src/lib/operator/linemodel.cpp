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
#include <RedshiftLibrary/statistics/pdfz.h>
#include <RedshiftLibrary/operator/pdfMargZLogResult.h>
#include <RedshiftLibrary/linemodel/templatesfitstore.h>

#include <RedshiftLibrary/spectrum/io/fitswriter.h>
#include <RedshiftLibrary/log/log.h>

#include <boost/numeric/conversion/bounds.hpp>
#include "boost/format.hpp"
#include <boost/chrono/thread_clock.hpp>

#include <RedshiftLibrary/processflow/datastore.h>

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
                                  const Float64 &opt_velocityfitmin,
                                  const Float64 &opt_velocityfitmax)
{
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
            Log.LogInfo( "Linemodel: FastFitLargeGrid auto disabled: raw %d redshifts will be calculated", sortedRedshifts.size());
        }else{
            Log.LogInfo( "Linemodel: FastFitLargeGrid enabled: %d redshifts will be calculated on the large grid (%d initially)", largeGridRedshifts.size(), sortedRedshifts.size());
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
                                         opt_rigidity);
    //CTemplateCatalog orthoTplCatalog = tplOrtho.getOrthogonalTplCatalog();
    CTemplatesOrthoStore orthoTplStore = tplOrtho.getOrthogonalTplStore();
    Int32 ctlgIdx = 0; //only one ortho config for now
    std::shared_ptr<CTemplateCatalog> orthoTplCatalog = orthoTplStore.getTplCatalog(ctlgIdx);


    Int32 nResults = sortedRedshifts.size();
    auto result = std::shared_ptr<CLineModelResult>( new CLineModelResult() );
    result->ChiSquare.resize( nResults );
    result->Redshifts.resize( nResults );
    result->Redshifts = sortedRedshifts;
    result->Status.resize( nResults );
    result->restRayList = restRayList;
    result->LineModelSolutions.resize( nResults );

    CLineModelElementList model( spectrum,
                                 spectrumContinuum,
                                 *orthoTplCatalog,
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
    Log.LogInfo( "Linemodel: initialized");

    //setup velocity fitting
    bool enableVelocityFitting = true;
    Float64 velfitMin = opt_velocityfitmin;
    Float64 velfitMax = opt_velocityfitmax;
    if(opt_velocityFitting != "yes"){
        enableVelocityFitting = false;
    }else{
        if(opt_lineWidthType=="combined"){

            Float64 infVelInstr = model.GetVelocityInfFromInstrumentResolution();
            if(velfitMin<infVelInstr){
                velfitMin=infVelInstr;
                Log.LogInfo( "Linemodel: velocity fitting min bound auto set from resolution");
            }
        }

        Log.LogInfo( "Linemodel: velocity fitting bounds: min=%.1f - max=%.1f", velfitMin, velfitMax);
    }

    //fit continuum
    bool enableFitContinuumPrecomputed = true;
    if(enableFitContinuumPrecomputed && opt_continuumcomponent == "tplfit")
    {
        boost::chrono::thread_clock::time_point start_tplfitprecompute = boost::chrono::thread_clock::now();
        Float64 redshiftStep = 1e-3;
        Float64 minRedshift = sortedRedshifts[0];
        Float64 maxRedshift = sortedRedshifts[sortedRedshifts.size()-1];
        Log.LogInfo( "Linemodel: continuum tpl fitting: step=%.5f, min=%.1f, max=%.1f", redshiftStep, minRedshift, maxRedshift);

        CTemplatesFitStore* tplfitStore = new CTemplatesFitStore(minRedshift, maxRedshift, redshiftStep);
        std::vector<Float64> redshiftsTplFit = tplfitStore->GetRedshiftList();
        Log.LogInfo( "Linemodel: continuum tpl redshift list n=%d", redshiftsTplFit.size());
        std::vector<std::shared_ptr<CChisquareResult>> chisquareResultsAllTpl;
        std::vector<std::string> chisquareResultsTplName;

        /*
        COperatorChiSquareLogLambda* chiSquareOperator;
        chiSquareOperator = new COperatorChiSquareLogLambda(opt_calibrationPath); //todo, delete this operator when done...
        //*/

        //*
        COperatorChiSquare2* chiSquareOperator;
        chiSquareOperator = new COperatorChiSquare2(opt_calibrationPath); //todo, delete this operator when done...
        //*/

        std::string opt_interp = "precomputedfinegrid"; // "lin"; //
        Int32 opt_dustFit = 1;
        Int32 opt_extinction = 0;
        Log.LogInfo( "linemodel: precomputing-fitContinuum_dustfit = %d", opt_dustFit );
        Log.LogInfo( "linemodel: precomputing-fitContinuum_igm = %d", opt_extinction );

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
                    Log.LogInfo( "Linemodel failed to compute chi square value for tpl=%s", tpl.GetName().c_str());
                }else{
                    chisquareResultsAllTpl.push_back(chisquareResult);
                    chisquareResultsTplName.push_back(tpl.GetName());
                }
            }
        }
        delete chiSquareOperator;

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
        Log.LogInfo( "Linemodel: tplfit-precompute done in %.4e sec", duration_tplfit_seconds);
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
    Log.LogInfo( "Linemodel: set abs lines limit to %f (ex: -1 means disabled)", absLinesLimit);

    Log.LogInfo( "Linemodel: processing");

    //Set model parameters to FIRST-PASS
    model.setPassMode(1);
    Log.LogInfo( "Linemodel: first-pass mode");

    //WARNING: HACK, first pass with continuum from spectrum.
    //model.SetContinuumComponent("fromspectrum");
    //
    Int32 indexLargeGrid = 0;
    boost::chrono::thread_clock::time_point start_mainloop = boost::chrono::thread_clock::now();
    for (Int32 i=0;i<nResults;i++)
    {
        if(enableFastFitLargeGrid==0 || i==0 || result->Redshifts[i] == largeGridRedshifts[indexLargeGrid])
        {
            ModelFit( model, lambdaRange, result->Redshifts[i], result->ChiSquare[i], result->LineModelSolutions[i], contreest_iterations, false);
            Log.LogDebug( "Z interval %d: Chi2 = %f", i, result->ChiSquare[i] );
            indexLargeGrid++;
            //Log.LogInfo( "\nLineModel Infos: large grid step %d", i);
        }else{
            result->ChiSquare[i] = result->ChiSquare[i-1] + 1e-2;
            result->LineModelSolutions[i] = result->LineModelSolutions[i-1];
        }
    }
    //WARNING: HACK, first pass with continuum from spectrum.
    //model.SetContinuumComponent(opt_continuumcomponent);
    //model.InitFitContinuum();
    //
    boost::chrono::thread_clock::time_point stop_mainloop = boost::chrono::thread_clock::now();
    Float64 duration_mainloop = boost::chrono::duration_cast<boost::chrono::microseconds>(stop_mainloop - start_mainloop).count();
    Float64 duration_firstpass_seconds = duration_mainloop/1e6;
    Log.LogInfo( "Linemodel: first-pass done in %.4e sec", duration_firstpass_seconds);
    Log.LogInfo("<proc-lm-fistpass><%d>", (Int32)duration_firstpass_seconds);

    // extrema
    Int32 extremumCount = opt_extremacount;
    Log.LogDebug( "opt_extremacount = %d", opt_extremacount );
    TPointList extremumList;
    TFloat64Range redshiftsRange(result->Redshifts[0], result->Redshifts[result->Redshifts.size()-1]);
    Log.LogDebug( "redshiftsRange.GetBegin() = %f, redshiftsRange.GetEnd() = %f",
		redshiftsRange.GetBegin(), redshiftsRange.GetEnd() );

    if(result->Redshifts.size() == 1)
    {
        extremumList.push_back(SPoint( result->Redshifts[0], result->ChiSquare[0] ));
        Log.LogInfo("Line Model, only 1 redshift calculated, only 1 extremum");
    }else
    {
        CExtremum extremum( redshiftsRange, extremumCount, true, 2);
        extremum.Find( result->Redshifts, result->ChiSquare, extremumList );
        if( extremumList.size() == 0 )
        {
            Log.LogError("Line Model, Extremum find method failed");
            return result;
        }

        Log.LogInfo( "Linemodel: found %d extrema", extremumList.size());
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
    Float64 extensionradius = 0.005;
    //TPointList extremumListExtended;
    //TBoolList isLocalExtrema;
    for( Int32 i=0; i<extremumList.size(); i++ )
    {
        Log.LogInfo("Linemodel: Raw extr #%d, z_e.X=%f, m_e.Y=%f", i, extremumList[i].X, extremumList[i].Y);
        Float64 x = extremumList[i].X;
        Float64 left_border = max(redshiftsRange.GetBegin(), x-extensionradius);
        Float64 right_border=min(redshiftsRange.GetEnd(), x+extensionradius);

        for (Int32 i=0;i<result->Redshifts.size();i++)
        {
            if(result->Redshifts[i] >= left_border && result->Redshifts[i] <= right_border){
                result->ExtremaExtendedRedshifts.push_back(result->Redshifts[i]);
                //SPoint pt;
                //pt.X = result->Redshifts[i];
                //pt.Y = result->ChiSquare[i];
                //extremumListExtended.push_back(pt);
                //Bool isExtrema = false;
                //if( x == pt.X){
                //    isExtrema = true;
                //}
                //isLocalExtrema.push_back(isExtrema);
            }
        }
    }
    //*/
    //TPointList extremumListExtended = extremumList;
   //todo: remove duplicate redshifts from the extended extrema list


    //Set model parameters to SECOND-PASS
    model.setPassMode(2);
    Log.LogInfo( "Linemodel: second-pass mode");

    std::vector<Float64> extrema_velocityEL;
    std::vector<Float64> extrema_velocityAL;
    extrema_velocityEL.resize(extremumList.size());
    extrema_velocityAL.resize(extremumList.size());

    TPointList extremumList2;
    extremumList2.resize(extremumList.size());

    //recompute the fine grid results around the extrema
    if(enableFastFitLargeGrid==1 || enableVelocityFitting)
    {
        for( Int32 i=0; i<extremumList.size(); i++ )
        {
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
                Log.LogInfo( "Problem. could not find extrema solution index...");
                continue;
            }


            // reestimate the model (eventually with continuum reestimation) on the extrema selected
            if(opt_continuumreest == "always"){
                contreest_iterations = 1;
            }else{
                contreest_iterations  = 0;
            }


            ModelFit( model, lambdaRange, result->Redshifts[idx], result->ChiSquare[idx], result->LineModelSolutions[idx], contreest_iterations, false);
            m = result->ChiSquare[idx];
            if(enableVelocityFitting){
                Bool enableManualStepVelocityFit = true;
                Bool enableLMVelocityFit = false;
                Bool enableLBFGSVelocityFit = false;
                if(enableLMVelocityFit){
                    //fit the emission and absorption width using the linemodel lmfit strategy
                    model.SetFittingMethod("lmfit");
                    model.SetElementIndexesDisabledAuto();
                    Float64 meritTmp;
                    ModelFit( model, lambdaRange, result->Redshifts[idx], meritTmp, result->LineModelSolutions[idx], contreest_iterations, false);


                    model.SetFittingMethod(opt_fittingmethod);
                    model.ResetElementIndexesDisabled();
                    Int32 velocityHasBeenReset = model.ApplyVelocityBound(velfitMin, velfitMax);
                    enableManualStepVelocityFit = velocityHasBeenReset;
                }

                if(enableLBFGSVelocityFit){
                    //fit the emission and absorption width using the linemodel lbfgsfit strategy
                    model.SetFittingMethod("lbfgsfit");
                    //model.SetElementIndexesDisabledAuto();
                    Float64 meritTmp;
                    ModelFit( model, lambdaRange, result->Redshifts[idx], meritTmp, result->LineModelSolutions[idx], contreest_iterations, false);


                    model.SetFittingMethod(opt_fittingmethod);
                    model.ResetElementIndexesDisabled();
                    Int32 velocityHasBeenReset = model.ApplyVelocityBound(velfitMin, velfitMax);
                    enableManualStepVelocityFit = velocityHasBeenReset;
                }

                if(enableManualStepVelocityFit){
                    //fit the emission and absorption width by minimizing the linemodel merit with linemodel "hybrid" fitting method
                    model.SetFittingMethod("hybrid");
                    Float64 vInfLim = velfitMin;
                    Float64 vSupLim = velfitMax;
                    Float64 vStep = 20.0;
                    Int32 nSteps = (int)((vSupLim-vInfLim)/vStep);


                    Float64 dzInfLim = -4e-4;
                    Float64 dzSupLim = 4e-4;
                    Float64 dzStep = 2e-4;
                    Int32 nDzSteps = (int)((dzSupLim-dzInfLim)/dzStep);

                    for(Int32 iLineType = 0; iLineType<2; iLineType++)
                    {
                        if(iLineType==0)
                        {
                            Log.LogInfo( "\nLineModel Infos: manualStep velocity fit ABSORPTION, for z = %.4f", result->Redshifts[idx]);

                        }else{
                            Log.LogInfo( "\nLineModel Infos: manualStep velocity fit EMISSION, for z = %.4f", result->Redshifts[idx]);

                        }

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
                                    model.SetVelocityAbsorption(vTest);
                                }else
                                {
                                    model.SetVelocityEmission(vTest);
                                }
                                Float64 meritv;
                                ModelFit( model, lambdaRange, result->Redshifts[idx]+dzTest, meritv, result->LineModelSolutions[idx], contreest_iterations, false);
                                //meritv = model.getLeastSquareMeritUnderElements();
                                //todo: eventually use the merit under the elements with a BIC estimator taking the n samples in the support for each velocity solution...

                                //Log.LogInfo("LineModel Solution: testing velocity: merit=%.1f for velocity = %.1f", meritv, vTest);
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
                            Log.LogInfo( "LineModel Solution: best Velocity found = %.1f", vOptim);
                            result->ChiSquare[idx] =  meritMin;
                            if(iLineType==0)
                            {
                                model.SetVelocityAbsorption(vOptim);
                            }else
                            {
                                model.SetVelocityEmission(vOptim);
                            }

                        }
                    }
                    model.SetFittingMethod(opt_fittingmethod);
                }
            }
            extrema_velocityEL[i]=model.GetVelocityEmission();
            extrema_velocityAL[i]=model.GetVelocityAbsorption();

            //finally compute the redshifts on the extended ExtremaExtendedRedshifts values
            Float64 left_border = max(redshiftsRange.GetBegin(), z-extensionradius);
            Float64 right_border=min(redshiftsRange.GetEnd(), z+extensionradius);
            //model.SetFittingMethod("nofit");
            extremumList2[i].Y = DBL_MAX;
            extremumList2[i].X = result->Redshifts[idx];
            Int32 idx2 = idx;
            for (Int32 iz=0;iz<result->Redshifts.size();iz++)
            {
                if(result->Redshifts[iz] >= left_border && result->Redshifts[iz] <= right_border){
                    ModelFit( model, lambdaRange, result->Redshifts[iz], result->ChiSquare[iz], result->LineModelSolutions[iz], contreest_iterations, false);
                    if(result->ChiSquare[iz]< extremumList2[i].Y)
                    {
                        extremumList2[i].X = result->Redshifts[iz];
                        extremumList2[i].Y = result->ChiSquare[iz];
                        idx2 = iz;
                    }
                }
            }
            //model.SetFittingMethod(opt_fittingmethod);

            Log.LogInfo("Linemodel: Recomputed extr #%d, idx=%d, z_e.X=%f, m_e.Y=%f", i, idx2, extremumList2[i].X, extremumList2[i].Y);

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
    }
    if( extremumListOrdered.size() == 0 )
    {
        Log.LogError("Line Model, Extremum Ordering failed");
    }

    Log.LogInfo("Line Model, Store extrema results");
    // store extrema results
    result->ResizeExtremaResults(extremumCount);

    Int32 start = spectrum.GetSpectralAxis().GetIndexAtWaveLength(lambdaRange.GetBegin());
    Int32 end = spectrum.GetSpectralAxis().GetIndexAtWaveLength(lambdaRange.GetEnd());
    Int32 nsamples = end - start + 1;
    Int32 savedModels = 0;
    m_savedModelSpectrumResults.clear();
    m_savedModelFittingResults.clear();
    m_savedModelRulesResults.clear();

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
        Log.LogInfo("Linemodel: Saving extr #%d, idx=%d, z=%f, m=%f", i, idx, result->Redshifts[idx], result->ChiSquare[idx]);


        // reestimate the model (eventually with continuum reestimation) on the extrema selected
        if(opt_continuumreest == "always" || opt_continuumreest == "onlyextrema"){
            contreest_iterations = 10; //4
        }else{
            contreest_iterations  = 0;
        }


        if(enableVelocityFitting){
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

        ModelFit( model, lambdaRange, result->Redshifts[idx], result->ChiSquare[idx], result->LineModelSolutions[idx], contreest_iterations, true);
        if(m!=result->ChiSquare[idx])
        {
            Log.LogInfo("Linemodel: m (%f for idx=%d) !=chi2 (%f) ", m, idx, result->ChiSquare[idx]);
        }
        m = result->ChiSquare[idx];//result->ChiSquare[idx];

        //save the model result
        static Int32 maxModelSave = std::min(m_maxModelSaveCount, extremumCount);
        if( savedModels<maxModelSave /*&& isLocalExtrema[i]*/)
        {
            // CModelSpectrumResult
            std::shared_ptr<CModelSpectrumResult>  resultspcmodel = std::shared_ptr<CModelSpectrumResult>( new CModelSpectrumResult(model.GetModelSpectrum()) );
            //std::shared_ptr<CModelSpectrumResult>  resultspcmodel = std::shared_ptr<CModelSpectrumResult>( new CModelSpectrumResult(model.GetObservedSpectrumWithLinesRemoved()) );
            m_savedModelSpectrumResults.push_back(resultspcmodel);

            // CModelFittingResult
            std::shared_ptr<CModelFittingResult>  resultfitmodel = std::shared_ptr<CModelFittingResult>( new CModelFittingResult(result->LineModelSolutions[idx], result->Redshifts[idx], result->ChiSquare[idx], result->restRayList, model.GetVelocityEmission(), model.GetVelocityAbsorption()) );
            m_savedModelFittingResults.push_back(resultfitmodel);

            // CModelRulesResult
            std::shared_ptr<CModelRulesResult>  resultrulesmodel = std::shared_ptr<CModelRulesResult>( new CModelRulesResult( model.GetModelRulesLog() ));
            m_savedModelRulesResults.push_back(resultrulesmodel);

            Int32 saveNLinemodelContinua = 1;
            if( savedModels < saveNLinemodelContinua && contreest_iterations>0)
            {
                // Save the reestimated continuum, only the first extrema
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

                std::string nameBaselineStr = (boost::format("linemodel_continuum_extrema_%1%") % savedModels).str();
                dataStore.StoreScopedGlobalResult(nameBaselineStr.c_str(), baselineResult);
            }

            savedModels++;
        }


        result->Extrema[i] = z;
        result->ExtremaMerit[i] = m;

        result->ExtremaLastPass[i] = z; //refined extremum is initialized here.

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
                Log.LogError("Linemodel: Deltaz computation failed");
            }
        }
        result->DeltaZ[i] = dz;

        //store the model norm
        result->mTransposeM[i] = model.EstimateMTransposeM(lambdaRange);

        //result->IsLocalExtrema[i]=isLocalExtrema[i];

        static Float64 cutThres = 3.0;
        Int32 nValidLines = result->GetNLinesOverCutThreshold(i, cutThres, cutThres);
        result->Posterior[i] = nValidLines;//m/Float64(1+nValidLines);
        Float64 cumulStrongELSNR = model.getCumulSNRStrongEL(); //getStrongerMultipleELAmpCoeff();
        result->StrongELSNR[i] = cumulStrongELSNR;

        result->LogArea[i] = -DBL_MAX;
        result->LogAreaCorrectedExtrema[i] = -1.0;


        Int32 nddl = model.GetNElements(); //get the total number of elements in the model
        nddl = result->LineModelSolutions[idx].nDDL; //override nddl by the actual number of elements in the fitted model

        //result->bic[i] = m + nddl*log(nsamples); //BIC
        Float64 aic = m + 2*nddl; //AIC
        result->bic[i] = aic;
        //result->bic[i] = aic + (2*nddl*(nddl+1) )/(nsamples-nddl-1);  //AICc, better when nsamples small

        //compute continuum indexes
        CContinuumIndexes continuumIndexes;
        result->ContinuumIndexes[i] = continuumIndexes.getIndexes( spectrum, z );

        //save the outsideLinesMask
        result->OutsideLinesMask[i] = model.getOutsideLinesMask();

        //save the continuum tpl fitting results
        result->FittedTplName[i] = model.getFitContinuum_tplName();
        result->FittedTplAmplitude[i] = model.getFitContinuum_tplAmplitude();
        result->FittedTplDustCoeff[i] = model.getFitContinuum_tplIsmDustCoeff();
        result->FittedTplMeiksinIdx[i] = model.getFitContinuum_tplIgmMeiksinIdx();

        //save the tplcorr results
        result->FittedTplcorrTplName[i] = model.getTplCorr_bestTplName();
    }

    //ComputeArea2(*result);

    /* ------------------------  COMPUTE POSTMARG PDF  --------------------------  */
    auto postmargZResult = std::shared_ptr<CPdfMargZLogResult>(new CPdfMargZLogResult());
    CPdfz pdfz;
    Float64 cstLog = model.getLikelihood_cstLog(lambdaRange);
    TFloat64List logProba;
    Float64 logEvidence;
    Int32 retPdfz = pdfz.Compute(result->ChiSquare, result->Redshifts, cstLog, logProba, logEvidence);
    if(retPdfz!=0)
    {
        Log.LogError("Linemodel: Pdfz computation failed");
    }else{
        postmargZResult->countTPL = result->Redshifts.size(); // assumed 1 model per z
        postmargZResult->Redshifts.resize(result->Redshifts.size());
        postmargZResult->valProbaLog.resize(result->Redshifts.size());
        for ( UInt32 k=0; k<result->Redshifts.size(); k++)
        {
            postmargZResult->Redshifts[k] = result->Redshifts[k] ;
            postmargZResult->valProbaLog[k] = logProba[k];
        }
        dataStore.StoreGlobalResult( "zPDF/logposterior.logMargP_Z_data", postmargZResult); //need to store this pdf with this exact same name so that zqual can load it. see zqual.cpp/ExtractFeaturesPDF
    }

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
       Log.LogInfo("Linemodel - Last Pass: begin");
       //
       //do the last pass on the 1st extremum range
       //
       Float64 halfRange = 1e-3;
       Float64 lastPassStep = 1e-4;
       //find the best extremum redhift
       Float64 bestRedshift=-1;
       Float64 bestMerit=DBL_MAX;
       Int32 bestIndex=-1;
       for(Int32 k=0; k<result->Extrema.size(); k++)
       {
           Log.LogInfo("Linemodel - Last Pass: k = %d", k);
           Log.LogInfo("Linemodel - Last Pass: result->Extrema[k] = %.5f", result->Extrema[k]);
           Log.LogInfo("Linemodel - Last Pass: result->ExtremaMerit[k] = %.5f", result->ExtremaMerit[k]);

            if(bestMerit>result->ExtremaMerit[k])
            {
                bestMerit = result->ExtremaMerit[k];
                bestRedshift = result->Extrema[k];
                bestIndex = k;
            }
       }
       Log.LogInfo("Linemodel - Last Pass: around extrema z = %.5f", bestRedshift);
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
       Log.LogInfo("Linemodel - Last Pass: range zmin=%.5f, zmax=%.5f", lastPassRedshifts[0], lastPassRedshifts[lastPassRedshifts.size()-1]);


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
        Float64 refinedExtremum = lastPassResult->Extrema[0];
        Log.LogInfo("Linemodel - Last Pass: found refined z=%.5f", refinedExtremum);

        result->ExtremaLastPass[bestIndex] = refinedExtremum;
    }else{
        Log.LogInfo("Linemodel - Last Pass: failed to do linemodel, rigidity=%s", opt_rigidity.c_str());
    }

    return result;
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
    for( Int32 i=0; i<results.Extrema.size(); i++ )
    {
        //find iz, izmin and izmax
        Int32 izmin= -1;
        Int32 iz= -1;
        Int32 izmax= -1;
        for( Int32 i2=0; i2<results.Redshifts.size(); i2++ )
        {
            if(iz == -1 && (results.Extrema[i]) <= results.Redshifts[i2]){
                iz = i2;
                if(i==0){
                    iz0=iz;
                }
            }
            if(izmin == -1 && (results.Extrema[i] - winsize/2.0) <= results.Redshifts[i2]){
                izmin = i2;
            }
            if(izmax == -1 && (results.Extrema[i] + winsize/2.0) <= results.Redshifts[i2]){
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
        Float64 gaussWidth = FitBayesWidth( spcSpectralAxis, spcFluxAxis, results.Extrema[i], izmin, izmax);
        Float64 gaussAmp = spcFluxAxis[iz];

        Float64 intsize = 0.001;
        Float64 area=0.0;
        for( Int32 i2=izmin; i2<izmax; i2++ )
        {
            Float64 x = spcSpectralAxis[i2];
            Float64 Yi = gaussAmp * exp (-1.*(x-results.Extrema[i])*(x-results.Extrema[i])/(2*gaussWidth*gaussWidth));
            area += Yi;
        }

        //Float64 area = gaussAmp*gaussWidth*sqrt(2.0*3.141592654);
        results.LogArea[i] = area;
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
    for( Int32 indz=0; indz<results.Extrema.size(); indz++ )
    {
        //find iz, izmin and izmax
        Int32 izmin= -1;
        Int32 iz= -1;
        Int32 izmax= -1;
        for( Int32 i2=0; i2<results.Redshifts.size(); i2++ )
        {
            if(iz == -1 && (results.Extrema[indz]) <= results.Redshifts[i2]){
                iz = i2;
                if(indz==0){
                    iz0=iz;
                }
            }
            if(izmin == -1 && (results.Extrema[indz] - winsize/2.0) <= results.Redshifts[i2]){
                izmin = i2;
            }
            if(izmax == -1 && (results.Extrema[indz] + winsize/2.0) <= results.Redshifts[i2]){
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

        double x0 = results.Extrema[indz];
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
            Log.LogInfo("Extrema: %g", results.Extrema[indz]);
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

        results.LogArea[indz] = logarea;
        results.SigmaZ[indz] = sigma;
        results.LogAreaCorrectedExtrema[indz] = zcorr;
    }
}

/**
 * \brief Calls model.fit() and sets the chisquare to the return value of that method.
 **/
Void COperatorLineModel::ModelFit(CLineModelElementList& model, const TFloat64Range& lambdaRange, Float64 redshift,
                  Float64& chiSquare, CLineModelResult::SLineModelSolution& modelSolution, Int32 contreest_iterations, bool enableLogging)
{
    chiSquare = boost::numeric::bounds<float>::highest();
    Float64 fit = model.fit( redshift, lambdaRange, modelSolution, contreest_iterations, enableLogging );
    chiSquare = fit;
    Log.LogDebug( "ModelFit: Chi2 = %f", fit );
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

