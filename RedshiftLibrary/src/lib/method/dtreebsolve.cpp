#include <RedshiftLibrary/method/dtreebsolve.h>
#include <RedshiftLibrary/method/dtreebsolveresult.h>

#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/operator/correlation.h>
#include <RedshiftLibrary/operator/chisquare.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/processflow/datastore.h>

#include <RedshiftLibrary/operator/correlation.h>
#include <RedshiftLibrary/operator/chicorr.h>
#include <RedshiftLibrary/method/blindsolveresult.h>

#include <RedshiftLibrary/method/blindsolve.h>
#include <RedshiftLibrary/method/blindsolveresult.h>
#include <RedshiftLibrary/operator/chisquare.h>
#include <RedshiftLibrary/operator/chisquare2.h>

#include <RedshiftLibrary/operator/peakdetection.h>
#include <RedshiftLibrary/operator/peakdetectionresult.h>
#include <RedshiftLibrary/operator/raydetection.h>
#include <RedshiftLibrary/operator/raydetectionresult.h>
#include <RedshiftLibrary/operator/raymatching.h>
#include <RedshiftLibrary/operator/raymatchingresult.h>

#include <RedshiftLibrary/method/chisquare2solve.h>
#include <RedshiftLibrary/method/correlationsolve.h>
#include <RedshiftLibrary/method/linematchingsolve.h>

#include <RedshiftLibrary/operator/linemodel.h>


#include <float.h>

using namespace NSEpic;
using namespace std;


CMethodDTreeBSolve::CMethodDTreeBSolve( std::string calibrationPath )
{
    m_calibrationPath = calibrationPath;
}

CMethodDTreeBSolve::~CMethodDTreeBSolve()
{

}

const std::string CMethodDTreeBSolve::GetDescription()
{
    std::string desc;

    desc = "Method Amazed0_2:\n";

    desc.append("\tparam: linemodel.linetypefilter = {""no"", ""E"", ""A""}\n");
    desc.append("\tparam: linemodel.lineforcefilter = {""no"", ""S""}\n");
    desc.append("\tparam: linemodel.fittingmethod = {""hybrid"", ""individual""}\n");
    desc.append("\tparam: linemodel.continuumcomponent = {""fromspectrum"", ""nocontinuum"", ""zero""}\n");
    desc.append("\tparam: linemodel.linewidthtype = {""instrumentdriven"", ""velocitydriven"",  ""combined"",  ""nispvsspsf201707"", ""fixed""}\n");
    desc.append("\tparam: linemodel.instrumentresolution = <float value>\n");
    desc.append("\tparam: linemodel.velocityemission = <float value>\n");
    desc.append("\tparam: linemodel.velocityabsorption = <float value>\n");

    desc.append("\tparam: linemodel.rules = {""all"", ""no""}\n");
    desc.append("\tparam: linemodel.continuumreestimation = {""no"", ""onlyextrema"", ""always""}\n");
    desc.append("\tparam: linemodel.extremacount = <float value>\n");
    desc.append("\tparam: linemodel.firstpass.largegridstep = <float value>, deactivated if negative or zero\n");

    desc.append("\tparam: chisquare.overlapthreshold = <float value>\n");
    desc.append("\tparam: chisquare.redshiftsupport = {""full"", ""extremaextended""}\n");
    desc.append("\tparam: chisquare.interpolation = {""precomputedfinegrid"", ""lin""}\n");
    desc.append("\tparam: chisquare.extinction = {""no"", ""yes""}\n");

    return desc;

}

std::shared_ptr<CDTreeBSolveResult> CMethodDTreeBSolve::Compute( CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                                        const CTemplateCatalog& tplCatalog, const TStringList &tplCategoryList, const CRayCatalog &restRayCatalog,
                                                        const TFloat64Range& lambdaRange, const TFloat64List &redshifts,
                                                        const Float64 radius)
{
    Bool storeResult = false;
    m_radius = radius;
    CDataStore::CAutoScope resultScope( resultStore, "dtreeBsolve" );

    storeResult = Solve(resultStore, spc, spcWithoutCont,
                                       tplCatalog, tplCategoryList, restRayCatalog,
                                       lambdaRange, redshifts );

    //storeResult = true;
    if( storeResult )
    {
        return std::shared_ptr<CDTreeBSolveResult>( new CDTreeBSolveResult() );
    }

    return NULL;
}

Bool CMethodDTreeBSolve::Solve(CDataStore &dataStore, const CSpectrum &spc, const CSpectrum &spcWithoutCont, const CTemplateCatalog &tplCatalog, const TStringList &tplCategoryList, const CRayCatalog &restRayCatalog, const TFloat64Range &lambdaRange, const TFloat64List &redshifts)
{
    CSpectrum _spcContinuum = spc;
    CSpectrumFluxAxis spcfluxAxis = _spcContinuum.GetFluxAxis();
    spcfluxAxis.Subtract(spcWithoutCont.GetFluxAxis());
    CSpectrumFluxAxis& sfluxAxisPtr = _spcContinuum.GetFluxAxis();
    sfluxAxisPtr = spcfluxAxis;

    std::string opt_linetypefilter;
    dataStore.GetScopedParam( "linemodel.linetypefilter", opt_linetypefilter, "no" );
    std::string opt_lineforcefilter;
    dataStore.GetScopedParam( "linemodel.lineforcefilter", opt_lineforcefilter, "no" );
    std::string opt_fittingmethod;
    dataStore.GetScopedParam( "linemodel.fittingmethod", opt_fittingmethod, "hybrid" );
    std::string opt_continuumcomponent;
    //dataStore.GetScopedParam( "linemodel.continuumcomponent", opt_continuumcomponent, "nocontinuum" );
    dataStore.GetScopedParam( "linemodel.continuumcomponent", opt_continuumcomponent, "fromspectrum" );
    std::string opt_lineWidthType;
    dataStore.GetScopedParam( "linemodel.linewidthtype", opt_lineWidthType, "combined" );
    Float64 opt_resolution;
    dataStore.GetScopedParam( "linemodel.instrumentresolution", opt_resolution, 2350.0 );
    Float64 opt_velocity_emission;
    dataStore.GetScopedParam( "linemodel.velocityemission", opt_velocity_emission, 100.0 );
    Float64 opt_velocity_absorption;
    dataStore.GetScopedParam( "linemodel.velocityabsorption", opt_velocity_absorption, 300.0 );
    std::string opt_velocityfit;
    dataStore.GetScopedParam( "linemodel.velocityfit", opt_velocityfit, "no" );
    std::string opt_continuumreest;
    dataStore.GetScopedParam( "linemodel.continuumreestimation", opt_continuumreest, "no" );
    std::string opt_rules;
    dataStore.GetScopedParam( "linemodel.rules", opt_rules, "all" );
    Float64 opt_extremacount;
    dataStore.GetScopedParam( "linemodel.extremacount", opt_extremacount, 10.0 );
    Float64 opt_twosteplargegridstep;
    dataStore.GetScopedParam( "linemodel.firstpass.largegridstep", opt_twosteplargegridstep, 0.001 );


    //_///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Compute Linemodel
    COperatorLineModel linemodel;
    auto result = dynamic_pointer_cast<CLineModelResult>(linemodel.Compute(dataStore,
                                                                           spc,
                                                                           _spcContinuum,
                                                                           tplCatalog,
                                                                           tplCategoryList,
                                                                           m_calibrationPath,
                                                                           restRayCatalog,
                                                                           opt_linetypefilter,
                                                                           opt_lineforcefilter,
                                                                           lambdaRange,
                                                                           redshifts,
                                                                           opt_extremacount,
                                                                           opt_fittingmethod,
                                                                           opt_continuumcomponent,
                                                                           opt_lineWidthType,
                                                                           opt_resolution,
                                                                           opt_velocity_emission,
                                                                           opt_velocity_absorption,
                                                                           opt_continuumreest,
                                                                           opt_rules,
                                                                           opt_velocityfit,
                                                                           opt_twosteplargegridstep) );

    /*
    static Float64 cutThres = 2.0;
    static Int32 bestSolutionIdx = 0;
    Int32 nValidLines = result->GetNLinesOverCutThreshold(bestSolutionIdx, cutThres, cutThres);
    Float64 bestExtremaMerit = result->GetExtremaMerit(0);
    Log.LogInfo( "Linemodelsolve : bestExtremaMerit, %f", bestExtremaMerit);
    Float64 thres = 0.001;
    Int32 idxNextValid = 1;
    for(Int32 idnext=1; idnext<result->Redshifts.size(); idnext++){
       if( std::abs(result->Redshifts[idnext]-result->Redshifts[0])> thres){
           idxNextValid = idnext;
           break;
       }
    }
    Float64 nextExtremaMerit = result->GetExtremaMerit(idxNextValid);
    Log.LogInfo( "Linemodelsolve : nextExtremaMerit, %f", nextExtremaMerit);
    //    if(nValidLines<2 || (bestExtremaMerit - nextExtremaMerit) > -50.0 ){
    //        result=0;
    //        Log.LogInfo( "Linemodelsolve : result set to 0" );
    //    }
    Log.LogInfo( "Linemodelsolve : for best solution, %d valid lines found", nValidLines);
    //*/

    if( !result )
    {
        //Log.LogInfo( "Failed to compute linemodel");
        return false;
    }else{
        // Store results
        dataStore.StoreScopedGlobalResult( "linemodel", result );
        //save linemodel fitting and spectrum-model results
        linemodel.storeGlobalModelResults(dataStore);
    }


    //*
    //_///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // compute chisquares
    if( result->ExtremaResult.Extrema.size() == 0 )
    {
        return false;
    }

    Float64 overlapThreshold;
    dataStore.GetScopedParam( "chisquare.overlapthreshold", overlapThreshold, 1.0 );
    std::string opt_redshiftsupport;
    dataStore.GetScopedParam( "chisquare.redshiftsupport", opt_redshiftsupport, "full" );
    // Compute merit function
//    TFloat64List extremumRedshifts( result->Extrema.size() );
//    for( Int32 i=0; i<result->Extrema.size(); i++ )
//    {
//        extremumRedshifts[i] = result->Extrema[i];
//    }
    std::string opt_interp;
    dataStore.GetScopedParam( "chisquare.interpolation", opt_interp, "precomputedfinegrid" );
    std::string opt_extinction;
    dataStore.GetScopedParam( "chisquare.extinction", opt_extinction, "no" );
    std::string opt_dustFit;
    dataStore.GetScopedParam( "chisquare.dustfit", opt_dustFit, "no" );

    TFloat64List redshiftsChi2;
    if(opt_redshiftsupport == "full"){
        redshiftsChi2 = redshifts;
    }else if(opt_redshiftsupport == "extremaextended"){
        //transform std::vector<TFloat64List> into a TFloat64List
        for(Int32 i = 0; i <result->ExtremaResult.ExtremaExtendedRedshifts.size(); i++){ 
            for(Int32 j = 0; j <result->ExtremaResult.ExtremaExtendedRedshifts[i].size(); j++){
                redshiftsChi2.push_back(result->ExtremaResult.ExtremaExtendedRedshifts[i][j]);
            }
        }
        //redshiftsChi2= result->ExtremaResult.ExtremaExtendedRedshifts;
    }else{
        redshiftsChi2= result->ExtremaResult.Extrema;
    }

    std::string spcComponent = "nocontinuum";

    // prepare the unused masks
    std::vector<CMask> maskList;

    CMethodChisquare2Solve chiSolve(m_calibrationPath);
    //achtung: override overlap threshold for chi2nocontinuum:
    Float64 overlapThresholdNC=overlapThreshold;
    overlapThresholdNC = 1.0;
    auto chisolveResultnc = chiSolve.Compute( dataStore, spc, spcWithoutCont,
                                                                        tplCatalog, tplCategoryList,
                                                                        lambdaRange, redshiftsChi2, overlapThresholdNC, maskList, spcComponent, m_radius, opt_interp, opt_extinction, opt_dustFit);
    if( chisolveResultnc ) {
        dataStore.StoreScopedGlobalResult( "redshiftresult", chisolveResultnc );
    }

    spcComponent = "continuum";
    Int32 enableFastContinuumFitLargeGrid = 1;
    std::vector<Float64> redshiftsChi2Continuum;
    //*
    if(enableFastContinuumFitLargeGrid==1){
        //calculate on a wider grid, defined by a minimum step
        Float64 dz_chi2c_thres = 1e-3;
        std::vector<Int32> removed_inds;
        Int32 lastKeptInd = 0;
        for(Int32 i=1; i<redshiftsChi2.size()-1; i++){
            if( abs(redshiftsChi2[i]-redshiftsChi2[lastKeptInd])<dz_chi2c_thres && abs(redshiftsChi2[i+1]-redshiftsChi2[i])<dz_chi2c_thres ){
                removed_inds.push_back(i);
            }else{
                lastKeptInd=i;
            }
        }
        Int32 rmInd = 0;
        for(Int32 i=1; i<redshiftsChi2.size()-1; i++){
            if(removed_inds[rmInd]==i){
                rmInd++;
            }else{
                redshiftsChi2Continuum.push_back(redshiftsChi2[i]);
            }
        }
    }else{
        redshiftsChi2Continuum=redshiftsChi2;
    }

    //*/
    auto chisolveResultcontinuum = chiSolve.Compute( dataStore, spc, spcWithoutCont,
                                                                        tplCatalog, tplCategoryList,
                                                                        lambdaRange, redshiftsChi2Continuum, overlapThreshold, maskList, spcComponent, m_radius, opt_interp, opt_extinction, opt_dustFit);
    //_///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //*/

    // Calculate the Combination //////////////////////////////////////////////////
    GetCombinedRedshift(dataStore);

    // /////////////////////////////////////////////////////////////////////////////


    return true;
}


Bool CMethodDTreeBSolve::GetCombinedRedshift(CDataStore& store)
{
    std::shared_ptr<CChisquareResult> resultCombined = std::shared_ptr<CChisquareResult>( new CChisquareResult() );

    std::string scope = "dtreeBsolve.linemodel";
    auto results = std::dynamic_pointer_cast<const CLineModelResult>( store.GetGlobalResult(scope.c_str()).lock() );

    //*
    //***********************************************************
    //retrieve chi2nc values
    std::vector<Float64> w_chi2nc;
    TFloat64List chi2nc;
    TFloat64List znc;
    Float64 minchi2nc= DBL_MAX;

    chi2nc = GetBestRedshiftChi2List(store, "chisquare_nocontinuum", minchi2nc, znc);

    //Log.LogInfo( "dtreeBsolve : znc size=%d", znc.size());
    for( Int32 i=0; i<chi2nc.size(); i++ )
    {
        w_chi2nc.push_back(chi2nc[i]/minchi2nc);
    }

    //*
    //***********************************************************
    //retrieve chi2 continuum values
    std::vector<Float64> w_chi2continuum;
    TFloat64List chi2continuum;
    TFloat64List chi2continuum_calcGrid;

    TFloat64List zcontinuum_calcGrid;
    Float64 minchi2continuum= DBL_MAX;

    chi2continuum_calcGrid = GetBestRedshiftChi2List(store, "chisquare_continuum", minchi2continuum, zcontinuum_calcGrid);

    //todo: interpolate the continuum results on the continuum fit grid
    Int32 enableFastContinuumFitLargeGrid = 1;
    if(enableFastContinuumFitLargeGrid){
        Int32 iCalcGrid=0;
        for( Int32 i=0; i<znc.size(); i++ )
        {
            if((iCalcGrid<zcontinuum_calcGrid.size() && znc[i]>=zcontinuum_calcGrid[iCalcGrid]) || iCalcGrid==0){
                chi2continuum.push_back(chi2continuum_calcGrid[iCalcGrid]);
                iCalcGrid++;
            }else if(iCalcGrid>=zcontinuum_calcGrid.size()){
                chi2continuum.push_back(chi2continuum_calcGrid[zcontinuum_calcGrid.size()-1]);
            }else{
                //interpolate between iCalcGrid-1 and iCalcGrid
                Float64 a = (chi2continuum_calcGrid[iCalcGrid]-chi2continuum_calcGrid[iCalcGrid-1])/(zcontinuum_calcGrid[iCalcGrid]-zcontinuum_calcGrid[iCalcGrid-1]);
                Float64 continuumVal = chi2continuum_calcGrid[iCalcGrid-1]+a*(znc[i]-zcontinuum_calcGrid[iCalcGrid-1]);
                chi2continuum.push_back(continuumVal);
            }
        }
    }else{
        chi2continuum = chi2continuum_calcGrid;
    }
    for( Int32 i=0; i<chi2continuum.size(); i++ )
    {
        w_chi2continuum.push_back(chi2continuum[i]/minchi2continuum);
    }


    //*
    //***********************************************************
    //retrieve linemodel values
    TFloat64List chi2lm;

    for( Int32 i=0; i<znc.size(); i++ )
    {
        for( Int32 iall=0; iall<results->Redshifts.size(); iall++ )
        {
            if(results->Redshifts[iall] == znc[i]){
                chi2lm.push_back(results->ChiSquare[iall]);
                break;
            }
        }
    }
    //Log.LogInfo( "dtreeBsolve : chi2lm size=%d", chi2lm.size());

    //save the combined chisquare result
    resultCombined->ChiSquare.resize( znc.size() );
    resultCombined->Redshifts.resize( znc.size() );
    resultCombined->Overlap.resize( znc.size() );

    // compute the average log likelihood
    for( Int32 i=0; i<znc.size(); i++ )
    {
        chi2lm[i] = chi2lm[i]/(float)results->nSpcSamples;
        chi2nc[i] = chi2nc[i]/(float)results->nSpcSamples;
        chi2continuum[i] = chi2continuum[i]/(float)results->nSpcSamples;
    }

    //*
    //***********************************************************
    // Estimate the combined chisquare result
    Int32 combineOption = 0; //1 = bayes, 0 = linear

    //Set the coefficients
    Float64 lmCoeff = 1.0;
    Float64 chi2ncCoeff = 1.0;
    Float64 chi2cCoeff = 1.0;
    if(combineOption==1) //bayesian combination
    {
        lmCoeff = 3.981e-6;
        chi2ncCoeff = 0.00316;
        chi2cCoeff = 3.981e-9;

        //logcoeff
        lmCoeff = -2.0*log(lmCoeff);
        chi2ncCoeff = -2.0*log(chi2ncCoeff);
        chi2cCoeff = -2.0*log(chi2cCoeff);

        Log.LogInfo( "dtreeBsolve : lmCoeff=%f", lmCoeff);
        Log.LogInfo( "dtreeBsolve : chi2ncCoeff=%f", chi2ncCoeff);
        Log.LogInfo( "dtreeBsolve : chi2cCoeff=%f", chi2cCoeff);
    }
    else //linear combination
    {
        //*
        // coeffs for PFS simu dec 2015.
        lmCoeff = 75.0;
        chi2ncCoeff = 100.0;
        chi2cCoeff = 50.0;
        //*/
        /*
        // coeffs for VVDS deep, udeep and VUDS
        lmCoeff = 12.0;
        chi2ncCoeff = 100.0;
        chi2cCoeff = 8.0;
        //*/
    }

    //Rescale if necessary
    Float64 minChi2WithCoeffs = DBL_MAX; //overall chi2 minimum for rescaling
    if(combineOption==1) //bayesian combination
    {
        std::string methodForMin = "";
        //find the rescaling value
        for( Int32 i=0; i<znc.size(); i++ )
        {
            if(chi2lm[i]+lmCoeff < minChi2WithCoeffs)
            {
                minChi2WithCoeffs = chi2lm[i]+lmCoeff;
                methodForMin = "lm";
            }
            if(chi2nc[i]+chi2ncCoeff < minChi2WithCoeffs)
            {
                minChi2WithCoeffs = chi2nc[i]+chi2ncCoeff;
                methodForMin = "chi2nc";
            }
            if(chi2continuum[i]+chi2cCoeff < minChi2WithCoeffs)
            {
                minChi2WithCoeffs = chi2continuum[i]+chi2cCoeff;
                methodForMin = "chi2c";
            }
        }
        Log.LogInfo( "dtreeBsolve : minChi2WithCoeffs=%f", minChi2WithCoeffs);
        Log.LogInfo( "dtreeBsolve : methodForMin=%s", methodForMin.c_str());
    }

    //Do the actual combination
    for( Int32 i=0; i<znc.size(); i++ )
    {
        Float64 post = DBL_MAX;
        if(combineOption==1) //bayesian combination
        {
            Float64 valf = 1.0;
            valf = chi2lm[i]+lmCoeff-minChi2WithCoeffs;
            Float64 lmLikelihood = exp(-valf/2.0);
            valf = chi2nc[i]+chi2ncCoeff-minChi2WithCoeffs;
            Float64 chi2ncLikelihood = exp(-valf/2.0);
            valf = chi2continuum[i]+chi2cCoeff-minChi2WithCoeffs;
            Float64 chi2continuumLikelihood = exp(-valf/2.0);

            post = -(lmLikelihood + chi2continuumLikelihood + chi2ncLikelihood );
        }
        else //linear combination
        {
            post = lmCoeff*chi2lm[i] + chi2ncCoeff*chi2nc[i] + chi2cCoeff*chi2continuum[i];
        }

        resultCombined->ChiSquare[i] = post;
        resultCombined->Redshifts[i] = znc[i];
        resultCombined->Overlap[i] = -1.0;

    }

    if( resultCombined ) {
        store.StoreScopedGlobalResult( "resultdtreeBCombined", resultCombined );
    }

    return true;
    //*/

}


TFloat64List CMethodDTreeBSolve::GetBestRedshiftChi2List( CDataStore& store, std::string scopeStr,  Float64& minmerit, TFloat64List& zList)
{
    std::string scope = "dtreeBsolve.chisquare2solve.";
    scope.append(scopeStr.c_str());

    TOperatorResultMap meritResults = store.GetPerTemplateResult(scope.c_str());

    TFloat64List meritList;
    //init meritresults
    TOperatorResultMap::const_iterator it0 = meritResults.begin();
    {
        auto meritResult = std::dynamic_pointer_cast<const CChisquareResult> ( (*it0).second );
        for( Int32 i=0; i<meritResult->ChiSquare.size(); i++ )
        {
            meritList.push_back(DBL_MAX);
            zList.push_back(meritResult->Redshifts[i]);
        }
    }

    //find best merit for each tpl
    minmerit = DBL_MAX;
    for( TOperatorResultMap::const_iterator it = meritResults.begin(); it != meritResults.end(); it++ )
    {
        auto meritResult = std::dynamic_pointer_cast<const CChisquareResult>((*it).second);
        for( Int32 i=0; i<meritResult->ChiSquare.size(); i++ )
        {
            if( meritList[i] > meritResult->ChiSquare[i]){
                meritList[i] = meritResult->ChiSquare[i];
            }
            if( minmerit > meritResult->ChiSquare[i]){
                minmerit = meritResult->ChiSquare[i];
            }
        }
    }

    return meritList;

}
