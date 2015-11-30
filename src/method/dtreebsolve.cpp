#include <epic/redshift/method/dtreebsolve.h>
#include <epic/redshift/method/dtreebsolveresult.h>

#include <epic/core/log/log.h>
#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/template/catalog.h>
#include <epic/redshift/operator/correlation.h>
#include <epic/redshift/operator/chisquare.h>
#include <epic/redshift/extremum/extremum.h>
#include <epic/redshift/processflow/datastore.h>

#include <epic/redshift/operator/correlation.h>
#include <epic/redshift/operator/chicorr.h>
#include <epic/redshift/method/blindsolveresult.h>

#include <epic/redshift/method/blindsolve.h>
#include <epic/redshift/method/blindsolveresult.h>
#include <epic/redshift/operator/chisquare.h>
#include <epic/redshift/operator/chisquare2.h>

#include <epic/redshift/operator/peakdetection.h>
#include <epic/redshift/operator/peakdetectionresult.h>
#include <epic/redshift/operator/raydetection.h>
#include <epic/redshift/operator/raydetectionresult.h>
#include <epic/redshift/operator/raymatching.h>
#include <epic/redshift/operator/raymatchingresult.h>

#include <epic/redshift/method/chisquare2solve.h>
#include <epic/redshift/method/fullsolve.h>
#include <epic/redshift/method/correlationsolve.h>
#include <epic/redshift/method/linematchingsolve.h>

#include <epic/redshift/operator/linemodel.h>

using namespace NSEpic;
using namespace std;


COperatorDTreeBSolve::COperatorDTreeBSolve()
{

}

COperatorDTreeBSolve::~COperatorDTreeBSolve()
{

}

const std::string COperatorDTreeBSolve::GetDescription()
{
    std::string desc;

    desc = "Method decisionaltreeb:\n";

    desc.append("\tparam: linemodel.linetypefilter = {""no"", ""E"", ""A""}\n");
    desc.append("\tparam: linemodel.lineforcefilter = {""no"", ""S""}\n");
    desc.append("\tparam: linemodel.fittingmethod = {""hybrid"", ""individual""}\n");
    desc.append("\tparam: linemodel.continuumcomponent = {""fromspectrum"", ""nocontinuum"", ""zero""}\n");
    desc.append("\tparam: linemodel.linewidthtype = {""psfinstrumentdriven"", ""zdriven"", ""fixedvelocity"", ""fixed""}\n");
    desc.append("\tparam: linemodel.instrumentresolution = <float value>\n");
    desc.append("\tparam: linemodel.velocityemission = <float value>\n");
    desc.append("\tparam: linemodel.velocityabsorption = <float value>\n");

    desc.append("\tparam: linemodel.continuumreestimation = {""no"", ""onlyextrema"", ""always""}\n");
    desc.append("\tparam: linemodel.extremacount = <float value>\n");

    desc.append("\tparam: chisquare.overlapthreshold = <float value>\n");
    desc.append("\tparam: chisquare.redshiftsupport = {""full"", ""extremaextended""}\n");

    return desc;

}

std::shared_ptr<const CDTreeBSolveResult> COperatorDTreeBSolve::Compute( CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                                        const CTemplateCatalog& tplCatalog, const TStringList &tplCategoryList, const CRayCatalog &restRayCatalog,
                                                        const TFloat64Range& lambdaRange, const TFloat64List &redshifts)
{
    Bool storeResult = false;

    CDataStore::CAutoScope resultScope( resultStore, "dtreeBsolve" );

    storeResult = Solve(resultStore, spc, spcWithoutCont,
                                       tplCatalog, tplCategoryList, restRayCatalog,
                                       lambdaRange, redshifts );

    //storeResult = true;
    if( storeResult )
    {
        return std::shared_ptr<const CDTreeBSolveResult>( new CDTreeBSolveResult() );
    }

    return NULL;
}

Bool COperatorDTreeBSolve::Solve(CDataStore &dataStore, const CSpectrum &spc, const CSpectrum &spcWithoutCont, const CTemplateCatalog &tplCatalog, const TStringList &tplCategoryList, const CRayCatalog &restRayCatalog, const TFloat64Range &lambdaRange, const TFloat64List &redshifts)
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
    dataStore.GetScopedParam( "linemodel.linewidthtype", opt_lineWidthType, "fixedvelocity" );
    Float64 opt_resolution;
    dataStore.GetScopedParam( "linemodel.instrumentresolution", opt_resolution, 3000.0 );
    Float64 opt_velocity_emission;
    dataStore.GetScopedParam( "linemodel.velocityemission", opt_velocity_emission, 100.0 );
    Float64 opt_velocity_absorption;
    dataStore.GetScopedParam( "linemodel.velocityabsorption", opt_velocity_absorption, 300.0 );
    std::string opt_continuumreest;
    dataStore.GetScopedParam( "linemodel.continuumreestimation", opt_continuumreest, "no" );
    Float64 opt_extremacount;
    dataStore.GetScopedParam( "linemodel.extremacount", opt_extremacount, 10.0 );


    // Compute merit function
    COperatorLineModel linemodel;
    auto result = dynamic_pointer_cast<CLineModelResult>(linemodel.Compute(dataStore, spc, _spcContinuum, restRayCatalog, opt_linetypefilter, opt_lineforcefilter, lambdaRange, redshifts, opt_extremacount, opt_fittingmethod, opt_continuumcomponent, opt_lineWidthType, opt_resolution, opt_velocity_emission, opt_velocity_absorption, opt_continuumreest) );

    /*
    static Float64 cutThres = 2.0;
    static Int32 bestSolutionIdx = 0;
    Int32 nValidLines = result->GetNLinesOverCutThreshold(bestSolutionIdx, cutThres);
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
    }

    //*
    //_///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // compute chisquare
    if( result->Extrema.size() == 0 )
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
    TFloat64List redshiftsChi2;
    if(opt_redshiftsupport == "full"){
        redshiftsChi2 = redshifts;
    }else if(opt_redshiftsupport == "extremaextended"){
        redshiftsChi2= result->ExtremaExtendedRedshifts;
    }else{
        redshiftsChi2= result->Extrema;
    }

    std::string spcComponent = "nocontinuum";
    CMethodChisquare2Solve chiSolve;
    auto chisolveResultnc = chiSolve.Compute( dataStore, spc, spcWithoutCont,
                                                                        tplCatalog, tplCategoryList,
                                                                        lambdaRange, redshiftsChi2, overlapThreshold, spcComponent);
    if( chisolveResultnc ) {
        dataStore.StoreScopedGlobalResult( "redshiftresult", chisolveResultnc );
    }

    spcComponent = "continuum";
    auto chisolveResultcontinuum = chiSolve.Compute( dataStore, spc, spcWithoutCont,
                                                                        tplCatalog, tplCategoryList,
                                                                        lambdaRange, redshiftsChi2, overlapThreshold, spcComponent);

    //_///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //*/
    return true;
}


