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

IMPLEMENT_MANAGED_OBJECT( COperatorDTreeBSolve )

COperatorDTreeBSolve::COperatorDTreeBSolve()
{

}

COperatorDTreeBSolve::~COperatorDTreeBSolve()
{

}

const std::string COperatorDTreeBSolve::GetDescription()
{
    std::string desc;

    desc = "Method DecisionalTreeB:\n";

    desc.append("\tparam: linemodel.linewidthtype = {""psfinstrumentdriven"", ""zdriven"", ""fixed""}\n");
    desc.append("\tparam: linemodel.continuumreestimation = {""no"", ""onlyextrema"", ""always""}\n");
    desc.append("\tparam: linemodel.extremacount = <float value>\n");


    return desc;

}

const CDTreeBSolveResult *COperatorDTreeBSolve::Compute( CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
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
        CDTreeBSolveResult*  SolveResult = new CDTreeBSolveResult();
        return SolveResult;
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

    std::string opt_lineWidthType;
    dataStore.GetScopedParam( "linemodel.linewidthtype", opt_lineWidthType, "psfinstrumentdriven" );
    std::string opt_continuumreest;
    dataStore.GetScopedParam( "linemodel.continuumreestimation", opt_continuumreest, "no" );
    Float64 opt_extremacount;
    dataStore.GetScopedParam( "linemodel.extremacount", opt_extremacount, 20.0 );


    // Compute merit function
    COperatorLineModel linemodel;
    CRef<CLineModelResult>  result = (CLineModelResult*)linemodel.Compute(dataStore, spc, _spcContinuum, restRayCatalog, lambdaRange, redshifts, opt_extremacount, opt_lineWidthType, opt_continuumreest);

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
        dataStore.StoreScopedGlobalResult( "linemodel", *result );
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
    // Compute merit function
    TFloat64List extremumRedshifts( result->Extrema.size() );
    for( Int32 i=0; i<result->Extrema.size(); i++ )
    {
        extremumRedshifts[i] = result->Extrema[i];
    }

    std::string spcComponent = "all";
    CMethodChisquare2Solve chiSolve;
    CConstRef<CChisquare2SolveResult> chisolveResult = chiSolve.Compute( dataStore, spc, spcWithoutCont,
                                                                        tplCatalog, tplCategoryList,
                                                                        lambdaRange, extremumRedshifts, overlapThreshold, spcComponent);
    if( chisolveResult ) {
        dataStore.StoreScopedGlobalResult( "redshiftresult", *chisolveResult );
    }
    //_///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //*/
    return true;
}


