#include <epic/redshift/method/dtreebsolve.h>
#include <epic/redshift/method/dtreebsolveresult.h>

#include <epic/core/log/log.h>
#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/template/catalog.h>
#include <epic/redshift/operator/correlation.h>
#include <epic/redshift/operator/chisquare.h>
#include <epic/redshift/extremum/extremum.h>
#include <epic/redshift/operator/resultstore.h>

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

const CDTreeBSolveResult *COperatorDTreeBSolve::Compute(COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                                        const CTemplateCatalog& tplCatalog, const TTemplateCategoryList& tplCategoryList, const CRayCatalog &restRayCatalog,
                                                        const TFloat64Range& lambdaRange, const TFloat64List &redshifts)
{
    Bool storeResult = false;

    COperatorResultStore::CAutoScope resultScope( resultStore, "dtreeBsolve" );

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

Bool COperatorDTreeBSolve::Solve(COperatorResultStore &resultStore, const CSpectrum &spc, const CSpectrum &spcWithoutCont, const CTemplateCatalog &tplCatalog, const TTemplateCategoryList &tplCategoryList, const CRayCatalog &restRayCatalog, const TFloat64Range &lambdaRange, const TFloat64List &redshifts)
{
    std::string scopeStr = "linemodel";

    CSpectrum _spcContinuum = spc;
    CSpectrumFluxAxis spcfluxAxis = _spcContinuum.GetFluxAxis();
    spcfluxAxis.Subtract(spcWithoutCont.GetFluxAxis());
    CSpectrumFluxAxis& sfluxAxisPtr = _spcContinuum.GetFluxAxis();
    sfluxAxisPtr = spcfluxAxis;

    Int32 widthType = CLineModelElement::nWidthType_PSFInstrumentDriven;
    //Int32 widthType = CLineModelElement::nWidthType_ZDriven;
    //Int32 widthType = CLineModelElement::nWidthType_Fixed;

    // Compute merit function
    COperatorLineModel linemodel;
    CRef<CLineModelResult>  result = (CLineModelResult*)linemodel.Compute(resultStore, spc, _spcContinuum, restRayCatalog, lambdaRange, redshifts, widthType);

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
        resultStore.StoreGlobalResult( scopeStr.c_str(), *result );
    }

    //*
    //_///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // compute chisquare
    if( result->Extrema.size() == 0 )
    {
        return false;
    }
    Float64 overlapThreshold = 1.0;
    // Compute merit function
    TFloat64List extremumRedshifts( result->Extrema.size() );
    for( Int32 i=0; i<result->Extrema.size(); i++ )
    {
        extremumRedshifts[i] = result->Extrema[i];
    }

    CMethodChisquare2Solve chiSolve;
    CConstRef<CChisquare2SolveResult> chisolveResult = chiSolve.Compute( resultStore, spc, spcWithoutCont,
                                                                        tplCatalog, tplCategoryList,
                                                                        lambdaRange, extremumRedshifts, overlapThreshold );
    if( chisolveResult ) {
        resultStore.StoreGlobalResult( "redshiftresult", *chisolveResult );
    }
    //_///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //*/
    return true;
}


