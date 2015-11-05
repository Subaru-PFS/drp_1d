#include <epic/redshift/method/linemodelsolve.h>

#include <epic/core/log/log.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/template/catalog.h>
#include <epic/redshift/operator/linemodel.h>
#include <epic/redshift/extremum/extremum.h>
#include <epic/redshift/operator/resultstore.h>

using namespace NSEpic;
using namespace std;

IMPLEMENT_MANAGED_OBJECT( CLineModelSolve )

CLineModelSolve::CLineModelSolve()
{

}

CLineModelSolve::~CLineModelSolve()
{

}


const CLineModelSolveResult* CLineModelSolve::Compute(  COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CRayCatalog& restraycatalog,
                                                        const TFloat64Range& lambdaRange, const TFloat64List& redshifts)
{
    Bool storeResult = false;
    COperatorResultStore::CAutoScope resultScope( resultStore, "linemodelsolve" );

    Solve( resultStore, spc, spcWithoutCont, restraycatalog, lambdaRange, redshifts);
    storeResult = true;


    if( storeResult )
    {
        CLineModelSolveResult*  result = new CLineModelSolveResult();
        return result;
    }

    return NULL;
}

Bool CLineModelSolve::Solve( COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CRayCatalog& restraycatalog,
                             const TFloat64Range& lambdaRange, const TFloat64List& redshifts)
{
    std::string scopeStr = "linemodel";

    //Int32 spcType = CLineModelSolveResult::nType_noContinuum;
    Int32 spcType = CLineModelSolveResult::nType_raw;
    //Int32 spcType = CLineModelSolveResult::nType_continuumOnly;
    Int32 _spctype = spcType;

    CSpectrum _spc;
    CSpectrum _spcContinuum = spc;
    CSpectrumFluxAxis spcfluxAxis = _spcContinuum.GetFluxAxis();
    spcfluxAxis.Subtract(spcWithoutCont.GetFluxAxis());
    CSpectrumFluxAxis& sfluxAxisPtr = _spcContinuum.GetFluxAxis();
    sfluxAxisPtr = spcfluxAxis;

    if(_spctype == CLineModelSolveResult::nType_continuumOnly){
        // use continuum only
        _spc = _spcContinuum;
        //scopeStr = "linemodel_continuum";
    }else if(_spctype == CLineModelSolveResult::nType_raw){
        // use full spectrum
        _spc = spc;
        //scopeStr = "linemodel";
    }else if(_spctype == CLineModelSolveResult::nType_noContinuum){
        // use spectrum without continuum
        _spc = spc;
        CSpectrumFluxAxis spcfluxAxis = spcWithoutCont.GetFluxAxis();
        CSpectrumFluxAxis& sfluxAxisPtr = _spc.GetFluxAxis();
        sfluxAxisPtr = spcfluxAxis;
        //scopeStr = "linemodel_nocontinuum";
    }


    Int32 widthType = CLineModelElement::nWidthType_PSFInstrumentDriven;
    //Int32 widthType = CLineModelElement::nWidthType_ZDriven;
    //Int32 widthType = CLineModelElement::nWidthType_Fixed;

    // Compute merit function
    COperatorLineModel linemodel;
    CRef<CLineModelResult>  result = (CLineModelResult*)linemodel.Compute(resultStore, _spc, _spcContinuum, restraycatalog, lambdaRange, redshifts, widthType);

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
//    if(nValidLines<2 /*|| (bestExtremaMerit - nextExtremaMerit) > -50.0*/ ){
//        result=0;
//        Log.LogInfo( "Linemodelsolve : result set to 0" );
//    }

    Log.LogInfo( "Linemodelsolve : for best solution, %d valid lines found", nValidLines);

    if( !result )
    {
        //Log.LogInfo( "Failed to compute linemodel");
        return false;
    }else{
        // Store results
        resultStore.StoreGlobalResult( scopeStr.c_str(), *result );
    }

    /*
    //_///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //**** ----- TEMP. Location for this add. data to be processed ----- ****
    if( result->Extrema.size() == 0 )
    {
        return false;
    }
    // Compute merit function
    TFloat64List extremumRedshifts( result->Extrema.size() );
    for( Int32 i=0; i<result->Extrema.size(); i++ )
    {
        extremumRedshifts[i] = result->Extrema[i];
    }

    //COperatorChiSquare2 meritChiSquare;
    //CRef<CCorrelationResult> chisquareResult = (CCorrelationResult*)meritChiSquare.Compute( spc, tpl, lambdaRange, extremumRedshifts, overlapThreshold );
    if( !chisquareResult )
    {
        return false;
    }

    // Store results
    resultStore.StorePerTemplateResult( tpl, "merit", *chisquareResult );
    //_///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/


    return true;
}
