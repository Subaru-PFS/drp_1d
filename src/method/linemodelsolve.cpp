#include <epic/redshift/method/linemodelsolve.h>

#include <epic/core/log/log.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/template/catalog.h>
#include <epic/redshift/operator/linemodel.h>
#include <epic/redshift/extremum/extremum.h>
#include <epic/redshift/processflow/datastore.h>

using namespace NSEpic;
using namespace std;

IMPLEMENT_MANAGED_OBJECT( CLineModelSolve )

CLineModelSolve::CLineModelSolve()
{

}

CLineModelSolve::~CLineModelSolve()
{

}

const std::string GetDescription(){
    std::string desc;

    desc = "Method LinemodelSolve:";
    desc.append("\tparam string: linewidthtype = {""psfinstrumentdriven"", ""linewidthtype"", ""linewidthtype""}");

    return desc;

}

const CLineModelSolveResult* CLineModelSolve::Compute(  CDataStore& dataStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CRayCatalog& restraycatalog,
                                                        const TFloat64Range& lambdaRange, const TFloat64List& redshifts)
{
    Bool storeResult = false;
    CDataStore::CAutoScope resultScope( dataStore, "linemodelsolve" );

    Solve( dataStore, spc, spcWithoutCont, restraycatalog, lambdaRange, redshifts);
    storeResult = true;


    if( storeResult )
    {
        CLineModelSolveResult*  result = new CLineModelSolveResult();
        return result;
    }

    return NULL;
}

Bool CLineModelSolve::Solve( CDataStore& dataStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CRayCatalog& restraycatalog,
                             const TFloat64Range& lambdaRange, const TFloat64List& redshifts)
{
    std::string scopeStr = "linemodel";

    CSpectrum _spc = spc;
    CSpectrum _spcContinuum = spc;
    CSpectrumFluxAxis spcfluxAxis = _spcContinuum.GetFluxAxis();
    spcfluxAxis.Subtract(spcWithoutCont.GetFluxAxis());
    CSpectrumFluxAxis& sfluxAxisPtr = _spcContinuum.GetFluxAxis();
    sfluxAxisPtr = spcfluxAxis;

    std::string widthType;
    dataStore.GetScopedParam( "linewidthtype", widthType, "psfinstrumentdriven" );

    // Compute merit function
    COperatorLineModel linemodel;
    CRef<CLineModelResult>  result = (CLineModelResult*)linemodel.Compute( dataStore, _spc, _spcContinuum, restraycatalog, lambdaRange, redshifts, widthType);

//    static Float64 cutThres = 2.0;
//    static Int32 bestSolutionIdx = 0;
//    Int32 nValidLines = result->GetNLinesOverCutThreshold(bestSolutionIdx, cutThres);
//    Float64 bestExtremaMerit = result->GetExtremaMerit(0);
//    Log.LogInfo( "Linemodelsolve : bestExtremaMerit, %f", bestExtremaMerit);
//    Float64 thres = 0.001;
//    Int32 idxNextValid = 1;
//    for(Int32 idnext=1; idnext<result->Redshifts.size(); idnext++){
//       if( std::abs(result->Redshifts[idnext]-result->Redshifts[0])> thres){
//           idxNextValid = idnext;
//           break;
//       }
//    }
//    Float64 nextExtremaMerit = result->GetExtremaMerit(idxNextValid);
//    Log.LogInfo( "Linemodelsolve : nextExtremaMerit, %f", nextExtremaMerit);
//    if(nValidLines<2 || (bestExtremaMerit - nextExtremaMerit) > -50.0 ){
//        result=0;
//        Log.LogInfo( "Linemodelsolve : result set to 0" );
//    }
//    Log.LogInfo( "Linemodelsolve : for best solution, %d valid lines found", nValidLines);

    if( !result )
    {
        //Log.LogInfo( "Failed to compute linemodel");
        return false;
    }else{
        // Store results
        dataStore.StoreScopedGlobalResult( scopeStr.c_str(), *result );
    }

    return true;
}
