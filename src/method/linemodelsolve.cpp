

#include <epic/redshift/method/linemodelsolve.h>

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

    std::string type;
    dataStore.GetScopedParam( "type", type, "nocontinuum" );

    CSpectrum _spc;
    CSpectrum _spcContinuum = spc;
    CSpectrumFluxAxis spcfluxAxis = _spcContinuum.GetFluxAxis();
    spcfluxAxis.Subtract(spcWithoutCont.GetFluxAxis());
    CSpectrumFluxAxis& sfluxAxisPtr = _spcContinuum.GetFluxAxis();
    sfluxAxisPtr = spcfluxAxis;

    if(type == "continuumonly"){
        // use continuum only
        _spc = _spcContinuum;
        //scopeStr = "linemodel_continuum";
    }else if(type == "raw"){
        // use full spectrum
        _spc = spc;
        //scopeStr = "linemodel";
    }else if(type == "nocontinuum" ){
        // use spectrum without continuum
        _spc = spc;
        CSpectrumFluxAxis spcfluxAxis = spcWithoutCont.GetFluxAxis();
        CSpectrumFluxAxis& sfluxAxisPtr = _spc.GetFluxAxis();
        sfluxAxisPtr = spcfluxAxis;
        //scopeStr = "linemodel_nocontinuum";
    }


    std::string widthType;
    dataStore.GetScopedParam( "widthtype", widthType, "psfinstrumentdriven" );

    // Compute merit function
    COperatorLineModel linemodel;
    CRef<CLineModelResult>  result = (CLineModelResult*)linemodel.Compute( _spc, _spcContinuum, restraycatalog, lambdaRange, redshifts, widthType);

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
