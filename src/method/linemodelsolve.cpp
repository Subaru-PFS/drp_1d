

#include <epic/redshift/method/linemodelsolve.h>

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
    CSpectrum _spc;

    std::string scopeStr = "linemodel";


    Int32 spcType = CLineModelSolveResult::nType_noContinuum;
    //Int32 spcType = CLineModelSolveResult::nType_raw;
    Int32 _spctype = spcType;

    if(_spctype == CLineModelSolveResult::nType_continuumOnly){
        // use continuum only
        _spc = spc;
        CSpectrumFluxAxis spcfluxAxis = _spc.GetFluxAxis();
        spcfluxAxis.Subtract(spcWithoutCont.GetFluxAxis());
        CSpectrumFluxAxis& sfluxAxisPtr = _spc.GetFluxAxis();
        sfluxAxisPtr = spcfluxAxis;
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
    //Int32 widthType = nWidthType_ZDriven;
    //Int32 widthType = nWidthType_Fixed;

    // Compute merit function
    COperatorLineModel linemodel;
    CRef<CLineModelResult>  result = (CLineModelResult*)linemodel.Compute( _spc, restraycatalog, lambdaRange, redshifts, widthType);

    if( !result )
    {
        //Log.LogInfo( "Failed to compute linemodel");
        return false;
    }else{
        // Store results
        resultStore.StoreGlobalResult( scopeStr.c_str(), *result );
    }


    return true;
}
