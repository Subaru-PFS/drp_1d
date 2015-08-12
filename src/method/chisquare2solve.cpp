

#include <epic/redshift/method/chisquare2solve.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/template/catalog.h>
#include <epic/redshift/operator/correlation.h>
#include <epic/redshift/operator/chisquare.h>
#include <epic/redshift/operator/chisquare2.h>
#include <epic/redshift/extremum/extremum.h>
#include <epic/redshift/operator/resultstore.h>

using namespace NSEpic;
using namespace std;

IMPLEMENT_MANAGED_OBJECT( CMethodChisquare2Solve )

CMethodChisquare2Solve::CMethodChisquare2Solve()
{

}

CMethodChisquare2Solve::~CMethodChisquare2Solve()
{

}


const CChisquare2SolveResult* CMethodChisquare2Solve::Compute(  COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                                        const CTemplateCatalog& tplCatalog, const TTemplateCategoryList& tplCategoryList,
                                                        const TFloat64Range& lambdaRange, const TFloat64List& redshifts, Float64 overlapThreshold )
{
    Bool storeResult = false;

    COperatorResultStore::CAutoScope resultScope( resultStore, "chisquare2solve" );

    for( UInt32 i=0; i<tplCategoryList.size(); i++ )
    {
        CTemplate::ECategory category = tplCategoryList[i];

        for( UInt32 j=0; j<tplCatalog.GetTemplateCount( category ); j++ )
        {
            const CTemplate& tpl = tplCatalog.GetTemplate( category, j );
            const CTemplate& tplWithoutCont = tplCatalog.GetTemplateWithoutContinuum( category, j );

            Solve( resultStore, spc, spcWithoutCont, tpl, tplWithoutCont, lambdaRange, redshifts, overlapThreshold, nType_raw );

            storeResult = true;
        }
    }


    if( storeResult )
    {
        CChisquare2SolveResult*  ChisquareSolveResult = new CChisquare2SolveResult();
        return ChisquareSolveResult;
    }

    return NULL;
}

Bool CMethodChisquare2Solve::Solve( COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CTemplate& tpl, const CTemplate& tplWithoutCont,
                               const TFloat64Range& lambdaRange, const TFloat64List& redshifts, Float64 overlapThreshold, Int32 spctype )
{
    CSpectrum _spc;
    CTemplate _tpl;
    std::string scopeStr = "chisquare";
    Int32 _ntype = 1;
    Int32 _spctype = spctype;
    Int32 _spctypetab[3] = {nType_raw, nType_noContinuum, nType_continuumOnly};


    //case: nType_all
    if(spctype = nType_all){
        _ntype = 3;
    }

    for( Int32 i=0; i<_ntype; i++){
        if(spctype = nType_all){
            _spctype = _spctypetab[i];
        }else{
            _spctype = spctype;
        }

        if(_spctype == nType_continuumOnly){
            // use continuum only
            _spc = spc;
            CSpectrumFluxAxis spcfluxAxis = _spc.GetFluxAxis();
            spcfluxAxis.Subtract(spcWithoutCont.GetFluxAxis());
            CSpectrumFluxAxis& sfluxAxisPtr = _spc.GetFluxAxis();
            sfluxAxisPtr = spcfluxAxis;
            _tpl = tpl;
            CSpectrumFluxAxis tplfluxAxis = _tpl.GetFluxAxis();
            tplfluxAxis.Subtract(tplWithoutCont.GetFluxAxis());
            CSpectrumFluxAxis& tfluxAxisPtr = _tpl.GetFluxAxis();
            tfluxAxisPtr = tplfluxAxis;
            scopeStr = "chisquare_continuum";
        }else if(_spctype == nType_raw){
            // use full spectrum
            _spc = spc;
            _tpl = tpl;
            scopeStr = "chisquare";

        }else if(_spctype == nType_noContinuum){
            // use spectrum without continuum
            _spc = spc;
            CSpectrumFluxAxis spcfluxAxis = spcWithoutCont.GetFluxAxis();
            CSpectrumFluxAxis& sfluxAxisPtr = _spc.GetFluxAxis();
            sfluxAxisPtr = spcfluxAxis;
            _tpl = tpl;
            CSpectrumFluxAxis tplfluxAxis = tplWithoutCont.GetFluxAxis();
            CSpectrumFluxAxis& tfluxAxisPtr = _tpl.GetFluxAxis();
            tfluxAxisPtr = tplfluxAxis;
            scopeStr = "chisquare_continuum";
            scopeStr = "chisquare_nocontinuum";
        }

        // Compute merit function
        COperatorChiSquare2 chiSquare;
        CRef<CChisquareResult>  chisquareResult = (CChisquareResult*)chiSquare.Compute( _spc, _tpl, lambdaRange, redshifts, overlapThreshold );
        if( !chisquareResult )
        {
            //Log.LogInfo( "Failed to compute chi square value");
            return false;
        }else{
            // Store results
            resultStore.StorePerTemplateResult( tpl, scopeStr.c_str(), *chisquareResult );
        }
    }

    return true;
}
