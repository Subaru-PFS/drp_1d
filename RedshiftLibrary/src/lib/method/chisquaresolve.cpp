
#include <RedshiftLibrary/method/chisquaresolve.h>

#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/operator/correlation.h>
#include <RedshiftLibrary/operator/chisquare.h>
#include <RedshiftLibrary/operator/chisquare2.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/processflow/datastore.h>

using namespace NSEpic;
using namespace std;


CMethodChisquareSolve::CMethodChisquareSolve(string calibrationPath)
{
    m_calibrationPath = calibrationPath;
}

CMethodChisquareSolve::~CMethodChisquareSolve()
{

}


std::shared_ptr<const CChisquareSolveResult>  CMethodChisquareSolve::Compute(  CDataStore& dataStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                                        const CTemplateCatalog& tplCatalog, const TStringList& tplCategoryList,
                                                        const TFloat64Range& lambdaRange, const TFloat64List& redshifts, Float64 overlapThreshold, std::string opt_interp )
{
    Bool storeResult = false;

    CDataStore::CAutoScope resultScope( dataStore, "chisquaresolve" );

    for( UInt32 i=0; i<tplCategoryList.size(); i++ )
    {
        std::string category = tplCategoryList[i];

        for( UInt32 j=0; j<tplCatalog.GetTemplateCount( category ); j++ )
        {
            const CTemplate& tpl = tplCatalog.GetTemplate( category, j );
            const CTemplate& tplWithoutCont = tplCatalog.GetTemplateWithoutContinuum( category, j );

            Solve( dataStore, spc, spcWithoutCont, tpl, tplWithoutCont, lambdaRange, redshifts, overlapThreshold, nType_full );

            storeResult = true;
        }
    }


    if( storeResult )
    {
        return std::shared_ptr<CChisquareSolveResult>( new CChisquareSolveResult() );
    }

    return NULL;
}

Bool CMethodChisquareSolve::Solve( CDataStore& dataStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CTemplate& tpl, const CTemplate& tplWithoutCont,
                               const TFloat64Range& lambdaRange, const TFloat64List& redshifts, Float64 overlapThreshold, Int32 spctype, std::string opt_interp )
{
    CSpectrum _spc;
    CTemplate _tpl;

    if(spctype == nType_continuumOnly){
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
    }else if(spctype == nType_full){
        // use full spectrum
        _spc = spc;
        _tpl = tpl;

    }else if(spctype == nType_noContinuum){
        // use spectrum without continuum
        _spc = spcWithoutCont;
        _tpl = tplWithoutCont;
    }
    // prepare the unused masks
    std::vector<CMask> maskList;

    // Compute merit function
    COperatorChiSquare2 chiSquare(m_calibrationPath);
    auto  chisquareResult = chiSquare.Compute( _spc, _tpl, lambdaRange, redshifts, overlapThreshold, maskList, opt_interp);
    if( !chisquareResult )
    {
        //Log.LogInfo( "Failed to compute chi square value");
        return false;
    }else{
        // Store results
        dataStore.StoreScopedPerTemplateResult( tpl, "chisquare", chisquareResult );
    }

    return true;
}
