#include <RedshiftLibrary/method/chisquaresolve.h>

#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/operator/correlation.h>
#include <RedshiftLibrary/operator/chisquare.h>
#include <RedshiftLibrary/operator/chisquare2.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/processflow/datastore.h>
#include <RedshiftLibrary/common/quicksort.h>
using namespace NSEpic;
using namespace std;


CMethodChisquareSolve::CMethodChisquareSolve(string calibrationPath)
{
    m_calibrationPath = calibrationPath;
}

CMethodChisquareSolve::~CMethodChisquareSolve()
{

}


std::shared_ptr<const CChisquareSolveResult> CMethodChisquareSolve::Compute(CDataStore& dataStore,
                                                                            const CSpectrum& spc,
                                                                            const CTemplateCatalog& tplCatalog,
                                                                            const TStringList& tplCategoryList,
                                                                            const TFloat64Range& lambdaRange,
                                                                            const TFloat64List& redshifts,
                                                                            Float64 overlapThreshold,
                                                                            const Float64 radius,
                                                                            std::string opt_interp)
{
    Bool storeResult = false;
    m_radius = radius;
    CDataStore::CAutoScope resultScope( dataStore, "chisquaresolve" );

    for( UInt32 i=0; i<tplCategoryList.size(); i++ )
    {
        std::string category = tplCategoryList[i];

        for( UInt32 j=0; j<tplCatalog.GetTemplateCount( category ); j++ )
        {
            const CTemplate& tpl = tplCatalog.GetTemplate( category, j );

            Solve(dataStore, spc, tpl, lambdaRange, redshifts, overlapThreshold, nType_full);

            storeResult = true;
        }
    }


    if( storeResult )
    {
        return std::shared_ptr<CChisquareSolveResult>( new CChisquareSolveResult() );
    }

    return NULL;
}

Bool CMethodChisquareSolve::Solve(CDataStore& dataStore,
                                  const CSpectrum& spc,
                                  const CTemplate& tpl,
                                  const TFloat64Range& lambdaRange,
                                  const TFloat64List& redshifts,
                                  Float64 overlapThreshold,
                                  Int32 spctype,
                                  std::string opt_interp)
{
    CSpectrum _spc = spc;
    CTemplate _tpl = tpl;

    _spc.SetType(spctype);
    _tpl.SetType(spctype);

    // prepare the unused masks
    std::vector<CMask> maskList;

    // Compute merit function
    COperatorChiSquare2 chiSquare(m_calibrationPath);
    //adding cast to be capable of reading redshift attribute
    auto  chisquareResult = std::dynamic_pointer_cast<CChisquareResult>(chiSquare.Compute( _spc, _tpl, lambdaRange, redshifts, overlapThreshold, maskList, opt_interp));
    
    chisquareResult->CallFindExtrema(m_radius);  
    
    if( !chisquareResult )
    {
        //Log.LogError( "Failed to compute chi square value");
        return false;
    }else{
        // Store results
        dataStore.StoreScopedPerTemplateResult( tpl, "chisquare", chisquareResult );
    }

    return true;
}
