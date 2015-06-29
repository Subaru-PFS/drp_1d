#include <epic/redshift/method/linematchingsolve.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/template/catalog.h>
#include <epic/redshift/operator/correlation.h>
#include <epic/redshift/operator/chisquare.h>
#include <epic/redshift/extremum/extremum.h>
#include <epic/redshift/operator/resultstore.h>

#include <epic/redshift/operator/peakdetection.h>
#include <epic/redshift/operator/raydetection.h>
#include <epic/redshift/operator/raymatching.h>

#include <epic/redshift/operator/peakdetectionresult.h>
#include <epic/redshift/operator/raydetectionresult.h>
#include <epic/redshift/operator/raymatchingresult.h>
#include <epic/redshift/method/linematchingsolveresult.h>

using namespace NSEpic;
using namespace std;

IMPLEMENT_MANAGED_OBJECT( COperatorLineMatchingSolve )

COperatorLineMatchingSolve::COperatorLineMatchingSolve()
{
    // Peak Detection
    m_winsize = 250.0;
    m_cut = 5.0;
    m_strongcut = 2.0;

    // Line Matching
    m_minMatchNum = 1;
    m_tol = 0.002;

    //eudlid overrides
    if(0)
    {
        m_winsize = 250.0;
        m_cut = 3.0;
    }

    //pfs TF overrides
    if(1)
    {
        m_winsize = 250.0;
        m_cut = 3.5;
    }
}

COperatorLineMatchingSolve::~COperatorLineMatchingSolve()
{

}

const CLineMatchingSolveResult* COperatorLineMatchingSolve::Compute(  COperatorResultStore& resultStore, const CSpectrum& spc,
                                                                      const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep, const CRayCatalog& restRayCatalog )
{
    Bool storeResult = false;

    COperatorResultStore::CAutoScope resultScope( resultStore, "linematchingsolve" );


    CPeakDetection peakDetection(m_winsize, m_cut);
    CConstRef<CPeakDetectionResult> peakDetectionResult = peakDetection.Compute( spc, lambdaRange);
    if( peakDetectionResult )
        resultStore.StoreGlobalResult( "peakdetection", *peakDetectionResult );

    CRayDetection rayDetection( m_cut, m_strongcut, m_winsize);
    CConstRef<CRayDetectionResult> rayDetectionResult = rayDetection.Compute( spc, lambdaRange, peakDetectionResult->PeakList, peakDetectionResult->EnlargedPeakList );

    if( rayDetectionResult ) {
        resultStore.StoreGlobalResult( "raycatalog", *rayDetectionResult );

        if(rayDetectionResult->RayCatalog.GetList().size()<1){
            //return NULL;
        }
    }

    // --- Match
    CRayMatching rayMatching;
    CRef<CRayMatchingResult> rayMatchingResult = rayMatching.Compute(rayDetectionResult->RayCatalog, restRayCatalog, redshiftsRange, m_minMatchNum, m_tol );

    // Store matching results
    if( rayMatchingResult )
        resultStore.StoreGlobalResult( "raymatching", *rayMatchingResult );


    storeResult = true; //always save a matching result
    if( storeResult )
    {
        CLineMatchingSolveResult*  SolveResult = new CLineMatchingSolveResult();
        return SolveResult;
    }

    return NULL;
}
