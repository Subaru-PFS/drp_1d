#include <epic/redshift/method/linematching2solve.h>

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
#include <epic/redshift/method/linematching2solveresult.h>

using namespace NSEpic;
using namespace std;

IMPLEMENT_MANAGED_OBJECT( COperatorLineMatching2Solve )

COperatorLineMatching2Solve::COperatorLineMatching2Solve()
{
    // Peak Detection
    m_winsize = 250.0;
    m_cut = 5.0;
    m_strongcut = 2.0;
    m_minsize = 3;
    m_maxsize = 70;
    m_enlargeRate = 2.0;


    // Line Matching
    m_minMatchNum = 1;
    m_tol = 0.002;

    //eudlid overrides
    if(0)
    {
        m_winsize = 250.0;
        m_cut = 3.0;
        m_maxsize = 120;
    }

    //pfs TF overrides
    if(1)
    {
        m_winsize = 250.0;
        m_cut = 150;
        m_maxsize = 120;
        m_enlargeRate = 1.0;
        m_tol = 0.0025;
    }
}

COperatorLineMatching2Solve::~COperatorLineMatching2Solve()
{

}

const CLineMatching2SolveResult* COperatorLineMatching2Solve::Compute(  COperatorResultStore& resultStore, const CSpectrum& spc,
                                                                      const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep, const CRayCatalog& restRayCatalog )
{
    Bool storeResult = false;

    COperatorResultStore::CAutoScope resultScope( resultStore, "linematchingsolve2" );


    CPeakDetection peakDetection(m_winsize, m_cut, 1, m_enlargeRate);
    CConstRef<CPeakDetectionResult> peakDetectionResult = peakDetection.Compute( spc, lambdaRange);
    if( peakDetectionResult )
        resultStore.StoreGlobalResult( "peakdetection", *peakDetectionResult );

    CRayDetection rayDetection( m_cut, m_strongcut, m_winsize, m_minsize, m_maxsize);
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
        CLineMatching2SolveResult*  SolveResult = new CLineMatching2SolveResult();
        return SolveResult;
    }

    return NULL;
}

