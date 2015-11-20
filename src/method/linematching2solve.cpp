#include <epic/core/debug/assert.h>
#include <epic/redshift/extremum/extremum.h>
#include <epic/redshift/method/linematching2solve.h>
#include <epic/redshift/method/linematching2solveresult.h>
#include <epic/redshift/operator/chisquare.h>
#include <epic/redshift/operator/correlation.h>
#include <epic/redshift/operator/peakdetection.h>
#include <epic/redshift/operator/peakdetectionresult.h>
#include <epic/redshift/operator/raydetection.h>
#include <epic/redshift/operator/raydetectionresult.h>
#include <epic/redshift/operator/raymatching.h>
#include <epic/redshift/operator/raymatchingresult.h>
#include <epic/redshift/processflow/datastore.h>
#include <epic/redshift/spectrum/template/catalog.h>

using namespace NSEpic;
using namespace std;

/**
 * Attribution constructor.
 */
COperatorLineMatching2Solve::COperatorLineMatching2Solve()
{    
    // Peak Detection
    m_winsize = 250.0;
    m_cut = 5.0;
    m_detectioncut = 5.0;
    m_detectionnoiseoffset = 0.0; // useful for simulation data
    m_strongcut = 2.0;
    m_minsize = 3;
    m_maxsize = 70;
    m_enlargeRate = 2.0;

    // Line Matching
    m_minMatchNum = 1;
    m_tol = 0.002;

    //EUCLID overrides
    if( 0 )
    {
        m_winsize = 250.0;
        m_cut = 3.0;
        m_detectioncut = 3.0;
        m_maxsize = 120;
    }

    //PFS Really just lines, TF overrides
    if( 0 )
    {
        // TF
        m_winsize = 250.0;
        m_cut = 1.0;
        m_detectioncut = 10.0;
        m_detectionnoiseoffset = 0.0001;
        m_maxsize = 120;
        m_enlargeRate = 1.0;
        m_tol = 0.0001;

        // + ErrF overrides
        //m_cut = 0.2;
    }

    //PFS LBGABS, TF overrides
    if( 0 )
    {
        // TF
        m_winsize = 500.0;
        m_cut = 1.0;
        m_detectioncut = 10.0;
        m_detectionnoiseoffset = 0.0001;
        m_maxsize = 120;
        m_enlargeRate = 1.0;
        m_tol = 0.0001;

        // + ErrF overrides
        //m_cut = 0.2;
    }

    //PFS Really just lines, F + ErrF overrides
    if( 1 )
    {   // F + ErrF
        m_winsize = 250.0;
        m_cut = 1.0;
        m_detectioncut = 1.0;
        m_maxsize = 120;
        m_enlargeRate = 2.0;
        m_tol = 0.0001;
    }
}

/**
 * Empty destructor.
 */
COperatorLineMatching2Solve::~COperatorLineMatching2Solve()
{

}

/**
 * Attempts to find a matching line using the input parameters.
 */
std::shared_ptr<const CLineMatching2SolveResult> COperatorLineMatching2Solve::Compute( CDataStore& resultStore, 
										       const CSpectrum& spc, 
										       const TFloat64Range& lambdaRange, 
										       const TFloat64Range& redshiftsRange, 
										       Float64 redshiftStep, 
										       const CRayCatalog& restRayCatalog )
{
    Bool storeResult = false;

    CDataStore::CAutoScope resultScope( resultStore, "linematching2solve" );
    //Int32 lineType = CRay::nType_Emission;
    Int32 lineType = CRay::nType_Absorption;

    CPeakDetection peakDetection(m_winsize, m_detectioncut, 1, m_enlargeRate, m_detectionnoiseoffset);
    CSpectrum _spc = spc;
    if( lineType == CRay::nType_Absorption )
      {
        _spc.InvertFlux();
      }

    auto peakDetectionResult = peakDetection.Compute( _spc, lambdaRange );
    if( peakDetectionResult )
      {
        resultStore.StoreScopedGlobalResult( "peakdetection", peakDetectionResult );
      }
    else
      {
        return NULL;
      }

    bool disableGaussianFitQualityCheck = true;
    CLineDetection lineDetection( lineType, m_cut, m_strongcut, m_winsize, m_minsize, m_maxsize, disableGaussianFitQualityCheck );
    auto lineDetectionResult = lineDetection.Compute( _spc, lambdaRange, peakDetectionResult->PeakList, peakDetectionResult->EnlargedPeakList );

    if( lineDetectionResult ) 
      {
        resultStore.StoreScopedGlobalResult( "raycatalog", lineDetectionResult );

        if( lineDetectionResult->RayCatalog.GetList().size()<1 )
	  {
            //return NULL;
	  }
      }

    // --- Match
    CRayMatching rayMatching;
    auto rayMatchingResult = rayMatching.Compute( lineDetectionResult->RayCatalog, restRayCatalog, redshiftsRange, m_minMatchNum, m_tol, lineType );


    if( rayMatchingResult )
      {
        //rayMatchingResult->FilterWithRules(_spc, lambdaRange, m_winsize);
        // Store matching results
        resultStore.StoreScopedGlobalResult( "raymatching", rayMatchingResult );
      }

    storeResult = true; //always save a matching result
    if( storeResult )
      {
        return std::shared_ptr<const CLineMatching2SolveResult>( new CLineMatching2SolveResult() );
      }

    return NULL;
}

