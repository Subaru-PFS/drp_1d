#include <epic/core/log/log.h>
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
  Log.LogDebug ( "COperatorLineMatchingSolve::COperatorLineMatching2Solve()" );

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
  
  if( 0 )
    { 
      Log.LogWarning ( "Using hardcoded 'PFS Really just lines, F + ErrF' overrides." );
      m_winsize = 250.0;
      m_cut = 1.0;
      m_detectioncut = 1.0;
      m_maxsize = 120;
      m_enlargeRate = 2.0;
      m_tol = 0.0001;
    }

  if( false )
    {
      Log.LogWarning ( "Using hardcoded 'PFS batch5' default parameters." );
      m_cut = 1.0;
    }
}

/**
 * Empty destructor.
 */
COperatorLineMatching2Solve::~COperatorLineMatching2Solve()
{
  Log.LogDebug ( "COperatorLineMatching2Solve::~COperatorLineMatching2Solve()" );
}

/**
 * Attempts to find a matching line using the input parameters.
 * This algorithm will create and run a CPeakDetection object, using the parameters as input.
 * If peaks are found, it will try to find lines in the peaks found, using a CLineDetection object (warning: this class is defined in raydetectionresult.h).
 * If lines are found, it will try to match catalogued templates to the found lines.
 */
std::shared_ptr<const CLineMatching2SolveResult> COperatorLineMatching2Solve::Compute( CDataStore& resultStore, 
										       const CSpectrum& spc, 
										       const TFloat64Range& lambdaRange, 
										       const TFloat64Range& redshiftsRange, 
										       Float64 redshiftStep, 
										       const CRayCatalog& restRayCatalog )
{
  Log.LogDebug ( "std::shared_ptr<const CLineMatching2SolveResult> COperatorLineMatching2Solve::Compute( CDataStore& resultStore, const CSpectrum& spc, const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep, const CRayCatalog& restRayCatalog )" );
  
  Bool storeResult = false;
  CDataStore::CAutoScope resultScope( resultStore, "linematching2solve" );

  Log.LogDebug ( "Attempting to load parameters from parameter JSON." );
  {
    resultStore.GetScopedParam( "linematching2.winsize", m_winsize, 250.0 );
    resultStore.GetScopedParam( "linematching2.cut", m_cut, 5.0 );
    resultStore.GetScopedParam( "linematching2.detectioncut", m_detectioncut, 5.0 );
    resultStore.GetScopedParam( "linematching2.detectionnoiseoffset", m_detectionnoiseoffset, 0.0 );
    resultStore.GetScopedParam( "linematching2.strongcut", m_strongcut, 2.0 );
    resultStore.GetScopedParam( "linematching2.minsize", m_minsize, 3.0 );
    resultStore.GetScopedParam( "linematching2.maxsize", m_maxsize, 70.0 );
    resultStore.GetScopedParam( "linematching2.enlargeRate", m_enlargeRate, 2.0 );
    resultStore.GetScopedParam( "linematching2.minMatchNum", m_minMatchNum, 1.0 );
    resultStore.GetScopedParam( "linematching2.tol", m_tol, 0.002 );
    Log.LogDebug ( "m_winsize = %f", m_winsize );
    Log.LogDebug ( "m_cut = %f", m_cut );
    Log.LogDebug ( "m_detectioncut = %f", m_detectioncut );
    Log.LogDebug ( "m_detectionnoiseoffset = %f", m_detectionnoiseoffset );
    Log.LogDebug ( "m_strongcut = %f", m_strongcut );
    Log.LogDebug ( "m_minsize = %f", m_minsize );
    Log.LogDebug ( "m_maxsize = %f", m_maxsize );
    Log.LogDebug ( "m_enlargeRate = %f", m_enlargeRate );
    Log.LogDebug ( "m_minMatchNum = %d", m_minMatchNum );
    Log.LogDebug ( "m_tol = %f", m_tol );
  }

  Int32 lineType = CRay::nType_Emission;
  //Int32 lineType = CRay::nType_Absorption;

  CPeakDetection peakDetection( m_winsize, m_detectioncut, 1, m_enlargeRate, m_detectionnoiseoffset );
  CSpectrum _spc = spc;
  if( lineType == CRay::nType_Absorption )
    {
      Log.LogDebug ( "lineType == CRay::nType_Absorption" );
      _spc.InvertFlux();
    }

  if ( lineType == CRay::nType_Emission )
    {
      Log.LogDebug ( "lineType == CRay::nType_Emission" );
    }

  auto peakDetectionResult = peakDetection.Compute( _spc, lambdaRange );
  if( peakDetectionResult )
    {
      Log.LogDebug ( "Storing %d peaks from PeakList in the result store.", peakDetectionResult->PeakList.size() );
      resultStore.StoreScopedGlobalResult( "peakdetection", peakDetectionResult );
    }
  else
    {
      Log.LogInfo ( "No peak detected - returning NULL." );
      return NULL;
    }

  // Since we detected at least one peak, try to detect lines related that correspond to those peaks.
  bool disableGaussianFitQualityCheck = false;
  if ( disableGaussianFitQualityCheck )
    {
      Log.LogDebug ( "Fit quality check disabled." );
    }
  else
    {
      Log.LogDebug ( "Fit quality check enabled." );
    }

  CLineDetection lineDetection( lineType, m_cut, m_strongcut, m_winsize, m_minsize, m_maxsize, disableGaussianFitQualityCheck );
  auto lineDetectionResult = lineDetection.Compute( _spc, lambdaRange, peakDetectionResult->PeakList, peakDetectionResult->EnlargedPeakList );
  
  if( lineDetectionResult ) 
    {
      Log.LogDebug ( "Storing %d lines from lineDetection in the result store.", lineDetectionResult->RayCatalog.GetList().size() );
      resultStore.StoreScopedGlobalResult( "raycatalog", lineDetectionResult );
      
      if( lineDetectionResult->RayCatalog.GetList().size()<1 )
	{
	  Log.LogError ( "Line detection returned no entries in RayCatalog." );
	  //return NULL;
	}
    }

    // Since we know at least one peak that corresponds to a line, let's try to match to a catalogued template.
    CRayMatching rayMatching;
    auto rayMatchingResult = rayMatching.Compute( lineDetectionResult->RayCatalog, restRayCatalog, redshiftsRange, m_minMatchNum, m_tol, lineType );

    if( rayMatchingResult )
      {
	rayMatchingResult->FilterWithRules( _spc, lambdaRange, m_winsize );
	Log.LogDebug ( "CRayMatching yielded %d sets of solutions and %d sets of filtered solutions.", rayMatchingResult->SolutionSetList.size(), rayMatchingResult->FilteredSolutionSetList.size() );
        // Store matching results
        resultStore.StoreScopedGlobalResult( "raymatching", rayMatchingResult );
      }

    storeResult = true; //always save a matching result
    if( storeResult )
      {
	Log.LogDebug ( "return std::shared_ptr<const CLineMatching2SolveResult>( new CLineMatching2SolveResult() );" );
        return std::shared_ptr<const CLineMatching2SolveResult>( new CLineMatching2SolveResult() );
      }

    Log.LogError ( "No result found - returning NULL." );
    return NULL;
}

const std::string COperatorLineMatching2Solve::GetDescription()
{
    std::string desc;
    desc = "Method linematching2:\n";
    desc.append("\tparam: linematching2.winsize = <float value in Angstrom>\n");
    desc.append("\tparam: linematching2.cut = <float value>\n");
    desc.append("\tparam: linematching2.detectioncut = <float value>\n");
    desc.append("\tparam: linematching2.detectionnoiseoffset = <float value>\n");
    desc.append("\tparam: linematching2.strongcut = <float value>\n");
    desc.append("\tparam: linematching2.minsize = <float value in Angstrom>\n");
    desc.append("\tparam: linematching2.maxsize = <float value in Angstrom>\n");
    desc.append("\tparam: linematching2.enlargeRate = <float value>\n");
    desc.append("\tparam: linematching2.minMatchNum = <int value>\n");
    desc.append("\tparam: linematching2.tol = <float value>\n");

    return desc;
}
