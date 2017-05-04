#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/method/linematching2solve.h>
#include <RedshiftLibrary/method/linematching2solveresult.h>
#include <RedshiftLibrary/operator/chisquare.h>
#include <RedshiftLibrary/operator/correlation.h>
#include <RedshiftLibrary/operator/peakdetection.h>
#include <RedshiftLibrary/operator/peakdetectionresult.h>
#include <RedshiftLibrary/operator/raydetection.h>
#include <RedshiftLibrary/operator/raydetectionresult.h>
#include <RedshiftLibrary/operator/raymatching.h>
#include <RedshiftLibrary/operator/raymatchingresult.h>
#include <RedshiftLibrary/processflow/datastore.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>

using namespace NSEpic;
using namespace std;

/**
 * \brief This constructor will attribute values to this method's parameters with default values.
 */
COperatorLineMatching2Solve::COperatorLineMatching2Solve()
{    
  Log.LogDebug ( "COperatorLineMatchingSolve::COperatorLineMatching2Solve()" );

  // Peak Detection
  m_winsize = 250.0;
  m_cut = 5.0;
  m_detectioncut = 5.0;
  m_detectionnoiseoffset = 0.0; /**\var Useful for simulation data. */
  m_strongcut = 2.0;
  m_minsize = 3;
  m_maxsize = 70;
  m_enlargeRate = 2.0;
  
  // Line Matching
  m_disablegaussianfitqualitycheck = false;
  m_minMatchNum = 1;
  m_tol = 0.002;

  // Log
  m_bypassDebug = true;
}

/**
 * Empty destructor.
 */
COperatorLineMatching2Solve::~COperatorLineMatching2Solve()
{
  Log.LogDebug ( "COperatorLineMatching2Solve::~COperatorLineMatching2Solve()" );
}

/**
 * \brief Attempts to find a matching line using the input parameters.
 * This algorithm will create and run a CPeakDetection object, using the parameters as input.
 *
 * The input spectrum is copied to the local scope.
 *
 * If the lines' type are absorption lines, the spectrum has its fluxes inverted.
 *
 * This local spectrum is passed to a lineDetection object.
 *
 * If the dynamic option is false, the following algorithm will be run:
 *
 * - If peaks are found, it will try to find lines in the peaks found, using a CLineDetection object (please note that this class is defined in raydetectionresult.h).
 *
 * - If lines are found, it will try to match these lines to the catalogue of known lines and determine the redshift using this information.
 *
 * If the dynamic option is true, the line detection is run several times, with the parameters being varied withing physically meaningful values.
 * When either a threshold number of peaks is detected, or all parameters are exhaustively searched, the algorithm continues as normal.
 */
std::shared_ptr<const CLineMatching2SolveResult> COperatorLineMatching2Solve::Compute( CDataStore& resultStore, 
										       const CSpectrum& spc, 
										       const TFloat64Range& lambdaRange, 
										       const TFloat64Range& redshiftsRange, 
										       Float64 redshiftStep, 
										       const CRayCatalog& restRayCatalog )
{
  Log.LogDebug ( "std::shared_ptr<const CLineMatching2SolveResult> COperatorLineMatching2Solve::Compute( CDataStore& resultStore, const CSpectrum& spc, const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep, const CRayCatalog& restRayCatalog )" );

  Int32 disablegaussianfitqualitycheck_value;
  Int32 lineType = CRay::nType_Emission;
  std::string linetypeStr = "E";
  Bool storeResult = false;

  CDataStore::CAutoScope resultScope ( resultStore, "linematching2solve" );

  Log.LogDebug ( "Attempting to load parameters from parameter JSON." );
  {
    resultStore.GetScopedParam( "linematching2.cut", m_cut, 5.0 );
    resultStore.GetScopedParam( "linematching2.detectioncut", m_detectioncut, 5.0 );
    resultStore.GetScopedParam( "linematching2.detectionnoiseoffset", m_detectionnoiseoffset, 0.0 );
    resultStore.GetScopedParam( "linematching2.disablegaussianfitqualitycheck", m_disablegaussianfitqualitycheck, 0 );
    resultStore.GetScopedParam( "linematching2.dynamicLinematching", m_dynamicLinematching, 0 );
    resultStore.GetScopedParam( "linematching2.enlargeRate", m_enlargeRate, 2.0 );
    resultStore.GetScopedParam( "linematching2.linetype", linetypeStr, "Emission" );
    resultStore.GetScopedParam( "linematching2.maxsize", m_maxsize, 70.0 );
    resultStore.GetScopedParam( "linematching2.minMatchNum", m_minMatchNum, 1.0 );
    resultStore.GetScopedParam( "linematching2.minsize", m_minsize, 3.0 );
    resultStore.GetScopedParam( "linematching2.strongcut", m_strongcut, 2.0 );
    resultStore.GetScopedParam( "linematching2.tol", m_tol, 0.002 );
    resultStore.GetScopedParam( "linematching2.winsize", m_winsize, 250.0 );
    lineType = CRay::nType_All;
    if( linetypeStr == "Absorption" )
      {
	lineType = CRay::nType_Absorption;
      }
    if( linetypeStr == "Emission" )
      {
    lineType = CRay::nType_Emission;
      }
  }
  
  Log.LogDebug( "Final parameters read:" );
  {
    Log.LogDebug( "m_cut = %f", m_cut );
    if( m_dynamicLinematching )
      {
	Log.LogDebug( "m_dynamicLinematching is true (any value, except 0, on json)" );
      }
    else
      {
	Log.LogDebug( "m_dynamicLinematching is false (value 0 on json)" );
      }
    Log.LogDebug( "m_detectioncut = %f", m_detectioncut );
    Log.LogDebug( "m_detectionnoiseoffset = %f", m_detectionnoiseoffset );
    if ( m_disablegaussianfitqualitycheck )
      {
	Log.LogDebug( "m_disablegaussianfitqualitycheck is true (any value, except 0, in json)" );
      }
    else
      {
	Log.LogDebug( "m_disablegaussianfitqualitycheck is false (value 0 in json)" );
      }
    Log.LogDebug( "m_enlargeRate = %f", m_enlargeRate );
    Log.LogDebug( "linetype = %s", linetypeStr.c_str ( ) );
    Log.LogDebug( "m_minMatchNum = %d", m_minMatchNum );
    Log.LogDebug( "m_minsize = %f", m_minsize );
    Log.LogDebug( "m_maxsize = %f", m_maxsize );
    Log.LogDebug( "m_strongcut = %f", m_strongcut );
    Log.LogDebug( "m_tol = %f", m_tol );
    Log.LogDebug( "m_winsize = %f", m_winsize );
  }

  CPeakDetection peakDetection ( m_winsize, m_detectioncut, 1, m_enlargeRate, m_detectionnoiseoffset );
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

  // Since we detected at least one peak, try to detect lines related to those peaks.
  Int32 bestMatchingNumber = -1;
  Float64 bestRedshift = -1.0;
  Int32 currentNumberOfPeaks = -1;
  Float64 cutBest = m_cut;
  auto cutCurrent = m_cut;
  Float64 cutMaximum = 1e1;
  Float64 cutMinimum = 1e-1;
  Float64 cutStep = ( cutMaximum - cutMinimum ) / 5;
  Float64 minimumFwhhBest = m_minsize;
  Float64 minimumFwhhCurrent = m_minsize;
  Float64 minimumFwhhMinimum = 1e-3;
  Float64 minimumFwhhMaximum = 1e1;
  Float64 minimumFwhhStep = ( minimumFwhhMaximum - minimumFwhhMinimum ) / 5;
  Int32 minimumNumberOfPeaks = 200;
  Bool newValues = false;
  Int32 numberOfPeaksBest = -1;
  Bool numberOfPeaksBestBypass = false;
  Int32 previousNumberOfPeaks = -1;
  Float64 strongcutBest = m_strongcut;
  auto strongcutCurrent = m_strongcut;
  Float64 strongcutMaximum = 1e1;
  Float64 strongcutMinimum = 1e-2;
  Float64 strongcutStep = ( strongcutMaximum - strongcutMinimum ) / 5;
  Float64 winsizeBest = m_winsize;
  auto winsizeCurrent = m_winsize;
  Float64 winsizeMaximum = 1e4;
  Float64 winsizeMinimum = 1e3;
  Float64 winsizeStep = ( winsizeMaximum - winsizeMinimum ) / 5;
  
  Int32 cmptMax = 100;
  Int32 iCmpt = 0;
  while ( iCmpt<cmptMax )
  {
      if ( ! m_bypassDebug )
      {
          Log.LogDebug( "*************************" );
          Log.LogDebug( "cutCurrent == %f", cutCurrent );
          Log.LogDebug( "minimumFwhhCurrent == %f", minimumFwhhCurrent );
          Log.LogDebug( "strongcutCurrent == %f", strongcutCurrent );
          Log.LogDebug( "winsizeCurrent == %f", winsizeCurrent );
      }
      CLineDetection lineDetection( lineType, cutCurrent, strongcutCurrent, winsizeCurrent, minimumFwhhCurrent, m_maxsize, m_disablegaussianfitqualitycheck );
      auto lineDetectionResult = lineDetection.Compute( _spc, lambdaRange, peakDetectionResult->PeakList, peakDetectionResult->EnlargedPeakList );
      if ( lineDetectionResult || ! m_dynamicLinematching)
      {
          currentNumberOfPeaks = lineDetectionResult->RayCatalog.GetList().size();
          if ( ! m_bypassDebug )
              Log.LogDebug ( "Found %d peaks.", currentNumberOfPeaks );
          if( currentNumberOfPeaks>=minimumNumberOfPeaks || numberOfPeaksBestBypass || ! m_dynamicLinematching )
          {
              Log.LogDebug ( "Storing %d lines from lineDetection in the result store.", lineDetectionResult->RayCatalog.GetList().size() );
              resultStore.StoreScopedGlobalResult( "raycatalog", lineDetectionResult );
              // Since we know at least one peak that corresponds to a line, let's try to match to a catalogued template.
              CRayMatching rayMatching;
              auto rayMatchingResult = rayMatching.Compute( lineDetectionResult->RayCatalog, restRayCatalog, redshiftsRange, m_minMatchNum, m_tol, lineType );
              if( rayMatchingResult )
              {
                  rayMatchingResult->FilterWithRules( _spc, lambdaRange, m_winsize );
                  Log.LogDebug ( "CRayMatching yielded %d sets of solutions and %d sets of filtered solutions.", rayMatchingResult->SolutionSetList.size(), rayMatchingResult->FilteredSolutionSetList.size() );
                  // Store matching results
                  resultStore.StoreScopedGlobalResult( "raymatching", rayMatchingResult );
                  rayMatchingResult->GetBestRedshift ( bestRedshift, bestMatchingNumber );
                  Log.LogDebug ( "bestRedshift == %f, bestMatchingNumber == %d", bestRedshift, bestMatchingNumber );
                  if( bestRedshift != -1.0 )
                  {
                      Log.LogDebug ( "return std::shared_ptr<const CLineMatching2SolveResult>( new CLineMatching2SolveResult() );" );
                      return std::shared_ptr<const CLineMatching2SolveResult>( new CLineMatching2SolveResult() );
                  }else if(! m_dynamicLinematching)
                  {
                      return std::shared_ptr<const CLineMatching2SolveResult>( new CLineMatching2SolveResult() );
                  }
              } // rayMatchingResult
              else if(! m_dynamicLinematching)
              {
                  return std::shared_ptr<const CLineMatching2SolveResult>( new CLineMatching2SolveResult() );
              }
          } // minimumNumberOfPeaks
      } // lineDetectionResult
      // update parameters
      {
          if( currentNumberOfPeaks>=numberOfPeaksBest )
          {
              numberOfPeaksBest = currentNumberOfPeaks;
              cutBest = cutCurrent;
              minimumFwhhBest = minimumFwhhCurrent;
              strongcutBest = strongcutCurrent;
              winsizeBest = winsizeCurrent;
          }
          if( currentNumberOfPeaks<previousNumberOfPeaks && ! newValues )
          {
              Log.LogError ( "Dynamic logic failed - number of peaks is falling." );
              return NULL;
          }
          else
          {
              previousNumberOfPeaks = currentNumberOfPeaks;
              if( ! m_dynamicLinematching )
              {
                  Log.LogError ( "No result found - returning NULL." );
                  return NULL;
              }
              newValues = false;
              if( minimumFwhhCurrent>minimumFwhhMinimum )
              {
                  minimumFwhhCurrent = std::max ( minimumFwhhMinimum, minimumFwhhCurrent - minimumFwhhStep );
                  continue;
              }
              newValues = true;
              minimumFwhhCurrent = minimumFwhhMaximum;
              if( cutCurrent>=cutMinimum )
              {
                  cutCurrent -= cutStep;
                  continue;
              }
              cutCurrent = strongcutCurrent;
              if( strongcutCurrent>strongcutMinimum )
              {
                  strongcutCurrent = std::max ( strongcutMinimum, strongcutCurrent - strongcutStep );
                  continue;
              }
              strongcutCurrent = strongcutMaximum;
              if( winsizeCurrent<=winsizeMaximum )
              {
                  winsizeCurrent += winsizeStep;
                  continue;
              }
              winsizeCurrent = winsizeMinimum;
              // if all checks failed, dynamic system failed.
              Log.LogDebug ( "Dynamic logic: all parameters searched. Using best values." );
              {
                  cutCurrent = cutBest;
                  minimumFwhhCurrent = minimumFwhhBest;
                  strongcutCurrent = strongcutBest;
                  winsizeCurrent = winsizeBest;
                  numberOfPeaksBestBypass = true;
                  continue;
              }
          }
      }
      Log.LogDebug( "Dynamically retrying linematching." );
      iCmpt++;
  }; // while

  if( iCmpt==cmptMax )
  {
      Log.LogWarning( "Warning. Stopped the linematching dynamic cut loop..." );
  }
  return std::shared_ptr<const CLineMatching2SolveResult>( new CLineMatching2SolveResult() );
}

const std::string COperatorLineMatching2Solve::GetDescription()
{
    std::string desc;
    desc = "Method linematching2:\n";
    desc.append("\tparam: linematching2.cut = <float value>\n");
    desc.append("\tparam: linematching2.detectioncut = <float value>\n");
    desc.append("\tparam: linematching2.detectionnoiseoffset = <float value>\n");
    desc.append("\tparam: linematching2.disablegaussianfitqualitycheck = <integer value>\n");
    desc.append("\tparam: linematching2.dynamicLinematching = <integer value>\n");
    desc.append("\tparam: linematching2.enlargeRate = <float value>\n");
    desc.append("\tparam: linematching2.linetype = {""Emission"", ""Absorption"", ""All""}\n");
    desc.append("\tparam: linematching2.maxsize = <float value in Angstrom>\n");
    desc.append("\tparam: linematching2.minMatchNum = <int value>\n");
    desc.append("\tparam: linematching2.minsize = <float value in Angstrom>\n");
    desc.append("\tparam: linematching2.strongcut = <float value>\n");
    desc.append("\tparam: linematching2.tol = <float value>\n");
    desc.append("\tparam: linematching2.winsize = <float value in Angstrom>\n");

    return desc;
}
