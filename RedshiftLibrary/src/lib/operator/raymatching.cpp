#include "RedshiftLibrary/operator/raymatching.h"

#include <algorithm>    // std::sort
#include <math.h>
#include <stdio.h>
#include "RedshiftLibrary/log/log.h"


using namespace NSEpic;

/**
 * This class's initial implementation duplicates the VIPGI method to match lines found with catalog lines, as in EZELmatch.py, Line 264
 */
CRayMatching::CRayMatching()
{

}

/**
 * Empty destructor.
 */
CRayMatching::~CRayMatching()
{

}

/**
 * Executes the algorithm for matching the input sets of lines against the catalogue of atomic transitions.
 * Parameters:
 * - detectedRayCatalog, the set of lines detected in the spectrum.
 * - restRayCatalog, the set of reference lines.
 * - redshiftRange.
 * - nThreshold.
 * - tol.
 * - typeFilter.
 * - detectedForceFilter.
 * - restForceFilter.
 * Returns a smart pointer of a CRayMatchingResult, containing the set of matches.
 */
std::shared_ptr<CRayMatchingResult> CRayMatching::Compute( const CRayCatalog& detectedRayCatalog, const CRayCatalog& restRayCatalog, const TFloat64Range& redshiftRange, Int32 nThreshold, Float64 tol , Int32 typeFilter, Int32 detectedForceFilter, Int32 restForceFilter )
{
  CRayCatalog::TRayVector detectedRayList = detectedRayCatalog.GetFilteredList( typeFilter, detectedForceFilter );
  CRayCatalog::TRayVector restRayList = restRayCatalog.GetFilteredList( typeFilter, restForceFilter );
  
  CRayMatchingResult::TSolutionSetList solutions;
  Int32 nDetectedRay = detectedRayList.size();
  if( nDetectedRay == 1 )
    {
      Log.LogDebug ( "CRayMatching: Only 1 ray detected. n=%d", nDetectedRay );
      // if there is only one detected line, then there are N=#restlines possible redshifts
      for( UInt32 iRestRay=0; iRestRay<restRayList.size(); iRestRay++ )
        {
	  CRayMatchingResult::TSolutionSet solution;
	  Float64 redShift = (detectedRayList[0].GetPosition()-restRayList[iRestRay].GetPosition())/restRayList[iRestRay].GetPosition();
	  if( redShift > 0 ) // we don't care about blueshift.
	    {
	      solution.push_back( CRayMatchingResult::SSolution( detectedRayList[0], restRayList[iRestRay], redShift ) );
	      solutions.push_back( solution );
	    }
        }
    }
  else
    {
      Log.LogDebug ( "CRayMatching: Rays detected. n=%d", nDetectedRay );

      for( UInt32 iDetectedRay=0; iDetectedRay<detectedRayList.size(); iDetectedRay++ )
      {
          for( UInt32 iRestRay=0; iRestRay<restRayList.size(); iRestRay++ )
          {
              // for each detected line / rest line couple, enumerate how many other couples fit within the tolerance
              CRayMatchingResult::TSolutionSet solution;
              Float64 redShift=(detectedRayList[iDetectedRay].GetPosition()-restRayList[iRestRay].GetPosition())/restRayList[iRestRay].GetPosition();
              if( redShift < 0 ) // we don't care about blueshifts.
              {
                  continue;
              }
              solution.push_back( CRayMatchingResult::SSolution( detectedRayList[iDetectedRay], restRayList[iRestRay], redShift) );

              for( UInt32 iDetectedRay2=0; iDetectedRay2<detectedRayList.size(); iDetectedRay2++ )
              {
                  if( iDetectedRay!=iDetectedRay2 )
                  {
                      for( UInt32 iRestRay2=0; iRestRay2<restRayList.size(); iRestRay2++ )
                      {
                          Float64 redShift2=(detectedRayList[iDetectedRay2].GetPosition()-restRayList[iRestRay2].GetPosition())/restRayList[iRestRay2].GetPosition();
                          if( redShift2 < 0 ) // we don't care about blueshifts.
                          {
                              continue;
                          }
                          Float64 redshiftTolerance = tol*(1+(redShift+redShift2)*0.5);
                          if( fabs( (redShift-redShift2) )<=redshiftTolerance )
                          {
                              Bool found = false;
                              //avoid repeated solution sets
                              for( UInt32 iSet=0; iSet<solution.size(); iSet++ )
                              {
                                  if( solution[iSet].DetectedRay.GetPosition()==detectedRayList[iDetectedRay2].GetPosition() )
                                  {
                                      found = true;
                                      break;
                                  }
                              }
                              if( !found )
                              {
                                  solution.push_back( CRayMatchingResult::SSolution( detectedRayList[iDetectedRay2], restRayList[iRestRay2], redShift2) );
                              }
                          }
                      } // for
                  }
              } // for
              if( solution.size()>0 )
              {
                  solutions.push_back( solution );
              }
          }
      }
  }
  Log.LogDebug ( "CRayMatching: non unique solutions found n=%d", solutions.size() );

  // delete duplicated solutions
  Int32 nSolFoundLTThres = 0;
  Int32 nSolFoundOutsideRedshiftRange = 0;
  Int32 nSolFoundDuplicated = 0;
  CRayMatchingResult::TSolutionSetList newSolutions;
  for( UInt32 iSol=0; iSol<solutions.size(); iSol++ )
  {
      CRayMatchingResult::TSolutionSet currentSet = solutions[iSol];
      sort( currentSet.begin(), currentSet.end() );
      Bool found = false;
      for( UInt32 iNewSol=0; iNewSol<newSolutions.size(); iNewSol++ )
      {
          if( AreSolutionSetsEqual( newSolutions[iNewSol],currentSet ) )
          {
              found = true;
          }
      }
      if( !found ) // if solution is new, let's check it is within the range
      {
          if( currentSet.size()>=nThreshold )
          {
              Float64 redshiftMean = 0.f;
              for( UInt32 i=0; i<currentSet.size(); i++ )
              {
                  redshiftMean += currentSet[i].Redshift;
              }
              redshiftMean /= currentSet.size();
              if( redshiftMean>redshiftRange.GetBegin() && redshiftMean<redshiftRange.GetEnd() )
              {
                  newSolutions.push_back( currentSet );
              }else{
                  nSolFoundOutsideRedshiftRange++;
              }
          }else{
              nSolFoundLTThres++;
          }
      }else{
          nSolFoundDuplicated++;
      }
  }
  Log.LogDebug ( "CRayMatching: Cropped (Duplicate) solutions found n=%d", nSolFoundDuplicated );
  Log.LogDebug ( "CRayMatching: Cropped (Outside zrange [ %f ; %f ]) solutions found n=%d",
                 redshiftRange.GetBegin(),
                 redshiftRange.GetEnd(),
                 nSolFoundOutsideRedshiftRange );
  Log.LogDebug ( "CRayMatching: Cropped (Lower than thres.) solutions found n=%d", nSolFoundLTThres );
  Log.LogDebug ( "CRayMatching: unique solutions found n=%d", newSolutions.size() );
  
  if( newSolutions.size()>0 ) 
    {
      auto result = std::shared_ptr<CRayMatchingResult>( new CRayMatchingResult() );
      result->SolutionSetList = newSolutions;
      result->m_RestCatalog = restRayCatalog;
      result->m_DetectedCatalog = restRayCatalog;

      return result;
    }
  return NULL;
}

/**
 * Given 2 TSolutionSets, returns true if they are equivalent (have the same line values), false otherwise.
 * The algorithm only works reliably if the inputs are sorted.
 */
Bool CRayMatching::AreSolutionSetsEqual( const CRayMatchingResult::TSolutionSet& s1, const CRayMatchingResult::TSolutionSet& s2 )
{
  if( s1.size() != s2.size() ) 
    {
      return false;
    }
  Bool diffFound = false;
  for( UInt32 iSet=0; iSet<s1.size(); iSet++ )
    {
      if( s1[iSet].DetectedRay!=s2[iSet].DetectedRay || s1[iSet].RestRay!=s2[iSet].RestRay || s1[iSet].Redshift!=s2[iSet].Redshift )
	{
	  diffFound = true;
	  return false;
        }
    }

  if( !diffFound )
    {
      return true;
    }
  return false;
}
