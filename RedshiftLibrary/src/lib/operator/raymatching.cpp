// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
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
      for( Int32 iRestRay=0; iRestRay<restRayList.size(); iRestRay++ )
      {
        CRayMatchingResult::TSolutionSet solution;
        Float64 redShift = (detectedRayList[0].GetPosition()-restRayList[iRestRay].GetPosition())/restRayList[iRestRay].GetPosition();
        if( redShift > 0 ) // we don't care about blueshift.
          {
            solution.push_back( CRayMatchingResult::SSolution( detectedRayList[0], std::move(restRayList[iRestRay]), redShift ) );
            solutions.push_back( solution );
          }
      }
    }
  else
    {
      Log.LogDebug ( "CRayMatching: Rays detected. n=%d", nDetectedRay );

      for( Int32 iDetectedRay=0; iDetectedRay<detectedRayList.size(); iDetectedRay++ )
      {
          for( Int32 iRestRay=0; iRestRay<restRayList.size(); iRestRay++ )
          {
              // for each detected line / rest line couple, enumerate how many other couples fit within the tolerance
              CRayMatchingResult::TSolutionSet solution;
              Float64 redShift=(detectedRayList[iDetectedRay].GetPosition()-restRayList[iRestRay].GetPosition())/restRayList[iRestRay].GetPosition();
              if( redShift < 0 ) // we don't care about blueshifts.
              {
                  continue;
              }
              solution.push_back( CRayMatchingResult::SSolution( detectedRayList[iDetectedRay], restRayList[iRestRay], redShift) );

              for( Int32 iDetectedRay2=0; iDetectedRay2<detectedRayList.size(); iDetectedRay2++ )
              {
                  if( iDetectedRay!=iDetectedRay2 )
                  {
                      for( Int32 iRestRay2=0; iRestRay2<restRayList.size(); iRestRay2++ )
                      {
                          Float64 redShift2=(detectedRayList[iDetectedRay2].GetPosition()-restRayList[iRestRay2].GetPosition())/restRayList[iRestRay2].GetPosition();
                          if( redShift2 < 0 ) // we don't care about blueshifts.
                          {
                              continue;
                          }
                          Float64 redshiftTolerance = tol*(1+(redShift+redShift2)*0.5);
                          if( fabs( (redShift-redShift2) )<=redshiftTolerance )
                          {
                              bool found = false;
                              //avoid repeated solution sets
                              for( Int32 iSet=0; iSet<solution.size(); iSet++ )
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
  for( Int32 iSol=0; iSol<solutions.size(); iSol++ )
  {
      CRayMatchingResult::TSolutionSet currentSet = solutions[iSol];
      sort( currentSet.begin(), currentSet.end() );
      bool found = false;
      for( Int32 iNewSol=0; iNewSol<newSolutions.size(); iNewSol++ )
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
              for( Int32 i=0; i<currentSet.size(); i++ )
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
      auto result = std::make_shared<CRayMatchingResult>();
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
bool CRayMatching::AreSolutionSetsEqual( const CRayMatchingResult::TSolutionSet& s1, const CRayMatchingResult::TSolutionSet& s2 )
{
  if( s1.size() != s2.size() ) 
    {
      return false;
    }
  bool diffFound = false;
  for( Int32 iSet=0; iSet<s1.size(); iSet++ )
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
