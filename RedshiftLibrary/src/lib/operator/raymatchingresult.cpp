#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/operator/raymatchingresult.h"
#include "RedshiftLibrary/ray/rules.h"
#include <stdio.h>
#include <math.h>
#include <algorithm>    // std::sort

using namespace NSEpic;

/**
 * Attribution constructor. Will only set the debugging bypass option.
 */
CRayMatchingResult::CRayMatchingResult()
{
  m_bypassDebug = true;
}

/**
 * Empty destructor.
 */
CRayMatchingResult::~CRayMatchingResult()
{

}


/**
 * @brief CRayMatchingResult::GetBestRedshift
 *
 * select the best redshift given a set of rules:
 * - the solution with the most strong lines
 * - (TODO) if no solution with strong lines, the missing strong lines' absence should be excused by:
 *      * high noise in the theoretical position of the strong lines
 *      * lambda range not able to cover the strong lines
 * @return
 */
Bool CRayMatchingResult::GetBestRedshift( Float64& Redshift, Int32& MatchingNumber ) const
{
  Int32 thresMatchingNumber = 2; //minimum matching number for a solution
  TSolutionSetList selectedResults = GetSolutionsListOverNumber( thresMatchingNumber-1 );
  
  if( selectedResults.size()>0 )
    {
      int iStrongMax = -1;
      int nStrongMax = -1;
      for( UInt32 iSol=0; iSol<selectedResults.size(); iSol++ )
	{
	  int currentNStrongRestLines = getNStrongRestLines( selectedResults[iSol] );
	  if( currentNStrongRestLines>nStrongMax )
	    {
	      iStrongMax = iSol;
	      nStrongMax = currentNStrongRestLines;
	    }
	  else 
	    {
	      if( currentNStrongRestLines==nStrongMax && selectedResults[iSol].size()>selectedResults[iStrongMax].size() )
		{
		  iStrongMax = iSol;
		  nStrongMax = currentNStrongRestLines;
		}
	    }
	}
      
      if( iStrongMax>-1 )
	{
	  Redshift = GetMeanRedshiftSolution( selectedResults[iStrongMax] );
	  MatchingNumber = selectedResults[iStrongMax].size();
	}
    }
  else
    {
      return false;
    }
   return true;
}

/**
 * @brief CRayMatchingResult::GetBestMatchNumRedshift
 *
 * get the best redshift in the sense of the highest matching number
 *
 * @return
 */
Bool CRayMatchingResult::GetBestMatchNumRedshift( Float64& Redshift, Int32& MatchingNumber ) const
{
  MatchingNumber = GetMaxMatchingNumber();
  TSolutionSetList selectedResults = GetSolutionsListOverNumber( MatchingNumber-1 );
  
  if( selectedResults.size()>0 )
    {
      Redshift = GetMeanRedshiftSolution( selectedResults[0] );
    }
  else
    {
      return false;
    }
  return true;
}

/**
 * Returns the set of results that have more than "number" entries.
 */
CRayMatchingResult::TSolutionSetList CRayMatchingResult::GetSolutionsListOverNumber( Int32 number ) const
{
  //select results by matching number
  TSolutionSetList selectedResults;
  for( UInt32 iSol=0; iSol<SolutionSetList.size(); iSol++ )
    {
      TSolutionSet currentSet = SolutionSetList[iSol];
      if( currentSet.size()>number )
	{
	  selectedResults.push_back( currentSet );
	}
    }
  return selectedResults;
}

/**
 * Returns a list of the mean redshifts for solutions that have at least "number" lines.
 */
TFloat64List CRayMatchingResult::GetAverageRedshiftListOverNumber( Int32 number ) const
{
  TFloat64List selectedRedshift;
  CRayMatchingResult::TSolutionSetList selectedResults = GetSolutionsListOverNumber( number );
  for( UInt32 j=0; j<selectedResults.size(); j++ )
    {
      Float64 z = GetMeanRedshiftSolution( selectedResults[j] );
      selectedRedshift.push_back( z );
    }
  return selectedRedshift;
}

/**
 * Returns a custom-computed round value for the mean redshift for each solution that has over "number" lines.
 */
TFloat64List CRayMatchingResult::GetRoundedRedshiftCandidatesOverNumber( Int32 number, Float64 step ) const
{
  TFloat64List selectedRedshift = GetAverageRedshiftListOverNumber( number );
  TFloat64List roundedRedshift;
  for( UInt32 j=0; j<selectedRedshift.size(); j++ )
    {
      Float64 zround = Float64( int( selectedRedshift[j]/step+0.5f )*step );
      roundedRedshift.push_back( zround );
    }
  return roundedRedshift;
}

/**
 * Retuns a list of redshifts that are "step"pedly within the "rangeWidth" interval around each mean redshift for solutions with over "number" lines.
 */
TFloat64List CRayMatchingResult::GetExtendedRedshiftCandidatesOverNumber( Int32 number, Float64 step, Float64 rangeWidth ) const
{
  TFloat64List roundedRedshift = GetRoundedRedshiftCandidatesOverNumber( number, step );
  TFloat64List extendedRedshifts;
  Int32 halfk = rangeWidth/step/2.0;
  for( UInt32 j=0; j<roundedRedshift.size(); j++ )
    {
      for( Int32 k=-halfk; k<halfk; k++ )
	{
	  Float64 z = roundedRedshift[j]+k*step;
	  extendedRedshifts.push_back( z );
	}
    }
  std::sort( extendedRedshifts.begin(), extendedRedshifts.end() );
  extendedRedshifts.erase( std::unique( extendedRedshifts.begin(), extendedRedshifts.end() ), extendedRedshifts.end() );
  return extendedRedshifts;
}

/**
 * If "index" is valid, return the mean of the redshifts in the solution set at "index".
 */
Float64 CRayMatchingResult::GetMeanRedshiftSolutionByIndex( Int32 index ) const
{
  if( index>SolutionSetList.size()-1 )
    {
      return -1.0;
    }

  Float64 redshiftMean = 0.0;
  const TSolutionSet& currentSet = SolutionSetList[index];
  for( UInt32 i=0; i<currentSet.size(); i++ )
    {
      redshiftMean += currentSet[i].Redshift;
    }
  redshiftMean /= currentSet.size();
  return redshiftMean;
}

/**
 * Returns the mean redshift for the solution set "s".
 */
Float64 CRayMatchingResult::GetMeanRedshiftSolution( const TSolutionSet& s ) const
{
  TSolutionSet currentSet = s;
  Float64 redshiftMean=0.0;
  for( UInt32 i=0; i<currentSet.size(); i++ )
    {
      redshiftMean += currentSet[i].Redshift;
    }
  redshiftMean /= (float)currentSet.size();

  return redshiftMean;
}

/**
 * Returns the number of lines that are near (tolerance = 0.11 Angstroms) to strong emission lines.
 */
Int32 CRayMatchingResult::getNStrongRestLines( const TSolutionSet& s ) const
{
  CRayCatalog::TRayVector strongRestRayList = m_RestCatalog.GetFilteredList( CRay::nType_Emission, CRay::nForce_Strong );
  Int32 ncatalog = strongRestRayList.size();

  TSolutionSet currentSet = s;
  Int32 nStrong=0;
  Float64 tol = 0.11;
  for( UInt32 i=0; i<currentSet.size(); i++ )
    {
      Int32 found  = 0;
      for( UInt32 c=0; c<ncatalog; c++ )
	{
	  if( fabs( currentSet[i].RestRay.GetPosition()-strongRestRayList[c].GetPosition() )<tol )
	    {
	      found = 1;
	      break;
	    }
	}
      if( found==1 )
	{
	  nStrong++;
	}
    }

  return nStrong;
}

/**
 * Returns the largest number of lines in a single solution for the current set of solutions.
 */
Int32 CRayMatchingResult::GetMaxMatchingNumber() const
{
  if( SolutionSetList.size()<1 )
    {
      return -1.0;
    }

  Int32 maxNumber = 0;
  for( UInt32 i=0; i<SolutionSetList.size(); i++ )
    {
      TSolutionSet currentSet = SolutionSetList[i];
      if( maxNumber<currentSet.size() )
	{
	  maxNumber = currentSet.size();
	}
    }
  return maxNumber;
}

/**
 * Attempts to apply all rules to each solution. If a rule is applicable, that solution is removed from the current set of solutions.
 */
void CRayMatchingResult::FilterWithRules( CSpectrum spc, TFloat64Range lambdaRange, Float64 winsize )
{
  Log.LogDebug ( "void CRayMatchingResult::FilterWithRules( CSpectrum spc, TFloat64Range lambdaRange, Float64 winsize )" );
  if( SolutionSetList.size()<1 )
    {
      Log.LogDebug ( "SolutionSetList.size()<1, returning" );
      return;
    }

  CRules rules( spc, m_DetectedCatalog, m_RestCatalog, lambdaRange, winsize );
  TSolutionSetList _solutionSetList;
  Log.LogDebug( "There are %d solutions to test.", SolutionSetList.size( ) );
  for( UInt32 i=0; i<SolutionSetList.size() ; i++ )
    {
      TSolutionSet currentSet = SolutionSetList[i];
      Float64 z = GetMeanRedshiftSolution( currentSet ); 
      Int32 ruleId = rules.check( z, currentSet );

      if( ruleId<=0 )
	{
	  if ( ! m_bypassDebug )
	    Log.LogDebug( "Solution %d passed!", i );
	  _solutionSetList.push_back( currentSet );
	}
      else
	{
	  if ( ! m_bypassDebug )
	    Log.LogDebug( "Solution %d failed rule %d", i, ruleId );
	  FilteredSolutionSetList.push_back( currentSet );
	  FilterTypeList.push_back( ruleId );
	}
    }
  SolutionSetList = _solutionSetList;
}
