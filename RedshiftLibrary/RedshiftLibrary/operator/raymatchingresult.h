#ifndef _REDSHIFT_OPERATOR_RAYMATCHINGRESULT_
#define _REDSHIFT_OPERATOR_RAYMATCHINGRESULT_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/ray/catalog.h>
#include <RedshiftLibrary/spectrum/spectrum.h>

#include <vector>

namespace NSEpic
{
  /**
   * \ingroup Redshift
   * Holds the data corresponding to a match between two sets of lines. Commonly, holds the result of matching detected lines with template lines.
   */
  class CRayMatchingResult : public COperatorResult
  {
  public:
  
    struct SSolution
    {
      CRay DetectedRay;
      CRay RestRay;
      Float64 Redshift;

      SSolution( CRay detectedRay, CRay restRay, Float64 redshift )
      {
	DetectedRay = detectedRay;
	RestRay = restRay;
	Redshift = redshift;
      }

      bool operator<( const SSolution& str ) const
      {
	if( DetectedRay.GetPosition()==str.DetectedRay.GetPosition() )
	  {
	    return ( RestRay.GetPosition()<str.RestRay.GetPosition() );
	  }
	else
	  {
	    return ( DetectedRay.GetPosition()<str.DetectedRay.GetPosition() );
	  }
      }  
    };

    typedef std::vector<SSolution> TSolutionSet; // a set of (detected line,rest line) couples for a given redshift
    typedef std::vector<TSolutionSet> TSolutionSetList; // a list of possible redshift solutions

    CRayMatchingResult();
    virtual ~CRayMatchingResult();

    void Save(std::ostream& stream ) const;
    void SaveLine(std::ostream& stream ) const;
    
    void SaveSolutionSetToStream( std::ostream& stream, TSolutionSetList selectedResults, Int32 type) const;

    Bool GetBestRedshift( Float64& Redshift, Int32& MatchingNumber ) const;
    Bool GetBestMatchNumRedshift( Float64& Redshift, Int32& MatchingNumber ) const;

    Int32 getNStrongRestLines( const TSolutionSet& s ) const;

    Int32 GetMaxMatchingNumber() const;
    Float64 GetMeanRedshiftSolution( const TSolutionSet& s ) const;
    Float64 GetMeanRedshiftSolutionByIndex( Int32 index ) const;
    TSolutionSetList GetSolutionsListOverNumber( Int32 number ) const;
    TFloat64List GetAverageRedshiftListOverNumber( Int32 number ) const;

    TFloat64List GetRoundedRedshiftCandidatesOverNumber( Int32 number, Float64 step ) const;
    TFloat64List GetExtendedRedshiftCandidatesOverNumber( Int32 number, Float64 step, Float64 rangeWidth ) const;

    void FilterWithRules( CSpectrum spc, TFloat64Range lambdaRange, Float64 winsize );

    TSolutionSetList SolutionSetList;
    TSolutionSetList FilteredSolutionSetList;
    std::vector<int> FilterTypeList;

    CRayCatalog m_RestCatalog;
    CRayCatalog m_DetectedCatalog;
  private:
    Bool m_bypassDebug;
  };
}

#endif // _REDSHIFT_OPERATOR_RAYMATCHINGRESULT_
