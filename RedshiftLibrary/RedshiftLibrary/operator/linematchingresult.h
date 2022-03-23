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
#ifndef _REDSHIFT_OPERATOR_RAYMATCHINGRESULT_
#define _REDSHIFT_OPERATOR_RAYMATCHINGRESULT_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/processflow/result.h"
#include "RedshiftLibrary/line/catalog.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

#include <vector>

namespace NSEpic
{
  /**
   * \ingroup Redshift
   * Holds the data corresponding to a match between two sets of lines. Commonly, holds the result of matching detected lines with template lines.
   */
  class CLineMatchingResult : public COperatorResult
  {
  public:
  
    struct SSolution
    {
      CLine DetectedLine;
      CLine RestLine;
      Float64 Redshift;

      SSolution( CLine detectedLine, CLine restLine, Float64 redshift )
      {
        DetectedLine = std::move(detectedLine);
        RestLine = std::move(restLine);
        Redshift = redshift;
      }

      bool operator<( const SSolution& str ) const
      {
	if( DetectedLine.GetPosition()==str.DetectedLine.GetPosition() )
	  {
	    return ( RestLine.GetPosition()<str.RestLine.GetPosition() );
	  }
	else
	  {
	    return ( DetectedLine.GetPosition()<str.DetectedLine.GetPosition() );
	  }
      }  
    };

    typedef std::vector<SSolution> TSolutionSet; // a set of (detected line,rest line) couples for a given redshift
    typedef std::vector<TSolutionSet> TSolutionSetList; // a list of possible redshift solutions

    CLineMatchingResult();
    virtual ~CLineMatchingResult();

    
    void SaveSolutionSetToStream( std::ostream& stream, TSolutionSetList selectedResults, Int32 type) const;

    bool GetBestRedshift( Float64& Redshift, Int32& MatchingNumber ) const;
    bool GetBestMatchNumRedshift( Float64& Redshift, Int32& MatchingNumber ) const;

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

    CLineCatalog m_RestCatalog;
    CLineCatalog m_DetectedCatalog;
  private:
    bool m_bypassDebug;
  };
}

#endif // _REDSHIFT_OPERATOR_RAYMATCHINGRESULT_
