#ifndef _REDSHIFT_OPERATOR_RAYMATCHINGRESULT_
#define _REDSHIFT_OPERATOR_RAYMATCHINGRESULT_

#include <epic/redshift/operator/result.h>
#include <epic/core/common/datatypes.h>

#include <vector>

namespace NSEpic
{

class CRayMatchingResult : public COperatorResult
{

    DEFINE_MANAGED_OBJECT( CRayMatchingResult )

public:

    struct SSolution
    {
       Float64 DetectedRay;
       Float64 RestRay;
       Float64 Redshift;

       SSolution( Float64 detectedRay, Float64 restRay, Float64 redshift )
       {
           DetectedRay = detectedRay;
           RestRay = restRay;
           Redshift = redshift;
       }

       bool operator < (const SSolution& str) const
       {
           if(DetectedRay == str.DetectedRay){
               return (RestRay < str.RestRay);
           }else{
               return (DetectedRay < str.DetectedRay);
           }
       }
    };

    typedef std::vector<SSolution>      TSolutionSet; // a set of (detected line,rest line) couples for a given redshift
    typedef std::vector<TSolutionSet>   TSolutionSetList; // a list of possible redshift solutions

    CRayMatchingResult();
    virtual ~CRayMatchingResult();

    Void                        Save( std::ostream& stream ) const;

    Bool                        GetBestRedshift(Float64& Redshift, Int32& MatchingNumber) const;
    Int32                       GetMaxMatchingNumber() const;
    Float64                     GetMeanRedshiftSolution( const TSolutionSet& s) const;
    Float64                     GetMeanRedshiftSolutionByIndex(Int32 index) const;
    TSolutionSetList            GetSolutionsListOverNumber(Int32 number) const;
    TFloat64List                GetAverageRedshiftListOverNumber(Int32 number) const;

    TFloat64List                GetRoundedRedshiftCandidatesOverNumber(Int32 number, Float64 step) const;
    TFloat64List                GetRedshiftRangeCandidatesOverNumber(Int32 number) const;

    TSolutionSetList    SolutionSetList;

};


}

#endif
