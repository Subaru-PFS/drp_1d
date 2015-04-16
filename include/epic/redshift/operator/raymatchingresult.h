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

    Void Save( std::ostream& stream ) const;

    TSolutionSetList SolutionSetList;
};


}

#endif
