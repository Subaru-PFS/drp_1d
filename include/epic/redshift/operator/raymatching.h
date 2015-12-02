#ifndef _REDSHIFT_RAY_MATCHING_
#define _REDSHIFT_RAY_MATCHING_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>
#include <epic/redshift/operator/raymatchingresult.h>
#include <epic/redshift/ray/ray.h>
#include <epic/redshift/ray/catalog.h>

namespace NSEpic
{
  class CRayMatching
  {
  public:
    CRayMatching();
    virtual ~CRayMatching();

    std::shared_ptr<CRayMatchingResult> Compute( const CRayCatalog& restRayCatalog, const CRayCatalog& detectedRayCatalog, const TFloat64Range& redshiftRange, Int32 nThreshold = 5, Float64 tol = 0.002, Int32 typeFilter = CRay::nType_Emission, Int32 detectedForceFilter = -1, Int32 restRorceFilter = -1 );
    
  private:
    Bool AreSolutionSetsEqual( const CRayMatchingResult::TSolutionSet& s1, const CRayMatchingResult::TSolutionSet& s2);
  };
}

#endif //_REDSHIFT_RAY_MATCHING_
