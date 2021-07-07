#ifndef _REDSHIFT_OPERATOR_RAYMATCHING_
#define _REDSHIFT_OPERATOR_RAYMATCHING_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/operator/raymatchingresult.h"
#include "RedshiftLibrary/ray/ray.h"
#include "RedshiftLibrary/ray/catalog.h"

namespace NSEpic
{

  /**
   * \ingroup Redshift
   * Holds the algorithms for calculating and comparing matches, when given reference and detection catalogue of rays (spectral lines) or matches themselves (instances of CRayMatchingResult).
   */
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

#endif // _REDSHIFT_OPERATOR_RAYMATCHING_
