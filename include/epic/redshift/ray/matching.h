#ifndef _REDSHIFT_RAY_MATCHING_
#define _REDSHIFT_RAY_MATCHING_

#include <epic/core/common/managedobject.h>
#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>

#include <epic/redshift/operator/raymatchingresult.h>
#include <epic/redshift/ray/ray.h>
#include <epic/redshift/ray/catalog.h>


namespace NSEpic
{

class CRayMatching : public CManagedObject
{
    DEFINE_MANAGED_OBJECT( CRayMatching )

public:

    CRayMatching();
    virtual ~CRayMatching();

    CRayMatchingResult* Compute(const CRayCatalog& restRayCatalog, const CRayCatalog& detectedRayCatalog, const TFloat64Range& redshiftRange, Int32 nThreshold = 5, Float64 tol = 0.002 );
    /*
    Float64 GetMeanRedshiftSolutionByIndex(Int32 index);
    Float64 GetMeanRedshiftSolution(TRedshiftSolutionSet& s);
    Int32 GetMaxMatchingNumber();
    const TRedshiftSolutionSetList GetSolutionsListOverNumber(Int32 number) const;
    const TRedshiftSolutionSetList GetResults() const;
    Bool GetBestRedshift(Float64& Redshift, Int32& MatchingNumber);
    */

private:

    Bool AreSolutionSetsEqual( const CRayMatchingResult::TSolutionSet& s1, const CRayMatchingResult::TSolutionSet& s2);


    //Int32       m_N;            // number of matching lines needed

    //Float64     m_Tol;          // tolerance (Angstrom) for the line matching
    //TFloat64Range     m_RedshiftRange;         // redshift range

    //Float64     m_ZStep;        // not used ? this parameter was specified in the EZ documentation
    //Float64     m_Delta;        // not used ? this parameter was specified in the EZ documentation
    //TRedshiftSolutionSetList            m_Results;

};

}

#endif
