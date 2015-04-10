#ifndef _REDSHIFT_RAY_MATCHING_
#define _REDSHIFT_RAY_MATCHING_

#include <epic/core/common/managedobject.h>
#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>

#include <epic/redshift/ray/ray.h>
#include <epic/redshift/ray/catalog.h>


namespace NSEpic
{

    struct SRedshiftSolution
    {
        SRedshiftSolution( )
        {
            DetectedRay= 0.0;
            RestRay = 0.0;
            Redshift = 0.0;
        }

        SRedshiftSolution( Float64 detectedRay, Float64 restRay, Float64 redshift)
        {
            DetectedRay = detectedRay;
            RestRay = restRay;
            Redshift = redshift;
        }
        Float64 DetectedRay;
        Float64 RestRay;
        Float64 Redshift;

        bool operator < (const SRedshiftSolution& str) const
        {
            if(DetectedRay == str.DetectedRay){
                return (RestRay < str.RestRay);
            }else{
                return (DetectedRay < str.DetectedRay);
            }
        }
    };

    typedef std::vector<SRedshiftSolution>   TRedshiftSolutionSet; // a set of (detected line,rest line) couples for a given redshift
    typedef std::vector<TRedshiftSolutionSet>   TRedshiftSolutionSetList; // a list of possible redshift solutions


class CRayMatching : public CManagedObject
{
    DEFINE_MANAGED_OBJECT( CRayMatching )

public:

    CRayMatching();
    virtual ~CRayMatching();

    Bool Compute(const CRayCatalog& restRayCatalog, const CRayCatalog& detectedRayCatalog, const TFloat64Range& redshiftRange, Int32 nThreshold = 5, Float64 tol = 0.002 );
    Float64 GetMeanRedshiftSolutionByIndex(Int32 index);
    Float64 GetMeanRedshiftSolution(TRedshiftSolutionSet& s);
    Int32 GetMaxMatchingNumber();
    const TRedshiftSolutionSetList GetSolutionsListOverNumber(Int32 number) const;
    const TRedshiftSolutionSetList GetResults() const;


private:

    Bool AreSolutionSetsEqual(TRedshiftSolutionSet s1, TRedshiftSolutionSet s2);


    //Int32       m_N;            // number of matching lines needed

    //Float64     m_Tol;          // tolerance (Angstrom) for the line matching
    //TFloat64Range     m_RedshiftRange;         // redshift range

    //Float64     m_ZStep;        // not used ? this parameter was specified in the EZ documentation
    //Float64     m_Delta;        // not used ? this parameter was specified in the EZ documentation
    TRedshiftSolutionSetList            m_Results;

};

}

#endif
