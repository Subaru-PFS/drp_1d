#ifndef _REDSHIFT_OPERATOR_RAYDETECTION_
#define _REDSHIFT_OPERATOR_RAYDETECTION_

#include <epic/core/common/managedobject.h>
#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>

namespace NSEpic
{

class CSpectrum;
class CRayDetectionResult;

class CRayDetection : public CManagedObject
{
    DEFINE_MANAGED_OBJECT( CRayDetection )

public:

    CRayDetection();
    virtual ~CRayDetection();

    const CRayDetectionResult* Compute( const CSpectrum& spectrum, const TLambdaRange& lambdaRange, const TInt32RangeList& resPeaks, const TInt32RangeList& resPeaksEnlarged );

private:

    Float64 XMadFind( const Float64* x, Int32 n, Float64 median );

};

}

#endif
