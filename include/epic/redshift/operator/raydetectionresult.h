#ifndef _REDSHIFT_OPERATOR_RAYDETECTIONRESULT_
#define _REDSHIFT_OPERATOR_RAYDETECTIONRESULT_

#include <epic/redshift/operator/result.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/ray/catalog.h>

#include <vector>

namespace NSEpic
{

class CRayDetectionResult : public COperatorResult
{

    DEFINE_MANAGED_OBJECT( CRayDetectionResult )

public:

    CRayDetectionResult();
    virtual ~CRayDetectionResult();

    Void Save( const COperatorResultStore& store, std::ostream& stream ) const;

    CRayCatalog RayCatalog;
};


}

#endif
