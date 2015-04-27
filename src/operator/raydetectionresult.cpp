#include <epic/redshift/operator/raydetectionresult.h>

#include <epic/redshift/ray/ray.h>

using namespace NSEpic;

IMPLEMENT_MANAGED_OBJECT( CRayDetectionResult )

CRayDetectionResult::CRayDetectionResult()
{

}

CRayDetectionResult::~CRayDetectionResult()
{

}

Void CRayDetectionResult::Save( std::ostream& stream ) const
{
    CRayCatalog::TRayVector::const_iterator it;

    for( it = RayCatalog.GetList().begin(); it != RayCatalog.GetList().end(); ++it )
    {
        it->Save( stream );
        stream << std::endl;
    }

}
