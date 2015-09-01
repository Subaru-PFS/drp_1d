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

Void CRayDetectionResult::Save( const COperatorResultStore& store, std::ostream& stream ) const
{
    CRayCatalog::TRayVector::const_iterator it;

    for( it = RayCatalog.GetList().begin(); it != RayCatalog.GetList().end(); ++it )
    {
        it->Save( stream );
        stream << std::endl;
    }

    stream << std::endl;
    stream << "#Peak detection status:" <<std::endl;
    for(Int32 i=0; i<PeakListDetectionStatus.size(); i++)
    {
        stream << "#" << i << "\t" << PeakListDetectionStatus[i];
        stream << std::endl;
    }



}



Void CRayDetectionResult::SaveLine( const COperatorResultStore& store, std::ostream& stream ) const
{
    stream << "not implemented" << std::endl;
}
