#include <RedshiftLibrary/operator/raydetectionresult.h>
#include <RedshiftLibrary/ray/ray.h>

using namespace NSEpic;

/**
 * Empty constructor.
 */
CLineDetectionResult::CLineDetectionResult()
{

}

/**
 * Empty destructor.
 */
CLineDetectionResult::~CLineDetectionResult()
{

}

/**
 * Call "Save" on each entry of "RayCatalog". Then output comments with each entry of "PeakListDetectionStatus".
 */
void CLineDetectionResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    CRayCatalog::TRayVector::const_iterator it;

    Bool SaveDescriptionDone=false;
    for( it = RayCatalog.GetList().begin(); it != RayCatalog.GetList().end(); ++it )
    {
        if(!SaveDescriptionDone){
            it->SaveDescription( stream );
            stream << std::endl;
            SaveDescriptionDone = true;
        }
        it->Save( stream );
        stream << std::endl;
    }

    stream << std::endl;
    stream << "#Peak detection status:" <<std::endl;
    for( Int32 i=0; i<PeakListDetectionStatus.size(); i++ )
    {
        stream << "#" << i << "\t" << PeakListDetectionStatus[i];
        stream << std::endl;
    }
}

/**
 * Stub method.
 */
void CLineDetectionResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    stream << "not implemented" << std::endl;
}
