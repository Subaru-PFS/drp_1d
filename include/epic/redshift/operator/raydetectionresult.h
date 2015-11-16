#ifndef _REDSHIFT_OPERATOR_RAYDETECTIONRESULT_
#define _REDSHIFT_OPERATOR_RAYDETECTIONRESULT_

#include <epic/redshift/processflow/result.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/ray/catalog.h>

#include <vector>

namespace NSEpic
{

/**
 * \ingroup Redshift
 */
class CLineDetectionResult : public COperatorResult
{

    DEFINE_MANAGED_OBJECT( CLineDetectionResult )

public:

    CLineDetectionResult();
    virtual ~CLineDetectionResult();

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;

    CRayCatalog RayCatalog;
    std::vector<std::string> PeakListDetectionStatus;
};


}

#endif
