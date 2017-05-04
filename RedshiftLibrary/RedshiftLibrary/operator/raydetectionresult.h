#ifndef _REDSHIFT_OPERATOR_RAYDETECTIONRESULT_
#define _REDSHIFT_OPERATOR_RAYDETECTIONRESULT_

#include <epic/core/common/datatypes.h>
#include <epic/redshift/processflow/result.h>
#include <epic/redshift/ray/catalog.h>

#include <vector>

namespace NSEpic
{

/**
 * \ingroup Redshift
 * Responsible for outputing the results in a standardized format.
 */
class CLineDetectionResult : public COperatorResult
{

public:

    CLineDetectionResult();
    virtual ~CLineDetectionResult();

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;

    CRayCatalog RayCatalog;
    std::vector<std::string> PeakListDetectionStatus;
};

}

#endif // _REDSHIFT_OPERATOR_RAYDETECTIONRESULT_