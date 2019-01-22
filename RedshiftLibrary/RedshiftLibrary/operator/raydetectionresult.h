#ifndef _REDSHIFT_OPERATOR_RAYDETECTIONRESULT_
#define _REDSHIFT_OPERATOR_RAYDETECTIONRESULT_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/ray/catalog.h>

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

    void Save( const CDataStore& store, std::ostream& stream ) const;
    void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    inline Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const
    {
        return 1;
    }

    CRayCatalog RayCatalog;
    std::vector<std::string> PeakListDetectionStatus;
};

}

#endif // _REDSHIFT_OPERATOR_RAYDETECTIONRESULT_
