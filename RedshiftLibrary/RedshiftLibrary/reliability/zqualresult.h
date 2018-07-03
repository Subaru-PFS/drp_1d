#ifndef _REDSHIFT_RELIABILITY_ZQUALRESULT_
#define _REDSHIFT_RELIABILITY_ZQUALRESULT_

#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/operator/operator.h>

#include <RedshiftLibrary/processflow/resultstore.h>
#include <RedshiftLibrary/processflow/datastore.h>

namespace NSEpic
{
class CProcessFlowContext;

class CQualzResult : public COperatorResult
{

public:

	CQualzResult();
	virtual ~CQualzResult();

    Void Save( const CDataStore& store, std::ostream& stream ) const;
	Void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    inline Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const
    {
        return 1;
    }

    Bool GetPredictedLabel( const CDataStore& store, std::string& predLabel ) const;

};

}

#endif
