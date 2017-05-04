#ifndef _REDSHIFT_RELIABILITY_ZQUALRESULT_
#define _REDSHIFT_RELIABILITY_ZQUALRESULT_

#include <epic/redshift/processflow/result.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/operator/operator.h>

#include <epic/redshift/processflow/resultstore.h>
#include <epic/redshift/processflow/datastore.h>

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


    Bool GetPredictedLabel( const CDataStore& store, std::string& predLabel ) const;

};

}

#endif
