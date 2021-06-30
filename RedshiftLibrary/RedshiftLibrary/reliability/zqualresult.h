#ifndef _REDSHIFT_RELIABILITY_ZQUALRESULT_
#define _REDSHIFT_RELIABILITY_ZQUALRESULT_

#include "RedshiftLibrary/processflow/result.h"
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/operator/operator.h"

#include "RedshiftLibrary/processflow/resultstore.h"
#include "RedshiftLibrary/processflow/datastore.h"

namespace NSEpic
{
class CProcessFlowContext;

class CQualzResult : public COperatorResult
{

public:

	CQualzResult();
	virtual ~CQualzResult();

    void Save(std::ostream& stream ) const;
    void SaveLine(std::ostream& stream ) const;
   

    Bool GetPredictedLabel( const CDataStore& store, std::string& predLabel ) const;

};

}

#endif
