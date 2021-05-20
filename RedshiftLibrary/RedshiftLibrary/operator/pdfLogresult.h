#ifndef _REDSHIFT_OPERATOR_PDFLOGRESULT_
#define _REDSHIFT_OPERATOR_PDFLOGRESULT_

#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/operator/operator.h>

namespace NSEpic
{

class CPdfLogResult : public COperatorResult
{

public:

    CPdfLogResult();
    virtual ~CPdfLogResult();

    void SetSize( UInt32 n );

    TFloat64List           		Redshifts;
    TFloat64List           		valProbaLog;
    TFloat64List            		Overlap;
    COperator::TStatusList   Status;

};


}

#endif
