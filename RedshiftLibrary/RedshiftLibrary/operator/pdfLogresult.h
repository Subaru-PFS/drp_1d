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

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    Void Load( std::istream& stream );

    Void SetSize( UInt32 n );

    TFloat64List           		Redshifts;
    TFloat64List           		valProbaLog;
    TFloat64List            		Overlap;
    COperator::TStatusList   Status;

};


}

#endif
