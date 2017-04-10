#ifndef _REDSHIFT_OPERATOR_PDFMARGZLOGRESULT_
#define _REDSHIFT_OPERATOR_PDFMARGZLOGRESULT_

#include <epic/redshift/processflow/result.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/operator/operator.h>


using namespace std;
namespace NSEpic
{

class CPdfMargZLogResult : public COperatorResult
{

  public:

    CPdfMargZLogResult();
    virtual ~CPdfMargZLogResult();

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;

    TFloat64List          Redshifts;
    TFloat64List          valProbaLog;
    UInt32 					countTPL;

};


}

#endif
