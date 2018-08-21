#ifndef _REDSHIFT_OPERATOR_PDFMARGZLOGRESULT_
#define _REDSHIFT_OPERATOR_PDFMARGZLOGRESULT_

#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/operator/operator.h>


using namespace std;
namespace NSEpic
{

class CPdfMargZLogResult : public COperatorResult
{

  public:

    CPdfMargZLogResult();
    virtual ~CPdfMargZLogResult();

    void Save( const CDataStore& store, std::ostream& stream ) const;
    void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    inline Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const
    {
        return 1;
    }
    Int32 Load( std::string filePath );

    Int32 getIndex( Float64 z ) const;

    TFloat64List          Redshifts;
    TFloat64List          valProbaLog;
    Float64                 valEvidenceLog;
    UInt32 					countTPL;

};


}

#endif
