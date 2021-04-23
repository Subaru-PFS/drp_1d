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
    CPdfMargZLogResult() = default;
    ~CPdfMargZLogResult() = default;
    CPdfMargZLogResult(const TFloat64List & redshifts);
    
   
    Int32 Load( std::string filePath );

    Int32 getIndex( Float64 z ) const;

    TFloat64List          Redshifts;
    TFloat64List          valProbaLog;
    Float64               valEvidenceLog;
    UInt32 				  countTPL;

  void getCandidateData(const int& rank,const std::string& name, Float64& v) const;
  void getCandidateData(const int& rank,const std::string& name, Int32& v) const;
  void getCandidateData(const int& rank,const std::string& name, std::string& v) const;
  void getCandidateData(const int& rank,const std::string& name, double **data, int *size) const;

  void getData(const std::string& name, Int32& v) const;
  void getData(const std::string& name, Float64& v) const;
  void getData(const std::string& name, std::string& v) const;
  void getData(const std::string& name, double **data, int *size) const;

};


}

#endif
