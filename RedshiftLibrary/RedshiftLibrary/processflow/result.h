#ifndef _REDSHIFT_PROCESSFLOW_OPERATORRESULT_
#define _REDSHIFT_PROCESSFLOW_OPERATORRESULT_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/continuum/indexes.h>
#include <vector>
#include <ostream>
#include <map>

namespace NSEpic
{

class CDataStore;

/**
 * \ingroup Redshift
 */
class COperatorResult
{

public:

    COperatorResult() = default;
    virtual ~COperatorResult() = default;

    virtual void Save(std::ostream& stream) const = 0;
    virtual void SaveLine(std::ostream& stream ) const = 0;
    virtual void SaveJSON(std::ostream& stream) const;
    //virtual void Load( std::istream& stream ) = 0;

    void SaveFloat64(std::ostream& stream,Float64 data) const;
    void SaveTFloat64List(std::ostream& stream,std::string name,TFloat64List data) const;
    void SaveTFloat64ListOfList(std::ostream& stream,std::string name,std::vector<TFloat64List> data) const;
    void SaveInt32Vector(std::ostream& stream,std::string name,std::vector<Int32> data) const;
    void SaveStringVector(std::ostream& stream,std::string name,std::vector<std::string> data) const;
    void SaveStringVectorOfVector(std::ostream& stream,std::string name,std::vector<std::vector<std::string>> data) const;
    void SaveTContinuumIndexListVector(std::ostream& stream,std::string name,std::vector<CContinuumIndexes::TContinuumIndexList> data) const;

    //should pure virtual, let's implement them first here before changing all operator results ....
  virtual void getCandidateData(const int& rank,const std::string& name, Float64& v) const;
  virtual void getCandidateData(const int& rank,const std::string& name, Int32& v) const;
  virtual void getCandidateData(const int& rank,const std::string& name, std::string& v) const;
  virtual void getCandidateData(const int& rank,const std::string& name, double **data, int *size) const;
  virtual void getCandidateData(const int& rank,const std::string& name, std::string *data, int *size) const;
  virtual void getCandidateData(const int& rank,const std::string& name, int  **data, int *size) const;

  virtual void getData(const std::string& name, Int32& v) const;
  virtual void getData(const std::string& name, Float64& v) const;
  virtual void getData(const std::string& name, std::string& v) const;
  virtual void getData(const std::string& name, double **data, int *size) const;
  virtual void getData(const std::string& name, std::string *data, int *size) const;
  virtual void getData(const std::string& name, int **data, int *size) const;

protected:


};

typedef std::vector< std::shared_ptr<COperatorResult> >           TOperatorResultList;
typedef std::map< std::string, std::shared_ptr< const COperatorResult> > TOperatorResultMap;

}

#endif
