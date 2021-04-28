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
