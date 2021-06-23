#ifndef _REDSHIFT_PROCESSFLOW_OPERATORRESULT_
#define _REDSHIFT_PROCESSFLOW_OPERATORRESULT_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/continuum/indexes.h>
#include <vector>
#include <ostream>
#include <map>

namespace NSEpic
{

/**
 * \ingroup Redshift
 */
class COperatorResult
{

public:

    COperatorResult() = default;
    virtual ~COperatorResult() = default;

  const std::string& getType() const {return m_type;}
  virtual const std::string& getCandidateDatasetType(const std::string& dataset) const
  {
    throw GlobalException(UNKNOWN_ATTRIBUTE,"This operator result does not support this operation");
  }
  virtual std::shared_ptr<const COperatorResult> getCandidate(const int& rank,const std::string& dataset) const
  {
    throw GlobalException(UNKNOWN_ATTRIBUTE,"This operator result does not support this operation");
  }

  virtual bool HasCandidateDataset(const std::string& dataset) const {
    throw GlobalException(UNKNOWN_ATTRIBUTE,"This operator result does not support this operation");
  }

protected:

  std::string m_type="COperatorResult";
};

typedef std::vector< std::shared_ptr<COperatorResult> >           TOperatorResultList;
typedef std::map< std::string, std::shared_ptr< const COperatorResult> > TOperatorResultMap;

}

#endif
