#ifndef _REDSHIFT_PROCESSFLOW_RELIABILITYRESULT_
#define _REDSHIFT_PROCESSFLOW_RELIABILITYRESULT_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/method/solveresult.h>

#include <vector>
#include <ostream>
#include <map>

namespace NSEpic
{

class CDataStore;

/**
 * \ingroup Redshift
 */
class CReliabilityResult : public CSolveResult
{

public:

  CReliabilityResult();

  void getData(const std::string& name, std::string& v) const;
  virtual void Save(std::ostream& stream) const {};
  virtual void SaveLine(std::ostream& stream ) const {};

  std::string m_ReliabilityLabel="C6";

};

}

#endif
