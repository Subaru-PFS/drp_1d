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


  std::string m_ReliabilityLabel="C6";

};

}

#endif
