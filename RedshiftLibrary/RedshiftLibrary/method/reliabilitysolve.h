#ifndef _METHOD_RELIABILITY_
#define _METHOD_RELIABILITY_

#include "RedshiftLibrary/method/solveresult.h"
#include "RedshiftLibrary/method/solve.h"
#include "RedshiftLibrary/method/reliabilityresult.h"


namespace NSEpic
{

  class CReliabilitySolve:public CSolve
  {

  public:

    CReliabilitySolve(TScopeStack &scope,std::string objectType);

    std::shared_ptr<CSolveResult> compute(std::shared_ptr<const CInputContext> inputContext,
                                          std::shared_ptr<COperatorResultStore> resultStore,
                                          TScopeStack &scope);
    
    
  };
}

#endif
