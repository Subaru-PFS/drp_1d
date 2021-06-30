#ifndef _METHOD_CLASSIFICATION_
#define _METHOD_CLASSIFICATION_

#include "RedshiftLibrary/method/solveresult.h"
#include "RedshiftLibrary/method/solve.h"
#include "RedshiftLibrary/method/classificationresult.h"


namespace NSEpic
{

  class CClassificationSolve:public CSolve
  {

  public:

    CClassificationSolve(TScopeStack &scope,std::string objectType);

  private:
  
    std::shared_ptr<CSolveResult> compute(std::shared_ptr<const CInputContext> inputContext,
                                          std::shared_ptr<COperatorResultStore> resultStore,
                 TScopeStack &scope) override;
    
    std::string typeLabel = "U";//"G"/"S"/"Q"


  };
}

#endif
