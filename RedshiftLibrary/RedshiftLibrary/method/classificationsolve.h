#ifndef _METHOD_CLASSIFICATION_
#define _METHOD_CLASSIFICATION_

#include <RedshiftLibrary/method/solveresult.h>
#include <RedshiftLibrary/method/solve.h>
#include <RedshiftLibrary/method/classificationresult.h>


namespace NSEpic
{

  class CClassificationSolve:public CSolve
  {

  public:

    CClassificationSolve(TScopeStack &scope,std::string objectType);

    std::shared_ptr<CSolveResult> compute(const CInputContext &inputContext,
                 COperatorResultStore &resultStore,
                 TScopeStack &scope);
    
    std::string typeLabel = "U";//"G"/"S"/"Q"


  };
}

#endif
