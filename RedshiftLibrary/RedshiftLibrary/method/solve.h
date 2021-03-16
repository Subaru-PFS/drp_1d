#ifndef _SOLVE_H_
#define _SOLVE_H_

#include <RedshiftLibrary/method/solveresult.h>
#include <RedshiftLibrary/processflow/inputcontext.h>
#include <RedshiftLibrary/processflow/resultstore.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/processflow/autoscope.h>

namespace NSEpic
{


  class CSolve
  {
  protected:
    const TStringList  m_categoryList; 
    //    std::string name;
    //std::vector<ParameterDescriptor> parameters; pd={name,type,constraint(union type)}
    TFloat64Range m_lambdaRange;
    TFloat64List m_redshifts;
    void InitRanges(const CInputContext& inputContext);
    CAutoScope m_objectTypeScope;
  public:

    CSolve(TScopeStack &scope,std::string objectType);
    ~CSolve();

    //Proposition 1 bis
    virtual std::shared_ptr<CSolveResult> Compute(const CInputContext &inputContext,
                                                 COperatorResultStore &resultStore,
                                                  TScopeStack &scope)=0;
    /*    {
      beginCompute(inputContext,resultStore,scope);
      compute(inputContext,resultStore,scope);
      endCompute(inputContext,resultStore,scope);
    };
    */
    virtual void saveToResultStore(){}//=0;
    
    // this method should implement at least populateParameters

    virtual void checkOrInit(){}//=0; // here we retrieve parameters for parameterStore to put them directly in local variables or into operators, rayCatalog and/or tplCatalog can also be checked
   

    
  };

  
}

#endif
