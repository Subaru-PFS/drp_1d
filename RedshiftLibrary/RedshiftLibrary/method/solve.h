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
    void InitRanges(std::shared_ptr<const CInputContext> inputContext);
    CAutoScope m_objectTypeScope;
    const std::string m_name;
    const std::string m_objectType;
    std::string m_redshiftSampling;
  public:

    CSolve(std::string name,TScopeStack &scope,std::string objectType);
    virtual ~CSolve()=default;
    CSolve(CSolve const& other) = default;
    CSolve& operator=(CSolve const& other) = default;
    CSolve(CSolve&& other) = default;
    CSolve& operator=(CSolve&& other) = default;


    void Compute(std::shared_ptr<const CInputContext> inputContext,
                 std::shared_ptr<COperatorResultStore> resultStore,
                 TScopeStack &scope);

    virtual std::shared_ptr<CSolveResult> compute(std::shared_ptr<const CInputContext> inputContext,
                                                  std::shared_ptr<COperatorResultStore> resultStore,
                                                  TScopeStack &scope)=0;
    
    virtual void saveToResultStore(std::shared_ptr<CSolveResult>,std::shared_ptr<COperatorResultStore> resultStore);
    
    // this method should implement at least populateParameters

    virtual void checkOrInit(){}//=0; // here we retrieve parameters for parameterStore to put them directly in local variables or into operators, rayCatalog and/or tplCatalog can also be checked
   

    
  };

  
}

#endif
