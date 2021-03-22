#include <RedshiftLibrary/method/solve.h>
#include <RedshiftLibrary/processflow/inputcontext.h>
#include <RedshiftLibrary/processflow/parameterstore.h>

using namespace NSEpic;

CSolve::CSolve(std::string name,TScopeStack &scope,std::string objectType):
  m_categoryList({objectType}),
  m_objectTypeScope(scope,objectType),
  m_objectType(objectType),
  m_name(name)
{

}

void CSolve::InitRanges(const CInputContext& inputContext)
{
  if (m_objectType != "classification")// TODO this is temporary hack, we can put a flag, or overload the method or intermediary CSolve class
    {
      m_lambdaRange=inputContext.m_lambdaRange;
      TFloat64Range redshiftRange=inputContext.m_ParameterStore->GetScoped<TFloat64Range>("redshiftrange");
      m_redshiftSampling=inputContext.m_ParameterStore->GetScoped<std::string>("redshiftsampling");
      Float64 redshiftStep = inputContext.m_ParameterStore->GetScoped<Float64>( "redshiftstep" );

      if(m_redshiftSampling=="log")
        {
          m_redshifts = redshiftRange.SpreadOverLog( redshiftStep ); //experimental: spreadover a grid at delta/(1+z), unusable because PDF needs regular z-step
        }
      else
        {
          m_redshifts = redshiftRange.SpreadOver( redshiftStep );
        }
    }
}


void CSolve::Compute(const CInputContext &inputContext,
                             COperatorResultStore &resultStore,
                             TScopeStack &scope)
{
  //      beginCompute(inputContext,resultStore,scope);
  InitRanges(inputContext);
  std::shared_ptr<CSolveResult> result=std::shared_ptr<CSolveResult>(nullptr);  
  {
    CAutoScope autoscope(scope,m_name);
    result = compute(inputContext,resultStore,scope);
  }
  saveToResultStore(result,resultStore);
  //endCompute(inputContext,resultStore,scope);
}

void CSolve::saveToResultStore(std::shared_ptr<CSolveResult> result,COperatorResultStore &resultStore)
{
  resultStore.StoreScopedGlobalResult("result",result); 
}
