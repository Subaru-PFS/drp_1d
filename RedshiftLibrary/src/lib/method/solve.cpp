#include <RedshiftLibrary/method/solve.h>
#include <RedshiftLibrary/processflow/inputContext->h>
#include <RedshiftLibrary/processflow/parameterstore.h>

using namespace NSEpic;

CSolve::CSolve(std::string name,TScopeStack &scope,std::string objectType):
  m_categoryList({objectType}),
  m_objectTypeScope(scope,objectType),
  m_objectType(objectType),
  m_name(name)
{

}
void CSolve::GetRedshiftSampling(std::shared_ptr<const CInputContext>, TFloat64Range& redshiftRange, Float64& redshiftStep)
{
    //default is to read from the scoped paramStore
    redshiftRange=inputContext->m_ParameterStore->GetScoped<TFloat64Range>("redshiftrange");
    redshiftStep = inputContext->m_ParameterStore->GetScoped<Float64>( "redshiftstep" );
    return;
}
void CSolve::InitRanges(std::shared_ptr<const CInputContext> inputContext)
{
  if (m_objectType == "star" || m_objectType=="qso" || m_objectType=="galaxy")// TODO this is temporary hack, we can put a flag, or overload the method or intermediary CSolve class
    {
      m_lambdaRange=inputContext->m_lambdaRange;//non-clamped
      m_redshiftSampling=inputContext->m_ParameterStore->GetScoped<std::string>("redshiftsampling");

      TFloat64Range redshiftRange;
      Float64 redshiftStep;
      GetRedshiftSampling(inputContext, redshiftRange, redshiftStep);
      
      if(m_redshiftSampling=="log")
          m_redshifts = redshiftRange.SpreadOverLogZplusOne( redshiftStep ); //experimental: spreadover a grid at delta/(1+z), unusable because PDF needs regular z-step
      else  
          m_redshifts = redshiftRange.SpreadOver( redshiftStep );
    }
}


void CSolve::Compute(std::shared_ptr<const CInputContext> inputContext,
                     std::shared_ptr<COperatorResultStore> resultStore,
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

void CSolve::saveToResultStore(std::shared_ptr<CSolveResult> result,std::shared_ptr<COperatorResultStore> resultStore)
{
  resultStore->StoreScopedGlobalResult("result",result); 
}
