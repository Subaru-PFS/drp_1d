#include <RedshiftLibrary/method/solve.h>
#include <RedshiftLibrary/processflow/inputcontext.h>
#include <RedshiftLibrary/processflow/parameterstore.h>

using namespace NSEpic;

CSolve::CSolve(std::string objectType):
m_categoryList({objectType})
{
  //  m_categoryList.push_back(objectType);
}

CSolve::~CSolve()
{
}


void CSolve::InitRanges(const CInputContext& inputContext)
{
  m_lambdaRange=inputContext.m_lambdaRange;
  TFloat64Range redshiftRange=inputContext.m_ParameterStore->Get<TFloat64Range>("redshiftrange");
  std::string redshiftSampling=inputContext.m_ParameterStore->Get<std::string>("redshiftsampling");
  Float64 redshiftStep = inputContext.m_ParameterStore->Get<Float64>( "redshiftstep" );

  if(redshiftSampling=="log")
    {
      m_redshifts = redshiftRange.SpreadOverLog( redshiftStep ); //experimental: spreadover a grid at delta/(1+z), unusable because PDF needs regular z-step
    }
  else
    {
      m_redshifts = redshiftRange.SpreadOver( redshiftStep );
    }

  
}
