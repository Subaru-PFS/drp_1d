// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#include "RedshiftLibrary/method/solve.h"
#include "RedshiftLibrary/processflow/inputcontext.h"
#include "RedshiftLibrary/processflow/parameterstore.h"

using namespace NSEpic;

CSolve::CSolve(std::string name,TScopeStack &scope,std::string objectType):
  m_categoryList({objectType}),
  m_objectTypeScope(scope,objectType),
  m_objectType(objectType),
  m_name(name)
{

}
void CSolve::GetRedshiftSampling(std::shared_ptr<const CInputContext> inputContext, TFloat64Range& redshiftRange, Float64& redshiftStep) 
{   
    auto searchLogRebin = inputContext->m_logRebin.find(m_objectType);
    if(searchLogRebin!=inputContext->m_logRebin.end())
    {
        redshiftRange = searchLogRebin->second.zrange;
        redshiftStep = searchLogRebin->second.logGridStep;
        if(m_redshiftSampling=="lin"){
                m_redshiftSampling = "log";
                Log.LogWarning("m_redshift sampling value is forced to log since FFTprocessing is used");
        }
    }else{
      //default is to read from the scoped paramStore
      redshiftRange = inputContext->GetParameterStore()->GetScoped<TFloat64Range>("redshiftrange");
      redshiftStep = inputContext->GetParameterStore()->GetScoped<Float64>( "redshiftstep" );
    }
    return;
}
void CSolve::InitRanges(std::shared_ptr<const CInputContext> inputContext)
{
  if (m_objectType == "star" || m_objectType=="qso" ||
      m_objectType=="galaxy" || m_objectType=="linemeas")// TODO this is temporary hack, we can put a flag, or overload the method or intermediary CSolve class
    {
      m_lambdaRange=inputContext->m_lambdaRange;//non-clamped

      //m_redshiftSampling could be overwritten if fftprocessing is activated
      m_redshiftSampling=inputContext->GetParameterStore()->GetScoped<std::string>("redshiftsampling");

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
    saveToResultStore(result,resultStore);
  }

  //endCompute(inputContext,resultStore,scope);
}

void CSolve::saveToResultStore(std::shared_ptr<CSolveResult> result,std::shared_ptr<COperatorResultStore> resultStore) const
{
  resultStore->StoreScopedGlobalResult("solveResult",result); 
}
