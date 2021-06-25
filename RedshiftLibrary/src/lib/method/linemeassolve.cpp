#include <RedshiftLibrary/method/linemeassolve.h>
#include <RedshiftLibrary/log/log.h>
#include "RedshiftLibrary/processflow/parameterstore.h"
#include "RedshiftLibrary/linemodel/modelfittingresult.h"

namespace NSEpic
{

  CLineMeasSolve::CLineMeasSolve(TScopeStack &scope,string objectType,string calibrationPath):
  CSolve("linemeassolve",scope,objectType),
    m_calibrationPath(calibrationPath)
{
}
  
  CLineMeasSolve::~CLineMeasSolve()
  {

  }

  void CLineMeasSolve::GetRedshiftSampling(std::shared_ptr<const CInputContext> inputContext, TFloat64Range& redshiftRange, Float64& redshiftStep) 
{
    //default is to read from the scoped paramStore
  Float64 rangeCenter = inputContext->GetParameterStore()->GetScoped<Float64>( "redshiftref" );
  Float64 halfRange = inputContext->GetParameterStore()->GetScoped<Float64>( "dzhalf" );

  redshiftRange = TFloat64Range(rangeCenter-halfRange,rangeCenter+halfRange);
    redshiftStep = inputContext->GetParameterStore()->GetScoped<Float64>( "redshiftstep" );

    
}

  
  void CLineMeasSolve::Init()
  {

  }



  std::shared_ptr<CSolveResult> CLineMeasSolve::compute(std::shared_ptr<const CInputContext> inputContext,
                                                        std::shared_ptr<COperatorResultStore> resultStore,
                                                        TScopeStack &scope)
  {
    
    CLineModelSolution cms;
    {
      CAutoScope autoscope(scope,"linemodel");

      cms =m_linemodel.computeForLineMeas(inputContext,m_calibrationPath,m_redshifts);
    }
    const CRayCatalog& restraycatalog=*(inputContext->GetRayCatalog().get());
    CRayCatalog::TRayVector restRayList = restraycatalog.GetFilteredList(-1,-1); // TODO should be retrievable directly from inputContext, with approprate filters

    std::shared_ptr<CModelFittingResult> res = std::make_shared<CModelFittingResult>(cms,
                                                                                    cms.Redshift,
                                                                                    1,
                                                                                    restRayList
                                                                                    );
    
    resultStore->StoreScopedGlobalResult("linemeas",res);
    return std::make_shared<CLineMeasSolveResult>(CLineMeasSolveResult());
  }


}
