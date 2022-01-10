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
#include <RedshiftLibrary/method/linemeassolve.h>
#include <RedshiftLibrary/log/log.h>
#include "RedshiftLibrary/processflow/parameterstore.h"

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
    cms.fillRayIds();
    /*
    const CRayCatalog& restraycatalog=*(inputContext->GetRayCatalog("galaxy"));
    CRayCatalog::TRayVector restRayList = restraycatalog.GetFilteredList(-1,-1); // TODO should be retrievable directly from inputContext, with approprate filters

    std::shared_ptr<CModelFittingResult> res = std::make_shared<CModelFittingResult>(cms,
                                                                                    cms.Redshift,
                                                                                    1,
                                                                                    restRayList
                                                                                    );
    */
    std::shared_ptr<CLineModelSolution> res = std::make_shared<CLineModelSolution>(cms);
    resultStore->StoreScopedGlobalResult("linemeas",res);
    resultStore->StoreScopedGlobalResult("linemeas_parameters",res);
    resultStore->StoreScopedGlobalResult("linemeas_model",m_linemodel.getFittedModel());
    return std::make_shared<CLineMeasSolveResult>(CLineMeasSolveResult());
  }


}
