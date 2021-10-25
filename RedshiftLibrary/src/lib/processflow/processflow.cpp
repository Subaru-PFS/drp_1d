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
#include "RedshiftLibrary/processflow/processflow.h"
#include "RedshiftLibrary/processflow/autoscope.h"

#include "RedshiftLibrary/linemodel/calibrationconfig.h"

#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/debug/assert.h"

#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/method/templatefittingsolve.h"
#include "RedshiftLibrary/method/linematchingsolve.h"
#include "RedshiftLibrary/method/tplcombinationsolve.h"
#include "RedshiftLibrary/method/linemodelsolve.h"
#include "RedshiftLibrary/method/classificationsolve.h"
#include "RedshiftLibrary/method/reliabilitysolve.h"
#include "RedshiftLibrary/method/classificationresult.h"
#include "RedshiftLibrary/statistics/pdfcandidateszresult.h"
#include "RedshiftLibrary/method/linemeassolve.h"

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <cstdio>
#include <cfloat>

using namespace std;
using namespace boost;
using namespace NSEpic;
namespace bfs = boost::filesystem;


void CProcessFlow::Process( CProcessFlowContext& ctx )
{

  try
    {
      Log.LogInfo("=====================================================================");
      std::ostringstream oss;
      oss << "Processing Spectrum: " << ctx.GetSpectrum()->GetName();
      Log.LogInfo(oss.str());
      Log.LogInfo("=====================================================================");

    Float64       maxCount; 
    Float64       redshiftseparation = 2*0.005;
    if(ctx.GetParameterStore()->Has<Float64>( "extremaredshiftseparation"))
      redshiftseparation = ctx.GetParameterStore()->Get<Float64>( "extremaredshiftseparation");//todo: decide on the default values using latest analyses plots

    //retrieve the calibration dir path
    std::string calibrationDirPath = ctx.GetParameterStore()->Get<std::string>( "calibrationDir");

    //************************************
    // Stellar method

    bool enableStarFitting = ctx.GetParameterStore()->Get<bool>( "enablestarsolve");
    Log.LogInfo( "Stellar solve enabled : %d", enableStarFitting);

    if(enableStarFitting){
        Log.LogInfo("Processing stellar fitting");
        CMethodTemplateFittingSolve solve(ctx.m_ScopeStack,"star");
        solve.Compute(ctx.GetInputContext(),
                      ctx.GetResultStore(),
                      ctx.m_ScopeStack);
   
      }
   
    // Quasar method
    //    std::shared_ptr<CSolveResult> qsoResult;
    bool enableQsoFitting = ctx.GetParameterStore()->Get<bool>( "enableqsosolve");
    Log.LogInfo( "QSO solve enabled : %d", enableQsoFitting);

    if(enableQsoFitting){
        Log.LogInfo("Processing QSO fitting");
        std::string qso_method = ctx.GetParameterStore()->Get<std::string>( "qso.method");
        boost::algorithm::to_lower(qso_method);
        if(qso_method=="templatefittingsolve"){
          CMethodTemplateFittingSolve solve(ctx.m_ScopeStack,"qso");
          solve.Compute( ctx.GetInputContext(),
                         ctx.GetResultStore(),
                         ctx.m_ScopeStack);
        }else if(qso_method=="tplcombinationsolve"){
            CMethodTplcombinationSolve solve(ctx.m_ScopeStack,"qso");
            solve.Compute(ctx.GetInputContext(),
                          ctx.GetResultStore(),
                          ctx.m_ScopeStack);
                                      }
        
        else if(qso_method=="linemodelsolve"){
            Log.LogInfo("Linemodel qso fitting...");
            CLineModelSolve Solve(ctx.m_ScopeStack,"qso",calibrationDirPath);
            Solve.Compute( ctx.GetInputContext(),
                           ctx.GetResultStore(),
                           ctx.m_ScopeStack);
        }else{
            throw GlobalException(INTERNAL_ERROR,"Problem found while parsing the qso method parameter");
        }
    }

    // Galaxy method

    if(ctx.GetParameterStore()->Get<bool>( "enablegalaxysolve"))
      {
        std::string methodName = ctx.GetParameterStore()->Get<std::string>( "galaxy.method");
        boost::algorithm::to_lower(methodName);

        std::string galaxy_method_pdf_reldir = "zPDF";
        if(methodName  == "linemodelsolve" ){

          CLineModelSolve Solve(ctx.m_ScopeStack,"galaxy",calibrationDirPath);
          Solve.Compute( ctx.GetInputContext(),
                         ctx.GetResultStore(),
                         ctx.m_ScopeStack);
          /*
            }else if(methodName  == "zweimodelsolve" ){

            CZweiModelSolve Solve(calibrationDirPath);
            mResult = Solve.Compute( ctx.GetDataStore(),
            ctx.GetSpectrum(),
            ctx.GetTemplateCatalog(),
            ctx.GetGalaxyCategoryList(),
            ctx.GetRayCatalog(),
            ctx.GetInputContext().m_lambdaRange,
            redshifts );
          */
        }else if(methodName  == "templatefittingsolve" ){
          CMethodTemplateFittingSolve solve(ctx.m_ScopeStack,"galaxy");
          solve.Compute( ctx.GetInputContext(),
                         ctx.GetResultStore(),
                         ctx.m_ScopeStack);

        }else if(methodName  == "tplcombinationsolve" ){
         CMethodTplcombinationSolve solve(ctx.m_ScopeStack,"galaxy");
          solve.Compute( ctx.GetInputContext(),
                         ctx.GetResultStore(),
                         ctx.m_ScopeStack);
        }
        /*
          }else if(methodName  == "linematching" ){
          COperatorLineMatchingSolve Solve;
          mResult = Solve.Compute(ctx.GetDataStore(),
          ctx.GetSpectrum(),
          lambdaRange,
          redshiftRange,
          redshiftStep,
          ctx.GetRayCatalog() );

        */
        else{
          throw GlobalException(INTERNAL_ERROR,"Problem found while parsing the method parameter");
        }
      }
    
    
    //estimate star/galaxy/qso classification
    Log.LogInfo("===============================================");
    if(!ctx.GetParameterStore()->Get<bool>( "enablelinemeassolve"))
      {
        {
          CClassificationSolve classifier(ctx.m_ScopeStack,"classification");
          classifier.Compute(ctx.GetInputContext(),
                             ctx.GetResultStore(),
                             ctx.m_ScopeStack);
        }
        {
          CReliabilitySolve reliability(ctx.m_ScopeStack,"reliability");
          reliability.Compute(ctx.GetInputContext(),
                              ctx.GetResultStore(),
                              ctx.m_ScopeStack);
        }
      }
    if(ctx.GetParameterStore()->Get<bool>( "enablelinemeassolve"))
      {
        CLineMeasSolve solve(ctx.m_ScopeStack,"linemeas");
        solve.Compute(ctx.GetInputContext(),
                          ctx.GetResultStore(),
                          ctx.m_ScopeStack);

      }
 }
    catch(GlobalException const&e)
    {
      throw e;
    }    
    catch(ParameterException const&e)
    {
      throw e;
    }    
  catch(std::exception const&e)
    {
      throw GlobalException(EXTERNAL_LIB_ERROR,Formatter()<<"ProcessFlow encountered an external lib error :"<<e.what());
    }    
}

