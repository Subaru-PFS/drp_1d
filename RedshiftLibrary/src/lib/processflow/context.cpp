#include <RedshiftLibrary/processflow/context.h>
#include <RedshiftLibrary/processflow/inputcontext.h>

#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/debug/assert.h>

#include <boost/filesystem.hpp>

namespace bfs = boost::filesystem;

using namespace NSEpic;

static void NewHandler(const char* reason,
                       const char* file,
                       int line,
                       int gsl_errno){

    Log.LogError(" gsl: %s:%d: ERROR: %s (Errtype: %s)",file, line, reason, gsl_strerror(gsl_errno));
    throw std::runtime_error("GSL Error");
    return ;
}

CProcessFlowContext::CProcessFlowContext()
{
    gsl_set_error_handler(NewHandler);
}

CProcessFlowContext::~CProcessFlowContext()
{

}
void CProcessFlowContext::Init(std::shared_ptr<CSpectrum> spectrum,
                               std::shared_ptr<CTemplateCatalog> templateCatalog,
                               std::shared_ptr<const CRayCatalog> rayCatalog,
                               const std::string& paramsJSONString)
{
  Log.LogInfo("Processing context initialization");

  std::shared_ptr<CParameterStore> parameterStore = std::shared_ptr<CParameterStore>(new CParameterStore(m_ScopeStack));
  parameterStore->FromString(paramsJSONString);

//  CInputContext *ic = new CInputContext(spectrum,templateCatalog,rayCatalog,parameterStore) ; 
  m_inputContext = std::shared_ptr<CInputContext>(new CInputContext(spectrum,templateCatalog,rayCatalog,parameterStore));

  m_ResultStore = std::shared_ptr<COperatorResultStore>( new COperatorResultStore(m_ScopeStack) );

  TFloat64Range lambdaRange = parameterStore->Get<TFloat64Range>("lambdarange");
  spectrum->GetSpectralAxis().ClampLambdaRange( lambdaRange, m_inputContext->m_lambdaRange );
  Log.LogInfo( "Processing spc: (CLambdaRange: %f-%f:%f)",
               m_inputContext->m_lambdaRange.GetBegin(),
               m_inputContext->m_lambdaRange.GetEnd(),
               spectrum->GetResolution());

  //************************************
  const Float64 lmin = GetInputContext()->m_lambdaRange.GetBegin();
  const Float64 lmax = GetInputContext()->m_lambdaRange.GetEnd();

  std::string enableInputSpcCorrectStr = parameterStore->Get<std::string>( "autocorrectinput");
  Bool enableInputSpcCorrect = enableInputSpcCorrectStr == "yes";
  if(enableInputSpcCorrect)
    {
      //Check if the Spectrum is valid on the lambdarange
      //correctInputSpectrum(ctx.GetInputContext()->m_lambdaRange);

      if( spectrum->correctSpectrum( lmin,lmax ))
        Log.LogInfo( "Successfully corrected noise on wavelength range (%.1f ; %.1f)",  lmin, lmax );
      }

   if( !spectrum->IsFluxValid( lmin, lmax ) ){
      Log.LogError("Failed to validate spectrum flux on wavelength range (%.1f ; %.1f)",
                   lmin, lmax );
      throw std::runtime_error("Failed to validate spectrum flux");
    }else{
      Log.LogDetail( "Successfully validated spectrum flux, on wavelength range (%.1f ; %.1f)", lmin, lmax );
    }
	//Check if the noise is valid in the lambdarange
    if( !spectrum->IsNoiseValid( lmin, lmax ) ){
      Log.LogError("Failed to validate noise on wavelength range (%.1f ; %.1f)",
                   lmin, lmax );
      throw std::runtime_error("Failed to validate noise from spectrum");
    }else{
      Log.LogDetail( "Successfully validated noise on wavelength range (%.1f ; %.1f)", lmin, lmax );
    }
  
}



