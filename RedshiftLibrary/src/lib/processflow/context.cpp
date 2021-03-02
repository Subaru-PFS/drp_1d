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
                               std::shared_ptr<const CTemplateCatalog> templateCatalog,
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
}



/*
CParameterStore& CProcessFlowContext::GetParameterStore()
{
    return *m_ParameterStore;
}

COperatorResultStore& CProcessFlowContext::GetResultStore()
{
    return *m_ResultStore;
}


CDataStore& CProcessFlowContext::GetDataStore()
{
    return *m_DataStore;
}
*/
bool CProcessFlowContext::correctSpectrum(Float64 LambdaMin,  Float64 LambdaMax)
{
    return m_Spectrum->correctSpectrum( LambdaMin, LambdaMax );
}
/*
const CSpectrum& CProcessFlowContext::GetSpectrum() const
{
    return *m_Spectrum;
}

const CTemplateCatalog& CProcessFlowContext::GetTemplateCatalog() const
{
    return *m_TemplateCatalog;
}
const TStringList&  CProcessFlowContext::GetGalaxyCategoryList() const
{
    return m_filteredGalaxyTemplateCategoryList;
}
const TStringList&  CProcessFlowContext::GetStarCategoryList() const
{
    return m_filteredStarTemplateCategoryList;                                ; 
}
const TStringList&  CProcessFlowContext::GetQSOCategoryList() const
{
    return m_filteredQSOTemplateCategoryList;
}

const CRayCatalog& CProcessFlowContext::GetRayCatalog() const
{
    return *m_RayCatalog;
}
*/
