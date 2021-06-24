#include "RedshiftLibrary/processflow/context.h"
#include "RedshiftLibrary/processflow/resultstore.h"
#include "RedshiftLibrary/processflow/inputcontext.h"

#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/debug/assert.h"

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
                               std::shared_ptr<CRayCatalog> galaxy_rayCatalog,
                               std::shared_ptr<CRayCatalog> qso_rayCatalog,
                               const std::string& paramsJSONString)
{
  Log.LogInfo("Processing context initialization");

  std::shared_ptr<CParameterStore> parameterStore = std::make_shared<CParameterStore>(m_ScopeStack);
  parameterStore->FromString(paramsJSONString);

//  CInputContext *ic = new CInputContext(spectrum,templateCatalog,rayCatalog,parameterStore) ; 
  m_inputContext = std::make_shared<const CInputContext>(spectrum,templateCatalog,galaxy_rayCatalog,qso_rayCatalog,parameterStore);

  m_ResultStore = std::make_shared<COperatorResultStore>(m_ScopeStack);

}

void CProcessFlowContext::testResultStore() {
  m_ResultStore = std::make_shared<COperatorResultStore>(m_ScopeStack);
    m_ResultStore->test();
  }
