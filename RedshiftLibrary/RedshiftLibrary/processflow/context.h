#ifndef _REDSHIFT_PROCESSFLOW_CONTEXT_
#define _REDSHIFT_PROCESSFLOW_CONTEXT_

#include <RedshiftLibrary/processflow/datastore.h>
#include <RedshiftLibrary/processflow/inputcontext.h>
#include <RedshiftLibrary/ray/catalog.h>
#include <RedshiftLibrary/ray/ray.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/io/reader.h>
#include <RedshiftLibrary/linemodel/calibrationconfig.h> 

#include <gsl/gsl_errno.h>

#include <map>
#include <string>
#include <memory>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class CRayCatalog;
class CParameterStore;
class COperatorResultStore;
class CDataStore;
class CClassifierStore;
class CInputContext;
/**
 * \ingroup Redshift
 * Store all data concerning computation and processing of a given spectrum.
 */
class CProcessFlowContext
{

public:

    CProcessFlowContext();
    ~CProcessFlowContext();

  void Init(std::shared_ptr<CSpectrum> spectrum,
            std::shared_ptr<CTemplateCatalog> templateCatalog,
            std::shared_ptr<CRayCatalog> rayCatalog,
            const std::string& paramsJSONString
            );
    
  std::shared_ptr<const CSpectrum> GetSpectrum() const {return m_inputContext->GetSpectrum();}
  std::shared_ptr<const CTemplateCatalog> GetTemplateCatalog() const {return m_inputContext->GetTemplateCatalog();}
  std::shared_ptr<const CRayCatalog> GetRayCatalog() const {return m_inputContext->GetRayCatalog();}
  std::shared_ptr<const CParameterStore> GetParameterStore() const {return m_inputContext->GetParameterStore();}
  std::shared_ptr<const CInputContext> GetInputContext() const {return m_inputContext;}
  std::shared_ptr<COperatorResultStore> GetResultStore(){return m_ResultStore;}

  TScopeStack                     m_ScopeStack;
private:

    std::shared_ptr<COperatorResultStore>  m_ResultStore;

  std::shared_ptr<const CInputContext>  m_inputContext;
 
    //added below variables - to discuss if we only define them here (and no more in processflow)
    TFloat64Range m_spclambdaRange;
    TFloat64Range m_redshiftRange;
    TFloat64List  m_redshifts;
    std::string   m_redshiftSampling;
};

}

#endif
