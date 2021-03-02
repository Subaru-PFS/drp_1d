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
            std::shared_ptr<const CTemplateCatalog> templateCatalog,
            std::shared_ptr<const CRayCatalog> rayCatalog,
            const std::string& paramsJSONString
            );

    

    bool correctSpectrum(Float64 LambdaMin, Float64 LambdaMax);

  /*
    const CSpectrum&                GetSpectrum() const;
    const CTemplateCatalog&         GetTemplateCatalog() const;
    const CRayCatalog&              GetRayCatalog() const;
    CParameterStore&                GetParameterStore();
    COperatorResultStore&           GetResultStore();
  */
    
  std::shared_ptr<const CSpectrum> GetSpectrum(){return m_inputContext->m_Spectrum;}
  std::shared_ptr<const CTemplateCatalog> GetTemplateCatalog(){return m_inputContext->m_TemplateCatalog;}
  std::shared_ptr<const CRayCatalog> GetRayCatalog(){return m_inputContext->m_RayCatalog;}
  std::shared_ptr<CParameterStore> GetParameterStore(){return m_inputContext->m_ParameterStore;}
  const CInputContext& GetInputContext(){return *(m_inputContext.get());}
  COperatorResultStore& GetResultStore(){return *(m_ResultStore.get());}

  /*
  const CSpectrum& GetSpectrum(){return *(m_inputContext->m_Spectrum);}
  const CTemplateCatalog& GetTemplateCatalog(){return *(m_inputContext->m_TemplateCatalog);}
  const CRayCatalog& GetRayCatalog() {return *(m_inputContext->m_RayCatalog);}
  const CParameterStore& GetParameterStore(){return *(m_inputContext->m_ParameterStore);}
  std::shared_ptr<COperatorResultStore> GetResultStore(){return m_ResultStore;}
  */  
//std::shared_ptr<CDataStore> GetDataStore();
  //  TScopeStack &getScopeStack(){return m_ScopeStack;}
 TScopeStack                     m_ScopeStack;
private:

    std::shared_ptr<CSpectrum>                 m_Spectrum;

    std::shared_ptr<const CTemplateCatalog>    m_TemplateCatalog;
    std::shared_ptr<const CRayCatalog>         m_RayCatalog;

    std::shared_ptr<CParameterStore>           m_ParameterStore;
    std::shared_ptr<COperatorResultStore>      m_ResultStore;

    std::shared_ptr<CDataStore>                m_DataStore;
  std::shared_ptr<CInputContext>  m_inputContext;
 


};


}

#endif
