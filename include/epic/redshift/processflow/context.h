#ifndef _REDSHIFT_PROCESSFLOW_CONTEXT_
#define _REDSHIFT_PROCESSFLOW_CONTEXT_

#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/processflow/datastore.h>

#include <map>
#include <string>
#include <memory>

namespace NSEpic
{

class CSpectrum;
class CContinuum;
class CPeakStore;
class CRayCatalog;
class CTemplateCatalog;
class CRayCatalog;
class CParameterStore;
class COperatorResultStore;
class CDataStore;

/**
 * \ingroup Redshift
 * Store all data concerning computation and processing of a given spectrum.
 */
class CProcessFlowContext
{

public:

    CProcessFlowContext();
    ~CProcessFlowContext();

    bool Init( const char* spectrumPath, const char* noisePath,
               const char* tempalteCatalogPath, const char* rayCatalogPath,
               std::shared_ptr<CParameterStore> paramStore  );

    bool Init( const char* spectrumPath, const char* noisePath,
               std::shared_ptr<const CTemplateCatalog> templateCatalog,
               std::shared_ptr<const CRayCatalog> rayCatalog,
               std::shared_ptr<CParameterStore> paramStore  );

    const CSpectrum&                GetSpectrum() const;
    const CSpectrum&                GetSpectrumWithoutContinuum() const;
    const CTemplateCatalog&         GetTemplateCatalog() const;
    const CRayCatalog&              GetRayCatalog() const;

    CParameterStore&                GetParameterStore();
    COperatorResultStore&           GetResultStore();
    CDataStore&                     GetDataStore();

private:

    std::shared_ptr<CSpectrum>                 m_Spectrum;
    std::shared_ptr<CSpectrum>                 m_SpectrumWithoutContinuum;

    std::shared_ptr<const CTemplateCatalog>    m_TemplateCatalog;
    std::shared_ptr<const CRayCatalog>         m_RayCatalog;

    std::shared_ptr<CParameterStore>           m_ParameterStore;
    std::shared_ptr<COperatorResultStore>      m_ResultStore;

    std::shared_ptr<CDataStore>                m_DataStore;

};


}

#endif
