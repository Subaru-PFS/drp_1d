#ifndef _REDSHIFT_PROCESSFLOW_CONTEXT_
#define _REDSHIFT_PROCESSFLOW_CONTEXT_

#include <epic/core/common/ref.h>
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
 * Store all data concerning computation and processign of a given spectrum
 *
 */
class CProcessFlowContext : public CManagedObject
{

    DEFINE_MANAGED_OBJECT( CProcessFlowContext )

public:

    CProcessFlowContext();
    ~CProcessFlowContext();

    bool Init( const char* spectrumPath, const char* noisePath,
               const char* tempalteCatalogPath, const char* rayCatalogPath,
               CParameterStore& paramStore  );

    bool Init( const char* spectrumPath, const char* noisePath,
               const CTemplateCatalog& templateCatalog, const CRayCatalog& rayCatalog,
               CParameterStore& paramStore  );

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

    CRef<CTemplateCatalog>          m_TemplateCatalog;
    CRef<CRayCatalog>               m_RayCatalog;


    CRef<CParameterStore>           m_ParameterStore;
    CRef<COperatorResultStore>      m_ResultStore;

    CRef<CDataStore>                m_DataStore;



};


}

#endif
