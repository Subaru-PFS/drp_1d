#ifndef _REDSHIFT_PROCESSFLOW_CONTEXT_
#define _REDSHIFT_PROCESSFLOW_CONTEXT_

#include <RedshiftLibrary/processflow/datastore.h>
#include <RedshiftLibrary/ray/catalog.h>
#include <RedshiftLibrary/ray/ray.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/io/reader.h>

#include <RedshiftLibrary/reliability/zclassifierstore.h>

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

/**
 * \ingroup Redshift
 * Store all data concerning computation and processing of a given spectrum.
 */
class CProcessFlowContext
{

public:

    CProcessFlowContext();
    ~CProcessFlowContext();

    bool Init(std::shared_ptr<CSpectrum> spectrum,
              std::string processingID,
              std::shared_ptr<const CTemplateCatalog> templateCatalog,
              std::shared_ptr<const CRayCatalog> rayCatalog,
              std::shared_ptr<CParameterStore> paramStore,
              std::shared_ptr<CClassifierStore> zqualStore);

    const CSpectrum&                GetSpectrum() const;
    const CTemplateCatalog&         GetTemplateCatalog() const;
    const CRayCatalog&              GetRayCatalog() const;

    bool correctSpectrum(Float64 LambdaMin, Float64 LambdaMax);

    CParameterStore&                GetParameterStore();
    COperatorResultStore&           GetResultStore();
    CDataStore&                     GetDataStore();
    CClassifierStore&               GetClassifierStore();

private:

    void                                       InitSpectrum();
    void                                       InitIsmIgm(const std::string & CalibrationDirPath);

    std::shared_ptr<CSpectrum>                 m_Spectrum;

    std::shared_ptr<const CTemplateCatalog>    m_TemplateCatalog;
    std::shared_ptr<const CRayCatalog>         m_RayCatalog;

    std::shared_ptr<CParameterStore>           m_ParameterStore;
    std::shared_ptr<COperatorResultStore>      m_ResultStore;

    std::shared_ptr<CDataStore>                m_DataStore;
    std::shared_ptr<CClassifierStore>          m_ClassifierStore;

};


}

#endif
