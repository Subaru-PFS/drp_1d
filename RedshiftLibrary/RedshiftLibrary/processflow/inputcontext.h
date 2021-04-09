#ifndef _INPUT_CONTEXT_H
#define _INPUT_CONTEXT_H

#include <memory>
#include <RedshiftLibrary/common/range.h>
namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class CRayCatalog;
class CParameterStore;

class CInputContext
{
 public:
  CInputContext(std::shared_ptr<CSpectrum> spc,
                std::shared_ptr<CTemplateCatalog> tmplCatalog,
                std::shared_ptr<const CRayCatalog> rayCatalog,
                std::shared_ptr<CParameterStore> paramStore);
  ~CInputContext(){}

  std::shared_ptr<CSpectrum> m_Spectrum;

  std::shared_ptr<CTemplateCatalog> m_TemplateCatalog;
  std::shared_ptr<const CRayCatalog> m_RayCatalog;

  std::shared_ptr<CParameterStore> m_ParameterStore;

  TFloat64Range m_lambdaRange;

private:

    void                                       InitSpectrum();
  //This should be a method of CTemplateCatalog, withou calibrationDirPath
    void                                       InitIsmIgm(const std::string & CalibrationDirPath);

};

} 
#endif
