#ifndef _INPUT_CONTEXT_H
#define _INPUT_CONTEXT_H

#include <memory>
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/log/log.h"
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
                std::shared_ptr<CRayCatalog> gal_rayCatalog,
                std::shared_ptr<CRayCatalog> qso_rayCatalog,
                std::shared_ptr<CParameterStore> paramStore);

  // const getters
  std::shared_ptr<const CSpectrum> GetSpectrum() const {return m_Spectrum;}
  std::shared_ptr<const CSpectrum>  GetRebinnedSpectrum() const {return m_rebinnedSpectrum;}
  std::shared_ptr<const CTemplateCatalog> GetTemplateCatalog() const {return m_TemplateCatalog;}
  std::shared_ptr<const CRayCatalog> GetRayCatalog(const std::string &objectType) const;
  std::shared_ptr<const CParameterStore> GetParameterStore() const {return m_ParameterStore;}
  void RebinInputs();

  // mutable getters
  std::shared_ptr<CSpectrum>  GetSpectrum() {return m_Spectrum;}
  std::shared_ptr<CSpectrum>  GetRebinnedSpectrum() {return m_rebinnedSpectrum;}
  std::shared_ptr<CTemplateCatalog>  GetTemplateCatalog() {return m_TemplateCatalog;}
  std::shared_ptr<CRayCatalog>  GetRayCatalog(const std::string &objectType);
  std::shared_ptr<CParameterStore> GetParameterStore() {return m_ParameterStore;}

  void SetRebinnedSpectrum(std::shared_ptr<CSpectrum> rebinnedSpc){m_rebinnedSpectrum = rebinnedSpc;}
  TFloat64Range   m_lambdaRange;    
  TFloat64Range   m_redshiftRangeFFT;//computed with logRebinning 
  Float64         m_redshiftStepFFT; 
  Bool            m_use_LogLambaSpectrum = 0;

private:

  std::shared_ptr<CSpectrum> m_Spectrum;
  std::shared_ptr<CSpectrum> m_rebinnedSpectrum;
  std::shared_ptr<CTemplateCatalog> m_TemplateCatalog;
  std::shared_ptr<CRayCatalog> m_gal_RayCatalog;
  std::shared_ptr<CRayCatalog> m_qso_RayCatalog;
  std::shared_ptr<CParameterStore> m_ParameterStore;

  void RebinInputWrapper();
  void validateSpectrum(std::shared_ptr<CSpectrum> spectrum, 
                        TFloat64Range lambdaRange, 
                        Bool enableInputSpcCorrect);
};

inline
std::shared_ptr<const CRayCatalog> CInputContext::GetRayCatalog(const std::string &objectType) const 
{
  return const_cast<CInputContext*>(this)->GetRayCatalog(objectType); 
}

inline
std::shared_ptr<CRayCatalog> CInputContext::GetRayCatalog(const std::string &objectType)  
{
  if (objectType=="galaxy")
    return m_gal_RayCatalog;
  else if (objectType=="qso")
    return m_qso_RayCatalog;
  else {
      Log.LogError("CInputContext::GetRayCatalog: invalid object type");
      throw std::runtime_error("CInputContext::GetRayCatalog: invalid object type");
  }
}

} 
#endif
