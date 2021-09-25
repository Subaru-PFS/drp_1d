// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#ifndef _INPUT_CONTEXT_H
#define _INPUT_CONTEXT_H

#include <memory>
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/spectrum/logrebinning.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/linemodel/templatesortho.h"
namespace NSEpic
{
class CSpectrumLogRebinning;//forward declaration
class CSpectrum;
class CTemplateCatalog;
class CRayCatalog;
class CParameterStore;
class CPhotBandCatalog;

class CInputContext
{
 public:
  CInputContext(std::shared_ptr<CSpectrum> spc,
                std::shared_ptr<CTemplateCatalog> tmplCatalog,
                std::shared_ptr<CRayCatalog> gal_rayCatalog,
                std::shared_ptr<CRayCatalog> qso_rayCatalog,
                std::shared_ptr<CPhotBandCatalog> photBandCatalog,
                std::shared_ptr<CParameterStore> paramStore);

  // const getters
  std::shared_ptr<const CSpectrum> GetSpectrum() const {return m_Spectrum;}
  std::shared_ptr<const CSpectrum>  GetRebinnedSpectrum() const {return m_rebinnedSpectrum;}
  std::shared_ptr<const CTemplateCatalog> GetTemplateCatalog() const {m_TemplateCatalog->resetCatalogState();return m_TemplateCatalog;}
  std::shared_ptr<const CRayCatalog> GetRayCatalog(const std::string &objectType) const;
  std::shared_ptr<const CPhotBandCatalog> GetPhotBandCatalog() const {return m_photBandCatalog;}
  std::shared_ptr<const CParameterStore> GetParameterStore() const {return m_ParameterStore;}

  // mutable getters
  std::shared_ptr<CSpectrum>  GetSpectrum() {return m_Spectrum;}
  std::shared_ptr<CSpectrum>  GetRebinnedSpectrum() {return m_rebinnedSpectrum;}
  std::shared_ptr<CTemplateCatalog>  GetTemplateCatalog() {m_TemplateCatalog->resetCatalogState();return m_TemplateCatalog;}
  std::shared_ptr<CRayCatalog>  GetRayCatalog(const std::string &objectType);
  std::shared_ptr<CPhotBandCatalog> GetPhotBandCatalog() {return m_photBandCatalog;}
  std::shared_ptr<CParameterStore> GetParameterStore() {return m_ParameterStore;}

  void SetRebinnedSpectrum(std::shared_ptr<CSpectrum> rebinnedSpc){m_rebinnedSpectrum = rebinnedSpc;}
  TFloat64Range   m_lambdaRange;    
  Bool            m_use_LogLambaSpectrum = 0;
  Float64         m_logGridStep;
  typedef struct{
    TFloat64Range zrange;
    }SRebinResults;
  std::map<std::string, SRebinResults> m_logRebin;
  
  std::vector<std::string> m_categories{"galaxy", "qso", "star"};
private:

  std::shared_ptr<CSpectrum> m_Spectrum;
  std::shared_ptr<CSpectrum> m_rebinnedSpectrum;
  std::shared_ptr<CTemplateCatalog> m_TemplateCatalog;
  std::shared_ptr<CRayCatalog> m_gal_RayCatalog;
  std::shared_ptr<CRayCatalog> m_qso_RayCatalog;
  std::shared_ptr<CParameterStore> m_ParameterStore;
  std::shared_ptr<CPhotBandCatalog> m_photBandCatalog;

  void OrthogonalizeTemplates();
  void RebinInputs();
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
      throw GlobalException(INTERNAL_ERROR,"CInputContext::GetRayCatalog: invalid object type");
  }
}

} 
#endif
