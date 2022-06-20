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

#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/linemodel/templatesortho.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/spectrum/logrebinning.h"
#include <memory>
namespace NSEpic {
class CSpectrumLogRebinning; // forward declaration
class CSpectrum;
class CTemplateCatalog;
class CLineCatalog;
class CParameterStore;
class CPhotBandCatalog;
class CLineCatalogsTplShape;
class CLine;
typedef std::vector<CLine> TLineVector;

class CInputContext {
public:
  CInputContext(std::shared_ptr<CParameterStore> paramStore);

  // const getters
  std::shared_ptr<const CSpectrum> GetSpectrum() const { return m_Spectrum; }
  std::shared_ptr<const CSpectrum> GetRebinnedSpectrum() const {
    return m_rebinnedSpectrum;
  }
  std::shared_ptr<const CTemplateCatalog> GetTemplateCatalog() const {
    return m_TemplateCatalog;
  }
  std::shared_ptr<const CLineCatalogsTplShape>
  GetTemplateRatioCatalog(const std::string &objectType) const;
  std::shared_ptr<const CLineCatalog>
  GetLineCatalog(const std::string &objectType,
                 const std::string &method) const;
  std::shared_ptr<const CPhotBandCatalog> GetPhotBandCatalog() const {
    return m_photBandCatalog;
  }
  std::shared_ptr<const CParameterStore> GetParameterStore() const {
    return m_ParameterStore;
  }

  // mutable getters
  const std::shared_ptr<CSpectrum> &GetSpectrum() { return m_Spectrum; }
  const std::shared_ptr<CSpectrum> &GetRebinnedSpectrum() {
    return m_rebinnedSpectrum;
  }
  std::shared_ptr<CTemplateCatalog> GetTemplateCatalog() {
    return m_TemplateCatalog;
  }
  const std::shared_ptr<CLineCatalogsTplShape> &
  GetTemplateRatioCatalog(const std::string &objectType);
  const std::shared_ptr<CLineCatalog> &
  GetLineCatalog(const std::string &objectType, const std::string &method);
  const CLineCatalog::TLineVector
  GetFilteredLineVector(const std::string &objectType,
                        const std::string &method, const std::string &type,
                        const std::string &force);
  const std::shared_ptr<CPhotBandCatalog> &GetPhotBandCatalog() {
    return m_photBandCatalog;
  }
  const std::shared_ptr<CParameterStore> &GetParameterStore() {
    return m_ParameterStore;
  }

  void SetRebinnedSpectrum(const std::shared_ptr<CSpectrum> &rebinnedSpc) {
    m_rebinnedSpectrum = rebinnedSpc;
  }
  TFloat64Range m_lambdaRange;
  bool m_use_LogLambaSpectrum = 0;
  Float64 m_logGridStep;
  typedef struct {
    TFloat64Range zrange;
  } SRebinResults;
  std::map<std::string, SRebinResults> m_logRebin;

  // fill this from api, based on a new parameter in parameters.json we call
  // object_types, it's a list consisting of the object types to process, it
  // replaces enable<object_type>solve
  TStringList
      m_categories; //{"galaxy", "qso", "star"}; rename this to object_type

  void setLineCatalog(const std::string &objectType, const std::string &method,
                      const std::shared_ptr<CLineCatalog> &catalog);
  void setLineRatioCatalogCatalog(
      const std::string &objectType,
      const std::shared_ptr<CLineCatalogsTplShape> &catalog);
  void
  setTemplateCatalog(const std::shared_ptr<CTemplateCatalog> &templateCatalog) {
    m_TemplateCatalog = templateCatalog;
  }
  void
  setPhotBandCatalog(const std::shared_ptr<CPhotBandCatalog> &photBandCatalog) {
    m_photBandCatalog = photBandCatalog;
  }
  void setSpectrum(const std::shared_ptr<CSpectrum> &spectrum) {
    m_Spectrum = spectrum;
  }
  void
  setfluxCorrectionMeiksin(const std::shared_ptr<CSpectrumFluxCorrectionMeiksin>
                               &igmcorrectionMeiksin) {
    m_igmcorrectionMeiksin = igmcorrectionMeiksin;
  }
  void setfluxCorrectionCalzetti(
      const std::shared_ptr<CSpectrumFluxCorrectionCalzetti>
          &ismcorrectionCalzetti) {
    m_ismcorrectionCalzetti = ismcorrectionCalzetti;
  }
  void Init();

  void resetSpectrumSpecific();

private:
  std::shared_ptr<CSpectrum> m_Spectrum;
  std::shared_ptr<CSpectrum> m_rebinnedSpectrum;
  std::shared_ptr<CTemplateCatalog> m_TemplateCatalog;
  std::map<std::string, std::map<std::string, std::shared_ptr<CLineCatalog>>>
      m_lineCatalogs;
  std::map<std::string, std::shared_ptr<CLineCatalogsTplShape>>
      m_lineRatioCatalogCatalogs;
  std::shared_ptr<CSpectrumFluxCorrectionMeiksin> m_igmcorrectionMeiksin;
  std::shared_ptr<CSpectrumFluxCorrectionCalzetti> m_ismcorrectionCalzetti;
  std::shared_ptr<CParameterStore> m_ParameterStore;
  std::shared_ptr<CPhotBandCatalog> m_photBandCatalog;
  void OrthogonalizeTemplates();
  void RebinInputs();
};

inline std::shared_ptr<const CLineCatalog>
CInputContext::GetLineCatalog(const std::string &objectType,
                              const std::string &method) const {
  return const_cast<CInputContext *>(this)->GetLineCatalog(objectType, method);
}

inline const std::shared_ptr<CLineCatalog> &
CInputContext::GetLineCatalog(const std::string &objectType,
                              const std::string &method) {
  //  if (std::findm_categories.find(objectType))
  // throw
  // GlobalException(ErrorCode::INTERNAL_ERROR,"CInputContext::GetLineCatalog:
  // invalid object type");
  return m_lineCatalogs[objectType][method];
}

inline const CLineCatalog::TLineVector CInputContext::GetFilteredLineVector(
    const std::string &objectType, const std::string &method,
    const std::string &type, const std::string &force) {
  //  if (std::findm_categories.find(objectType))
  // throw
  // GlobalException(ErrorCode::INTERNAL_ERROR,"CInputContext::GetLineCatalog:
  // invalid object type");
  return m_lineCatalogs[objectType][method]->GetFilteredList(type, force);
}

inline std::shared_ptr<const CLineCatalogsTplShape>
CInputContext::GetTemplateRatioCatalog(const std::string &objectType) const {
  return const_cast<CInputContext *>(this)->GetTemplateRatioCatalog(objectType);
}

inline const std::shared_ptr<CLineCatalogsTplShape> &
CInputContext::GetTemplateRatioCatalog(const std::string &objectType) {
  //  if (std::findm_categories.find(objectType))
  // throw
  // GlobalException(ErrorCode::INTERNAL_ERROR,"CInputContext::GetTemplateRatioCatalog:
  // invalid object type");
  return m_lineRatioCatalogCatalogs[objectType];
}

} // namespace NSEpic
#endif
