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
#ifndef _REDSHIFT_PROCESSFLOW_CONTEXT_
#define _REDSHIFT_PROCESSFLOW_CONTEXT_

#include <map>
#include <memory>
#include <string>

#include <gsl/gsl_errno.h>

#include "RedshiftLibrary/line/catalog.h"
#include "RedshiftLibrary/line/line.h"
#include "RedshiftLibrary/processflow/inputcontext.h"
#include "RedshiftLibrary/processflow/scopestack.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"
#include "RedshiftLibrary/spectrum/template/template.h"

#define Context (CProcessFlowContext::GetInstance())

namespace NSEpic {

class CSpectrum;
class CTemplateCatalog;
class CLineCatalog;
class CParameterStore;
class COperatorResultStore;
class CClassifierStore;
class CInputContext;
class CLine;
class CLineCatalogsTplRatio;
/**
 * \ingroup Redshift
 * Store all data concerning computation and processing of a given spectrum.
 */
class CProcessFlowContext : public CSingleton<CProcessFlowContext> {

public:
  void addSpectrum(const std::shared_ptr<CSpectrum> &spectrum) {
    m_inputContext->addSpectrum(spectrum);
  }
  void
  setTemplateCatalog(const std::shared_ptr<CTemplateCatalog> &templateCatalog) {
    m_inputContext->setTemplateCatalog(templateCatalog);
  }
  void
  setPhotBandCatalog(const std::shared_ptr<CPhotBandCatalog> &photBandCatalog) {
    m_inputContext->setPhotBandCatalog(photBandCatalog);
  }
  void setLineCatalog(const std::string &spectrumModel,
                      const std::string &method,
                      const std::shared_ptr<CLineCatalog> &catalog) {
    m_inputContext->setLineCatalog(spectrumModel, method, catalog);
  }
  void setLineRatioCatalogCatalog(
      const std::string &spectrumModel,
      const std::shared_ptr<CLineCatalogsTplRatio> &catalog) {
    m_inputContext->setLineRatioCatalogCatalog(spectrumModel, catalog);
  }
  void
  setfluxCorrectionMeiksin(const std::shared_ptr<CSpectrumFluxCorrectionMeiksin>
                               &igmcorrectionMeiksin) {
    m_inputContext->setfluxCorrectionMeiksin(igmcorrectionMeiksin);
  };
  void setfluxCorrectionCalzetti(
      const std::shared_ptr<CSpectrumFluxCorrectionCalzetti>
          &ismcorrectionCalzetti) {
    m_inputContext->setfluxCorrectionCalzetti(ismcorrectionCalzetti);
  };

  void Init();
  void reset();

  std::shared_ptr<const CParameterStore>
  LoadParameterStore(const std::string &paramsJSONString);

  std::shared_ptr<const CSpectrum> GetSpectrum(bool rebinned = false) const {
    return m_inputContext->GetSpectrum(rebinned);
  }
  std::shared_ptr<const CSpectrum> GetRebinnedSpectrum() const {
    return m_inputContext->GetRebinnedSpectrum();
  }
  std::shared_ptr<const CTemplateCatalog> GetTemplateCatalog() const {
    return m_inputContext->GetTemplateCatalog();
  }

  std::shared_ptr<const CLineCatalog>
  GetLineCatalog(const std::string &spectrumModel,
                 const std::string &method) const {
    return m_inputContext->GetLineCatalog(spectrumModel, method);
  }
  std::shared_ptr<const CParameterStore> GetParameterStore() const {
    return m_inputContext->GetParameterStore();
  }
  std::shared_ptr<const CInputContext> GetInputContext() const {
    return m_inputContext;
  }
  const std::shared_ptr<COperatorResultStore> &GetResultStore() {
    return m_ResultStore;
  }

  std::shared_ptr<const TFloat64Range> GetLambdaRange(int i = 0) const {
    return m_inputContext->getLambdaRange(i);
  }

  std::shared_ptr<const TFloat64Range> GetClampedLambdaRange(bool rebinned,
                                                             int i = 0) const {
    return m_inputContext->getClampedLambdaRange(rebinned, i);
  }

  const std::string &GetCurrentCategory() const {
    if (!m_ScopeStack->has_type(ScopeType::SPECTRUMMODEL))
      THROWG(SCOPESTACK_ERROR, Formatter() << "cannot get current "
                                           << ScopeType::SPECTRUMMODEL
                                           << " since not in scope");
    return m_ScopeStack->get_type_value(ScopeType::SPECTRUMMODEL);
  }

  const std::string &GetCurrentStage() const {
    if (!m_ScopeStack->has_type(ScopeType::STAGE))
      THROWG(SCOPESTACK_ERROR, Formatter()
                                   << "cannot get current " << ScopeType::STAGE
                                   << " since not in scope");
    return m_ScopeStack->get_type_value(ScopeType::STAGE);
  }

  const std::string &GetCurrentMethod() const {
    if (!m_ScopeStack->has_type(ScopeType::METHOD))
      THROWG(SCOPESTACK_ERROR, Formatter()
                                   << "cannot get current " << ScopeType::METHOD
                                   << " since not in scope");
    return m_ScopeStack->get_type_value(ScopeType::METHOD);
  }

  CLineMap getCLineMap();
  std::shared_ptr<CLineCatalogsTplRatio> GetTplRatioCatalog();
  std::shared_ptr<const CPhotBandCatalog> GetPhotBandCatalog();

  const std::vector<std::shared_ptr<const TFloat64Range>> &
  getClampedLambdaRanges(bool rebinned = false) const {
    if (rebinned)
      return m_inputContext->getRebinnedClampedLambdaRanges();
    else
      return m_inputContext->getClampedLambdaRanges();
  }

  const std::vector<std::shared_ptr<const TFloat64Range>> &
  getRebinnedClampedLambdaRanges() const {
    return m_inputContext->getRebinnedClampedLambdaRanges();
  }

  const std::vector<std::shared_ptr<const CSpectrum>> &
  getSpectra(bool rebinned = false) const {
    if (rebinned)
      return m_inputContext->getRebinnedSpectra();
    else
      return m_inputContext->getSpectra();
  }

  const std::vector<std::shared_ptr<const CSpectrum>> &
  getRebinnedSpectra() const {
    return m_inputContext->getRebinnedSpectra();
  }

  std::shared_ptr<CScopeStack> m_ScopeStack;

private:
  friend class CSingleton<CProcessFlowContext>;

  CProcessFlowContext();
  ~CProcessFlowContext() = default;
  std::shared_ptr<COperatorResultStore> m_ResultStore;
  std::shared_ptr<CParameterStore> m_parameterStore;
  std::shared_ptr<CInputContext> m_inputContext;

  // added below variables - to discuss if we only define them here (and no more
  // in processflow)

  TFloat64Range m_redshiftRange;
  TFloat64List m_redshifts;
  std::string m_redshiftSampling;
};

} // namespace NSEpic

#endif
