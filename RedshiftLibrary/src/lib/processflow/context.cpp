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
#include "RedshiftLibrary/processflow/context.h"
#include "RedshiftLibrary/processflow/autoscope.h"
#include "RedshiftLibrary/processflow/inputcontext.h"
#include "RedshiftLibrary/processflow/resultstore.h"

#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/log/log.h"

#include <boost/filesystem.hpp>

namespace bfs = boost::filesystem;

using namespace NSEpic;

static void NewHandler(const char *reason, const char *file, int line,
                       int gsl_errno) {
  THROWG(EXTERNAL_LIB_ERROR,
         Formatter() << "GSL Error : "
                     << " gsl: " << file << ":" << line << ": ERROR:" << reason
                     << " (Errtype: " << gsl_strerror(gsl_errno) << ")");
  return;
}

CProcessFlowContext::CProcessFlowContext() {
  gsl_set_error_handler(NewHandler);
  m_parameterStore = std::make_shared<CParameterStore>(m_ScopeStack);
  m_ResultStore = std::make_shared<COperatorResultStore>(m_ScopeStack);
  m_inputContext = std::make_shared<CInputContext>(m_parameterStore);
}

std::shared_ptr<const CParameterStore>
CProcessFlowContext::LoadParameterStore(const std::string &paramsJSONString) {
  m_parameterStore->FromString(paramsJSONString);
  return m_parameterStore;
}

void CProcessFlowContext::Init() {
  Log.LogInfo("Processing context initialization");
  m_inputContext->Init();
}

void CProcessFlowContext::reset() {
  m_ResultStore->reset();
  m_inputContext->resetSpectrumSpecific();
}

CLineMap CProcessFlowContext::getCLineMap() {
  CAutoScope autoscope(m_ScopeStack, "lineModel");
  return m_inputContext->GetFilteredLineMap(
      GetCurrentCategory(), GetCurrentMethod(),
      m_parameterStore->GetScoped<std::string>("lineTypeFilter"),
      m_parameterStore->GetScoped<std::string>("lineForceFilter"));
}

std::shared_ptr<CLineCatalogsTplRatio>
CProcessFlowContext::GetTplRatioCatalog() {
  return m_inputContext->GetTemplateRatioCatalog(m_ScopeStack[0]);
}

std::shared_ptr<const CPhotBandCatalog>
CProcessFlowContext::GetPhotBandCatalog() {
  return m_inputContext->GetPhotBandCatalog();
}
