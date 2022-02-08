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


#include "RedshiftLibrary/processflow/inputcontext.h"
#include "RedshiftLibrary/ray/catalog.h"
#include "RedshiftLibrary/ray/ray.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/linemodel/calibrationconfig.h"

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

  void setSpectrum(std::shared_ptr<CSpectrum> spectrum){ m_inputContext->setSpectrum(spectrum);}
  void setTemplateCatalog(std::shared_ptr<CTemplateCatalog> templateCatalog){ m_inputContext->setTemplateCatalog(templateCatalog);}
  void setPhotBandCatalog(std::shared_ptr<CPhotBandCatalog> photBandCatalog){ m_inputContext->setPhotBandCatalog(photBandCatalog);}
  void setLineCatalog(const std::string& objectType,std::shared_ptr<CRayCatalog> catalog) { m_inputContext->setLineCatalog(objectType,catalog);}
  void setLineRatioCatalogCatalog(const std::string& objectType,std::shared_ptr<CRayCatalogsTplShape> catalog) { m_inputContext->setLineRatioCatalogCatalog(objectType,catalog);}
  void Init();
  std::shared_ptr<const CParameterStore> LoadParameterStore(const std::string& paramsJSONString);
  std::shared_ptr<const CSpectrum> GetRebinnedSpectrum() const {return m_inputContext->GetRebinnedSpectrum();}  
  std::shared_ptr<const CSpectrum> GetSpectrum() const {return m_inputContext->GetSpectrum();}
  std::shared_ptr<const CTemplateCatalog> GetTemplateCatalog() const {return m_inputContext->GetTemplateCatalog();}
  std::shared_ptr<const CRayCatalog> GetRayCatalog(const std::string &objectType) const {return m_inputContext->GetRayCatalog(objectType);}
  std::shared_ptr<const CParameterStore> GetParameterStore() const {return m_inputContext->GetParameterStore();}
  std::shared_ptr<const CInputContext> GetInputContext() const {return m_inputContext;}
  std::shared_ptr<COperatorResultStore> GetResultStore(){return m_ResultStore;}

  TScopeStack                     m_ScopeStack;

  void testResultStore(); 
private:

    std::shared_ptr<COperatorResultStore> m_ResultStore;
    std::shared_ptr<CParameterStore>      m_parameterStore;
    std::shared_ptr<CInputContext>  m_inputContext;
 
    //added below variables - to discuss if we only define them here (and no more in processflow)
    TFloat64Range m_spclambdaRange;
    TFloat64Range m_redshiftRange;
    TFloat64List  m_redshifts;
    std::string   m_redshiftSampling;
};

}

#endif
