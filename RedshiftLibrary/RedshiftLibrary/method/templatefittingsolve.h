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
#ifndef _REDSHIFT_METHOD_TEMPLATEFITTINGSOLVE_
#define _REDSHIFT_METHOD_TEMPLATEFITTINGSOLVE_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/method/solve.h"
#include "RedshiftLibrary/processflow/resultstore.h"
#include "RedshiftLibrary/processflow/inputcontext.h"
#include "RedshiftLibrary/method/templatefittingsolveresult.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "RedshiftLibrary/operator/templatefittingBase.h"
#include "RedshiftLibrary/operator/pdfz.h"
#include "RedshiftLibrary/operator/pdfMargZLogResult.h"

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;

/**
 * \ingroup Redshift
 */
  class CMethodTemplateFittingSolve : public CSolve
{

 public:

    enum EType
    {
             nType_raw = 1,
             nType_continuumOnly = 2,
             nType_noContinuum = 3,
             nType_all = 4,
    };


  CMethodTemplateFittingSolve(TScopeStack &scope,std::string objectType);

private:

  std::shared_ptr<CSolveResult> compute(std::shared_ptr<const CInputContext> inputContext,
                                        std::shared_ptr<COperatorResultStore> resultStore,
                                        TScopeStack &scope) override;


  bool Solve(std::shared_ptr<COperatorResultStore> resultStore,
               const CSpectrum& spc,
               const std::shared_ptr<const CTemplate> & tpl,
               Float64 overlapThreshold,
               std::vector<CMask> maskList,
               EType spctype=nType_raw,
               std::string opt_interp="lin",
               bool opt_extinction=false,
               bool opt_dustFitting=false);

  ChisquareArray BuildChisquareArray(std::shared_ptr<const COperatorResultStore> store, const std::string & scopeStr) const;

  std::shared_ptr<const ExtremaResult>  SaveExtremaResult(  shared_ptr<const COperatorResultStore> store,                        
                                                            const std::string & scopeStr,
                                                            const TCandidateZbyRank & ranked_zCandidates,
                                                            const CTemplateCatalog& tplCatalog,
                                                            const TStringList& tplCategoryList,
                                                            Float64 overlapThreshold,
                                                            std::string opt_interp);

    void StoreExtremaResults(std::shared_ptr<COperatorResultStore> dataStore,
                             std::shared_ptr<const ExtremaResult> & ExtremaResult) const ;
    
    std::shared_ptr<COperatorTemplateFittingBase> m_templateFittingOperator;

    std::string m_opt_pdfcombination;
    Float64 m_redshiftSeparation;
    Int64 m_opt_maxCandidate;
    bool m_opt_enableSaveIntermediateTemplateFittingResults=false;

};


}

#endif
