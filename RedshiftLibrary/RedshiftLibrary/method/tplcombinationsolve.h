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
#ifndef _REDSHIFT_METHOD_TPLCOMBINATIONSOLVE_
#define _REDSHIFT_METHOD_TPLCOMBINATIONSOLVE_

#include "RedshiftLibrary/common/datatypes.h"

#include "RedshiftLibrary/method/objectSolve.h"
#include "RedshiftLibrary/operator/pdfz.h"
#include "RedshiftLibrary/operator/tplcombination.h"

#include "RedshiftLibrary/operator/tplCombinationExtremaResult.h"
namespace NSEpic {

class CSpectrum;
class CTemplateCatalog;
class CResultStore;

/**
 * \ingroup Redshift
 */
class CTplCombinationSolve : public CObjectSolve {

public:
  enum class EType {
    raw,
    continuumOnly,
    noContinuum,
    all,
  };

  CTplCombinationSolve();

  std::shared_ptr<CSolveResult> compute() override;

private:
  bool Solve(std::shared_ptr<COperatorResultStore> resultStore,
             const CSpectrum &spc, const CTemplateCatalog &tplCatalog,
             const TFloat64Range &lambdaRange, const TFloat64List &redshifts,
             Float64 overlapThreshold, const std::vector<CMask> &maskList,
             EType spctype = EType::raw, const std::string &opt_interp = "lin",
             bool opt_extinction = false, bool opt_dustFitting = false);

  ChisquareArray
  BuildChisquareArray(std::shared_ptr<COperatorResultStore> store,
                      const std::string &scopeStr) const;
  void StoreExtremaResults(
      std::shared_ptr<COperatorResultStore> resultStore,
      std::shared_ptr<const TplCombinationExtremaResult> &extremaResult) const;
  std::shared_ptr<const TplCombinationExtremaResult> buildExtremaResults(
      std::shared_ptr<const COperatorResultStore> store,
      const std::string &scopeStr, const TCandidateZbyRank &ranked_zCandidates,
      const CSpectrum &spc, const CTemplateCatalog &tplCatalog,
      const TFloat64Range &lambdaRange, Float64 overlapThreshold);
  std::string getResultName(CSpectrum::EType _spctype);
  void checkTemplates(const TTemplateConstRefList &tplList) const;
  COperatorTplcombination m_tplcombinationOperator;

  std::string m_opt_pdfcombination;
  Float64 m_redshiftSeparation;
  Int32 m_opt_maxCandidate;
  std::string m_opt_saveintermediateresults;
};

} // namespace NSEpic

#endif
