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
#ifndef _REDSHIFT_METHOD_LINEMODELSOLVE_
#define _REDSHIFT_METHOD_LINEMODELSOLVE_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/method/twopasssolve.h"
#include "RedshiftLibrary/operator/linemodel.h"

#include "RedshiftLibrary/processflow/inputcontext.h"
#include "RedshiftLibrary/processflow/resultstore.h"

namespace lineModelSolve_test {
class continuumChi2CorrectlySet_test;
}

namespace NSEpic {

class CSpectrum;
class CTemplateCatalog;

/**
 * \ingroup Redshift
 */
class CLineModelSolve : public CTwoPassSolve {
public:
  CLineModelSolve();

  bool
  PopulateParameters(std::shared_ptr<const CParameterStore> parameterStore);

  std::shared_ptr<CSolveResult> compute() override;

  void Solve();

private:
  friend class lineModelSolve_test::continuumChi2CorrectlySet_test;
  COperatorPdfz initializePdfz(Int32 maxPeakPerWindow, Int32 peakSeparation,
                               Int32 cutThreshold, Int32 extremaCount) const;
  ChisquareArray
  BuildContinuumChisquareArray(const Float64 cstLog,
                               const TZGridListParams &spZgridParams) const;

  ChisquareArray
  BuildChisquareArray(const std::shared_ptr<const CLineModelResult> &result,
                      const TZGridListParams &zgridParams = {},
                      const TCandidateZbyRank &parentZCand = {}) const;

  void GetZpriorsOptions(bool &zPriorStrongLinePresence,
                         bool &zPriorHaStrongestLine, bool &zPriorNLineSNR,
                         Float64 &opt_nlines_snr_penalization_factor,
                         bool &zPriorEuclidNHa) const;

  TFloat64List
  BuildZpriors(const std::shared_ptr<const CLineModelResult> &result,
               Int32 kTplRatio = -1) const;

  void storeExtremaResults(
      std::shared_ptr<COperatorResultStore> dataStore,
      std::shared_ptr<const LineModelExtremaResult> ExtremaResult) const;

  void fillChisquareArrayForTplRatio(
      const std::shared_ptr<const CLineModelResult> &result,
      ChisquareArray &chisquarearray) const;

  void initSkipSecondPass() override;
  void initTwoPassZStepFactor() override;

  COperatorLineModel m_linemodel;

  std::string m_opt_lineratiotype;
  std::string m_opt_continuumreest;
  TContinuumComponent m_opt_continuumcomponent;

  std::string m_opt_pdfcombination;
  Int32 m_opt_extremacount;
  Int32 m_opt_extremacountB;
  Int32 m_opt_maxCandidate;

  Float64 m_opt_stronglinesprior;
  Float64 m_opt_haPrior;
  Float64 m_opt_euclidNHaEmittersPriorStrength;

  Float64 m_opt_secondpass_halfwindowsize;
  Float64 m_opt_candidatesLogprobaCutThreshold;

  bool m_useloglambdasampling;
};

} // namespace NSEpic

#endif
