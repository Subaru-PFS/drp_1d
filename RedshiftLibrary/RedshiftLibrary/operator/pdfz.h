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
#ifndef _REDSHIFT_OPERATOR_PDFZ_
#define _REDSHIFT_OPERATOR_PDFZ_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/operator/operator.h"

#include "RedshiftLibrary/operator/pdfMargZLogResult.h"
#include "RedshiftLibrary/statistics/pdfcandidatesz.h"
#include "RedshiftLibrary/statistics/pdfcandidateszresult.h"

#include <string>

namespace NSEpic {

struct ChisquareArray {
  TFloat64List redshifts;
  std::vector<TFloat64List> chisquares;
  std::vector<TFloat64List> zpriors;
  TFloat64List modelpriors;
  Float64 cstLog = 0.;
};

/**
 * \ingroup Redshift
 * Pdfz
 */
class COperatorPdfz : public COperator {

public:
  COperatorPdfz(
      const std::string &opt_combine,
      Float64 peakSeparation = 0.0, // no minimal seaparation
      Float64 meritcut = 0.0,       // no cut
      Int32 maxCandidate = 10,      // max number of candidate at the end
      const std::string &Id_prefix = "EXT", bool allow_extrema_at_border = true,
      Int32 maxPeakCount_per_window =
          0, // <=0 will be set to maxCandidate (default to one window)
      const std::vector<TFloat64List> &candidatesRedshifts =
          std::vector<TFloat64List>(1),
      const TCandidateZbyRank &parentCand = TCandidateZbyRank(1));

  std::shared_ptr<CPdfCandidateszResult<TCandidateZ>>
  Compute(const ChisquareArray &chisquares, bool integ = true);

  void CombinePDF(const ChisquareArray &chisquares);

  bool checkPdfSum() const;

  std::shared_ptr<CPdfMargZLogResult> m_postmargZResult;

  // static member function to do calculation on pdfs
  ///////////////////////////////////////////////////

  static Int32 getIndex(const TFloat64List &redshifts, Float64 z);

  static void ComputePdf(const TFloat64List &merits,
                         const TFloat64List &redshifts, const Float64 cstLog,
                         const TFloat64List &zPrior, TFloat64List &logPdf,
                         Float64 &logEvidence);

  static Float64 getSumTrapez(const TRedshiftList &redshifts,
                              const TFloat64List &valprobalog);

  static Float64 logSumExpTrick(const TFloat64List &valproba,
                                const TFloat64List &redshifts);

  const std::shared_ptr<const CPdfMargZLogResult>
  compressFirstpassPDF(Int32 ratio);

private:
  void ComputeEvidenceAll(const TFloat64List &LogEvidencesWPriorM,
                          Float64 &MaxiLogEvidence);

  void ComputeAllPdfs(const ChisquareArray &chisquarearray,
                      std::vector<TFloat64List> &logProbaList,
                      TFloat64List &LogEvidencesWPriorM,
                      TFloat64List &logPriorModel, Float64 &MaxiLogEvidence);

  void Marginalize(const ChisquareArray &chisquarearray);
  void BestProba(const ChisquareArray &chisquarearray);
  void BestChi2(const ChisquareArray &chisquarearray);

  void isPdfValid() const;

  TCandidateZbyID searchMaxPDFcandidates() const;

  const std::string m_opt_combine;
  TCandidateZRangebyID m_candidatesZRanges;
  const Int32 m_maxPeakCount_per_window;
  const Int32 m_maxCandidate;
  const Float64 m_peakSeparation;
  const bool m_allow_extrema_at_border;
  const Float64 m_meritcut;
  const std::string m_Id_prefix; // =  "EXT"; // for "extrema"

  const TCandidateZbyID m_parentCandidates;
};

} // namespace NSEpic

#endif
