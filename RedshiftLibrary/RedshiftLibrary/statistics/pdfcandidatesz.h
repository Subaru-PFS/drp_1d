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
#ifndef _REDSHIFT_STATISTICS_PDFCANDIDATESZ_
#define _REDSHIFT_STATISTICS_PDFCANDIDATESZ_

#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/operator/operator.h"
#include "RedshiftLibrary/processflow/result.h"

namespace Statistics_pdfcandidatesz { // boost_test_suite
// all boost_auto_test_case that use private method
class SortByValSumProbaInt_test;
class getCandidateSumTrapez_test;
class getCandidateRobustGaussFit_test;
class getCandidateGaussFit_test;
} // namespace Statistics_pdfcandidatesz
namespace NSEpic {

#include "RedshiftLibrary/statistics/pdfcandidatesz.i"

typedef std::map<std::string, std::shared_ptr<TCandidateZ>> TCandidateZbyID;
typedef std::vector<std::pair<std::string, std::shared_ptr<TCandidateZ>>>
    TCandidateZbyRank;
typedef std::map<std::string, TFloat64Range> TCandidateZRangebyID;
typedef std::pair<std::string, std::shared_ptr<TCandidateZ>>
    pair_Id_TCandidateZ;
template <typename T> class CPdfCandidateszResult;

class CPdfCandidatesZ : public COperator {

public:
  CPdfCandidatesZ(const TCandidateZbyID &candidates);
  CPdfCandidatesZ(const TFloat64List &redshifts);

  std::shared_ptr<CPdfCandidateszResult<TCandidateZ>>
  Compute(TRedshiftList const &PdfRedshifts, TFloat64List const &PdfProbaLog);
  void computeCandidatesDeltaz(const TRedshiftList &PdfRedshifts,
                               const TFloat64List &PdfProbaLog);
  TStringList SetIntegrationRanges(const TFloat64Range PdfZRange,
                                   TCandidateZRangebyID &ranges);
  static Float64 integrationWindowHalfWidth(Float64 DeltaZ);

  Int32 m_optMethod = 0; // 0: direct integration, 1:gaussian fit
  Float64 m_dzDefault =
      1e-3; // default value in case deltaz couldnt be calculted, should be
            // instrument dependant (parameter ?)

  TCandidateZbyID m_candidates;

private:
  friend class Statistics_pdfcandidatesz::SortByValSumProbaInt_test;
  friend class Statistics_pdfcandidatesz::getCandidateSumTrapez_test;
  friend class Statistics_pdfcandidatesz::getCandidateRobustGaussFit_test;
  friend class Statistics_pdfcandidatesz::getCandidateGaussFit_test;

  void getCandidateSumTrapez(const TRedshiftList &redshifts,
                             const TFloat64List &valprobalog,
                             const TFloat64Range &zrange,
                             std::shared_ptr<TCandidateZ> &candidate) const;

  bool
  getCandidateRobustGaussFit(const TRedshiftList &redshifts,
                             const TFloat64List &valprobalog,
                             const TFloat64Range &zrange,
                             std::shared_ptr<TCandidateZ> &candidate) const;

  void getCandidateGaussFit(const TRedshiftList &redshifts,
                            const TFloat64List &valprobalog,
                            const TFloat64Range &zrange,
                            std::shared_ptr<TCandidateZ> &candidate) const;

  void SortByValSumProbaInt(TCandidateZbyRank &ranked_candidates) const;
};

} // namespace NSEpic

#endif
