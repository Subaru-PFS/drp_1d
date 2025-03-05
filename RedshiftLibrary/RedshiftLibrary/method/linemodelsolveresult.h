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
#ifndef _REDSHIFT_METHOD_LINEMODELSOLVERESULT_
#define _REDSHIFT_METHOD_LINEMODELSOLVERESULT_

#include <memory>
#include <vector>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/line/catalog.h"
#include "RedshiftLibrary/linemodel/linemodelextremaresult.h"
#include "RedshiftLibrary/method/solveresult.h"

namespace NSEpic {

/**
 * \ingroup Redshift
 */
class CLineModelSolveResult : public CPdfSolveResult {

public:
  CLineModelSolveResult(
      const std::shared_ptr<const TCandidateZ> &BestExtremumResult,
      const std::string &opt_pdfcombination, Float64 evidence,
      Float64 continuumEvidence = NAN)
      : CPdfSolveResult("CLineModelSolveResult", BestExtremumResult,
                        opt_pdfcombination, evidence),
        m_continuumEvidence(continuumEvidence) {
    if (!std::isnan(m_continuumEvidence)) {
      setSwitchedToFromSpectrum(true);
    }
  };
  Float64 getContinuumEvidence() const override { return m_continuumEvidence; };
  bool getSwitchedToFromSpectrum() const override {
    return m_switchedToFromSpectrum;
  };
  void setSwitchedToFromSpectrum(const bool switched) {
    m_switchedToFromSpectrum = switched;
  };

private:
  std::string tplratioName = undefStr;
  std::string tplContinuumName = undefStr;
  Float64 sigma = NAN;
  Float64 snrHa = NAN;
  Float64 lfHa = NAN;
  Float64 snrOII = NAN;
  Float64 lfOII = NAN;
  Float64 LyaWidthCoeff = NAN;
  Float64 LyaAlpha = NAN;
  Float64 LyaDelta = NAN;
  Int32 LyaIgm = undefIdx;
  bool m_switchedToFromSpectrum = false;
  Float64 m_continuumEvidence =
      NAN; // For cases where switches in fromSpectrum, keeps evidence of
           // continuum for usage in classification
  Float64 minContinuumReducedChi2 = 0.;
};

} // namespace NSEpic

#endif
