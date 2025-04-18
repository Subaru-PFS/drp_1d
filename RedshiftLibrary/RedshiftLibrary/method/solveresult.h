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
#ifndef _SOLVE_RESULT_
#define _SOLVE_RESULT_

#include <vector>

#include "RedshiftLibrary/processflow/result.h"

namespace NSEpic {

class TCandidateZ;

class CSolveResult : public COperatorResult {
public:
  CSolveResult(const std::string &type) : COperatorResult(type){};
  virtual ~CSolveResult() = 0;
  CSolveResult(const CSolveResult &) = default;
  CSolveResult(CSolveResult &&) = default;
  CSolveResult &operator=(const CSolveResult &) = default;
  CSolveResult &operator=(CSolveResult &&) = default;
};

class CPdfSolveResult : public CSolveResult {

public:
  CPdfSolveResult(const std::string &type,
                  const std::shared_ptr<const TCandidateZ> &ExtremaResult,
                  const std::string &opt_pdfcombination, Float64 evidence);
  virtual ~CPdfSolveResult() = default;
  CPdfSolveResult(CPdfSolveResult const &other) = default;
  CPdfSolveResult &operator=(CPdfSolveResult const &other) = default;
  CPdfSolveResult(CPdfSolveResult &&other) = default;
  CPdfSolveResult &operator=(CPdfSolveResult &&other) = default;
  Int32 m_bestRedshiftMethod = 2; // 0:best chi2 or proba, 2: best marg proba

  inline Float64 getMerit() const { return m_merit; }
  inline Float64 getEvidence() const { return m_evidence; }
  virtual Float64 getContinuumEvidence() const { return m_evidence; };
  virtual bool getSwitchedToFromSpectrum() const { return false; };

protected:
  Float64 m_redshift;
  Float64 m_merit;
  Float64 m_evidence = -INFINITY;
};
} // namespace NSEpic

#endif
