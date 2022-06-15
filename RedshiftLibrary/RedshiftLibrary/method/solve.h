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
#ifndef _SOLVE_H_
#define _SOLVE_H_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/method/solveresult.h"
#include "RedshiftLibrary/processflow/autoscope.h"
#include "RedshiftLibrary/processflow/context.h"
#include "RedshiftLibrary/processflow/resultstore.h"

namespace NSEpic {

class CSolve {
protected:
  virtual std::shared_ptr<CSolveResult>
  compute(std::shared_ptr<const CInputContext> inputContext,
          std::shared_ptr<COperatorResultStore> resultStore,
          TScopeStack &scope) = 0;

  virtual void InitRanges(std::shared_ptr<const CInputContext> inputContext);
  virtual void GetRedshiftSampling(std::shared_ptr<const CInputContext>,
                                   TFloat64Range &redshiftRange,
                                   Float64 &redshiftStep);

  virtual void
  saveToResultStore(std::shared_ptr<CSolveResult>,
                    std::shared_ptr<COperatorResultStore> resultStore) const;

  // this method should implement at least populateParameters
  virtual void checkOrInit() {
  } //=0; // here we retrieve parameters for parameterStore to put them directly
    // in local variables or into operators, lineCatalog and/or tplCatalog can
    // also be checked

  const TStringList m_categoryList;
  TFloat64Range m_lambdaRange;
  TFloat64List m_redshifts;
  CAutoScope m_objectTypeScope;
  const std::string m_name;
  const std::string m_objectType;
  std::string m_redshiftSampling;

public:
  CSolve(std::string name, TScopeStack &scope, std::string objectType);
  virtual ~CSolve() = default;
  CSolve(CSolve const &other) = default;
  CSolve &operator=(CSolve const &other) = default;
  CSolve(CSolve &&other) = default;
  CSolve &operator=(CSolve &&other) = default;

  void Compute();
};

} // namespace NSEpic

#endif
