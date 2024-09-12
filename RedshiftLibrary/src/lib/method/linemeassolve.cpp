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
#include "RedshiftLibrary/method/linemeassolve.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/processflow/parameterstore.h"

namespace NSEpic {

CLineMeasSolve::CLineMeasSolve() : CObjectSolve("lineMeasSolve") {}

void CLineMeasSolve::GetRedshiftSampling(const CInputContext &inputContext,
                                         TFloat64Range &redshiftRange,
                                         Float64 &redshiftStep) {
  // default is to read from the scoped paramStore
  Float64 rangeCenter = inputContext.GetParameterStore()->GetScopedAt<Float64>(
      "redshiftref", ScopeType::SPECTRUMMODEL);
  Float64 halfRange = inputContext.GetParameterStore()->GetScopedAt<Float64>(
      "lineMeasDzHalf", ScopeType::SPECTRUMMODEL);
  redshiftStep = inputContext.GetParameterStore()->GetScopedAt<Float64>(
      "lineMeasRedshiftStep", ScopeType::SPECTRUMMODEL);

  halfRange = std::ceil(halfRange / redshiftStep) * redshiftStep;

  Float64 half_r = halfRange;
  Float64 half_l = halfRange;
  if (m_redshiftSampling == "log") {
    half_r = (exp(halfRange) - 1.0) * (1. + rangeCenter);
    half_l = (1.0 - exp(-halfRange)) * (1. + rangeCenter);
  }
  redshiftRange = TFloat64Range(rangeCenter - half_l, rangeCenter + half_r);
}

std::shared_ptr<CSolveResult> CLineMeasSolve::compute() {
  auto const &inputContext = Context.GetInputContext();
  auto const &resultStore = Context.GetResultStore();

  Float64 opt_nsigmasupport =
      inputContext->GetParameterStore()->GetScoped<Float64>(
          "lineModel.nSigmaSupport"); // try with 16 (-> parameters.json)

  m_linemodel.Init(m_redshifts, m_redshiftStep, m_redshiftSampling);

  CLineModelSolution bestModelSolution;
  Float64 bestz = NAN;
  {

    bestModelSolution =
        m_linemodel.computeForLineMeas(inputContext, m_redshifts, bestz);
  }

  std::shared_ptr<CModelSpectrumResult> modelspc =
      std::make_shared<CModelSpectrumResult>();
  modelspc->addModel(
      m_linemodel.getFittedModelWithoutcontinuum(bestModelSolution), "");
  std::shared_ptr<const CLineModelSolution> res =
      std::make_shared<CLineModelSolution>(std::move(bestModelSolution));
  resultStore->StoreScopedGlobalResult("linemeas", res);
  resultStore->StoreScopedGlobalResult("linemeas_parameters", res);
  resultStore->StoreScopedGlobalResult("linemeas_model", modelspc);
  return std::make_shared<CLineMeasSolveResult>(CLineMeasSolveResult());
}

} // namespace NSEpic
