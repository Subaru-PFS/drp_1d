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
#include "RedshiftLibrary/processflow/parameterstore.h"
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/method/linemeassolve.h>

namespace NSEpic {

CLineMeasSolve::CLineMeasSolve(TScopeStack &scope, string objectType)
    : CObjectSolve("LineMeasSolve", scope, objectType) {}

void CLineMeasSolve::GetRedshiftSampling(
    std::shared_ptr<const CInputContext> inputContext,
    TFloat64Range &redshiftRange, Float64 &redshiftStep) {
  // default is to read from the scoped paramStore
  Float64 rangeCenter =
      inputContext->GetParameterStore()->GetScoped<Float64>("redshiftref");
  Float64 halfRange =
      inputContext->GetParameterStore()->GetScoped<Float64>("linemeas_dzhalf");

  redshiftRange =
      TFloat64Range(rangeCenter - halfRange, rangeCenter + halfRange);
  redshiftStep = inputContext->GetParameterStore()->GetScoped<Float64>(
      "linemeas_redshiftstep");
}

void CLineMeasSolve::Init() {}

std::shared_ptr<CSolveResult>
CLineMeasSolve::compute(std::shared_ptr<const CInputContext> inputContext,
                        std::shared_ptr<COperatorResultStore> resultStore,
                        TScopeStack &scope) {

  const CSpectrum &spc = *(inputContext->GetSpectrum());
  const CLineCatalog &restlinecatalog =
      *(inputContext->GetLineCatalog(m_objectType, m_name));
  // We keep only emission lines, absorption lines are not handled yet (need to
  // manage continuum appropriately)
  const CLineCatalog::TLineVector restLineList =
      restlinecatalog.GetFilteredList(CLine::nType_Emission, -1);
  Log.LogDebug("restLineList.size() = %d", restLineList.size());

  Float64 opt_nsigmasupport =
      inputContext->GetParameterStore()->GetScoped<Float64>(
          "linemodel.nsigmasupport"); // try with 16 (-> parameters.json)
  const std::string &opt_continuumcomponent =
      "nocontinuum"; // params->GetScoped<std::string>("continuumcomponent");

  bool opt_enableImproveBalmerFit =
      inputContext->GetParameterStore()->GetScoped<bool>(
          "linemodel.improveBalmerFit");

  m_linemodel.Init(spc, m_redshifts, std::move(restLineList), m_categoryList,
                   opt_continuumcomponent, opt_nsigmasupport,
                   opt_enableImproveBalmerFit);

  CLineModelSolution bestModelSolution;
  Float64 bestz = NAN;
  {
    CAutoScope autoscope(scope, "linemodel");

    bestModelSolution =
        m_linemodel.computeForLineMeas(inputContext, m_redshifts, bestz);
  }

  std::shared_ptr<const CModelSpectrumResult> modelspc =
      std::make_shared<const CModelSpectrumResult>(
          m_linemodel.getFittedModelWithoutcontinuum(bestModelSolution));
  std::shared_ptr<const CLineModelSolution> res =
      std::make_shared<CLineModelSolution>(std::move(bestModelSolution));
  resultStore->StoreScopedGlobalResult("linemeas", res);
  resultStore->StoreScopedGlobalResult("linemeas_parameters", res);
  resultStore->StoreScopedGlobalResult("linemeas_model", modelspc);
  return std::make_shared<CLineMeasSolveResult>(CLineMeasSolveResult());
}

} // namespace NSEpic
