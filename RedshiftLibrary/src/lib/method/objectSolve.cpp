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

#include "RedshiftLibrary/method/objectSolve.h"
#include "RedshiftLibrary/common/zgridparam.h"

using namespace NSEpic;

void CObjectSolve::InitRanges(const CInputContext &inputContext) {
  m_lambdaRange = TFloat64Range(*inputContext.getLambdaRange()); // non-clamped

  // m_zLogSampling could be overwritten if fftprocessing is activated
  //  Warning for LineMeas :  we consider linemeas use the same redshiftsampling
  //  as the objectMethod if linemeas is called in standalone, redshiftsampling
  //  must be present in parameters
  // TODO change here
  m_zLogSampling = inputContext.GetParameterStore()->GetScopedAt<std::string>(
                       "redshiftSampling", ScopeType::SPECTRUMMODEL) == "log";

  TFloat64Range redshiftRange;

  GetRedshiftSampling(inputContext, redshiftRange, m_redshiftStep);

  Log.LogInfo(Formatter() << "Init redshift range with " << redshiftRange
                          << " and " << m_redshiftStep);
  createRedshiftGrid(inputContext, redshiftRange);
}

// common for all except linemodelsolve
void CObjectSolve::createRedshiftGrid(const CInputContext &inputContext,
                                      const TFloat64Range &redshiftRange) {
  CZGridParam zp(redshiftRange, m_redshiftStep);
  m_redshifts = zp.getZGrid(m_zLogSampling);
}

void CObjectSolve::GetRedshiftSampling(const CInputContext &inputContext,
                                       TFloat64Range &redshiftRange,
                                       Float64 &redshiftStep) {
  auto searchLogRebin = inputContext.m_logRebin.find(m_category);
  if (searchLogRebin != inputContext.m_logRebin.end()) {
    redshiftRange = searchLogRebin->second.zrange;
    redshiftStep = inputContext.getLogGridStep();
    if (!m_zLogSampling)
      THROWG(ErrorCode::BAD_PARAMETER_VALUE,
             Formatter() << "CSolve::" << __func__
                         << ": redshiftsampling param should be set to log "
                            "since FFTprocessing is used");

  } else {
    // default is to read from the scoped paramStore
    redshiftRange =
        inputContext.GetParameterStore()->GetScopedAt<TFloat64Range>(
            "redshiftRange", ScopeType::SPECTRUMMODEL);
    redshiftStep = inputContext.GetParameterStore()->GetScopedAt<Float64>(
        "redshiftStep", ScopeType::SPECTRUMMODEL);
  }
  return;
}
