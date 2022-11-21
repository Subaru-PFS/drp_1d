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
#include "RedshiftLibrary/linemodel/svdlcfitter.h"
#include "RedshiftLibrary/line/linetags.h"
#include "RedshiftLibrary/processflow/context.h"
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlin.h>
using namespace NSEpic;
using namespace std;

CSvdlcFitter::CSvdlcFitter(CLineModelElementList &elements,
                           std::shared_ptr<const CSpectrum> inputSpectrum,
                           std::shared_ptr<const TLambdaRange> lambdaRange,
                           std::shared_ptr<CSpectrumModel> spectrumModel,
                           const CLineCatalog::TLineVector &restLineList,
                           std::shared_ptr<CContinuumManager> continuumManager,
                           Int32 polyOrder)
    : CAbstractFitter(elements, inputSpectrum, lambdaRange, spectrumModel,
                      restLineList),
      m_fitc_polyOrder(polyOrder), m_continuumManager(continuumManager),
      m_spectralAxis(inputSpectrum->GetSpectralAxis())

{}

// fit the amplitude of all elements AND continuum amplitude together
// with linear solver: gsl_multifit_wlinear

void CSvdlcFitter::fit(Float64 redshift) {

  // 1. fit only the current continuum
  // prepare continuum on the observed grid
  // Log.LogDebug("    model: fitting svdlc, with continuum-tpl=%s",
  //	       m_fitContinuum_tplName.c_str());

  // re-interpolate the continuum on the grid

  m_continuumManager->reinterpolateContinuumResetAmp();
  TInt32List validEltsIdx = m_Elements.GetModelValidElementsIndexes();
  TFloat64List ampsfitted;
  TFloat64List errorsfitted;
  Float64 chi2_cl = INFINITY;

  fitAmplitudesLinesAndContinuumLinSolve(validEltsIdx, m_spectralAxis,
                                         ampsfitted, errorsfitted, chi2_cl,
                                         redshift);

  m_continuumManager->setFitContinuumFromFittedAmps(ampsfitted, validEltsIdx);
  m_model->initModelWithContinuum();
}

/**
 * @brief CLineModelFitting::fitAmplitudesLinesAndContinuumLinSolve
 * @param EltsIdx:  elements to be fitted
 * @param lambdaRange
 * @param spectralAxis
 * @param fluxAxis: must contain the observed spectrum
 * @param continuumfluxAxis: must contain the continuum to be amp-fitted
 * @param ampsfitted: output
 * @param errorsfitted: output
 * @param polyOrder: order of the polynom to be fitted along with the
 * continuum
 * (-1 = disabled)
 * @return
 *
 * WARNING: not sure about fitting abs. lines with this method...
 */
Int32 CSvdlcFitter::fitAmplitudesLinesAndContinuumLinSolve(
    const TInt32List &EltsIdx, const CSpectrumSpectralAxis &spectralAxis,
    TFloat64List &ampsfitted, TFloat64List &errorsfitted, Float64 &chisquare,
    Float64 redshift) {
  // boost::chrono::thread_clock::time_point start_prep =
  // boost::chrono::thread_clock::now();
  const CSpectrumFluxAxis &fluxAxis = m_model->getSpcFluxAxis();
  const CSpectrumFluxAxis &continuumfluxAxis = m_model->getContinuumFluxAxis();
  const auto &ErrorNoContinuum = m_inputSpc.GetFluxAxis().GetError();

  if (EltsIdx.size() < 1)
    return -1;
  Int32 nddl = EltsIdx.size() + 1 +
               std ::max(m_fitc_polyOrder + 1,
                         0); // number of param to be fitted=nlines+continuum

  Int32 imin = -1, imax = -1;
  bool b = m_lambdaRange.getClosedIntervalIndices(
      spectralAxis.GetSamplesVector(), imin, imax);
  if (!b)
    return -1;

  ampsfitted.resize(nddl, 0.0);
  errorsfitted.resize(nddl, 1e12);

  Int32 n = imax - imin + 1;
  if (n < nddl)
    return -1;

  for (Int32 eltIdx : EltsIdx)
    m_Elements.SetElementAmplitude(eltIdx, 1.0, 0.0);

  // Linear fit
  double chisq;
  gsl_matrix *X = gsl_matrix_alloc(n, nddl);
  gsl_matrix *cov = gsl_matrix_alloc(nddl, nddl);
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_vector *w = gsl_vector_alloc(n);
  gsl_vector *c = gsl_vector_alloc(nddl);

  // Normalize
  Float64 normFactor = 1.0 / fluxAxis.computeMaxAbsValue(imin, imax);

  // Prepare the fit data
  fillMatrix(imin, imax, redshift, EltsIdx, spectralAxis, continuumfluxAxis, X);

  for (Int32 i = 0, idx = imin; idx <= imax; ++i, ++idx) {
    Float64 yi = fluxAxis[idx] * normFactor;
    Float64 ei = ErrorNoContinuum[idx] * normFactor;
    Float64 ci = continuumfluxAxis[idx];
    gsl_vector_set(y, i, yi);
    gsl_vector_set(w, i, 1.0 / (ei * ei));
  }

  //
  // boost::chrono::thread_clock::time_point stop_prep =
  // boost::chrono::thread_clock::now();
  // Float64 duration_prep =
  // boost::chrono::duration_cast<boost::chrono::microseconds>(stop_prep -
  // start_prep).count();
  // boost::chrono::thread_clock::time_point start_fit =
  // boost::chrono::thread_clock::now();

  gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n, nddl);
  gsl_multifit_wlinear(X, w, y, c, cov, &chisq, work);
  gsl_multifit_linear_free(work);

  //
  // boost::chrono::thread_clock::time_point stop_fit =
  // boost::chrono::thread_clock::now(); Float64 duration_fit =
  // boost::chrono::duration_cast<boost::chrono::microseconds>(stop_fit -
  // start_fit).count(); Log.LogInfo("LineModel linear fit: prep = %.3f - fit
  // =
  // %.3f", duration_prep, duration_fit);

#define C(i) (gsl_vector_get(c, (i)))
#define COV(i, j) (gsl_matrix_get(cov, (i), (j)))
  Log.LogDebug("# best fit: Y = %g X1 + %g X2 ...", C(0), C(1));
  Log.LogDebug("# covariance matrix:");
  Log.LogDebug("[ %+.5e, %+.5e", COV(0, 0), COV(0, 1));
  Log.LogDebug("  %+.5e, %+.5e ]", COV(1, 0), COV(1, 1));
  Log.LogDebug("# chisq = %g, chisq/n = %g", chisq, chisq / n);

  for (Int32 iddl = 0; iddl < nddl; iddl++) {
    Float64 a = gsl_vector_get(c, iddl) / normFactor;
    Log.LogDetail("# Found amplitude %d: %+.5e", iddl, a);
  }

  Int32 sameSign = 1;
  Float64 a0 = gsl_vector_get(c, 0) / normFactor;
  for (Int32 iddl = 1; iddl < EltsIdx.size(); iddl++) {
    Float64 a = gsl_vector_get(c, iddl) / normFactor;
    sameSign &= std::signbit(a0 * a);
  }

  Log.LogDetail("# Found n=%d amplitudes with sameSign=%d", EltsIdx.size(),
                sameSign);

  for (Int32 iddl = 0; iddl < nddl; iddl++) {
    Float64 cova = gsl_matrix_get(cov, iddl, iddl);
    errorsfitted[iddl] = sqrt(cova) / normFactor;
    ampsfitted[iddl] = gsl_vector_get(c, iddl) / normFactor;
    if (iddl < EltsIdx.size())
      m_Elements.SetElementAmplitude(EltsIdx[iddl], ampsfitted[iddl],
                                     errorsfitted[iddl]);
  }

  if (m_fitc_polyOrder >= 0) {
    for (Int32 kCoeff = 0; kCoeff < m_fitc_polyOrder + 1; kCoeff++) {
      Float64 p = gsl_vector_get(c, EltsIdx.size() + 1 + kCoeff) / normFactor;
      Log.LogDetail("# Found p%d poly amplitude = %+.5e", kCoeff, p);
    }
  }

  Log.LogDetail("# Returning (L+C) n=%d amplitudes", ampsfitted.size());

  gsl_matrix_free(X);
  gsl_vector_free(y);
  gsl_vector_free(w);
  gsl_vector_free(c);
  gsl_matrix_free(cov);

  chisquare = chisq;
  return sameSign;
}

void CSvdlcFitter::fillMatrix(Int32 imin, Int32 imax, Float64 redshift,
                              const TInt32List &EltsIdx,
                              const CSpectrumSpectralAxis &spectralAxis,
                              const CSpectrumFluxAxis &continuumfluxAxis,
                              gsl_matrix *X) const {
  Int32 n = imax - imin + 1;

  for (Int32 i = 0, idx = imin; idx <= imax; ++i, ++idx) {
    Float64 xi = spectralAxis[idx];
    Float64 ci = continuumfluxAxis[idx];
    for (Int32 iddl = 0, e = EltsIdx.size(); iddl < e; iddl++) {
      Float64 fval =
          m_Elements[EltsIdx[iddl]]->getModelAtLambda(xi, redshift, ci);
      gsl_matrix_set(X, i, iddl, fval);
    }

    gsl_matrix_set(X, i, EltsIdx.size(), ci);

    for (Int32 kCoeff = 0; kCoeff < m_fitc_polyOrder + 1; kCoeff++) {
      Float64 vect = std::pow(xi, kCoeff);
      gsl_matrix_set(X, i, EltsIdx.size() + 1 + kCoeff, vect);
    }
  }
  return;
}