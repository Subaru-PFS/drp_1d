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
                           std::shared_ptr<CContinuumManager> continuumManager,
                           Int32 polyOrder)
    : CAbstractFitter(elements, inputSpectrum, lambdaRange, spectrumModel),
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

  m_continuumManager->reinterpolateContinuum();
  TInt32List validEltsIdx = m_Elements.GetModelValidElementsIndexes();
  TFloat64List ampsfitted;
  TFloat64List errorsfitted;
  Float64 chi2_cl = INFINITY;

  fitAmplitudesLinesAndContinuumLinSolve(
      validEltsIdx, m_spectralAxis, m_model->getSpcFluxAxis(),
      m_model->getContinuumFluxAxis(), ampsfitted, errorsfitted, chi2_cl,
      redshift, m_fitc_polyOrder);

  m_continuumManager->setFitContinuumFromFittedAmps(ampsfitted, validEltsIdx);
  m_model->initModelWithContinuum();
  m_model->substractContToFlux();
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
    const CSpectrumFluxAxis &fluxAxis,
    const CSpectrumFluxAxis &continuumfluxAxis, TFloat64List &ampsfitted,
    TFloat64List &errorsfitted, Float64 &chisquare, Float64 redshift,
    Int32 polyOrder) {
  // boost::chrono::thread_clock::time_point start_prep =
  // boost::chrono::thread_clock::now();

  Int32 idx = 0;

  Int32 nddl =
      EltsIdx.size() + 1; // number of param to be fitted=nlines+continuum
  if (nddl < 2) {
    return -1;
  }

  if (polyOrder >= 0) {
    nddl += polyOrder + 1;
  }

  TInt32List xInds;
  for (Int32 i = 0; i < spectralAxis.GetSamplesCount(); i++) {
    if (spectralAxis[i] >= m_lambdaRange.GetBegin() &&
        spectralAxis[i] <= m_lambdaRange.GetEnd()) {
      xInds.push_back(i);
    }
  }
  if (xInds.size() < 1) {
    return -1;
  }

  for (Int32 iddl = 0; iddl < EltsIdx.size(); iddl++) {
    m_Elements.SetElementAmplitude(EltsIdx[iddl], 1.0, 0.0);
  }

  const Float64 *spectral = spectralAxis.GetSamples();
  const Float64 *flux = fluxAxis.GetSamples();
  const auto &ErrorNoContinuum = m_inputSpc.GetFluxAxis().GetError();

  // Linear fit
  int i, n;
  Float64 fval;
  double chisq;
  gsl_matrix *X, *cov;
  gsl_vector *y, *w, *c;

  n = xInds.size();
  if (n < nddl) {
    ampsfitted.resize(nddl);
    errorsfitted.resize(nddl);
    for (Int32 iddl = 0; iddl < nddl; iddl++) {
      ampsfitted[iddl] = 0.0;
      errorsfitted[iddl] = 1e12; // some high number
    }
    return -1;
  }

  X = gsl_matrix_alloc(n, nddl);
  y = gsl_vector_alloc(n);
  w = gsl_vector_alloc(n);
  c = gsl_vector_alloc(nddl);
  cov = gsl_matrix_alloc(nddl, nddl);

  // Normalize
  Float64 maxabsval = DBL_MIN;
  for (i = 0; i < n; i++) {
    idx = xInds[i];
    if (maxabsval < std::abs(flux[idx])) {
      maxabsval = std::abs(flux[idx]);
    }
  }
  Float64 normFactor = 1.0 / maxabsval;
  Log.LogDetail("normFactor = '%.3e'\n", normFactor);

  // Prepare the fit data
  for (i = 0; i < n; i++) {
    double xi, yi, ei, ci;
    idx = xInds[i];
    xi = spectral[idx];
    yi = flux[idx] * normFactor;
    ei = ErrorNoContinuum[idx] * normFactor;
    ci = continuumfluxAxis[idx];

    for (Int32 iddl = 0; iddl < EltsIdx.size(); iddl++) {
      fval = m_Elements[EltsIdx[iddl]]->getModelAtLambda(xi, redshift, ci);
      gsl_matrix_set(X, i, iddl, fval);

      Log.LogDebug("fval = '%.3e'", fval);
    }

    gsl_matrix_set(X, i, EltsIdx.size(), ci);

    if (polyOrder >= 0) {
      for (Int32 kCoeff = 0; kCoeff < polyOrder + 1; kCoeff++) {
        Float64 vect = 1.0;
        for (Int32 kt = 0; kt < kCoeff; kt++) {
          vect *= xi;
        }
        gsl_matrix_set(X, i, EltsIdx.size() + 1 + kCoeff, vect);
      }
    }

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
  Log.LogDebug("[");
  Log.LogDebug("  %+.5e, %+.5e", COV(0, 0), COV(0, 1));
  Log.LogDebug("  %+.5e, %+.5e", COV(1, 0), COV(1, 1));

  //        Log.LogDebug("[ %+.5e, %+.5e, %+.5e  \n", COV(0,0), COV(0,1),
  //        COV(0,2)); Log.LogDebug("  %+.5e, %+.5e, %+.5e  \n", COV(1,0),
  //        COV(1,1), COV(1,2)); Log.LogDebug("  %+.5e, %+.5e, %+.5e ]\n",
  //        COV(2,0), COV(2,1), COV(2,2));

  Log.LogDebug("]");
  Log.LogDebug("# chisq = %g", chisq);
  Log.LogDebug("# chisq/n = %g", chisq / n);

  for (Int32 iddl = 0; iddl < nddl; iddl++) {
    Float64 a = gsl_vector_get(c, iddl) / normFactor;
    Log.LogDetail("# Found amplitude %d: %+.5e", iddl, a);
  }

  Int32 sameSign = 1;
  Float64 a0 = gsl_vector_get(c, 0) / normFactor;
  for (Int32 iddl = 1; iddl < EltsIdx.size(); iddl++) {
    Float64 a = gsl_vector_get(c, iddl) / normFactor;
    Float64 product = a0 * a;
    if (product < 0) {
      sameSign = 0;
    }
  }

  Log.LogDetail("# Found n=%d amplitudes with sameSign=%d", EltsIdx.size(),
                sameSign);

  ampsfitted.resize(nddl);
  errorsfitted.resize(nddl);
  for (Int32 iddl = 0; iddl < nddl; iddl++) {
    Float64 a = gsl_vector_get(c, iddl) / normFactor;
    Float64 cova = gsl_matrix_get(cov, iddl, iddl);
    Float64 sigma = sqrt(cova) / normFactor;
    if (iddl < EltsIdx.size()) {
      m_Elements.SetElementAmplitude(EltsIdx[iddl], a, sigma);
    }
    ampsfitted[iddl] = (a);
    errorsfitted[iddl] = (sigma);
  }

  if (polyOrder >= 0) {
    for (Int32 kCoeff = 0; kCoeff < polyOrder + 1; kCoeff++) {
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
