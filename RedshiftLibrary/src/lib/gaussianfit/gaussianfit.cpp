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
#include "RedshiftLibrary/gaussianfit/gaussianfit.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_vector.h>
#include <math.h>

using namespace NSEpic;

CGaussianFit::CGaussianFit()
    : m_AbsTol(0.0), m_Amplitude(0.0), m_AmplitudeErr(0.0), m_C(0.0),
      m_CErr(0.0), m_Mu(0.0), m_MuErr(0.0), m_RelTol(1e-12), m_coeff0(0.0),
      m_PolyOrder(2) {}

/**
 * Sets amplitude to m_Amplitude, position to m_Mu and width to m_C.
 */
void CGaussianFit::GetResults(Float64 &amplitude, Float64 &position,
                              Float64 &width) const {
  amplitude = m_Amplitude;
  position = m_Mu;
  width = m_C;
}

/**
 * Sets coeff0 to m_coeff0.
 */
void CGaussianFit::GetResultsPolyCoeff0(Float64 &coeff0) const {
  coeff0 = m_coeff0;
}

/**
 * Sets amplitude to m_AmplitudeErr, position to m_MuErr and width to m_CErr.
 */
void CGaussianFit::GetResultsError(Float64 &amplitude, Float64 &position,
                                   Float64 &width) const {
  amplitude = m_AmplitudeErr;
  position = m_MuErr;
  width = m_CErr;
}
/**
 * Calculate a first approximate value for the gaussian fit, and set the
 * referenced parameters to these values.
 */
void CGaussianFit::ComputeFirstGuess(const CSpectrum &spectrum,
                                     const TInt32Range &studyRange,
                                     Int32 polyOrder, Float64 &peakValue,
                                     Float64 &peakPos, Float64 &gaussAmp) {
  Int32 n = studyRange.GetLength();
  const Float64 *x =
      spectrum.GetSpectralAxis().GetSamples() + studyRange.GetBegin();
  const Float64 *y =
      spectrum.GetFluxAxis().GetSamples() + studyRange.GetBegin();

  Int32 i;
  // Copy flux axis to v
  Float64 *v = (Float64 *)calloc(n, sizeof(Float64));
  for (i = 0; i < n; i++) {
    v[i] = y[i];
  }

  // Sort v
  gsl_sort(v, 1, n);
  // Extract median
  Float64 y_median = gsl_stats_median_from_sorted_data(v, 1, n);

  for (i = 0; i < n; i++) {
    v[i] = y[i] - y_median;
  }

  Float64 max = gsl_stats_max(v, 1, n);
  Float64 min = gsl_stats_min(v, 1, n);

  // if (fabs(max) > fabs(min)) // TODO, WARNING, Gaussian fit forced to
  // positive amplitudes, aschmitt, 20150827
  if (max > min) {
    // Peak value
    peakValue = max;
    // Peak position: mu
    peakPos = x[(Int32)gsl_stats_max_index(v, 1, n)];
  } else {
    // Peak value
    peakValue = min;
    // Peak position: mu
    peakPos = x[(Int32)gsl_stats_min_index(v, 1, n)];
  }

  // Gaussian amplitude
  Float64 std = 0;
  for (i = 0; i < n; i++) {
    std += (y[i] - std) * (y[i] - std);
  }
  std /= (n - 1);

  Int32 count = 0;
  for (i = 0; i < n; i++) {
    if (y[i] > y_median + 3. * std) {
      count++;
    }
  }

  if (count > 2) {
    gaussAmp = (x[count] - x[0]) / 6.;
  } else {
    gaussAmp = (x[n - 1] - x[1]) / 6.;
  }

  free(v);
}

/**
 *
 */
CGaussianFit::EStatus CGaussianFit::Compute(const CSpectrum &spectrum,
                                            const TInt32Range &studyRange) {
  Int32 np = (3 + m_PolyOrder + 1);
  Int32 n = studyRange.GetLength();
  if (n < 1 || n < np) {
    return nStatus_IllegalInput;
  }

  // Create suitable first guess
  // To BE improved, trivial version now
  Float64 firstGuessData[np];
  gsl_vector_view firstGuessView = gsl_vector_view_array(firstGuessData, np);

  Float64 peakValue;
  Float64 peakPos;
  Float64 gaussAmp;

  ComputeFirstGuess(spectrum, studyRange, m_PolyOrder, peakValue, peakPos,
                    gaussAmp);

  firstGuessData[0] = peakValue;
  firstGuessData[1] = peakPos;
  firstGuessData[2] = gaussAmp;

  SUserData userData;
  userData.spectrum = &spectrum;
  userData.studyRange = &studyRange;
  userData.polyOrder = m_PolyOrder;

  // Use polinomial part in the form P(x) = a0 + a1*(x-mu) + ... + an*(x-mu)^n
  // mu is the center of the gaussian
  gsl_multifit_function_fdf multifitFunction;
  multifitFunction.f = &GaussF;
  multifitFunction.df = &GaussDF;
  multifitFunction.fdf = &GaussFDF;
  multifitFunction.n = n;
  multifitFunction.p = np;
  multifitFunction.params = &userData;

  // Use a derivative solver Levenberg-Marquardt
  gsl_multifit_fdfsolver *multifitSolver =
      gsl_multifit_fdfsolver_alloc(gsl_multifit_fdfsolver_lmsder, n, np);
  gsl_multifit_fdfsolver_set(multifitSolver, &multifitFunction,
                             &firstGuessView.vector);

  // Iterate
  int status = 0;

  if (0) {
    int iter = 0;
    do {
      iter++;
      status = gsl_multifit_fdfsolver_iterate(multifitSolver);

      if (status)
        break;

      status = gsl_multifit_test_delta(multifitSolver->dx, multifitSolver->x,
                                       m_AbsTol, m_RelTol);
    } while (status == GSL_CONTINUE && iter < 500);
  } else {
    /* solve the system with a maximum of 500 iterations */
    int info;
    gsl_multifit_fdfsolver_driver(multifitSolver, 500, m_RelTol, m_RelTol, 0.0,
                                  &info);
  }

  // Set values and errors
  gsl_matrix *covarMatrix = gsl_matrix_alloc(np, np);
  // gsl_multifit_covar( multifitSolver->J, 0.0, covarMatrix ); //WARNING: fit
  // broken since using GSL2.1 instead of GSL1.16 = OLD version
  gsl_matrix *J =
      gsl_matrix_alloc(n, np); // WARNING: fit broken since using GSL2.1 instead
                               // of GSL1.16  = replacement version
  gsl_multifit_fdfsolver_jac(multifitSolver,
                             J); // WARNING: fit broken since using GSL2.1
                                 // instead of GSL1.16  = replacement version
  gsl_multifit_covar(
      J, 0.0, covarMatrix); // WARNING: fit broken since using GSL2.1 instead of
                            // GSL1.16     = replacement version

  // Float64 chi = gsl_blas_dnrm2( multifitSolver->f );  //WARNING: fit broken
  // since using GSL2.1 instead of GSL1.16 = OLD version
  gsl_vector *res_f;
  res_f = gsl_multifit_fdfsolver_residual(
      multifitSolver); // WARNING: fit broken since using GSL2.1 instead of
                       // GSL1.16     = replacement version
  Float64 chi =
      gsl_blas_dnrm2(res_f); // WARNING: fit broken since using GSL2.1 instead
                             // of GSL1.16     = replacement version

  Float64 dof = n - np;
  Float64 c = GSL_MAX_DBL(1, chi / sqrt(dof));

  Float64 *output = (Float64 *)calloc(n, sizeof(Float64));
  Float64 *outputError = (Float64 *)calloc(n, sizeof(Float64));

  for (Int32 i = 0; i < np; i++) {
    output[i] = gsl_vector_get(multifitSolver->x, i);
    outputError[i] = c * sqrt(gsl_matrix_get(covarMatrix, i, i));
  }

  // Set amplitude > 0 by default
  output[2] = fabs(output[2]);

  m_Amplitude = output[0];
  m_AmplitudeErr = outputError[0];
  m_Mu = output[1];
  m_MuErr = outputError[1];
  m_C = output[2] / sqrt(2);
  // WARNING: this coefficient (sqrt(2)) has been added to compensate for the
  // gaussian expression used in GaussF and GaussDF
  //  so that the gaussian function denominator ( originally:
  //  c0*exp(-1*(x-mu)**2/c2**2) ) becomes c0*exp(-1*(x-mu)**2/ (2*c2**2))
  m_CErr = outputError[2] / sqrt(2);
  m_coeff0 = output[3];

  gsl_multifit_fdfsolver_free(multifitSolver);
  gsl_matrix_free(covarMatrix);
  gsl_matrix_free(J);

  free(output);
  free(outputError);
  return getReturnCode(status);
}

CGaussianFit::EStatus CGaussianFit::getReturnCode(int status) const {

  if (status == GSL_ETOLX)
    return nStatus_Success; // GSL_ETOLX, the change in the position
                            // vector falls below machine precision
  if (status == GSL_CONTINUE)
    return nStatus_IterationHasNotConverged;

  if (status != GSL_SUCCESS)
    return nStatus_FailToReachTolerance;

  return nStatus_Success;
}
/**
 * Set of functions used to model a fit of Gaussian + polynomial term P(mu)
 * = c3
 *
 * mu = c1
 * Y[i] = c0*exp(-1*(x-mu)**2/c2**2) + c3 + c4*(x-mu) + c5*(x-mu)**2 + ... +
 * cn*(x-mu)**(n-3)
 *
 */
int CGaussianFit::GaussF(const gsl_vector *param, void *data, gsl_vector *f) {
  SUserData *userData = (SUserData *)data;
  Int32 n = userData->studyRange->GetLength();
  const Float64 *x = userData->spectrum->GetSpectralAxis().GetSamples() +
                     userData->studyRange->GetBegin();
  const Float64 *y = userData->spectrum->GetFluxAxis().GetSamples() +
                     userData->studyRange->GetBegin();
  // const Float64* err = userData->spectrum->GetFluxAxis().GetError() +
  // userData->studyRange->GetBegin();
  TFloat64List err(n);
  if (true) // WARNING: Hardcoded disable the use of err vect. = 1.0
  {
    for (Int32 i = 0; i < n; i++) {
      err[i] = 1.0;
    }
  }

  Float64 A = gsl_vector_get(param, 0);
  Float64 mu = gsl_vector_get(param, 1);
  Float64 c = gsl_vector_get(param, 2);

  // Polynomial order
  int order = (int)userData->polyOrder;

  for (Int32 i = 0; i < n; i++) {
    /* Polynomial term */
    Float64 Pi = 0;
    for (Int32 k = 0; k <= order; k++) {
      Pi += gsl_vector_get(param, 3 + k) * pow(x[i] - mu, k);
    }

    // Add gaussian term to polynomial term
    Float64 Yi = A * exp(-1. * (x[i] - mu) * (x[i] - mu) / (c * c)) + Pi;
    gsl_vector_set(f, i, (Yi - y[i]) / err[i]);
  }

  return GSL_SUCCESS;
}

/**
 * Jacobian function definition
 */
int CGaussianFit::GaussDF(const gsl_vector *param, void *data, gsl_matrix *J) {
  SUserData *userData = (SUserData *)data;
  Int32 n = userData->studyRange->GetLength();
  const Float64 *x = userData->spectrum->GetSpectralAxis().GetSamples() +
                     userData->studyRange->GetBegin();
  // Float64* err = userData->spectrum->GetFluxAxis().GetError() +
  // userData->studyRange->GetBegin();
  TFloat64List err(n);
  if (true) // WARNING: Hardcoded disable the use of err vect. = 1.0
  {
    for (Int32 i = 0; i < n; i++) {
      err[i] = 1.0;
    }
  }

  Float64 A = gsl_vector_get(param, 0);
  Float64 mu = gsl_vector_get(param, 1);
  Float64 c = gsl_vector_get(param, 2);

  Float64 A_d = 2 * A / (c * c);

  // Polynomial order
  int order = (int)userData->polyOrder;

  // Number of parameters
  for (Int32 i = 0; i < n; i++) {
    // Jacobian matrix J(i,j) = dfi / dxj,
    // where fi = (Yi - yi)/err[i],
    //       Yi = A * exp(-1*(xi-mu)**2/c**2) + P(n, x-mu)

    // Exponential term */
    Float64 e = exp(-1.0 * (x[i] - mu) * (x[i] - mu) / (c * c)) / err[i];

    // Gaussian term
    gsl_matrix_set(J, i, 0, e);
    Float64 P_mu = 0;
    for (Int32 k = 1; k <= order; k++) {
      P_mu +=
          (-1. * (k * pow(x[i] - mu, k - 1) * gsl_vector_get(param, 3 + k)));
    }

    gsl_matrix_set(J, i, 1, A_d * (x[i] - mu) * e + P_mu / err[i]);
    gsl_matrix_set(J, i, 2, A_d * (x[i] - mu) * (x[i] - mu) * e / c);

    // Polynomial term
    for (Int32 k = 0; k <= order; k++)
      gsl_matrix_set(J, i, 3 + k, pow(x[i] - mu, k) / err[i]);
  }
  return GSL_SUCCESS;
}

/**
 * Calls GaussF, then GaussDF and returns GSL_SUCCESS.
 */
int CGaussianFit::GaussFDF(const gsl_vector *param, void *data, gsl_vector *f,
                           gsl_matrix *J) {
  GaussF(param, data, f);
  GaussDF(param, data, J);

  return GSL_SUCCESS;
}
