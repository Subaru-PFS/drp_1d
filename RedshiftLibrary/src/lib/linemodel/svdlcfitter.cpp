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
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlin.h>

#include "RedshiftLibrary/common/size.h"
#include "RedshiftLibrary/line/linetags.h"
#include "RedshiftLibrary/linemodel/svdlcfitter.h"
#include "RedshiftLibrary/processflow/context.h"

using namespace NSEpic;
using namespace std;

CSvdlcFitter::CSvdlcFitter(
    const std::shared_ptr<CLMEltListVector> &elementsVector,
    const CCSpectrumVectorPtr &inputSpcs,
    const CTLambdaRangePtrVector &lambdaRanges,
    const CSpcModelVectorPtr &spectrumModels, const CLineMap &restLineList,
    const CSpectraGlobalIndex &spcIndex,
    std::shared_ptr<CContinuumManager> continuumManager, Int32 polyOrder,
    bool enableAmplitudeOffset, bool enableLambdaOffsetsFit)
    : CAbstractFitter(elementsVector, inputSpcs, lambdaRanges, spectrumModels,
                      restLineList, spcIndex, enableAmplitudeOffset,
                      enableLambdaOffsetsFit),
      m_fitc_polyOrder(polyOrder), m_continuumManager(continuumManager),
      m_spectralAxis(getSpectrum().GetSpectralAxis())

{}

// fit the amplitude of all elements AND continuum amplitude together
// with linear solver: gsl_multifit_wlinear

void CSvdlcFitter::doFit(Float64 redshift) {

  // 1. fit only the current continuum
  // prepare continuum on the observed grid

  // re-interpolate the continuum on the grid
  m_continuumManager->reinterpolateContinuumResetAmp();

  m_spectraIndex.setAtBegining(); // TODO dummy impl

  TInt32List validEltsIdx = m_ElementsVector->getValidElementIndices();

  std::string fitGroupTag = "svdlc";
  for (auto const &param : m_ElementsVector->getElementParam())
    param->SetFittingGroupInfo(fitGroupTag);

  TFloat64List ampsfitted;
  TFloat64List errorsfitted;
  Float64 chi2_cl = INFINITY;

  fitAmplitudesLinesAndContinuumLinSolve(validEltsIdx, m_spectralAxis,
                                         ampsfitted, errorsfitted, chi2_cl,
                                         redshift);
  // TODO multiobs loop here
  m_continuumManager->setFitContinuumFromFittedAmps(ampsfitted, validEltsIdx);
  for (auto &spcIndex : m_spectraIndex)
    getModel().initModelWithContinuum();
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
void CSvdlcFitter::fitAmplitudesLinesAndContinuumLinSolve(
    const TInt32List &EltsIdx, const CSpectrumSpectralAxis &spectralAxis,
    TFloat64List &ampsfitted, TFloat64List &errorsfitted, Float64 &chisquare,
    Float64 redshift) {
  m_spectraIndex.setAtBegining(); // TODO dummy impl
  const CSpectrumFluxAxis &fluxAxis = getModel().getSpcFluxAxis();
  const CSpectrumFluxAxis &continuumfluxAxis =
      getModel().getContinuumFluxAxis();
  const auto &ErrorNoContinuum = getSpectrum().GetFluxAxis().GetError();

  if (EltsIdx.size() < 1)
    THROWG(ErrorCode::EMPTY_LIST, Formatter()
                                      << "Input elements list is empty");
  Int32 nddl_ini =
      EltsIdx.size() + 1 +
      std ::max(m_fitc_polyOrder + 1,
                0); // number of param to be fitted=nlines+continuum

  Int32 imin = -1, imax = -1;
  getLambdaRange().getClosedIntervalIndices(spectralAxis.GetSamplesVector(),
                                            imin, imax);
  ampsfitted.assign(nddl_ini, NAN);
  errorsfitted.assign(nddl_ini, NAN);

  Int32 n = imax - imin + 1;
  if (n < nddl_ini)
    THROWG(ErrorCode::LESS_OBSERVED_SAMPLES_THAN_AMPLITUDES_TO_FIT,
           Formatter() << " SVD ill ranked:"
                       << " number of samples = " << n
                       << ", number of parameters to fit = " << nddl_ini);

  for (Int32 eltIdx : EltsIdx)
    m_ElementsVector->SetElementAmplitude(eltIdx, 1.0, 0.0);

  // Linear fit
  double chisq;
  gsl_matrix *X = gsl_matrix_alloc(n, nddl_ini);
  gsl_vector *y = gsl_vector_alloc(n);
  gsl_vector *w = gsl_vector_alloc(n);

  // Normalize
  Float64 normFactor = 1.0 / fluxAxis.computeMaxAbsValue(imin, imax);

  // Prepare the fit data
  fillMatrix(imin, imax, redshift, EltsIdx, spectralAxis, continuumfluxAxis, X);

  for (Int32 i = 0, idx = imin; idx <= imax; ++i, ++idx) {
    Float64 const yi = fluxAxis[idx] * normFactor;
    Float64 const ei = ErrorNoContinuum[idx] * normFactor;
    Float64 const ci = continuumfluxAxis[idx];
    gsl_vector_set(y, i, yi);
    gsl_vector_set(w, i, 1.0 / (ei * ei));
  }

  // remove eventually null columns
  TInt32List EltsIdxToFit;
  TInt32List valid_col_indices;
  X = cleanMatrix(EltsIdx, X, EltsIdxToFit, valid_col_indices);
  if (EltsIdxToFit.empty()) {
    gsl_vector_free(y);
    gsl_vector_free(w);
    return;
  }
  Int32 nddl = X->size2; // nb of columns
  gsl_vector *c = gsl_vector_alloc(nddl);
  gsl_matrix *cov = gsl_matrix_alloc(nddl, nddl);

  gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n, nddl);
  gsl_multifit_wlinear(X, w, y, c, cov, &chisq, work);
  gsl_multifit_linear_free(work);

#define C(i) (gsl_vector_get(c, (i)))
#define COV(i, j) (gsl_matrix_get(cov, (i), (j)))
  Log.LogDebug(Formatter() << "# best fit: Y = " << C(0) << " X1 + " << C(1)
                           << " X2 ...");
  Log.LogDebug(Formatter() << "# covariance matrix:");
  Log.LogDebug(Formatter() << "[ " << COV(0, 0) << ", " << COV(0, 1));
  Log.LogDebug(Formatter() << "  " << COV(1, 0) << ", " << COV(1, 1) << " ]");
  Log.LogDebug(Formatter() << "# chisq = " << chisq
                           << ", chisq/n = " << chisq / n);

  for (Int32 iddl = 0; iddl < nddl; iddl++) {
    Float64 a = gsl_vector_get(c, iddl) / normFactor;
    Log.LogDetail(Formatter() << "# Found amplitude " << iddl << ": " << a);
  }

  Int32 sameSign = 1;
  Float64 a0 = gsl_vector_get(c, 0) / normFactor;
  for (Int32 iddl = 1; iddl < ssize(EltsIdxToFit); iddl++) {
    Float64 const a = gsl_vector_get(c, iddl) / normFactor;
    sameSign &= std::signbit(a0 * a);
  }

  Log.LogDetail(Formatter() << "# Found n=" << EltsIdx.size()
                            << " amplitudes with sameSign=" << sameSign);

  for (Int32 iddl = 0; iddl < nddl; iddl++) {
    Int32 const valid_col = valid_col_indices[iddl];
    Float64 const var = gsl_matrix_get(cov, iddl, iddl);
    errorsfitted[valid_col] = sqrt(var) / normFactor;
    ampsfitted[valid_col] = gsl_vector_get(c, iddl) / normFactor;
    if (iddl < ssize(EltsIdxToFit))
      m_ElementsVector->SetElementAmplitude(
          EltsIdxToFit[iddl], ampsfitted[valid_col], errorsfitted[valid_col]);
  }

  if (m_fitc_polyOrder >= 0) {
    for (Int32 kCoeff = 0; kCoeff < m_fitc_polyOrder + 1; kCoeff++) {
      Float64 const p =
          gsl_vector_get(c, EltsIdx.size() + 1 + kCoeff) / normFactor;
      Log.LogDetail(Formatter()
                    << "# Found p" << kCoeff << " poly amplitude = " << p);
    }
  }

  Log.LogDetail(Formatter() << "# Returning (L+C) n=" << ampsfitted.size()
                            << " amplitudes");

  gsl_matrix_free(X);
  gsl_vector_free(y);
  gsl_vector_free(w);
  gsl_vector_free(c);
  gsl_matrix_free(cov);

  chisquare = chisq;
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
          getElementList()[EltsIdx[iddl]]->getModelAtLambda(xi, redshift, ci);
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

gsl_matrix *CSvdlcFitter::cleanMatrix(const TInt32List &EltsIdx,
                                      gsl_matrix *Xini,
                                      TInt32List &EltsIdxToFit,
                                      TInt32List &valid_col_indices) const {
  Int32 nline = Xini->size1;
  Int32 ncol_ini = Xini->size2;
  gsl_matrix *X = Xini;
  for (Int32 iddl = 0; iddl < ssize(EltsIdx); ++iddl) {
    auto const col_view = gsl_matrix_const_column(Xini, iddl);
    if (gsl_blas_dnrm2(&col_view.vector) > DBL_MIN) {
      valid_col_indices.push_back(iddl);
      EltsIdxToFit.push_back(EltsIdx[iddl]);
    } else {
      // set the amplitude to NAN
      Int32 elt_idx = EltsIdx[iddl];
      m_ElementsVector->SetElementAmplitude(elt_idx, NAN, NAN);
      m_ElementsVector->getElementParam()[elt_idx]->m_nullLineProfiles = true;
      Flag.warning(WarningCode::NULL_LINES_PROFILE,
                   Formatter() << "Null lines profile"
                               << " of elt " << elt_idx);
    }
  }

  if (valid_col_indices.empty()) {
    Flag.warning(WarningCode::NULL_LINES_PROFILE,
                 Formatter() << "SVD reduced to continuum only since all lines "
                                "profiles are null");
  }

  // TODO handle the null continuum case

  for (Int32 icol = EltsIdx.size(); icol < nline; ++icol)
    valid_col_indices.push_back(icol);

  if (EltsIdxToFit.size() == EltsIdx.size())
    return X;

  // replace matrix with a new one containig a copy of valid columns
  Int32 ncol = ncol_ini - (EltsIdx.size() - EltsIdxToFit.size());
  X = gsl_matrix_alloc(nline, ncol);
  // copy the valid columns
  for (Int32 icol = 0; icol != ssize(EltsIdxToFit); ++icol) {
    Int32 const valid_col_idx = valid_col_indices[icol];
    auto col_view = gsl_matrix_column(Xini, valid_col_idx);
    gsl_matrix_set_col(X, icol, &col_view.vector);
  }

  // copy the other columns (continuum + polynomial)
  Int32 ncol_to_copy = ncol_ini - EltsIdx.size();
  auto Xini_view =
      gsl_matrix_const_submatrix(Xini, 0, EltsIdx.size(), nline, ncol_to_copy);
  auto X_view =
      gsl_matrix_submatrix(X, 0, EltsIdxToFit.size(), nline, ncol_to_copy);
  gsl_matrix_memcpy(&X_view.matrix, &Xini_view.matrix);
  gsl_matrix_free(Xini);

  return X;
}
