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
#include "RedshiftLibrary/linemodel/svdfitter.h"
#include "RedshiftLibrary/common/flag.h"
#include "RedshiftLibrary/line/linetags.h"
#include "RedshiftLibrary/processflow/context.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_vector.h>

using namespace NSEpic;
using namespace std;

// set all the amplitudes to 1.0
void CSvdFitter::doFit(Float64 redshift) {
  TInt32List validEltsIdx = m_Elements.GetModelValidElementsIndexes();

  if (validEltsIdx.empty())
    return;

  std::string fitGroupTag = "svd";
  for (Int32 idx : validEltsIdx)
    m_Elements[idx]->SetFittingGroupInfo(fitGroupTag);

  if (m_enableAmplitudeOffsets)
    m_Elements.resetAmplitudeOffset();

  fitAmplitudesLinSolveAndLambdaOffset(validEltsIdx, m_enableLambdaOffsetsFit,
                                       redshift);
}

/**
 * \brief Use GSL to fit linearly the elements listed in argument EltsIdx.
 * If size of argument EltsIdx is less than 1 return -1.
 **/
bool CSvdFitter::fitAmplitudesLinSolve(const TInt32List &EltsIdx,
                                       TFloat64List &ampsfitted,
                                       TFloat64List &errorsfitted,
                                       Float64 redshift,
                                       const TInt32List &IdxToFit) {

  const CSpectrumSpectralAxis &spectralAxis = m_inputSpc.GetSpectralAxis();
  const CSpectrumFluxAxis &fluxAxis = m_model->getSpcFluxAxisNoContinuum();
  const CSpectrumFluxAxis &continuumfluxAxis = m_model->getContinuumFluxAxis();

  bool useAmpOffset = m_enableAmplitudeOffsets;

  if (EltsIdx.size() < 1)
    THROWG(INTERNAL_ERROR, "empty Line element list to fit");

  TInt32List xInds = m_Elements.getSupportIndexes(EltsIdx);
  Int32 n = xInds.size();
  if (n < 1)
    THROWG(INTERNAL_ERROR, "no observed samples for the line Element to fit");

  // list of elements to fit
  TInt32List EltsIdxToFit;
  if (IdxToFit.empty())
    EltsIdxToFit = EltsIdx;
  else {
    EltsIdxToFit.reserve(IdxToFit.size());
    for (auto idx : IdxToFit)
      EltsIdxToFit.push_back(EltsIdx[idx]);
  }
  Int32 nddl = EltsIdxToFit.size();

  if (useAmpOffset)
    nddl += TPolynomCoeffs::degree + 1;

  if (n < nddl) {
    Flag.warning(WarningCode::LINEARFIT_RANK_DEFICIENT,
                 Formatter() << __func__ << " SVD ill ranked:"
                             << " number of samples = " << n
                             << ", number of parameters to fit = " << nddl);
    ampsfitted.assign(EltsIdx.size(), 0.0);
    errorsfitted.assign(EltsIdx.size(), INFINITY);
    for (Int32 iddl = 0; iddl < EltsIdx.size(); iddl++)
      m_Elements.SetElementAmplitude(EltsIdx[iddl], 0., INFINITY);
    if (useAmpOffset) {
      for (Int32 iddl = 0; iddl < EltsIdx.size(); ++iddl)
        m_Elements[EltsIdx[iddl]]->SetPolynomCoeffs({0., 0., 0.});
    }
    return true;
  }

  for (Int32 iddl = 0; iddl < EltsIdxToFit.size(); iddl++)
    m_Elements.SetElementAmplitude(EltsIdxToFit[iddl], 1.0, 0.0);

  const auto &ErrorNoContinuum = m_inputSpc.GetErrorAxis();

  // Linear fit
  Float64 fval;
  double chisq;
  gsl_matrix *X, *cov;
  gsl_vector *y, *w, *c;

  X = gsl_matrix_alloc(n, nddl);
  y = gsl_vector_alloc(n);
  w = gsl_vector_alloc(n);
  c = gsl_vector_alloc(nddl);
  cov = gsl_matrix_alloc(nddl, nddl);

  // Normalize
  Float64 maxabsval = DBL_MIN;
  for (auto idx : xInds) {
    if (maxabsval < std::abs(fluxAxis[idx]))
      maxabsval = std::abs(fluxAxis[idx]);
  }
  Float64 normFactor = 1.0 / maxabsval;
  Log.LogDetail("normFactor = '%.3e'\n", normFactor);

  // Prepare the fit data
  for (Int32 i = 0; i < n; i++) {
    double xi, yi, ei;
    Int32 idx = xInds[i];
    xi = spectralAxis[idx];
    yi = fluxAxis[idx] * normFactor;
    ei = ErrorNoContinuum[idx] * normFactor;

    gsl_vector_set(y, i, yi);
    gsl_vector_set(w, i, 1.0 / (ei * ei));

    for (Int32 iddl = 0; iddl < EltsIdxToFit.size(); iddl++) {
      fval = m_Elements[EltsIdxToFit[iddl]]->getModelAtLambda(
          xi, redshift, continuumfluxAxis[idx]);
      gsl_matrix_set(X, i, iddl, fval);
      Log.LogDebug("fval = '%.3e'", fval);
    }

    if (useAmpOffset) {
      gsl_matrix_set(X, i, EltsIdxToFit.size(), 1.0);
      if (TPolynomCoeffs::degree == 0)
        continue;
      gsl_matrix_set(X, i, EltsIdxToFit.size() + 1, xi);
      if (TPolynomCoeffs::degree == 1)
        continue;
      gsl_matrix_set(X, i, EltsIdxToFit.size() + 2, xi * xi);
    }
  }

  {
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n, nddl);
    gsl_multifit_wlinear(X, w, y, c, cov, &chisq, work);
    gsl_multifit_linear_free(work);
  }

  for (Int32 iddl = 0; iddl < nddl; iddl++) {
    Float64 a = gsl_vector_get(c, iddl) / normFactor;
    Log.LogDetail("# Found amplitude %d: %+.5e", iddl, a);
  }

  bool allPositive = true;

  Log.LogDetail("# Found negative amplitudes with sameSign");
  ampsfitted.resize(EltsIdxToFit.size());
  errorsfitted.resize(EltsIdxToFit.size());
  for (Int32 iddl = 0; iddl < EltsIdxToFit.size(); iddl++) {
    Float64 a = gsl_vector_get(c, iddl) / normFactor;
    if (a < 0)
      allPositive = false;
    Float64 cova = gsl_matrix_get(cov, iddl, iddl);
    Float64 sigma = sqrt(cova) / normFactor;
    m_Elements.SetElementAmplitude(EltsIdxToFit[iddl], a, sigma);
    ampsfitted[iddl] = (a);
    errorsfitted[iddl] = (sigma);
  }

  if (useAmpOffset) {
    Float64 x0 = gsl_vector_get(c, EltsIdxToFit.size()) / normFactor;
    Float64 x1 = 0.0;
    Float64 x2 = 0.0;
    if (TPolynomCoeffs::degree > 0)
      x1 = gsl_vector_get(c, EltsIdxToFit.size() + 1) / normFactor;
    if (TPolynomCoeffs::degree > 1)
      x2 = gsl_vector_get(c, EltsIdxToFit.size() + 2) / normFactor;
    // set the polynomial coeffs for all elements, even those not fitted and
    // fixed at zero
    for (Int32 iddl = 0; iddl < EltsIdx.size(); ++iddl)
      m_Elements[EltsIdx[iddl]]->SetPolynomCoeffs({x0, x1, x2});
  }

  gsl_matrix_free(X);
  gsl_vector_free(y);
  gsl_vector_free(w);
  gsl_vector_free(c);
  gsl_matrix_free(cov);

  return allPositive;
}

// poor-man positive constraint
// if some amplitudes are negative, force them to zero and refit the
// positive ones
void CSvdFitter::fitAmplitudesLinSolvePositive(const TInt32List &EltsIdx,
                                               Float64 redshift) {

  TFloat64List ampsfitted;
  TFloat64List errorsfitted;

  bool allPositive =
      fitAmplitudesLinSolve(EltsIdx, ampsfitted, errorsfitted, redshift);

  if (!allPositive) {
    TInt32List idx_positive;
    for (Int32 ifit = 0; ifit < EltsIdx.size(); ifit++) {
      if (ampsfitted[ifit] < 0) {
        m_Elements.SetElementAmplitude(EltsIdx[ifit], 0.0, errorsfitted[ifit]);
      } else {
        idx_positive.push_back(ifit);
      }
    }
    // refit the positive elements together
    if (!m_enableAmplitudeOffsets && idx_positive.size() == 1) {
      fitAmplitude(EltsIdx[idx_positive.front()], redshift, undefIdx);
    } else if (idx_positive.size() > 0) {
      bool allPositive2 = fitAmplitudesLinSolve(
          EltsIdx, ampsfitted, errorsfitted, redshift, idx_positive);

      if (!allPositive2) {
        for (Int32 irefit = 0; irefit < idx_positive.size(); ++irefit) {
          if (ampsfitted[irefit] > 0) {
            fitAmplitude(EltsIdx[idx_positive[irefit]], redshift, undefIdx);
          } else {
            m_Elements.SetElementAmplitude(EltsIdx[idx_positive[irefit]], 0.0,
                                           errorsfitted[irefit]);
          }
        }
      }
    }
  }
}

void CSvdFitter::fitAmplitudesLinSolveAndLambdaOffset(TInt32List EltsIdx,
                                                      bool enableOffsetFitting,
                                                      Float64 redshift) {

  const CSpectrumFluxAxis &fluxAxis = m_model->getSpcFluxAxisNoContinuum();

  bool atLeastOneOffsetToFit =
      HasLambdaOffsetFitting(EltsIdx, enableOffsetFitting);
  Int32 nSteps = GetLambdaOffsetSteps(atLeastOneOffsetToFit);

  Float64 bestMerit = DBL_MAX;
  Int32 idxBestMerit = -1;
  for (Int32 iO = 0; iO < nSteps; iO++) {
    // set offset value
    if (atLeastOneOffsetToFit)
      setLambdaOffset(EltsIdx, iO);

    // fit for this offset
    fitAmplitudesLinSolvePositive(EltsIdx, redshift);

    // check fitting
    if (!atLeastOneOffsetToFit)
      continue;

    Float64 sumFit = 0.0;
    m_model->refreshModelUnderElements(EltsIdx);

    // todo: replace lambdarange using elements limits for speed
    for (Int32 iE : EltsIdx)
      sumFit += m_model->getModelErrorUnderElement(iE, fluxAxis);

    if (sumFit < bestMerit) {
      bestMerit = sumFit;
      idxBestMerit = iO;
    }
  }

  if (idxBestMerit == -1 || !atLeastOneOffsetToFit)
    return;

  // set offset value
  if (atLeastOneOffsetToFit)
    setLambdaOffset(EltsIdx, idxBestMerit);
  // fit again for this offset
  fitAmplitudesLinSolvePositive(EltsIdx, redshift);
}