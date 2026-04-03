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
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_vector.h>

#include "RedshiftLibrary/common/flag.h"
#include "RedshiftLibrary/common/size.h"
#include "RedshiftLibrary/linemodel/svdfitter.h"
#include "RedshiftLibrary/processflow/context.h"

using namespace NSEpic;
using namespace std;

CSvdFitter::CSvdFitter(const std::shared_ptr<CLMEltListVector> &elementsVector,
                       const CCSpectrumVectorPtr &inputSpcs,
                       const CTLambdaRangePtrVector &lambdaRanges,
                       const CSpcModelVectorPtr &spectrumModels,
                       const CLineMap &restLineList,
                       const CSpectraGlobalIndex &spcIndex,
                       bool enableAmplitudeOffsets, bool enableLambdaOffsetsFit)
    : CAbstractFitter(elementsVector, inputSpcs, lambdaRanges, spectrumModels,
                      restLineList, spcIndex, enableAmplitudeOffsets,
                      enableLambdaOffsetsFit) {
  if (HasLambdaOffsetFitting() && m_enableLambdaOffsetsFit) {
    std::shared_ptr<const CParameterStore> const &ps =
        Context.GetParameterStore();
    m_LambdaOffsetMax = ps->GetScoped<Float64>("lbdaOffsetMax");
    m_LambdaOffsetMin = -m_LambdaOffsetMax;
    auto const &opt_fitting_method =
        ps->GetScoped<std::string>("fittingMethod");
    if (opt_fitting_method == "svd" || opt_fitting_method == "hybrid")
      m_LambdaOffsetStep = ps->GetScoped<Float64>("lbdaOffsetStep");
  }
};

void CSvdFitter::doFit(Float64 redshift) {
  // NB dummy multiobs implementation (functional for one obs only)
  m_spectraIndex.setAtBegining(); // temporary multiobs implementation
  TInt32List validEltsIdx = m_ElementsVector->getValidElementIndices();
  if (validEltsIdx.empty())
    return;

  std::string fitGroupTag = "svd";
  for (auto const &param : m_ElementsVector->getElementsParams())
    param->SetFittingGroupInfo(fitGroupTag);

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

  const CSpectrumSpectralAxis &spectralAxis = getSpectrum().GetSpectralAxis();
  const CSpectrumFluxAxis &fluxAxis = getModel().getSpcFluxAxisNoContinuum();
  const CSpectrumFluxAxis &continuumfluxAxis =
      getModel().getContinuumFluxAxis();

  if (EltsIdx.size() < 1)
    THROWG(ErrorCode::INTERNAL_ERROR, "empty Line element list to fit");

  TInt32List xInds =
      getElementList().getSupportIndexes(EltsIdx, m_enableAmplitudeOffsets);
  Int32 n = xInds.size();
  if (n < 1)
    THROWG(ErrorCode::INTERNAL_ERROR,
           "no observed samples for the line Element to fit");

  // list of elements to fit
  TInt32List EltsIdxToFit_ini;
  if (IdxToFit.empty())
    EltsIdxToFit_ini = EltsIdx;
  else {
    EltsIdxToFit_ini.reserve(IdxToFit.size());
    for (auto idx : IdxToFit)
      EltsIdxToFit_ini.push_back(EltsIdx[idx]);
  }
  Int32 nddl_ini = EltsIdxToFit_ini.size();

  Int32 ncol_polynome = 0;
  if (m_enableAmplitudeOffsets) {
    ncol_polynome = CPolynomCoeffs::degree + 1;
    nddl_ini += ncol_polynome;
  }

  if (n < nddl_ini) {
    Flag.warning(WarningCode::LESS_OBSERVED_SAMPLES_THAN_AMPLITUDES_TO_FIT,
                 Formatter() << "SVD aborted since ill ranked:"
                             << " number of samples = " << n
                             << ", number of parameters to fit = " << nddl_ini);
    for (Int32 iElt : EltsIdxToFit_ini) {
      m_ElementsVector->SetElementAmplitude(iElt, NAN, NAN);
      m_ElementsVector->getElementsParams()[iElt]->m_fitFailed = true;
    }
    if (m_enableAmplitudeOffsets) {
      for (Int32 iElt : EltsIdxToFit_ini)
        m_ElementsVector->getElementsParams()[iElt]->SetPolynomCoeffs(
            {NAN, NAN, NAN});
    }
    return true;
  }

  for (auto const iElt : EltsIdxToFit_ini) {
    m_ElementsVector->getElementsParams()[iElt]->SetAbsLinesLimit(INFINITY);
    m_ElementsVector->SetElementAmplitude(iElt, 1.0, 0.0);
    m_ElementsVector->getElementsParams()[iElt]->SetAbsLinesLimit(1.0);
  }

  // Linear fit
  Float64 fval;
  double chisq;
  gsl_matrix *Xini, *X, *cov;
  gsl_vector *y, *w, *c;

  Xini = gsl_matrix_alloc(n, nddl_ini);
  y = gsl_vector_alloc(n);
  w = gsl_vector_alloc(n);

  // Normalize
  Int32 maxabsval_idx = *std::max_element(
      xInds.cbegin(), xInds.cend(), [&fluxAxis](Int32 l, Int32 r) {
        return std::abs(fluxAxis[l]) < std::abs(fluxAxis[r]);
      });
  Float64 normFactor = 1.0 / std::abs(fluxAxis[maxabsval_idx]);

  Log.LogDetail(Formatter() << "normFactor = '" << std::fixed
                            << std::setprecision(3) << normFactor << "'");

  // handle continuum with ampOffset (polynomial under line) enabled
  std::function getContinuum = [&continuumfluxAxis, normFactor](Int32 idx) {
    return continuumfluxAxis[idx] * normFactor;
  };
  if (m_enableAmplitudeOffsets) {
    // in that case, we cannot perform a linear leastsquare on both the
    // absorption amplitude and the polynomial coeffs. we thus assumes the
    // continnum is constant under the line, and will correct for it after the
    // fit
    getContinuum = [](Int32 idx) { return 1; };
  }

  // Prepare the fit data
  for (Int32 i = 0; i < n; i++) {
    double xi, yi, wi;
    Int32 idx = xInds[i];
    xi = spectralAxis[idx];
    yi = fluxAxis[idx] * normFactor;
    wi = fluxAxis.GetWeight(idx, normFactor);

    gsl_vector_set(y, i, yi);
    gsl_vector_set(w, i, wi);

    for (Int32 iddl = 0; iddl < ssize(EltsIdxToFit_ini); iddl++) {
      fval = getElementList()[EltsIdxToFit_ini[iddl]]->getModelAtLambda(
          xi, redshift, getContinuum(idx));
      gsl_matrix_set(Xini, i, iddl, fval);
      Log.LogDebug(Formatter() << "fval = '" << fval << "'");
    }

    if (m_enableAmplitudeOffsets) {
      gsl_matrix_set(Xini, i, EltsIdxToFit_ini.size(), 1.0);
      if (CPolynomCoeffs::degree == 0)
        continue;
      gsl_matrix_set(Xini, i, EltsIdxToFit_ini.size() + 1, xi);
      if (CPolynomCoeffs::degree == 1)
        continue;
      gsl_matrix_set(Xini, i, EltsIdxToFit_ini.size() + 2, xi * xi);
    }
  }

  // detect null colunms (null profiles)
  // reduce the matrix accordingly, and set the null elements m_null_profiles
  // true
  TInt32List valid_col_indices;
  TInt32List EltsIdxToFit;
  Int32 nddl = nddl_ini;
  X = Xini;
  for (Int32 iddl = 0; iddl < ssize(EltsIdxToFit_ini); iddl++) {
    auto const col_view = gsl_matrix_const_column(Xini, iddl);
    if (gsl_blas_dnrm2(&col_view.vector) > DBL_MIN) {
      valid_col_indices.push_back(iddl);
      EltsIdxToFit.push_back(EltsIdxToFit_ini[iddl]);
    } else {
      // set the amplitude to NAN
      Int32 const elt_idx = EltsIdxToFit_ini[iddl];
      m_ElementsVector->SetElementAmplitude(elt_idx, NAN, NAN);
      m_ElementsVector->getElementsParams()[elt_idx]->m_nullLineProfiles = true;
      Flag.warning(WarningCode::NULL_LINES_PROFILE,
                   Formatter() << "Null lines profile"
                               << " of elt " << elt_idx);
    }
  }

  if (valid_col_indices.empty()) {
    Flag.warning(WarningCode::NULL_LINES_PROFILE,
                 Formatter() << "SVD aborted since all lines null");
    if (m_enableAmplitudeOffsets) {
      for (Int32 elt_idx : EltsIdxToFit_ini)
        m_ElementsVector->getElementsParams()[elt_idx]->SetPolynomCoeffs(
            {NAN, NAN, NAN});
    }
    gsl_matrix_free(Xini);
    gsl_vector_free(y);
    gsl_vector_free(w);
    return true;
  }

  if (valid_col_indices.size() != EltsIdxToFit_ini.size()) {
    nddl -= EltsIdxToFit_ini.size() - valid_col_indices.size();
    X = gsl_matrix_alloc(n, nddl);
    // copy the valid columns
    for (Int32 icol = 0; icol != ssize(valid_col_indices); ++icol) {
      Int32 const valid_col_idx = valid_col_indices[icol];
      auto col_view = gsl_matrix_column(Xini, valid_col_idx);
      gsl_matrix_set_col(X, icol, &col_view.vector);
    }
    // copy the polynomial columns
    if (m_enableAmplitudeOffsets) {
      auto Xini_view = gsl_matrix_const_submatrix(
          Xini, 0, nddl_ini - ncol_polynome, n, ncol_polynome);
      auto X_view =
          gsl_matrix_submatrix(X, 0, nddl - ncol_polynome, n, ncol_polynome);
      gsl_matrix_memcpy(&X_view.matrix, &Xini_view.matrix);
    }
    gsl_matrix_free(Xini);
  }

  c = gsl_vector_alloc(nddl);
  cov = gsl_matrix_alloc(nddl, nddl);

  // perform the fit
  {
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n, nddl);
    gsl_multifit_wlinear(X, w, y, c, cov, &chisq, work);
    gsl_multifit_linear_free(work);
  }

  auto const eltSize = ssize(EltsIdxToFit);
#define COV(i, j) (gsl_matrix_get(cov, i, j))

  // get the continuum polynomials solution and push to the LineElements
  if (m_enableAmplitudeOffsets) {
    Float64 x0 = gsl_vector_get(c, eltSize);
    Float64 x1 = 0.0;
    Float64 x2 = 0.0;
    if (CPolynomCoeffs::degree > 0)
      x1 = gsl_vector_get(c, eltSize + 1);
    if (CPolynomCoeffs::degree > 1)
      x2 = gsl_vector_get(c, eltSize + 2);
    // set the polynomial coeffs for all elements, even those not fitted and
    // fixed at zero
    auto const s = eltSize;
    Eigen::Matrix3d polyCoeffsCovar{
        {COV(s, s), COV(s, s + 1), COV(s, s + 2)},
        {COV(s + 1, s), COV(s + 1, s + 1), COV(s + 1, s + 2)},
        {COV(s + 2, s), COV(s + 2, s + 1), COV(s + 2, s + 2)}};
    CPolynomCoeffs polyCoeffs{x0, x1, x2, polyCoeffsCovar};
    for (Int32 iddl = 0; iddl < ssize(EltsIdx); ++iddl)
      m_ElementsVector->getElementsParams()[EltsIdx[iddl]]->SetPolynomCoeffs(
          polyCoeffs * (1 / normFactor));
  }

  bool allPositive = true;

  // get the amplitudes solution and push into the LineElements
  ampsfitted.assign(EltsIdxToFit_ini.size(), NAN);
  errorsfitted.assign(EltsIdxToFit_ini.size(), NAN);
  for (Int32 iddl = 0; iddl < eltSize; iddl++) {
    auto const &elt = m_ElementsVector->getElementList()[EltsIdxToFit[iddl]];
    auto const &elt_param =
        m_ElementsVector->getElementsParams()[EltsIdxToFit[iddl]];
    bool const rescaleAbsToFittedCont =
        m_enableAmplitudeOffsets && elt_param->IsAbsorption();
    Float64 a = elt_param->IsEmission() || rescaleAbsToFittedCont
                    ? gsl_vector_get(c, iddl) / normFactor
                    : gsl_vector_get(c, iddl);
    Log.LogDetail(Formatter() << "# Found amplitude " << iddl << ": " << a);
    Float64 const var = COV(iddl, iddl);
    Float64 std = elt_param->IsEmission() || rescaleAbsToFittedCont
                      ? sqrt(var) / normFactor
                      : sqrt(var);

    // rescale the absorption amplitudes to continuum + fitted polynomial
    if (rescaleAbsToFittedCont) {
      auto const &[cont, _] = elt->GetContinuumAtCenterProfile(
          spectralAxis, redshift, continuumfluxAxis, true);
      a /= cont;
      std /= cont;
    }
    if (a < 0)
      allPositive = false;
    m_ElementsVector->SetElementAmplitude(EltsIdxToFit[iddl], a, std);
    ampsfitted[valid_col_indices[iddl]] = a;
    errorsfitted[valid_col_indices[iddl]] = std;
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
  // NB dummy multiobs implementation (functional for one obs only)
  TFloat64List ampsfitted;
  TFloat64List errorsfitted;

  // filter out not fittable lines
  auto const ValidEltsIdx = m_ElementsVector->getValidElementIndices(EltsIdx);
  bool const allPositive =
      fitAmplitudesLinSolve(ValidEltsIdx, ampsfitted, errorsfitted, redshift);
  if (allPositive)
    return;

  TInt32List idx_positive;
  for (Int32 ifit = 0; ifit < ssize(ValidEltsIdx); ifit++) {
    if (ampsfitted[ifit] > 0)
      idx_positive.push_back(ifit);
    else if (isfinite(ampsfitted[ifit]) && ampsfitted[ifit] < 0) {
      m_ElementsVector->SetElementAmplitude(ValidEltsIdx[ifit], 0.0,
                                            errorsfitted[ifit]);
    }
  }
  if (idx_positive.empty())
    return;

  // refit the positive elements together
  if (!m_enableAmplitudeOffsets && idx_positive.size() == 1) {
    fitAmplitude(ValidEltsIdx[idx_positive.front()], redshift, undefIdx);
    m_spectraIndex.setAtBegining(); // temporary multiobs implementation
    return;
  }
  bool const allPositive2 = fitAmplitudesLinSolve(
      ValidEltsIdx, ampsfitted, errorsfitted, redshift, idx_positive);
  if (allPositive2) {
    m_spectraIndex.setAtBegining(); // temporary multiobs implementation
    return;
  }
  for (Int32 irefit = 0; irefit < ssize(idx_positive); ++irefit) {
    // set at zero negative amplitudes
    if (isfinite(ampsfitted[irefit]) && ampsfitted[irefit] < 0) {
      m_ElementsVector->SetElementAmplitude(ValidEltsIdx[idx_positive[irefit]],
                                            0.0, errorsfitted[irefit]);
    }
  }
  m_spectraIndex.setAtBegining(); // temporary multiobs implementation
}

void CSvdFitter::fitAmplitudesLinSolveAndLambdaOffset(TInt32List EltsIdx,
                                                      bool enableOffsetFitting,
                                                      Float64 redshift) {

  bool atLeastOneOffsetToFit =
      HasLineElementToOffset(EltsIdx, enableOffsetFitting);
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
    getModel().refreshModelUnderElements(EltsIdx);

    sumFit += getModelResidualRmsUnderElements(EltsIdx);

    if (sumFit < bestMerit) {
      bestMerit = sumFit;
      idxBestMerit = iO;
    }
  }

  if (idxBestMerit == -1 || !atLeastOneOffsetToFit)
    return;

  // set offset value
  setLambdaOffset(EltsIdx, idxBestMerit);
  // fit again for this offset
  fitAmplitudesLinSolvePositive(EltsIdx, redshift);
}
