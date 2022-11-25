#include "RedshiftLibrary/linemodel/abstractfitter.h"
#include "RedshiftLibrary/linemodel/hybridfitter.h"
#include "RedshiftLibrary/linemodel/individualfitter.h"
#include "RedshiftLibrary/linemodel/onesfitter.h"
#include "RedshiftLibrary/linemodel/randomfitter.h"
#include "RedshiftLibrary/linemodel/svdfitter.h"
#include "RedshiftLibrary/linemodel/svdlcfitter.h"
#include "RedshiftLibrary/processflow/context.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_vector.h>

using namespace NSEpic;

CAbstractFitter::CAbstractFitter(
    CLineModelElementList &elements,
    std::shared_ptr<const CSpectrum> inputSpectrum,
    std::shared_ptr<const TLambdaRange> lambdaRange,
    std::shared_ptr<CSpectrumModel> spectrumModel,
    const CLineCatalog::TLineVector &restLineList)
    : m_Elements(elements), m_inputSpc(*(inputSpectrum)),
      m_RestLineList(restLineList), m_lambdaRange(*(lambdaRange)),
      m_model(spectrumModel) {}

Int32 CAbstractFitter::fitAmplitudesLinSolveAndLambdaOffset(
    TInt32List EltsIdx, const CSpectrumSpectralAxis &spectralAxis,
    std::vector<Float64> &ampsfitted, std::vector<Float64> &errorsfitted,
    bool enableOffsetFitting, Float64 redshift) {

  const CSpectrumFluxAxis &fluxAxis = m_model->getSpcFluxAxisNoContinuum();
  const CSpectrumFluxAxis &continuumfluxAxis = m_model->getContinuumFluxAxis();

  Int32 ret = -1;

  bool atLeastOneOffsetToFit = false;
  if (enableOffsetFitting) {
    for (Int32 iE : EltsIdx)
      for (const auto &line : m_Elements[iE]->m_Lines)
        // check if the line is to be fitted
        if (line.GetOffsetFitEnabled()) {
          atLeastOneOffsetToFit = true;
          break;
        }
  }

  Int32 nSteps =
      atLeastOneOffsetToFit
          ? int((m_LambdaOffsetMax - m_LambdaOffsetMin) / m_LambdaOffsetStep +
                0.5)
          : 1;

  Float64 bestMerit = DBL_MAX;
  Int32 idxBestMerit = -1;
  for (Int32 iO = 0; iO < nSteps; iO++) {
    // set offset value
    if (atLeastOneOffsetToFit)
      setOffset(EltsIdx, iO);

    // fit for this offset
    ret = fitAmplitudesLinSolve(EltsIdx, spectralAxis, fluxAxis,
                                continuumfluxAxis, ampsfitted, errorsfitted,
                                redshift);

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
    return ret;

  // set offset value
  if (atLeastOneOffsetToFit)
    setOffset(EltsIdx, idxBestMerit);
  // fit again for this offset
  ret =
      fitAmplitudesLinSolve(EltsIdx, spectralAxis, fluxAxis, continuumfluxAxis,
                            ampsfitted, errorsfitted, redshift);

  return ret;
}

void CAbstractFitter::setOffset(const TInt32List &EltsIdx,
                                Int32 offsetCount) const {

  Float64 offset = m_LambdaOffsetMin + m_LambdaOffsetStep * offsetCount;
  for (Int32 iE : EltsIdx) {
    for (auto &line : m_Elements[iE]->m_Lines)
      if (line.GetOffsetFitEnabled())
        line.SetOffset(offset);
  }
  return;
}
/**
 * \brief Use GSL to fit linearly the elements listed in argument EltsIdx.
 * If size of argument EltsIdx is less than 1 return -1.
 **/
Int32 CAbstractFitter::fitAmplitudesLinSolve(
    const TInt32List &EltsIdx, const CSpectrumSpectralAxis &spectralAxis,
    const CSpectrumFluxAxis &fluxAxis,
    const CSpectrumFluxAxis &continuumfluxAxis, TFloat64List &ampsfitted,
    TFloat64List &errorsfitted, Float64 redshift) {

  bool useAmpOffset = m_enableAmplitudeOffsets;
  Int32 idxAmpOffset = -1;

  Int32 nddl = EltsIdx.size();
  if (nddl < 1)
    return -1;

  TInt32List xInds = m_Elements.getSupportIndexes(EltsIdx);
  if (xInds.size() < 1)
    return -1;

  if (useAmpOffset) {
    nddl += m_AmplitudeOffsetsDegree + 1;
    // find the amplitudeOffset Support that corresponds to these elts
    idxAmpOffset = m_Elements.getIndexAmpOffset(xInds[0]);
  }

  for (Int32 iddl = 0; iddl < EltsIdx.size(); iddl++)
    m_Elements.SetElementAmplitude(EltsIdx[iddl], 1.0, 0.0);

  const auto &ErrorNoContinuum = m_inputSpc.GetErrorAxis();

  // Linear fit
  Float64 fval;
  double chisq;
  gsl_matrix *X, *cov;
  gsl_vector *y, *w, *c;

  Int32 n = xInds.size();
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

    for (Int32 iddl = 0; iddl < EltsIdx.size(); iddl++) {
      fval = m_Elements[EltsIdx[iddl]]->getModelAtLambda(
          xi, redshift, continuumfluxAxis[idx]);
      gsl_matrix_set(X, i, iddl, fval);
      Log.LogDebug("fval = '%.3e'", fval);
    }

    if (useAmpOffset) {
      gsl_matrix_set(X, i, EltsIdx.size(), 1.0);
      if (m_AmplitudeOffsetsDegree == 0)
        continue;
      gsl_matrix_set(X, i, EltsIdx.size() + 1, xi);
      if (m_AmplitudeOffsetsDegree == 1)
        continue;
      gsl_matrix_set(X, i, EltsIdx.size() + 2, xi * xi);
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

  Int32 sameSign = 1;
  Float64 a0 = gsl_vector_get(c, 0) / normFactor;
  for (Int32 iddl = 1; iddl < EltsIdx.size(); iddl++) {
    Float64 a = gsl_vector_get(c, iddl) / normFactor;
    Float64 product = a0 * a;
    if (product < 0)
      sameSign = 0;
  }

  Log.LogDetail("# Found amplitudes with sameSign=%d", sameSign);
  if (sameSign) {
    for (Int32 iddl = 0; iddl < EltsIdx.size(); iddl++) {
      Float64 a = gsl_vector_get(c, iddl) / normFactor;
      Float64 cova = gsl_matrix_get(cov, iddl, iddl);
      Float64 sigma = sqrt(cova) / normFactor;
      m_Elements.SetElementAmplitude(EltsIdx[iddl], a, sigma);
    }
    // refreshModel();
  } else {
    ampsfitted.resize(EltsIdx.size());
    errorsfitted.resize(EltsIdx.size());
    for (Int32 iddl = 0; iddl < EltsIdx.size(); iddl++) {
      Float64 a = gsl_vector_get(c, iddl) / normFactor;
      Float64 cova = gsl_matrix_get(cov, iddl, iddl);
      Float64 sigma = sqrt(cova) / normFactor;
      m_Elements.SetElementAmplitude(EltsIdx[iddl], a, sigma);
      ampsfitted[iddl] = (a);
      errorsfitted[iddl] = (sigma);
    }
  }

  if (useAmpOffset) {
    Float64 x0 = gsl_vector_get(c, EltsIdx.size()) / normFactor;
    Float64 x1 = 0.0;
    Float64 x2 = 0.0;
    if (m_AmplitudeOffsetsDegree > 0)
      x1 = gsl_vector_get(c, EltsIdx.size() + 1) / normFactor;
    if (m_AmplitudeOffsetsDegree > 1)
      x2 = gsl_vector_get(c, EltsIdx.size() + 2) / normFactor;
    m_Elements.setAmplitudeOffsetsCoeffsAt(idxAmpOffset, {x0, x1, x2});
  }

  gsl_matrix_free(X);
  gsl_vector_free(y);
  gsl_vector_free(w);
  gsl_vector_free(c);
  gsl_matrix_free(cov);

  return sameSign;
}

std::unique_ptr<CAbstractFitter> CAbstractFitter::makeFitter(
    std::string fittingMethod, CLineModelElementList &elements,
    std::shared_ptr<const CSpectrum> inputSpectrum,
    std::shared_ptr<const TLambdaRange> lambdaRange,
    std::shared_ptr<CSpectrumModel> spectrumModel,
    const CLineCatalog::TLineVector &restLineList,
    std::shared_ptr<CContinuumManager> continuumManager) {
  if (fittingMethod == "hybrid")
    return std::unique_ptr<CHybridFitter>(new CHybridFitter(
        elements, inputSpectrum, lambdaRange, spectrumModel, restLineList));
  else if (fittingMethod == "svd")
    return std::unique_ptr<CSvdFitter>(new CSvdFitter(
        elements, inputSpectrum, lambdaRange, spectrumModel, restLineList));
  else if (fittingMethod == "svdlc")
    return std::unique_ptr<CSvdlcFitter>(
        new CSvdlcFitter(elements, inputSpectrum, lambdaRange, spectrumModel,
                         restLineList, continuumManager));
  else if (fittingMethod == "svdlcp2")
    return std::unique_ptr<CSvdlcFitter>(
        new CSvdlcFitter(elements, inputSpectrum, lambdaRange, spectrumModel,
                         restLineList, continuumManager, 2));

  else if (fittingMethod == "ones")
    return std::unique_ptr<COnesFitter>(new COnesFitter(
        elements, inputSpectrum, lambdaRange, spectrumModel, restLineList));
  else if (fittingMethod == "random")
    return std::unique_ptr<CRandomFitter>(new CRandomFitter(
        elements, inputSpectrum, lambdaRange, spectrumModel, restLineList));
  else if (fittingMethod == "individual")
    return std::unique_ptr<CIndividualFitter>(new CIndividualFitter(
        elements, inputSpectrum, lambdaRange, spectrumModel, restLineList));
  else
    THROWG(INTERNAL_ERROR, Formatter()
                               << "Unknown fitting method " << fittingMethod);
}
