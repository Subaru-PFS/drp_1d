#include "RedshiftLibrary/linemodel/abstractfitter.h"
#include "RedshiftLibrary/linemodel/hybridfitter.h"
#include "RedshiftLibrary/linemodel/individualfitter.h"
#include "RedshiftLibrary/linemodel/onesfitter.h"
#include "RedshiftLibrary/linemodel/randomfitter.h"
#include "RedshiftLibrary/linemodel/svdfitter.h"
#include "RedshiftLibrary/linemodel/svdlcfitter.h"
#include "RedshiftLibrary/processflow/context.h"

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