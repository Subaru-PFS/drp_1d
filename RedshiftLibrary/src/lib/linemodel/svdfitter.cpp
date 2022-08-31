#include "RedshiftLibrary/linemodel/svdfitter.h"
#include "RedshiftLibrary/line/linetags.h"
#include "RedshiftLibrary/processflow/context.h"

using namespace NSEpic;
using namespace std;

CSvdFitter::CSvdFitter(CLineModelElementList &elements,
                       std::shared_ptr<const CSpectrum> inputSpectrum,
                       std::shared_ptr<const TLambdaRange> lambdaRange,
                       std::shared_ptr<CSpectrumModel> spectrumModel)
    : CAbstractFitter(elements, inputSpectrum, lambdaRange, spectrumModel)

{}

// set all the amplitudes to 1.0
void CSvdFitter::fit(Float64 redshift) {
  TInt32List validEltsIdx = m_Elements.GetModelValidElementsIndexes();
  TFloat64List ampsfitted;
  TFloat64List errorsfitted;
  fitAmplitudesLinSolveAndLambdaOffset(
      validEltsIdx, m_inputSpc.GetSpectralAxis(),
      m_model->getSpcFluxAxisNoContinuum(), m_model->getContinuumFluxAxis(),
      ampsfitted, errorsfitted, m_enableLambdaOffsetsFit, redshift);
}
