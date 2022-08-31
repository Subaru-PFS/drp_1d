#include "RedshiftLibrary/linemodel/individualfitter.h"
#include "RedshiftLibrary/line/linetags.h"
#include "RedshiftLibrary/processflow/context.h"

using namespace NSEpic;
using namespace std;

CIndividualFitter::CIndividualFitter(
    CLineModelElementList &elements,
    std::shared_ptr<const CSpectrum> inputSpectrum,
    std::shared_ptr<const TLambdaRange> lambdaRange,
    std::shared_ptr<CSpectrumModel> spectrumModel)
    : CAbstractFitter(elements, inputSpectrum, lambdaRange, spectrumModel)

{}

// fit the amplitudes of each element independently, unless there is
// overlap
void CIndividualFitter::fit(Float64 redshift) {
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    m_Elements[iElts]->fitAmplitude(m_inputSpc.GetSpectralAxis(),
                                    m_model->getSpcFluxAxisNoContinuum(),
                                    m_model->getContinuumFluxAxis(), redshift);
  }
}
