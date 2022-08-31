#include "RedshiftLibrary/linemodel/onesfitter.h"
#include "RedshiftLibrary/line/linetags.h"
#include "RedshiftLibrary/processflow/context.h"

using namespace NSEpic;
using namespace std;

COnesFitter::COnesFitter(CLineModelElementList &elements,
                         std::shared_ptr<const CSpectrum> inputSpectrum,
                         std::shared_ptr<const TLambdaRange> lambdaRange,
                         std::shared_ptr<CSpectrumModel> spectrumModel)
    : CAbstractFitter(elements, inputSpectrum, lambdaRange, spectrumModel)

{}

// set all the amplitudes to 1.0
void COnesFitter::fit(Float64 redshift) {
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    m_Elements[iElts]->SetFittedAmplitude(1.0, 1.0);
  }
}
