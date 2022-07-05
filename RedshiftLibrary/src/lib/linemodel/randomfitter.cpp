#include "RedshiftLibrary/linemodel/randomfitter.h"
#include "RedshiftLibrary/line/linetags.h"
#include "RedshiftLibrary/processflow/context.h"

using namespace NSEpic;
using namespace std;

CRandomFitter::CRandomFitter(CLineModelElementList &elements,
                             std::shared_ptr<const CSpectrum> inputSpectrum,
                             std::shared_ptr<const TLambdaRange> lambdaRange,
                             std::shared_ptr<CSpectrumModel> spectrumModel)
    : CAbstractFitter(elements, inputSpectrum, lambdaRange, spectrumModel)

{}

void CRandomFitter::fit(Float64 redshift) {
  srand(time(0));
  Float64 randNumFloat = (Float64)rand() / (Float64)(RAND_MAX);

  Float64 coeffAmpEmission = pow(10.0, randNumFloat * 3.0 - 1.0);
  randNumFloat = (Float64)rand() / (Float64)(RAND_MAX);
  Float64 coeffAmpAbsorption = pow(10.0, randNumFloat * 1.0 - 1.0);
  Log.LogInfo("\nLineModel simulation: coeffAmpEmission = %.2f",
              coeffAmpEmission);
  Log.LogInfo("LineModel simulation: coeffAmpAbsorption = %.2f",
              coeffAmpAbsorption);
  // fit the model amplitudes individually
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    Float64 meanContinuum = getContinuumMeanUnderElement(iElts);
    Float64 err = 1e-22;
    Float64 amax = meanContinuum;
    if (m_Elements[iElts]->m_Lines[0].GetType() == CLine::nType_Absorption) {
      amax = meanContinuum * 0.5 * coeffAmpAbsorption;
    } else {
      amax = meanContinuum * coeffAmpEmission;
    }
    randNumFloat = (Float64)rand() / (Float64)(RAND_MAX);
    Float64 a = randNumFloat * amax;
    if (a < 0.0) {
      a = 0.0;
    }
    // get the max nominal amplitude
    Int32 nLines = m_Elements[iElts]->GetSize();
    Float64 maxNominalAmp = -1.0;
    for (Int32 j = 0; j < nLines; j++) {
      if (maxNominalAmp < m_Elements[iElts]->GetNominalAmplitude(j)) {
        maxNominalAmp = m_Elements[iElts]->GetNominalAmplitude(j);
      }
    }

    m_Elements.SetElementAmplitude(iElts, a / maxNominalAmp, err);
  }
}

Float64 CRandomFitter::getContinuumMeanUnderElement(Int32 eltId) const {
  Int32 n = 0;
  Float64 m = 0.0;
  // Float64 sumErr=0.0;

  TInt32RangeList support;
  Int32 iElts = eltId;
  {
    if (m_Elements[iElts]->IsOutsideLambdaRange()) {
      return 0.0;
    }
    TInt32RangeList s = m_Elements[iElts]->getSupport();
    for (Int32 iS = 0; iS < s.size(); iS++) {
      support.push_back(s[iS]);
    }
  }

  const auto &ContinuumFluxAxis = m_model->getContinuumFluxAxis();
  // const auto & ErrorNoContinuum = m_inputSpc->GetFluxAxis().GetError();
  for (Int32 iS = 0; iS < support.size(); iS++) {
    for (Int32 j = support[iS].GetBegin(); j < support[iS].GetEnd(); j++) {
      n++;
      // w = 1.0 / ErrorNoContinuum[j];
      // sumErr += w;
      // m += ContinuumFluxAxis[j] * w;
      m += ContinuumFluxAxis[j];
    }
  }

  return m / Float64(n);
}
