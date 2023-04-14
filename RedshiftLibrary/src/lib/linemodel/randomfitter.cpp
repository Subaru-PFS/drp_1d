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
#include "RedshiftLibrary/linemodel/randomfitter.h"
#include "RedshiftLibrary/line/linetags.h"
#include "RedshiftLibrary/processflow/context.h"

using namespace NSEpic;
using namespace std;

void CRandomFitter::doFit(Float64 redshift) {
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
    if (m_Elements[iElts]->GetElementType() == CLine::nType_Absorption) {
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
    for (Int32 j = support[iS].GetBegin(); j <= support[iS].GetEnd(); j++) {
      n++;
      // w = 1.0 / ErrorNoContinuum[j];
      // sumErr += w;
      // m += ContinuumFluxAxis[j] * w;
      m += ContinuumFluxAxis[j];
    }
  }

  return m / Float64(n);
}
