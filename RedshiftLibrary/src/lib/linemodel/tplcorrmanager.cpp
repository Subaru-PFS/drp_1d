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

#include "RedshiftLibrary/linemodel/tplcorrmanager.h"
#include "RedshiftLibrary/line/catalogsTplRatio.h"
#include "RedshiftLibrary/linemodel/spectrummodel.h"

using namespace NSEpic;

CTplCorrManager::CTplCorrManager(
    CLineModelElementList &elements, std::shared_ptr<CSpectrumModel> model,
    std::shared_ptr<const CSpectrum> inputSpc,
    std::shared_ptr<const TFloat64Range> lambdaRange,
    std::shared_ptr<CContinuumManager> continuumManager,
    const TLineVector &restLineList)
    : CTplratioManager(elements, model, inputSpc, lambdaRange, continuumManager,
                       restLineList) {}

Float64 CTplCorrManager::computeMerit(Int32 itratio) {

  m_model->refreshModel();
  Int32 s = m_RestLineList.size();
  TFloat64List Amplitudes(s);
  TFloat64List AmplitudesUncertainties(s); // noise sigma
  for (Int32 iRestLine = 0; iRestLine < s; iRestLine++) {
    Int32 subeIdx = undefIdx;
    Int32 eIdx = m_Elements.findElementIndex(iRestLine, subeIdx);

    if (eIdx == undefIdx || subeIdx == undefIdx ||
        m_Elements[eIdx]->IsOutsideLambdaRange(subeIdx)) {
      continue; // data already set to its default values
    }
    Float64 amp = m_Elements[eIdx]->GetFittedAmplitude(subeIdx);
    Amplitudes[iRestLine] = amp;
    Float64 ampError = m_Elements[eIdx]->GetFittedAmplitudeErrorSigma(subeIdx);
    AmplitudesUncertainties[iRestLine] = ampError;
  }
  TFloat64List correctedAmplitudes(Amplitudes.size());
  m_CatalogTplRatio->GetBestFit(m_RestLineList, Amplitudes,
                                AmplitudesUncertainties, correctedAmplitudes,
                                m_tplratioBestTplName);
  for (Int32 iRestLine = 0; iRestLine < m_RestLineList.size(); iRestLine++) {
    Int32 subeIdx = undefIdx;
    Int32 eIdx = m_Elements.findElementIndex(iRestLine, subeIdx);
    if (eIdx == undefIdx || subeIdx == undefIdx)
      continue;

    Float64 er = AmplitudesUncertainties[subeIdx]; // not modifying the
    // fitting error for now
    Float64 nominalAmp = m_Elements[eIdx]->GetNominalAmplitude(subeIdx);
    m_Elements[eIdx]->SetElementAmplitude(
        correctedAmplitudes[iRestLine] / nominalAmp, er);
  }
  m_model->refreshModel();
  return getLeastSquareMerit();
}

void CTplCorrManager::saveResults(Int32 itratio) {
  //  m_tplratioBestTplName = bestTplratioName; done in computeMerit directly
}
