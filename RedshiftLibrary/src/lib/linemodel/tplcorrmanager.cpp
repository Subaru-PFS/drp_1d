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
    const std::shared_ptr<CLMEltListVector> &elementsVector,
    const CSpcModelVectorPtr &models, const CCSpectrumVectorPtr &inputSpcs,
    const CTLambdaRangePtrVector &lambdaRanges,
    std::shared_ptr<CContinuumManager> continuumManager,
    const CLineMap &restLineList, const std::shared_ptr<Int32> &curObs)
    : CTplratioManager(elementsVector, models, inputSpcs, lambdaRanges,
                       continuumManager, restLineList, curObs) {}

Float64 CTplCorrManager::computeMerit(Int32 itratio) {
  *m_curObs = 0; // dummy implementation

  getModel().refreshModel();
  TFloat64List Amplitudes;
  TFloat64List AmplitudesUncertainties; // noise sigma
  std::vector<std::pair<Int32, Int32>> eIdxList;
  TInt32List validLinesIndex;
  for (auto [i_lineCatalog, line_it] =
           std::pair(Int32(0), m_RestLineList.cbegin());
       line_it != m_RestLineList.cend(); ++i_lineCatalog, ++line_it) {
    Int32 const line_id = line_it->first;
    auto const [elt_idx, elt_line_idx] =
        getElementList().findElementIndex(line_id);
    if (elt_idx == undefIdx)
      continue;
    auto const &elt_ptr = getElementList()[elt_idx];
    if (elt_ptr->IsOutsideLambdaRange(elt_line_idx))
      continue;
    eIdxList.push_back(std::pair(
        elt_idx, elt_line_idx)); // save elt_idx, line_idx for next loop
    validLinesIndex.push_back(i_lineCatalog);
    Amplitudes.push_back(elt_ptr->GetFittedAmplitude(elt_line_idx));
    AmplitudesUncertainties.push_back(
        elt_ptr->GetFittedAmplitudeErrorSigma(elt_line_idx));
  }
  TFloat64List correctedAmplitudes;
  m_CatalogTplRatio->GetBestFit(validLinesIndex, Amplitudes,
                                AmplitudesUncertainties, correctedAmplitudes,
                                m_tplratioBestTplName);

  for (Int32 iValidLine = 0; iValidLine != validLinesIndex.size();
       ++iValidLine) {
    auto const [elt_idx, line_idx] = eIdxList[iValidLine];
    auto const &elt_ptr = getElementList()[elt_idx];

    Float64 const er =
        AmplitudesUncertainties[iValidLine]; // not modifying the fitting error
                                             // for now
    Float64 const nominalAmp = elt_ptr->GetNominalAmplitude(line_idx);
    elt_ptr->SetElementAmplitude(correctedAmplitudes[iValidLine] / nominalAmp,
                                 er);
  }
  getModel().refreshModel();
  return getLeastSquareMerit();
}

void CTplCorrManager::saveResults(Int32 itratio) {
  //  m_tplratioBestTplName = bestTplratioName; done in computeMerit directly
}
