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

#include "RedshiftLibrary/linemodel/ratiotofreemanager.h"
#include "RedshiftLibrary/common/size.h"
#include "RedshiftLibrary/operator/linemodelresult.h"

using namespace NSEpic;

CRatioToFreeManager::CRatioToFreeManager(
    const std::shared_ptr<CLMEltListVector> &elementsVector,
    const CSpcModelVectorPtr &models, const CCSpectrumVectorPtr &inputSpcs,
    const CTLambdaRangePtrVector &lambdaRanges,
    std::shared_ptr<CContinuumManager> continuumManager,
    const CLineMap &restLineList, const CSpectraGlobalIndex &spcIndex)
    : CLineRatioManager(elementsVector, models, inputSpcs, lambdaRanges,
                        continuumManager, restLineList, spcIndex),
      CRulesManager(elementsVector, models, inputSpcs, lambdaRanges,
                    continuumManager, restLineList, spcIndex),
      CTplratioManager(elementsVector, models, inputSpcs, lambdaRanges,
                       continuumManager, restLineList, spcIndex){};
void CRatioToFreeManager::setPassMode(Int32 iPass) {
  CLineRatioManager::setPassMode(iPass);
  if (m_pass != 1 && m_pass != 2)
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "m_pass must be 1 or 2, got " << m_pass);
  m_pass = iPass;
}

int CRatioToFreeManager::prepareFit(Float64 redshift) {
  int prepare = 0;
  if (m_pass == 1)
    prepare = CTplratioManager::prepareFit(redshift);
  else if (m_pass == 2)
    prepare = CLineRatioManager::prepareFit(redshift);
  return prepare;
}

std::pair<Float64, Float64> CRatioToFreeManager::computeMerit(Int32 itratio) {
  std::pair<Float64, Float64> merit;
  if (m_pass == 1) {
    merit = CTplratioManager::computeMerit(itratio);
  } else if (m_pass == 2) {
    merit = CRulesManager::computeMerit(itratio);
  }
  return merit;
}

void CRatioToFreeManager::saveResults(Int32 itratio) {
  if (m_pass == 1)
    CTplratioManager::saveResults(itratio);
  else if (m_pass == 2)
    CLineRatioManager::saveResults(itratio);
};

void CRatioToFreeManager::setChiSquareRatioResult(
    const Int32 index_z, const std::shared_ptr<CLineModelResult> &lmResult) {
  if (m_pass == 1) {
    CTplratioManager::setChiSquareRatioResult(index_z, lmResult);
  } else if (m_pass == 2) {
    if (GetChisquareTplratio().size() < 1)
      return;

    if (index_z >= ssize(lmResult->Redshifts))
      THROWG(ErrorCode::INTERNAL_ERROR, "Invalid z index");
    auto const &nRatios = getTplratio_count();
    for (Int32 k = 0; k < nRatios; k++) {
      lmResult->ChiSquareTplratios[k][index_z] = GetChisquareTplratio()[0];
      lmResult->ScaleMargCorrectionTplratios[k][index_z] =
          GetScaleMargTplratio()[0];
      lmResult->StrongELPresentTplratios[k][index_z] =
          GetStrongELPresentTplratio()[0];
      lmResult->StrongHalphaELPresentTplratios[k][index_z] =
          getHaELPresentTplratio()[0];
      lmResult->NLinesAboveSNRTplratios[k][index_z] =
          GetNLinesAboveSNRTplratio()[0];
      lmResult->PriorLinesTplratios[k][index_z] = GetPriorLinesTplratio()[0];
    }
  }
}