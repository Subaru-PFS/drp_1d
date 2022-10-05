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
#include "RedshiftLibrary/linemodel/linemodelsolution.h"

using namespace NSEpic;

CLineModelSolution::CLineModelSolution() {
  this->m_type = "CLineModelSolution";
}

CLineModelSolution::CLineModelSolution(
    const CLineCatalog::TLineVector &restLineList) {
  this->m_type = "CLineModelSolution";
  // filling lineIds
  lineId.assign(restLineList.size(), -1);
  for (Int32 i = 0; i < restLineList.size(); i++)
    lineId[i] = restLineList[i].GetID();

  Int32 s = lineId.size();
  ElementId.assign(s, undefIdx);
  Amplitudes.assign(s, NAN);
  AmplitudesUncertainties.assign(s, NAN);
  FittingError.assign(s, NAN);
  LambdaObs.assign(s, NAN);
  Offset.assign(s, NAN);
  Velocity.assign(s, NAN);
  CenterContinuumFlux.assign(s, NAN);
  ContinuumError.assign(s, NAN);
  Sigmas.assign(s, NAN);
  Fluxs.assign(s, NAN);
  FluxErrors.assign(s, NAN);
  FluxDirectIntegration.assign(s, NAN);
  FluxDirectIntegrationError.assign(s, NAN);
  OutsideLambdaRange.assign(s, true);
  fittingGroupInfo.assign(s, "undefined");
  continuum_pCoeff0.assign(s, NAN);
  continuum_pCoeff1.assign(s, NAN);
  continuum_pCoeff2.assign(s, NAN);
}

bool CLineModelSolution::isLineValid(Int32 lineIdx) const {
  if (!lineId.size())
    THROWG(INTERNAL_ERROR, "lineModelSolution is empty");
  return !OutsideLambdaRange[lineIdx] & Amplitudes[lineIdx] > 0.0;
}