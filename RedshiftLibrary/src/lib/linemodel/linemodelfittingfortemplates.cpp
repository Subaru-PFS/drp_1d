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

#include "RedshiftLibrary/linemodel/linemodelfittingfortemplates.h"

using namespace NSEpic;
using namespace std;

CLineModelFittingForTemplates::CLineModelFittingForTemplates(
    const std::shared_ptr<const CTemplate> &template_,
    const TLambdaRange &lambdaRange,
    const std::shared_ptr<COperatorContinuumFitting> &continuumFittingOperator)
    : CLineModelFitting(1) {
  m_inputSpcs =
      std::make_shared<std::vector<std::shared_ptr<const CSpectrum>>>();

  m_inputSpcs->push_back(template_);
  m_lambdaRanges.push_back(std::make_shared<const TLambdaRange>(lambdaRange));
  initParameters();
  // override ortho specific parameters
  m_fittingmethod = "hybrid";
  m_enableAmplitudeOffsets = false;
  m_enableLbdaOffsets = false;

  auto lineRatioType = CLineRatioManager::EType::rules;
  initMembers(continuumFittingOperator, lineRatioType,
              ElementComposition::Default);
  // temporary options override to be removed when full tpl ortho is implemented
  setLineRatioManager(lineRatioType);

  dynamic_cast<CRulesManager *>(m_lineRatioManager.get())->setRulesOption("no");
  setContinuumComponent(TContinuumComponent("fromSpectrum"));
}