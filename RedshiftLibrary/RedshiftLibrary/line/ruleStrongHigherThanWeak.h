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
#ifndef _REDSHIFT_LINE_RULESTRONGHIGHERTHANWEAK_
#define _REDSHIFT_LINE_RULESTRONGHIGHERTHANWEAK_

#include <boost/format.hpp>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/line/rule.h"
#include "RedshiftLibrary/linemodel/elementlist.h"

class RuleStrongHigherThanWeak_fixture;

namespace NSEpic {

/**
 * \ingroup Redshift
 * Rule to limit lines according to their pairing.
 */
class CRuleStrongHigherThanWeak : public CRule {
public:
  CRuleStrongHigherThanWeak();

  bool Check(CLMEltListVector &LineModelElementList);
  void SetUp(bool EnabledArgument, ...);

private:
  friend class ::RuleStrongHigherThanWeak_fixture;
  CLine::EType m_LineType;
  void Correct(CLMEltListVector &LineModelElementList);
  Float64 FindHighestStrongLineAmp(Int32 linetype, Float64 &er,
                                   std::string &name,
                                   CLMEltListVector &LineModelElementList);
  std::pair<Int32, Int32>
  FindLowestStrongLineIndex(const CLMEltListVector &LineModelElementList);
  void constructLogMsg(const std::string &nameWeak,
                       const std::string &strongName, Float64 ampB,
                       Float64 maxB);
  Float64 maxAmplitude(Float64 ampStrong, Float64 erStrong);
  void correctLineModelElement(CLineModelElement &element,
                               const TLineModelElementParam_ptr &elt_param_ptr,
                               Float64 maxAmplitude,
                               const std::string &nameStrong);
};
} // namespace NSEpic

#endif
