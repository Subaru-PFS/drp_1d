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
#ifndef _REDSHIFT_RAY_RULESTRONGHIGHERTHANWEAK_
#define _REDSHIFT_RAY_RULESTRONGHIGHERTHANWEAK_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/linemodel/elementlist.h"
#include "RedshiftLibrary/ray/rule.h"
#include <boost/format.hpp>

namespace NSEpic {
/**
 * \ingroup Redshift
 * Rule to limit lines according to their pairing.
 */
class CRuleStrongHigherThanWeak : public CRule
{
  public:
    CRuleStrongHigherThanWeak();
    ~CRuleStrongHigherThanWeak();

    Bool Check(CLineModelElementList &LineModelElementList);
    void SetUp(Bool EnabledArgument, ...);

  private:
    Int32 m_LineType;
    void Correct(CLineModelElementList &LineModelElementList);
    Float64
    FindHighestStrongLineAmp(Int32 linetype, Float64 &er, std::string &name,
                             CLineModelElementList &LineModelElementList);
};
} // namespace NSEpic

#endif
