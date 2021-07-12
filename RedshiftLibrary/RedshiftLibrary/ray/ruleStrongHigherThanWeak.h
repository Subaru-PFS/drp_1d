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
