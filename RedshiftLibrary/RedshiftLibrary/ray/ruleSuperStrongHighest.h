#ifndef _REDSHIFT_RAY_RULESUPERSTRONG_
#define _REDSHIFT_RAY_RULESUPERSTRONG_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/linemodel/elementlist.h>
#include <RedshiftLibrary/ray/rule.h>
#include <boost/format.hpp>

namespace NSEpic {
/**
 * \ingroup Redshift
 * Rule to limit 'not super strong' lines to be lower than the super strong
 * ones.
 */
class CRuleSuperStrong : public CRule
{
  public:
    CRuleSuperStrong();
    Bool Check(CLineModelElementList &LineModelElementList);
    void SetUp(Bool EnabledArgument, ...);

  private:
    Int32 m_LineType = 0;
    TStringList m_SuperStrongTags;
    void Correct(CLineModelElementList &LineModelElementList);
    Float64
    FindHighestSuperStrongLineAmp(TStringList superstrongTags, Float64 &er,
                                  std::string &name,
                                  CLineModelElementList &LineModelElementList);
};
} // namespace NSEpic

#endif
