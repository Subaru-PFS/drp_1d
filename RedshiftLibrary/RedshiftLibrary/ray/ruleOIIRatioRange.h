#ifndef _REDSHIFT_RAY_RULEOIIRATIORANGE_
#define _REDSHIFT_RAY_RULEOIIRATIORANGE_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/linemodel/elementlist.h>
#include <RedshiftLibrary/ray/rule.h>
#include <boost/format.hpp>

namespace NSEpic {
/**
 * \ingroup Redshift
 */
class CRuleRatioRange : public CRule
{
  public:
    CRuleRatioRange();
    ~CRuleRatioRange();
    Bool Check(CLineModelElementList &LineModelElementList);
    void SetUp(Bool EnabledArgument, ...);

  private:
    Int32 m_LineType;
    std::string m_LineA;
    std::string m_LineB;
    Float64 m_Coefficient;
    void Correct(CLineModelElementList &LineModelElementList);
};
} // namespace NSEpic

#endif
