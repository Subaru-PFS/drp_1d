#ifndef _REDSHIFT_RAY_RULEOIIRATIORANGE_
#define _REDSHIFT_RAY_RULEOIIRATIORANGE_

#include <epic/core/common/datatypes.h>
#include <boost/format.hpp>
#include <epic/redshift/linemodel/elementlist.h>
#include <epic/redshift/ray/rule.h>

namespace NSEpic
{
  /**
   * \ingroup Redshift
   */
  class CRuleOIIRatioRange : public CRule
  {
  public:
    Bool Check( CLineModelElementList& LineModelElementList );
    void SetUp( Bool EnabledArgument, ... );
  private:
    Int32 m_LineType;
    std::string m_LineA;
    std::string m_LineB;
    Float64 m_Coefficient;
    void Correct( CLineModelElementList& LineModelElementList );
  };
}

#endif
