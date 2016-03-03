#ifndef _REDSHIFT_RAY_RULESTRONGHIGHERTHANWEAK_
#define _REDSHIFT_RAY_RULESTRONGHIGHERTHANWEAK_

#include <epic/core/common/datatypes.h>
#include <boost/format.hpp>
#include <epic/redshift/linemodel/elementlist.h>
#include <epic/redshift/ray/rule.h>

namespace NSEpic
{
  /**
   * \ingroup Redshift
   * Rule to limit lines according to their pairing.
   */
  class CRuleStrongHigherThanWeak : public CRule
  {
  public:
    Bool Check( CLineModelElementList& LineModelElementList );
    void SetUp( Bool EnabledArgument, ... );
  private:
    Int32 m_LineType;
    void Correct( CLineModelElementList& LineModelElementList );
    Float64 FindHighestStrongLineAmp( Int32 linetype , Float64 &er, CLineModelElementList& LineModelElementList );
  };
}

#endif
