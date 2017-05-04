#ifndef _REDSHIFT_RAY_RULESUPERSTRONG_
#define _REDSHIFT_RAY_RULESUPERSTRONG_

#include <epic/core/common/datatypes.h>
#include <boost/format.hpp>
#include <epic/redshift/linemodel/elementlist.h>
#include <epic/redshift/ray/rule.h>

namespace NSEpic
{
  /**
   * \ingroup Redshift
   * Rule to limit 'not super strong' lines to be lower than the super strong ones.
   */
  class CRuleSuperStrong : public CRule
  {
  public:
    Bool Check( CLineModelElementList& LineModelElementList );
    void SetUp( Bool EnabledArgument, ... );
  private:
    Int32 m_LineType;
    TStringList m_SuperStrongTags;
    void Correct( CLineModelElementList& LineModelElementList );
    Float64 FindHighestSuperStrongLineAmp( TStringList superstrongTags, Float64 &er, std::string &name, CLineModelElementList& LineModelElementList );
  };
}

#endif
