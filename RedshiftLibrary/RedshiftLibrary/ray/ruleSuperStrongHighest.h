#ifndef _REDSHIFT_RAY_RULESUPERSTRONG_
#define _REDSHIFT_RAY_RULESUPERSTRONG_

#include <RedshiftLibrary/common/datatypes.h>
#include <boost/format.hpp>
#include <RedshiftLibrary/linemodel/elementlist.h>
#include <RedshiftLibrary/ray/rule.h>

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
