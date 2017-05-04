#ifndef _REDSHIFT_RAY_RULE2SINGLELINESAMPLITUDE_
#define _REDSHIFT_RAY_RULE2SINGLELINESAMPLITUDE_

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
  class CRule2SingleLinesAmplitude : public CRule
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
