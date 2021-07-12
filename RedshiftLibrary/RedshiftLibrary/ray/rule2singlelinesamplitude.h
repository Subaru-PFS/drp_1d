#ifndef _REDSHIFT_RAY_RULE2SINGLELINESAMPLITUDE_
#define _REDSHIFT_RAY_RULE2SINGLELINESAMPLITUDE_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/linemodel/elementlist.h"
#include "RedshiftLibrary/ray/rule.h"
#include <boost/format.hpp>

namespace NSEpic {
/**
 * \ingroup Redshift
 * Rule to limit lines according to their pairing.
 */
class CRule2SingleLinesAmplitude : public CRule
{
  public:
    CRule2SingleLinesAmplitude();
    ~CRule2SingleLinesAmplitude();
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
