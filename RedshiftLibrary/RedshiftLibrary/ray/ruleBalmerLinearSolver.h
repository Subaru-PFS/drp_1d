#ifndef _REDSHIFT_RAY_RULEBALMERLINEARSOLVER_
#define _REDSHIFT_RAY_RULEBALMERLINEARSOLVER_

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
  class CRuleBalmerLinearSolver : public CRule
  {
  public:
    Bool Check( CLineModelElementList& LineModelElementList );
    void SetUp( Bool EnabledArgument, ... );
  private:
    void Correct( CLineModelElementList& LineModelElementList );
    TFloat64List BalmerModelLinSolve( std::vector<Float64> lambdax, std::vector<Float64> continuumx, std::vector<Float64> datax, std::vector<Float64> errdatax );
  };
}

#endif
