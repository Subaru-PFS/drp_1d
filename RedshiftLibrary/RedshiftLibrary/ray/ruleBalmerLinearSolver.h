#ifndef _REDSHIFT_RAY_RULEBALMERLINEARSOLVER_
#define _REDSHIFT_RAY_RULEBALMERLINEARSOLVER_

#include <RedshiftLibrary/common/datatypes.h>
#include <boost/format.hpp>
#include <RedshiftLibrary/linemodel/elementlist.h>
#include <RedshiftLibrary/ray/rule.h>

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
