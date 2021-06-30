#ifndef _REDSHIFT_LINEMODEL_MODELFITTINGRESULT_
#define _REDSHIFT_LINEMODEL_MODELFITTINGRESULT_


#include "RedshiftLibrary/processflow/result.h"
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/operator/operator.h"

#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/operator/linemodelresult.h"


namespace NSEpic
{

  /**
   * \ingroup Redshift
   */
#include "RedshiftLibrary/linemodel/modelfittingresult.i"
inline
const CLineModelSolution& CModelFittingResult::GetLineModelSolution() const
{
    return LineModelSolution;
}


}

#endif
