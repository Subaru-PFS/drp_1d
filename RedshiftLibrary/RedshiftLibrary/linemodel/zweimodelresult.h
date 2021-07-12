#ifndef _REDSHIFT_LINEMODEL_ZWEIMODELRESULT_
#define _REDSHIFT_LINEMODEL_ZWEIMODELRESULT_


#include "RedshiftLibrary/processflow/result.h"
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/operator/operator.h"

namespace NSEpic
{

  /**
   * \ingroup Redshift
   */
class CZweiModelResult : public COperatorResult
{

public:

    CZweiModelResult(std::vector<Float64> redshifts_s1, std::vector<Float64> redshifts_s2, std::vector<std::vector<Float64>> combined_merits);
    CZweiModelResult();
    virtual ~CZweiModelResult();


private:
    std::vector<Float64> m_redshifts_s1;
    std::vector<Float64> m_redshifts_s2;
    std::vector<std::vector<Float64>> m_combined_merits;

};


}

#endif
