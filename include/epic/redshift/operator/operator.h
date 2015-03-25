#ifndef _REDSHIFT_OPERATOR_OPERATOR_
#define _REDSHIFT_OPERATOR_OPERATOR_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>

#include <vector>

namespace NSEpic
{

class CSpectrum;
class CTemplate;
class CRedshifts;

class COperator
{


public:


    enum EStatus
    {
        nStatus_OK = 0,
        nStatus_DataError,
        nStatus_NoOverlap
    };

    typedef std::vector<EStatus> TStatusList;

    COperator();
    virtual ~COperator();

    virtual Bool Compute( const CSpectrum& spectrum, const CTemplate& tpl,
                          const TFloat64Range& lambdaRange, const CRedshifts& redshifts, Float64 overlapThreshold ) = 0;


    const TStatusList& GetStatus() const;
    const TFloat64List&         GetOverlap() const;
    const TFloat64List&         GetResults() const;

protected:

    TFloat64List            m_Result;
    TFloat64List            m_Overlap;
    std::vector<EStatus>    m_Status;

};


}

#endif
