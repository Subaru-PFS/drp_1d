#ifndef _REDSHIFT_OPERATOR_CORRELATION_
#define _REDSHIFT_OPERATOR_CORRELATION_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>
#include <epic/redshift/common/redshifts.h>
#include <epic/redshift/spectrum/fluxaxis.h>
#include <epic/redshift/spectrum/spectralaxis.h>
#include <epic/redshift/operator/operator.h>

#include <string>

#include <boost/date_time/posix_time/posix_time.hpp>

namespace NSEpic
{

class CSpectrum;
class CTemplate;

class COperatorCorrelation : public COperator
{

public:

    enum EStatus
    {
        nStatus_OK = 0,
        nStatus_AboveOverlapThreshold,
        nStatus_InvalidInputData,
        nStatus_Count
    };

    typedef std::vector<EStatus> TStatusList;

    COperatorCorrelation();
    ~COperatorCorrelation();

    Bool Compute( const CSpectrum& s1, const CTemplate& s2, const TFloat64Range& r, const CRedshifts& redhisfts, Float64 overlap );

    const TFloat64List& GetResults() const;
    const TStatusList& GetStatus() const;
    Float64 GetComputationDuration() const;

private:

    TFloat64List            m_Correlation;
    TFloat64List            m_Overlap;
    TStatusList             m_Status;
    Float64                 m_TotalDuration;
};


}

#endif
