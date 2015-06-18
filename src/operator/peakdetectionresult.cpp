#include <epic/redshift/operator/peakdetectionresult.h>

#include <epic/redshift/ray/ray.h>

using namespace NSEpic;

IMPLEMENT_MANAGED_OBJECT( CPeakDetectionResult )

CPeakDetectionResult::CPeakDetectionResult()
{

}

CPeakDetectionResult::~CPeakDetectionResult()
{

}

Void CPeakDetectionResult::Save( const COperatorResultStore& store, std::ostream& stream ) const
{
    stream << "not implemented" << std::endl;
}

Void CPeakDetectionResult::SaveLine( const COperatorResultStore& store, std::ostream& stream ) const
{
    stream << "not implemented" << std::endl;
}

