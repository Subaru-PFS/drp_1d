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
    stream << "#index\tpeak.begin index\tpeak.end index" << std::endl;


    UInt32 nPeaks = EnlargedPeakList.size();
    for( UInt32 j=0; j<nPeaks; j++ ){
        stream << j << "\t";
        stream << EnlargedPeakList[j].GetBegin() << "\t";
        stream << EnlargedPeakList[j].GetEnd() << "\t";
        stream << std::endl;

    }
}

Void CPeakDetectionResult::SaveLine( const COperatorResultStore& store, std::ostream& stream ) const
{
    stream << "not implemented" << std::endl;
}

