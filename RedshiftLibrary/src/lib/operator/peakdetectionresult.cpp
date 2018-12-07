#include <RedshiftLibrary/operator/peakdetectionresult.h>

#include <RedshiftLibrary/ray/ray.h>

using namespace NSEpic;

CPeakDetectionResult::CPeakDetectionResult()
{

}

CPeakDetectionResult::~CPeakDetectionResult()
{

}

void CPeakDetectionResult::Save( const CDataStore& store, std::ostream& stream ) const
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

void CPeakDetectionResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    stream << "not implemented" << std::endl;
}

