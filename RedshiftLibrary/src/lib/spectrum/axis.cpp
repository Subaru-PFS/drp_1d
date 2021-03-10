#include <RedshiftLibrary/spectrum/axis.h>
#include <RedshiftLibrary/log/log.h>
#include <numeric>
#include <algorithm>

using namespace NSEpic;
using namespace std;


CSpectrumAxis::CSpectrumAxis( UInt32 n, Float64 value) :
    m_Samples( n , value)
{
}

CSpectrumAxis::CSpectrumAxis( const Float64* samples, UInt32 n ) :
    m_Samples( n )
{
    for( UInt32 i=0; i<n; i++ )
    {
        m_Samples[i] = samples[i];
    }
}

CSpectrumAxis& CSpectrumAxis::operator*=(const Float64 op)
{
    for(Int32 i = 0; i < m_Samples.size(); i++){
        m_Samples[i] *= op;
    }
    return *this;
}

void CSpectrumAxis::SetSize( UInt32 s )
{
    m_Samples.resize( s );
}

UInt32 CSpectrumAxis::GetSamplesCount( )
{
    return m_Samples.size();
}

Int32 CSpectrumAxis::extractFrom(const CSpectrumAxis& other, Int32 startIdx, Int32 endIdx)
{
    m_Samples.resize(endIdx-startIdx +1);
    for(Int32 i = startIdx; i < endIdx + 1; i++){
        m_Samples[i - startIdx] = other.m_Samples[i];
    }
    return 0;
}
/*
    maskedAxis is the output axis after applying the mask on the current object
*/
Int32 CSpectrumAxis::MaskAxis(TFloat64List& mask, CSpectrumAxis& maskedAxis)const//mask is 0. or 1.
{
    if(mask.size()!= GetSamplesCount()){
        Log.LogError("CSpectrumAxis::MaskAxis: mask and axis sizes are not equal. Abort");
        throw runtime_error("CSpectrumAxis::MaskAxis: mask and axis sizes are not equal. Abort");
    }
    UInt32 sum = UInt32(std::accumulate(m_Samples.begin(), m_Samples.end(), 0));
    maskedAxis.m_Samples.clear();
    maskedAxis.m_Samples.reserve(sum);
    for(Int32 i = 0; i < GetSamplesCount(); i++){
        if(mask[i]==1.)
            maskedAxis.m_Samples.push_back(m_Samples[i]);
    }
    return 0;
}