#include <RedshiftLibrary/spectrum/axis.h>

#include <RedshiftLibrary/debug/assert.h>
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

Int32 CSpectrumAxis::MaskAxis(const CSpectrumAxis& other, TFloat64List& mask)//mask is 0. or 1.
{
    UInt32 sum = std::accumulate(other.m_Samples.begin(), other.m_Samples.end(), 1.);
    m_Samples.reserve(sum);
    for(Int32 i = 0; i < other.GetSamplesCount(); i++){
        if(mask[i]==1.)
            m_Samples.push_back(other.m_Samples[i]);
    }
    return 0;
}