#include <RedshiftLibrary/spectrum/axis.h>

#include <RedshiftLibrary/debug/assert.h>

#include <algorithm>

using namespace NSEpic;
using namespace std;

CSpectrumAxis::CSpectrumAxis()
{

}

CSpectrumAxis::CSpectrumAxis( UInt32 n ) :
    m_Samples( n )
{
    for( UInt32 i=0; i<n; i++ )
    {
        m_Samples[i] = 0.0;
    }
}

CSpectrumAxis::CSpectrumAxis( const Float64* samples, UInt32 n ) :
    m_Samples( n )
{
    for( UInt32 i=0; i<n; i++ )
    {
        m_Samples[i] = samples[i];
    }
}

CSpectrumAxis::~CSpectrumAxis()
{

}

CSpectrumAxis& CSpectrumAxis::operator=(const CSpectrumAxis& other)
{
    m_Samples = other.m_Samples;
    return *this;
}

void CSpectrumAxis::SetSize( UInt32 s )
{
    m_Samples.resize( s );
}


