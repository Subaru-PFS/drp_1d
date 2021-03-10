#include <RedshiftLibrary/spectrum/noiseaxis.h>
#include <RedshiftLibrary/common/mask.h>

#include <math.h>
#include <RedshiftLibrary/log/log.h>

using namespace NSEpic;
using namespace std;


CSpectrumNoiseAxis::CSpectrumNoiseAxis( UInt32 n ) :
    CSpectrumAxis( n, 1.0) 
{

}

CSpectrumNoiseAxis::CSpectrumNoiseAxis( const Float64* samples, UInt32 n ) :
    CSpectrumAxis( samples, n ) 
{

}



Bool CSpectrumNoiseAxis::Invert()
{
    Int32 N = GetSamplesCount();
    for( UInt32 i=0; i<N; i++ )
    {
        m_Samples[i] = 1/m_Samples[i];
    }
    return true;
}

//default to 1. instead of O.
void CSpectrumNoiseAxis::SetSize( UInt32 s, Float64 valueDef)
{
    m_Samples.assign( s, valueDef);
}