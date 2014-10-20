#include <epic/redshift/spectrum/axis.h>

#include <epic/core/debug/assert.h>
#include <epic/core/serializer/serializer.h>

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

Void CSpectrumAxis::SetSize( UInt32 s )
{
    m_Samples.resize( s );
}


Bool CSpectrumAxis::Serialize( CSerializer& ar )
{
    Int16 version = 1;

    if( ar.BeginScope( "Axis", version ) == version )
    {
        ar.Serialize( m_Samples, "Samples" );
        ar.EndScope();

        return true;
    }

    return false;
}
