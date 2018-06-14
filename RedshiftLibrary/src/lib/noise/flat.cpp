#include <RedshiftLibrary/noise/flat.h>

#include <RedshiftLibrary/spectrum/axis.h>
#include <RedshiftLibrary/spectrum/spectrum.h>

using namespace NSEpic;
using namespace std;

CNoiseFlat::CNoiseFlat() :
    m_StatErrorLevel( 1.0 )
{
}

CNoiseFlat::~CNoiseFlat()
{
}

Void CNoiseFlat::SetStatErrorLevel( Float64 level )
{
    if( level <= 0 )
        level = 1;

    m_StatErrorLevel = level;
}

Void CNoiseFlat::AddNoise( CSpectrum& s1 ) const
{
    Float64* statError = s1.GetFluxAxis().GetError();

    for( UInt32 i=0; i<s1.GetFluxAxis().GetSamplesCount(); i++ )
    {
        statError[i] = m_StatErrorLevel;
    }
}
