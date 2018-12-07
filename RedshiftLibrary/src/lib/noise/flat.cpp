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

void CNoiseFlat::SetStatErrorLevel( Float64 level )
{
    if( level <= 0 )
        level = 1;

    m_StatErrorLevel = level;
}

void CNoiseFlat::AddNoise( CSpectrum& s1 ) const
{
    TFloat64List& statError = s1.GetFluxAxis().GetError();
    statError.assign(s1.GetFluxAxis().GetSamplesCount(), m_StatErrorLevel);
}
