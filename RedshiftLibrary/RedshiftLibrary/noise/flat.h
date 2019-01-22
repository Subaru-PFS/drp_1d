#ifndef _REDSHIFT_NOISE_FLAT_
#define _REDSHIFT_NOISE_FLAT_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/noise/noise.h>

namespace NSEpic
{

class CSpectrum;

/**
 * \ingroup Redshift
 */
class CNoiseFlat : public CNoise
{

public:

    CNoiseFlat();
    ~CNoiseFlat();

    void SetStatErrorLevel( Float64 level );

    void AddNoise( CSpectrum& s1 ) const;

private:

    Float64 m_StatErrorLevel;

};


}

#endif
