#ifndef _REDSHIFT_NOISE_FLAT_
#define _REDSHIFT_NOISE_FLAT_

#include <epic/core/common/datatypes.h>
#include <epic/redshift/noise/noise.h>

namespace NSEpic
{

class CSpectrum;

/**
 * \ingroup Redshift
 */
class CNoiseFlat : public CNoise
{

    DEFINE_MANAGED_OBJECT( CNoiseFlat )

public:

    CNoiseFlat();
    ~CNoiseFlat();

    Void SetStatErrorLevel( Float64 level );

    Bool AddNoise( CSpectrum& s1 ) const; 

private:

    Float64 m_StatErrorLevel;

};


}

#endif
