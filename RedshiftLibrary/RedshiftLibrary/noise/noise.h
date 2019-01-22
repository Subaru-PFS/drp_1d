#ifndef _REDSHIFT_NOISE_NOISE_
#define _REDSHIFT_NOISE_NOISE_

#include <RedshiftLibrary/common/datatypes.h>

namespace NSEpic
{

class CSpectrum;

/**
 * \ingroup Redshift
 */
class CNoise
{

public:

    CNoise();
    virtual ~CNoise();

    virtual void AddNoise( CSpectrum& s1 ) const = 0;

private:


};


}

#endif
