#ifndef _REDSHIFT_NOISE_NOISE_
#define _REDSHIFT_NOISE_NOISE_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/managedobject.h>

namespace NSEpic
{

class CSpectrum;

class CNoise : public CManagedObject
{

public:

    CNoise();
    virtual ~CNoise();

    virtual Bool AddNoise( CSpectrum& s1 ) const = 0; 

private:


};


}

#endif
