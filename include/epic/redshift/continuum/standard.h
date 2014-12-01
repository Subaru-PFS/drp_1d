#ifndef _REDSHIFT_CONTINUUM_STANDARD_
#define _REDSHIFT_CONTINUUM_STANDARD_

#include <epic/redshift/continuum/continuum.h>

namespace NSEpic
{

class CContinuumStandard : public CContinuum
{

public:

    CContinuumStandard();
    ~CContinuumStandard();

    Bool Remove( CSpectrum& s );

private:



};


}

#endif