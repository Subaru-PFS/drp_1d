#ifndef _REDSHIFT_CONTINUUM_CONTINUUM_
#define _REDSHIFT_CONTINUUM_CONTINUUM_

#include <epic/core/common/managedobject.h>
#include <epic/core/common/datatypes.h>

namespace __NS__
{

class CSpectrum;
class CSpectrumFluxAxis;

class CContinuum : public CManagedObject
{

public:

    CContinuum();
    virtual ~CContinuum();

    virtual Bool RemoveContinuum( const CSpectrum& s, CSpectrumFluxAxis& noContinuumFluxAxis ) = 0;

private:


};


}

#endif
