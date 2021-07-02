#ifndef _REDSHIFT_CONTINUUM_CONTINUUM_
#define _REDSHIFT_CONTINUUM_CONTINUUM_

#include "RedshiftLibrary/common/datatypes.h"

namespace NSEpic
{

class CSpectrum;
class CSpectrumFluxAxis;

/**
 * \ingroup Redshift
 * Common ancestral class for continuum estimators.
 **/
class CContinuum
{

public:

    CContinuum();
    virtual ~CContinuum();

    virtual Bool RemoveContinuum( const CSpectrum& s, CSpectrumFluxAxis& noContinuumFluxAxis ) = 0;

private:


};


}

#endif
