#ifndef _REDSHIFT_SPECTRUM_NOISEAXIS_
#define _REDSHIFT_SPECTRUM_NOISEAXIS_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/spectrum/axis.h"

namespace NSEpic
{
/**
 * \ingroup Redshift
 */
class  CSpectrumNoiseAxis : public CSpectrumAxis
{

public:
    using CSpectrumAxis::CSpectrumAxis;

    CSpectrumNoiseAxis() = default;
    explicit  CSpectrumNoiseAxis( UInt32 n );
    CSpectrumNoiseAxis( const Float64* samples, UInt32 n );

    void                SetSize( UInt32 s, Float64 valueDef = 1.0);
    Bool                Invert();

private:

};

}

#endif
