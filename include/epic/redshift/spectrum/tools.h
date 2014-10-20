#ifndef _REDSHIFT_SPECTRUM_TOOLS_
#define _REDSHIFT_SPECTRUM_TOOLS_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>

namespace NSEpic
{

class CSpectrumAxis;
class CMask;

class CSpectrumTools
{

public:

    CSpectrumTools( );
    ~CSpectrumTools( );

    static Void Interpolate( const CSpectrumAxis& Xorg, const CSpectrumAxis& Yorg, Int32 offsetOrg, Int32 nOrg, const CSpectrumAxis& Xint, CSpectrumAxis& Yint, CMask& Wint );

private:

};


}

#endif
