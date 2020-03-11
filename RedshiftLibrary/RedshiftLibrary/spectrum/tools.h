#ifndef _REDSHIFT_SPECTRUM_TOOLS_
#define _REDSHIFT_SPECTRUM_TOOLS_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>

namespace NSEpic
{

class CSpectrumAxis;
class CMask;

/**
 * \ingroup Redshift
 */
class CSpectrumTools
{

public:

    CSpectrumTools( );
    ~CSpectrumTools( );

    static void Interpolate( const CSpectrumAxis& Xorg, const CSpectrumAxis& Yorg, Int32 offsetOrg, Int32 nOrg, const CSpectrumAxis& Xint, CSpectrumAxis& Yint, CMask& Wint );

private:

};


}

#endif
