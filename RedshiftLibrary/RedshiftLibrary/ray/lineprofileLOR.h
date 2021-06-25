#ifndef _REDSHIFT_LINE_PROFILE_LOR_
#define _REDSHIFT_LINE_PROFILE_LOR_
#include <string>
#include "RedshiftLibrary/log/log.h"
#include <math.h>
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/ray/lineprofile.h"
namespace NSEpic
{
    /**
     * \ingroup Redshift
     */
    class CLineProfileLOR: public CLineProfile
    {
        public:
            CLineProfileLOR(const Float64 nsigmasupport = 8.0);
            Float64 GetLineProfile(Float64 x, Float64 x0, const Float64 sigma) override; //override is optional, but a good practice with c++11
            Float64 GetLineFlux( Float64 A, const Float64 sigma) override;
            Float64 GetLineProfileDerivZ(Float64 x, Float64 lambda0, Float64 redshift, const Float64 sigma) override;
            Float64 GetLineProfileDerivSigma(Float64 x, Float64 x0, const Float64 sigma) override;
            Float64 GetNSigmaSupport() override;
    };
}
#endif