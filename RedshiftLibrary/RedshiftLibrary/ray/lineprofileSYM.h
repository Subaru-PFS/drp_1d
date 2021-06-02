#ifndef _REDSHIFT_LINE_PROFILE_SYM_
#define _REDSHIFT_LINE_PROFILE_SYM_
#include <string>
#include <math.h>
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/ray/lineprofile.h"
namespace NSEpic
{
    /**
     * \ingroup Redshift
     */
    class CLineProfileSYM : public CLineProfile
    {
        public:
            CLineProfileSYM(const Float64 nsigmasupport=8.0);
            Float64 GetLineProfile(Float64 x, Float64 x0, const Float64 sigma) override;
            Float64 GetLineFlux( Float64 A, const Float64 sigma ) override;
            Float64 GetLineProfileDerivZ(Float64 x, Float64 lambda0, Float64 redshift, const Float64 sigma) override;
            Float64 GetLineProfileDerivSigma(Float64 x, Float64 x0, const Float64 sigma) override;
    };
}
#endif