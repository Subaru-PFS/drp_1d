#ifndef _REDSHIFT_LINE_PROFILE_ASYMFIT_
#define _REDSHIFT_LINE_PROFILE_ASYMFIT_
#include <string>
#include <math.h>
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/ray/lineprofileASYM.h"

namespace NSEpic
{
    /**
     * \ingroup Redshift
     */
    class CLineProfileASYMFIT: public CLineProfileASYM
    {
        public:
            CLineProfileASYMFIT(const Float64 nsigmasupport = 8.0, 
                                TAsymParams params = {2., 2., 0.},
                                const std::string centeringMethod = "mean");
            Bool    isAsymFit() override;
            Bool    isAsymFixed() override;
            void    SetAsymParams(TAsymParams params);
            void resetAsymFitParams() override;

    };
}
#endif