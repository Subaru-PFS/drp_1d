#ifndef _REDSHIFT_LINE_PROFILE_ASYMFIT_
#define _REDSHIFT_LINE_PROFILE_ASYMFIT_
#include <string>
#include <math.h>
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/ray/lineprofileASYMFIXED.h"
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

            CLineProfileASYMFIT(const CLineProfileASYMFIT & other) = default; 
            CLineProfileASYMFIT(CLineProfileASYMFIT && other) = default; 
            CLineProfileASYMFIT& operator=(const CLineProfileASYMFIT& other) = default;  
            CLineProfileASYMFIT& operator=(CLineProfileASYMFIT&& other) = default; 
            ~CLineProfileASYMFIT()=default; 
    };
inline
CLineProfileASYMFIT::CLineProfileASYMFIT(const Float64 nsigmasupport, TAsymParams params, const std::string centeringMethod):
CLineProfileASYM(ASYMFIT, nsigmasupport, params, centeringMethod)
{

}

inline
void CLineProfileASYMFIT::SetAsymParams(TAsymParams params)
{
    m_asym_sigma_coeff = params.sigma;
    m_asym_alpha = params.alpha;
    m_asym_delta = params.delta;
}
inline
Bool CLineProfileASYMFIT::isAsymFit(){
    return 1;
}
inline
Bool CLineProfileASYMFIT::isAsymFixed(){
    return 0;
}
}
#endif