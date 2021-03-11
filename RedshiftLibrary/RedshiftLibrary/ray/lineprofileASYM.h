#ifndef _REDSHIFT_LINE_PROFILE_ASYM_
#define _REDSHIFT_LINE_PROFILE_ASYM_
#include <string>
#include <math.h>
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/ray/lineprofile.h"
namespace NSEpic
{
    /**
     * \ingroup Redshift
     */
    class CLineProfileASYM: public CLineProfile
    {
        public:
            CLineProfileASYM(const Float64 nsigmasupport = 8.0);
            ~CLineProfileASYM()=default; 
            Float64 GetLineProfile(Float64 x, Float64 x0, const Float64 sigma) override;
            Float64 GetLineFlux(Float64 A, const Float64 sigma) override;
            Float64 GetLineProfileDerivZ(Float64 x, Float64 lambda0, Float64 redshift, const Float64 sigma) override;
            Float64 GetLineProfileDerivSigma(Float64 x, Float64 x0, const Float64 sigma) override;
            Float64 GetNSigmaSupport() override;
            TFloat64List GetLineProfileVector() override;//equivalent to ::computeKernel

            CLineProfileASYM(const CLineProfileASYM & other) = default; 
            CLineProfileASYM(CLineProfileASYM && other) = default; 
            CLineProfileASYM& operator=(const CLineProfileASYM& other) = default;  
            CLineProfileASYM& operator=(CLineProfileASYM&& other) = default; 
        private:
            Float64 m_asym_sigma_coeff = 1.0;
            Float64 m_asym_alpha = 4.5;
    };
inline
CLineProfileASYM::CLineProfileASYM(const Float64 nsigmasupport):
CLineProfile(nsigmasupport, "ASYM")
{
}
inline
Float64 CLineProfileASYM::GetLineProfile(Float64 x, Float64 x0, Float64 sigma)
{
    Float64 xc = x-x0;
    Float64 val = 0.0;
    Float64 xsurc;

    sigma = sigma*m_asym_sigma_coeff;
    xsurc = xc/sigma;
    val = exp(-0.5*xsurc*xsurc)*(1.0+erf(m_asym_alpha/sqrt(2.0)*xsurc));
    return val;
}
inline
Float64 CLineProfileASYM::GetNSigmaSupport(){
    return m_nsigmasupport*m_asym_sigma_coeff;
}
inline
Float64 CLineProfileASYM::GetLineFlux( Float64 A , Float64 sigma){
    return A*sigma*sqrt(2*M_PI);
}
inline
Float64 CLineProfileASYM::GetLineProfileDerivZ(Float64 x, Float64 lambda0, Float64 redshift, Float64 sigma)
{
    Float64 xc = x-lambda0*(1+redshift);
    Float64 val=0.0;
    Float64 xsurc, xsurc2,  xcd;

    sigma = sigma*m_asym_sigma_coeff;
    xsurc = xc/sigma;
    xsurc2 = xsurc*xsurc;
    val = lambda0 /sigma * xsurc * exp(-0.5*xsurc2) *(1.0+erf(m_asym_alpha/sqrt(2.0)*xsurc)) -m_asym_alpha * lambda0 /sqrt(2*M_PI) /sigma * exp(-(1+m_asym_alpha*m_asym_alpha)/2 * xsurc2);
    return val;
}

inline
Float64 CLineProfileASYM::GetLineProfileDerivSigma(Float64 x, Float64 x0, Float64 sigma)
{
    Float64 val=0.0;
    Float64 xc = x-x0;
    Float64 xsurc, xsurc2, alpha2, valsym, valsymd, valasym, valasymd;

    sigma = sigma*m_asym_sigma_coeff;
    xsurc = xc/sigma;
    xsurc2 = xsurc*xsurc;

    valsym = exp(-0.5*xsurc2);
    valsymd = m_asym_sigma_coeff * xsurc2/sigma * exp(-0.5*xsurc2);

    valasym = (1.0+erf(m_asym_alpha/sqrt(2.0)*xsurc));
    valasymd = -m_asym_sigma_coeff*m_asym_alpha*sqrt(2)/sqrt(M_PI)*xsurc/sigma*exp(-0.5*xsurc2*m_asym_alpha*m_asym_alpha);
    val = valsym*valasymd+valsymd*valasym;
    return val;
}
inline
TFloat64List CLineProfileASYM::GetLineProfileVector(){
    TFloat64List v;
    return v;
}
}
#endif