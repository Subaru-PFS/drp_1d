#ifndef _REDSHIFT_LINE_PROFILE_ASYMFIT_
#define _REDSHIFT_LINE_PROFILE_ASYMFIT_
#include <string>
#include <math.h>
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/ray/lineprofile.h"
namespace NSEpic
{
    /**
     * \ingroup Redshift
     */
    class CLineProfileASYMFIT: public CLineProfile
    {
        public:
            CLineProfileASYMFIT(const Float64 nsigmasupport = 8.0);
            Float64 GetLineProfile(Float64 x, Float64 x0, Float64 sigma) override;
            Float64 GetLineFlux( Float64 A, Float64 sigma) override;
            Float64 GetLineProfileDerivZ(Float64 x, Float64 lambda0, Float64 redshift, Float64 sigma) override;
            Float64 GetLineProfileDerivSigma(Float64 x, Float64 x0, Float64 sigma) override;
            Float64 GetNSigmaSupport() override;
            TFloat64List GetLineProfileVector() override;//equivalent to ::computeKernel

            CLineProfileASYMFIT(const CLineProfileASYMFIT & other) = default; 
            CLineProfileASYMFIT(CLineProfileASYMFIT && other) = default; 
            CLineProfileASYMFIT& operator=(const CLineProfileASYMFIT& other) = default;  
            CLineProfileASYMFIT& operator=(CLineProfileASYMFIT&& other) = default; 

            ~CLineProfileASYMFIT()=default; 
        private:
            Float64 m_asymfit_delta = 0.0;
            Float64 m_asymfit_alpha = 2.0;
            Float64 m_asymfit_sigma_coeff = 2.0;
            // TAsymParams     m_asymParams = {NAN, NAN, NAN};
            /*
            m_SourceSizeDispersion = 0.1;

            m_asym_sigma_coeff = 1.0;
            m_asym_alpha = 4.5;

            m_symxl_sigma_coeff = 5.0;

            m_asym2_sigma_coeff = 2.0;
            m_asym2_alpha = 2.0;*/

    };
inline
CLineProfileASYMFIT::CLineProfileASYMFIT(const Float64 nsigmasupport):
CLineProfile(nsigmasupport, "ASYMFIT")
{}

inline
Float64 CLineProfileASYMFIT::GetLineProfile(Float64 x, Float64 x0, Float64 sigma)
{  
    Float64 xc = x-x0, xcd;
    Float64  delta, muz;
    Float64 val = 0.0;
    Float64 xsurc;

    sigma = sigma*m_asymfit_sigma_coeff;
    //correction in order to have the line shifted on the mean: from https://en.wikipedia.org/wiki/Skew_normal_distribution
    delta = m_asymfit_alpha/std::sqrt(1.+m_asymfit_alpha*m_asymfit_alpha);
    muz = delta*sqrt(2./M_PI);
    xc = xc + sigma*muz;

    /*
    //correction in order to have the line shifted on the mode: from https://en.wikipedia.org/wiki/Skew_normal_distribution
    delta = alpha/std::sqrt(1.+alpha*alpha);
    muz = delta*sqrt(2./M_PI);
    sigmaz = std::sqrt(1-muz*muz);
    gamma1 = ((4-M_PI)/2.0)*pow(delta*std::sqrt(2/M_PI), 3.)/pow(1-2*delta*delta/M_PI, 3./2.);
    m0 = muz - gamma1*sigmaz/2.0 - 0.5*exp(-2*M_PI/alpha);
    xc = xc + sigma*m0;
    //*/

    xcd = xc + m_asymfit_delta;
    xsurc = xcd/sigma;
    val = exp(-0.5*xsurc*xsurc)*(1.0+erf(m_asymfit_alpha/sqrt(2.0)*xsurc));
    return val;
}
inline
Float64 CLineProfileASYMFIT::GetNSigmaSupport(){
    return m_nsigmasupport*m_asymfit_sigma_coeff*2.5;
}
inline
Float64 CLineProfileASYMFIT::GetLineFlux( Float64 A , Float64 sigma){
    return A*sigma*m_asymfit_sigma_coeff*sqrt(2*M_PI);
}
inline
Float64 CLineProfileASYMFIT::GetLineProfileDerivZ(Float64 x, Float64 lambda0, Float64 redshift, Float64 sigma)
{
    Float64 xc = x-lambda0*(1+redshift), xcd;
    Float64 val=0.0;
    Float64 xsurc, xsurc2;
    Float64 delta, muz; 
    Float64 alpha2;

    sigma = sigma*m_asymfit_sigma_coeff;
    alpha2 = m_asymfit_alpha*m_asymfit_alpha;

    //correction in order to have the line shifted on the mean: from https://en.wikipedia.org/wiki/Skew_normal_distribution
    delta = m_asymfit_alpha/std::sqrt(1.+alpha2);
    muz = delta*sqrt(2./M_PI);
    xc = xc + sigma*muz;

    /*
    //correction in order to have the line shifted on the mode: from https://en.wikipedia.org/wiki/Skew_normal_distribution
    delta = alpha/std::sqrt(1.+alpha2);
    muz = delta*sqrt(2./M_PI);
    sigmaz = std::sqrt(1-muz*muz);
    gamma1 = ((4-M_PI)/2.0)*pow(delta*std::sqrt(2/M_PI), 3.)/pow(1-2*delta*delta/M_PI, 3./2.);
    m0 = muz - gamma1*sigmaz/2.0 - 0.5*exp(-2*M_PI/alpha);
    xc = xc + sigma*m0;
    //*/

    xcd = xc+m_asymfit_delta;
    xsurc = xcd/sigma;
    xsurc2 = xsurc*xsurc;
    val = lambda0 /sigma * xsurc * exp(-0.5*xsurc2) *(1.0+erf(m_asymfit_alpha/sqrt(2.0)*xsurc)) -m_asymfit_alpha * lambda0 /sqrt(2*M_PI) /sigma * exp(-(1+alpha2)/2 * xsurc2);

    return val;
}
inline
Float64 CLineProfileASYMFIT::GetLineProfileDerivSigma(Float64 x, Float64 x0, Float64 sigma)
{
    Float64 val=0.0;
    Float64 xc = x-x0;
    Float64 xsurc, xsurc2, alpha2, valsym, valsymd, valasym, valasymd, xcd;
    Float64 muz;

    sigma = sigma*m_asymfit_sigma_coeff;
    alpha2 = m_asymfit_alpha*m_asymfit_alpha;

    //correction in order to have the line shifted on the mean: from https://en.wikipedia.org/wiki/Skew_normal_distribution
    Float64 delta = m_asymfit_alpha/std::sqrt(1.+alpha2);
    muz = delta*sqrt(2./M_PI);
    xc = xc + sigma*muz;

    /*
    //correction in order to have the line shifted on the mode: from https://en.wikipedia.org/wiki/Skew_normal_distribution
    delta = alpha/std::sqrt(1.+alpha2);
    muz = delta*sqrt(2./M_PI);
    sigmaz = std::sqrt(1-muz*muz);
    gamma1 = ((4-M_PI)/2.0)*pow(delta*std::sqrt(2/M_PI), 3.)/pow(1-2*delta*delta/M_PI, 3./2.);
    m0 = muz - gamma1*sigmaz/2.0 - 0.5*exp(-2*M_PI/alpha);
    xc = xc + sigma*m0;
    //*/

    xcd = xc + m_asymfit_delta;
    xsurc = xcd/sigma;
    xsurc2 = xsurc * xsurc;
    valsym = exp(-0.5*xsurc2);
    valsymd = m_asymfit_sigma_coeff * (xsurc2 /sigma - muz*xsurc/sigma) * exp(-0.5*xsurc2); // for mean centering
    /* valsymd = coeff * (xsurc2 /sigma - m0*xsurc/sigma) * exp(-0.5*xsurc2); //for mode centering*/

    valasym = (1.0+erf(m_asymfit_alpha/sqrt(2.0)*xsurc));
    valasymd = m_asymfit_sigma_coeff  * m_asymfit_alpha*sqrt(2)/(sigma*sqrt(M_PI))*(muz - xsurc) * exp(-0.5*xsurc2*alpha2); // for mean centering
    /* valasymd = coeff  * alpha*sqrt(2)/(sigma*sqrt(M_PI))*(m0 - xsurc) * exp(-0.5*xsurc2*alpha2); // for mode centering */

    val = valsym * valasymd + valsymd * valasym;
    return val;
}
inline
TFloat64List CLineProfileASYMFIT::GetLineProfileVector(){
    TFloat64List v;
    return v;
}
}
#endif