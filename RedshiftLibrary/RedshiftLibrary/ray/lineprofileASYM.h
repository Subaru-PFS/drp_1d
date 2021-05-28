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
            /*enum TASYMType {
                BASE,
                FIT,
                FIXED
            };*/
            /*
            enum TCenteringMethod{
                NONE,
                MEAN,
                MODE
            }*/
            CLineProfileASYM(const Float64 nsigmasupport = 8.0, 
                            const TAsymParams params = {1., 4.5, 0.}, 
                            const std::string centeringMethod = "none");
            CLineProfileASYM(const TProfile pltype,
                            const Float64 nsigmasupport = 8.0, 
                            const TAsymParams params = {2., 2.5, 0.}, 
                            const std::string centeringMethod = "mean");//mainly called by asymfit
 
            Float64 GetLineProfile(Float64 x, Float64 x0, const Float64 sigma) override;
            Float64 GetLineFlux(Float64 A, const Float64 sigma) override;
            Float64 GetLineProfileDerivZ(Float64 x, Float64 x0, Float64 redshift, const Float64 sigma) override;
            Float64 GetLineProfileDerivSigma(Float64 x, Float64 x0, const Float64 sigma) override;
            Float64 GetNSigmaSupport() override;
            TFloat64List GetLineProfileVector() override;//equivalent to ::computeKernel

            Float64 GetXSurc(Float64 xc, Float64& sigma, Float64& xsurc);
            Float64 GetAsymDelta();

            virtual Bool    isAsymFixed() override;
            virtual Bool    isAsymFit()   override;

            ~CLineProfileASYM() = default;
            CLineProfileASYM(const CLineProfileASYM & other) = default; 
            CLineProfileASYM(CLineProfileASYM && other) = default; 
            CLineProfileASYM& operator=(const CLineProfileASYM& other) = default;  
            CLineProfileASYM& operator=(CLineProfileASYM&& other) = default; 
        protected: 
            Float64 m_asym_sigma_coeff = 1.0;//vs 2. for asymFit/Fixed
            Float64 m_asym_alpha = 4.5;
            Float64 m_asym_delta = 0.;
            std::string m_centeringMethod = "none";
            Float64 m_constSigma = 1;//vs 2.5 for AsymFit and AsymFixed
    };
inline
CLineProfileASYM::CLineProfileASYM(const Float64 nsigmasupport, TAsymParams params, const std::string centeringMethod):
    CLineProfile(nsigmasupport, ASYM),
    m_centeringMethod(centeringMethod),
    m_asym_sigma_coeff(params.sigma),
    m_asym_alpha(params.alpha),
    m_asym_delta(params.delta)
{
    if(m_centeringMethod == "mean" || m_centeringMethod == "mode")
        m_constSigma = 2.5; //not sure if should be also a param
}

inline
CLineProfileASYM::CLineProfileASYM(const TProfile pltype, const Float64 nsigmasupport, const TAsymParams params, const std::string centeringMethod):
    CLineProfile(nsigmasupport, pltype),
    m_centeringMethod(centeringMethod),
    m_asym_sigma_coeff(params.sigma),
    m_asym_alpha(params.alpha),
    m_asym_delta(params.delta)
{
    if(m_centeringMethod == "mean" || m_centeringMethod == "mode")
        m_constSigma = 2.5; //not sure if should be also a param
}

inline 
Float64 CLineProfileASYM::GetXSurc(Float64 xc, Float64& sigma, Float64& xsurc)
{
    sigma = sigma*m_asym_sigma_coeff;
    Float64 delta = m_asym_alpha/std::sqrt(1.+m_asym_alpha*m_asym_alpha);
    Float64 muz = delta*sqrt(2./M_PI);

    Float64 m0 = 0;
    if(m_centeringMethod == "none"){ //ASYM
        m0 = 0; //i.e., no centering;
        delta = 0.; //reset it cause no centering here

        if(m_asym_delta)//temporary check
            throw std::runtime_error("Problem in configuring ASYM lineprofile: asym_delta should be null");
    }else if(m_centeringMethod == "mean")//ASYMFIT, ASYMFIXED
    {//correction in order to have the line shifted on the mean: from https://en.wikipedia.org/wiki/Skew_normal_distribution
        m0 = muz;
    }else if(m_centeringMethod == "mode")
    {//correction in order to have the line shifted on the mode: from https://en.wikipedia.org/wiki/Skew_normal_distribution
        Float64 sigmaz = std::sqrt(1-muz*muz);
        Float64 gamma1 = ((4-M_PI)/2.0)*pow(delta*std::sqrt(2/M_PI), 3.)/pow(1-2*delta*delta/M_PI, 3./2.);
        m0 = muz - gamma1*sigmaz/2.0 - 0.5*exp(-2*M_PI/m_asym_alpha);
    }

    xc = xc + sigma*m0;

    Float64 xcd = xc + m_asym_delta;
    xsurc = xcd/sigma;
    return m0;
}

inline
Float64 CLineProfileASYM::GetLineProfile(Float64 x, Float64 x0, Float64 sigma)
{
    Float64 xsurc;
    GetXSurc(x-x0, sigma, xsurc);

    Float64 val = exp(-0.5*xsurc*xsurc)*(1.0+erf(m_asym_alpha/sqrt(2.0)*xsurc));
    return val;
}
inline
Float64 CLineProfileASYM::GetNSigmaSupport(){
    return m_nsigmasupport*m_asym_sigma_coeff*m_constSigma;
}
inline
Float64 CLineProfileASYM::GetLineFlux(Float64 A, Float64 sigma)
{
    return A*sigma*m_asym_sigma_coeff*sqrt(2*M_PI);
}
inline
Float64 CLineProfileASYM::GetLineProfileDerivZ(Float64 x, Float64 x0, Float64 redshift, Float64 sigma)
{
    Float64 xsurc;
    Float64 alpha2 = m_asym_alpha*m_asym_alpha;
    Float64 xc = x-x0*(1+redshift);

    GetXSurc(xc, sigma, xsurc);

    Float64 xsurc2 = xsurc*xsurc;
    Float64 val = x0 /sigma * xsurc * exp(-0.5*xsurc2) *(1.0+erf(m_asym_alpha/sqrt(2.0)*xsurc)) -m_asym_alpha * x0 /sqrt(2*M_PI) /sigma * exp(-(1+alpha2)/2 * xsurc2);
    return val;
}

inline
Float64 CLineProfileASYM::GetLineProfileDerivSigma(Float64 x, Float64 x0, Float64 sigma)
{
    Float64 alpha2 = m_asym_alpha*m_asym_alpha;

    Float64 xsurc;
    Float64 muz = GetXSurc(x-x0, sigma, xsurc);
    Float64 xsurc2 = xsurc*xsurc;

    Float64 valSym = exp(-0.5*xsurc2);
    Float64 valAsym = (1.0+erf(m_asym_alpha/sqrt(2.0)*xsurc));

    Float64 valSymD, valAsymD;
    if(m_centeringMethod == "none"){
        //temporary: double check muz == 0
        if(muz)
            throw std::runtime_error("Problem: mu is not null for ASYM profile!");
    }
    valSymD  = m_asym_sigma_coeff * (xsurc2/sigma - muz*xsurc/sigma)*exp(-0.5*xsurc2); 
    valAsymD = m_asym_sigma_coeff * m_asym_alpha*sqrt(2)/(sigma*sqrt(M_PI))*(muz - xsurc)*exp(-0.5*xsurc2*alpha2);

    Float64 val = valSym * valAsymD + valSymD * valAsym;
    return val;
}
inline
TFloat64List CLineProfileASYM::GetLineProfileVector(){
    TFloat64List v;
    return v;
}
inline
Bool CLineProfileASYM::isAsymFixed(){ //probably we need to clean this
    return 1;
}
inline
Bool CLineProfileASYM::isAsymFit(){ //probably we need to clean this
    return 0;
}
inline
Float64 CLineProfileASYM::GetAsymDelta(){
    return m_asym_delta;//default. Mainly used for asymfit/fixed
}
}
#endif