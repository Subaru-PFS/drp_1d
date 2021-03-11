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
            ~CLineProfileLOR()=default; 
            virtual Float64 GetLineProfile(Float64 x, Float64 x0, const Float64 sigma) override; //override is optional, but a good practice with c++11
            virtual Float64 GetLineFlux( Float64 A, const Float64 sigma) override;
            virtual Float64 GetLineProfileDerivZ(Float64 x, Float64 lambda0, Float64 redshift, const Float64 sigma) override;
            virtual Float64 GetLineProfileDerivSigma(Float64 x, Float64 x0, const Float64 sigma) override;
            virtual Float64 GetNSigmaSupport() override;
            virtual TFloat64List GetLineProfileVector() override;//equivalent to ::computeKernel

            CLineProfileLOR(const CLineProfileLOR & other) = default; 
            CLineProfileLOR(CLineProfileLOR && other) = default; 
            CLineProfileLOR& operator=(const CLineProfileLOR& other) = default;  
            CLineProfileLOR& operator=(CLineProfileLOR&& other) = default; 

    };
inline
CLineProfileLOR::CLineProfileLOR(const Float64 nsigmasupport):
CLineProfile(nsigmasupport,  "LOR")
{}
inline
Float64 CLineProfileLOR::GetLineProfile(Float64 x, Float64 x0, Float64 sigma)
{
    Float64 xc = x-x0;
    Float64 val = 0.0;
    Float64 xsurc;

    xsurc = xc/sigma;
    val = 1.0/(1+xsurc*xsurc);
    return val;
}

inline
Float64 CLineProfileLOR::GetNSigmaSupport(){
    return m_nsigmasupport*2.0;
}
inline
Float64 CLineProfileLOR::GetLineFlux( Float64 A , Float64 sigma){
    return A*sigma*M_PI;
}
inline
Float64 CLineProfileLOR::GetLineProfileDerivZ(Float64 x, Float64 lambda0, Float64 redshift, Float64 sigma)
{
      Log.LogError("Deriv for Z not IMPLEMENTED for profile LOR");
      throw std::runtime_error("Deriv for Z not IMPLEMENTED for profile LOR");
}
inline
Float64 CLineProfileLOR::GetLineProfileDerivSigma(Float64 x, Float64 x0, Float64 sigma)
{
      Log.LogError("No derivate sigma for LOR profile");
      throw std::runtime_error("No derivate sigma for LOR profile");
}
inline
TFloat64List CLineProfileLOR::GetLineProfileVector(){
    TFloat64List v;
    return v;
}
}
#endif