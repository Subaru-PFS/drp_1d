#ifndef _REDSHIFT_LINE_PROFILE_SYM_
#define _REDSHIFT_LINE_PROFILE_SYM_
#include <string>
#include <math.h>
#include "RedshiftLibrary/common/datatypes.h"
namespace NSEpic
{
    /**
     * \ingroup Redshift
     */
    class CLineProfileSYM : public CLineProfile
    {
        public:
            CLineProfileSYM(const Float64 nsigmasupport=8.0);
            ~CLineProfileSYM()=default; 
            Float64 GetLineProfile(Float64 x, Float64 x0, const Float64 sigma) override;
            Float64 GetLineFlux( Float64 A, const Float64 sigma ) override;
            Float64 GetLineProfileDerivZ(Float64 x, Float64 lambda0, Float64 redshift, const Float64 sigma) override;
            Float64 GetLineProfileDerivSigma(Float64 x, Float64 x0, const Float64 sigma) override;

            TFloat64List GetLineProfileVector() override;//equivalent to ::computeKernel

            CLineProfileSYM(const CLineProfileSYM & other) = default; 
            CLineProfileSYM(CLineProfileSYM && other) = default; 
            CLineProfileSYM& operator=(const CLineProfileSYM& other) = default;  
            CLineProfileSYM& operator=(CLineProfileSYM&& other) = default; 

    };
inline
CLineProfileSYM::CLineProfileSYM(const Float64 nsigmasupport):
CLineProfile(nsigmasupport, "SYM")
{}

inline
Float64 CLineProfileSYM::GetLineProfile(Float64 x, Float64 x0, Float64 sigma){
    Float64 xc = x-x0;
    Float64 val = 0.0;
    Float64 xsurc;

    xsurc = xc/sigma;
    val = exp(-0.5*xsurc*xsurc);
    return val;
}

inline
Float64 CLineProfileSYM::GetLineFlux( Float64 A , Float64 sigma){
    return A*sigma*sqrt(2*M_PI);
}
inline
Float64 CLineProfileSYM::GetLineProfileDerivZ(Float64 x, Float64 lambda0, Float64 redshift, Float64 sigma)
{
    Float64 xc = x-lambda0*(1+redshift);
    Float64 val=0.0;
    Float64 xsurc;
    
    xsurc = xc/sigma;
    val = lambda0 /sigma * xsurc * exp(-0.5*xsurc*xsurc);
    return val;
}
inline
Float64 CLineProfileSYM::GetLineProfileDerivSigma(Float64 x, Float64 x0, Float64 sigma)
{
    Float64 val=0.0;
    Float64 xc = x-x0;
    Float64 xsurc, xsurc2;

    xsurc = xc/sigma;
    xsurc2 = xsurc*xsurc;
    val = xsurc2/sigma * exp(-0.5*xsurc2);
    return val;
}
inline
TFloat64List CLineProfileSYM::GetLineProfileVector(){
    TFloat64List v;
    return v;
}
}
#endif