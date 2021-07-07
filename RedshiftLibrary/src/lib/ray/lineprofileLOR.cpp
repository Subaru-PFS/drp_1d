#include "RedshiftLibrary/ray/lineprofileLOR.h"

using namespace NSEpic;
using namespace std;

CLineProfileLOR::CLineProfileLOR(const Float64 nsigmasupport):
CLineProfile(nsigmasupport, LOR)
{}
 
Float64 CLineProfileLOR::GetLineProfile(Float64 x, Float64 x0, Float64 sigma)
{
    Float64 xc = x-x0;
    Float64 val = 0.0;
    Float64 xsurc;

    xsurc = xc/sigma;
    val = 1.0/(1+xsurc*xsurc);
    return val;
}
 
Float64 CLineProfileLOR::GetNSigmaSupport() const
{
    return m_nsigmasupport*2.0;
}
 
Float64 CLineProfileLOR::GetLineFlux( Float64 A , Float64 sigma)
{
    return A*sigma*M_PI;
}
 
Float64 CLineProfileLOR::GetLineProfileDerivZ(Float64 x, Float64 lambda0, Float64 redshift, Float64 sigma)
{
      Log.LogError("Deriv for Z not IMPLEMENTED for profile LOR");
      throw std::runtime_error("Deriv for Z not IMPLEMENTED for profile LOR");
}
 
Float64 CLineProfileLOR::GetLineProfileDerivSigma(Float64 x, Float64 x0, Float64 sigma)
{
      Log.LogError("No derivate sigma for LOR profile");
      throw std::runtime_error("No derivate sigma for LOR profile");
}