#include "RedshiftLibrary/ray/lineprofileSYM.h"

using namespace NSEpic;
using namespace std;
 
CLineProfileSYM::CLineProfileSYM(const Float64 nsigmasupport):
CLineProfile(nsigmasupport, SYM)
{}

Float64 CLineProfileSYM::GetLineProfile(Float64 x, Float64 x0, Float64 sigma)
{
    Float64 xc = x-x0;
    Float64 val = 0.0;
    Float64 xsurc;

    xsurc = xc/sigma;
    val = exp(-0.5*xsurc*xsurc);
    return val;
}

Float64 CLineProfileSYM::GetLineFlux( Float64 A , Float64 sigma)
{
    return A*sigma*sqrt(2*M_PI);
}

Float64 CLineProfileSYM::GetLineProfileDerivZ(Float64 x, Float64 lambda0, Float64 redshift, Float64 sigma)
{
    Float64 xc = x-lambda0*(1+redshift);
    Float64 val=0.0;
    Float64 xsurc;
    
    xsurc = xc/sigma;
    val = lambda0 /sigma * xsurc * exp(-0.5*xsurc*xsurc);
    return val;
}
 
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