#include "RedshiftLibrary/spectrum/LSF_NISPSIM_2016.h"
#include <RedshiftLibrary/log/log.h>

using namespace NSEpic;
using namespace std;
  
CLSFGaussianNISPSIM2016::CLSFGaussianNISPSIM2016():
    CLSF(GaussianNISPSIM2016, std::make_shared<CLineProfileSYM>())
{
    
}

Float64 CLSFGaussianNISPSIM2016::GetWidth(Float64 lambda) const
{
    Float64 instrumentSigma = (lambda*8.121e-4 + 7.4248)/2.35;
    return instrumentSigma;
}

bool CLSFGaussianNISPSIM2016::IsValid() const
{
    return true;//TODO
}