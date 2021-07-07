#include "RedshiftLibrary/spectrum/LSFConstantResolution.h"
#include "RedshiftLibrary/log/log.h"


using namespace NSEpic;
using namespace std;

CLSFGaussianConstantResolution::CLSFGaussianConstantResolution(const Float64 resolution):
    CLSF(GaussianConstantResolution, std::make_shared<CLineProfileSYM>()),
    m_Resolution(resolution)
{
    IsValid();
}

Float64 CLSFGaussianConstantResolution::GetWidth(Float64 lambda) const
{
    Float64 defaultSigma = lambda/m_Resolution*m_instrumentResolutionEmpiricalFactor;
    return defaultSigma;
}

bool CLSFGaussianConstantResolution::IsValid() const
{
    return (m_Resolution>1.);
}