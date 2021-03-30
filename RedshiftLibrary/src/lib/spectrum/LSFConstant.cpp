#include <RedshiftLibrary/spectrum/LSFConstant.h>
#include <RedshiftLibrary/log/log.h>


using namespace NSEpic;
using namespace std;

/**
 * Constructor.
 */
CLSFConstantGaussian::CLSFConstantGaussian(const Float64 sigma):
    m_sigma(sigma)
{

}

/**
 * Destructor.
 */
CLSFConstantGaussian::~CLSFConstantGaussian()
{

}

/**
 * Return the spectral resolution.
 */
Float64 CLSFConstantGaussian::GetSigma(Float64 lambda) const
{
    return m_sigma;
}

/**
 * Set the spectral resolution.
 */
void CLSFConstantGaussian::SetSigma(const Float64 sigma)
{
    m_sigma = sigma;
}

/**
 * Check validity of the LSF.
 */
bool CLSFConstantGaussian::IsValid() const
{
    return (m_sigma>0.0);
}