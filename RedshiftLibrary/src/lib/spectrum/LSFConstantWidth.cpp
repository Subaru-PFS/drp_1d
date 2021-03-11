#include <RedshiftLibrary/spectrum/LSFConstantWidth.h>
#include <RedshiftLibrary/log/log.h>


using namespace NSEpic;
using namespace std;

/**
 * Constructor.
 */
CLSFConstantGaussianWidth::CLSFConstantGaussianWidth(const Float64 width):
    m_width(width)
{

}

/**
 * Destructor.
 */
CLSFConstantGaussianWidth::~CLSFConstantGaussianWidth()
{

}

/**
 * Return the spectral resolution.
 */
Float64 CLSFConstantGaussianWidth::GetWidth(Float64 lambda) const
{
    return m_width;
}

/**
 * Set the spectral resolution.
 */
void CLSFConstantGaussianWidth::SetWidth(const Float64 width)
{
    m_width = width;
}

/**
 * Check validity of the LSF.
 */
bool CLSFConstantGaussianWidth::IsValid() const
{
    return (m_width>0.0);
}
