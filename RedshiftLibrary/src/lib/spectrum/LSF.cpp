#include <RedshiftLibrary/spectrum/LSF.h>
#include <RedshiftLibrary/log/log.h>


using namespace NSEpic;
using namespace std;

/**
 * Constructor.
 */
CLSF::CLSF(const Float64 sigma):
    m_sigma(sigma)
{

}

/**
 * Destructor.
 */
CLSF::~CLSF()	
{

}

/**
 * Return the spectral resolution.
 */
const Float64 CLSF::GetSigma() const
{
    return m_sigma;
}

/**
 * Set the spectral resolution.
 */
void CLSF::SetSigma(const Float64 sigma)
{
    m_sigma = sigma;
}

/**
 * Check validity of the LSF
 */
bool CLSF::IsValid() const
{
    return (m_sigma>0.0);
}
