#ifndef _REDSHIFT_SPECTRUM_LSFCONSTANT_
#define _REDSHIFT_SPECTRUM_LSFCONSTANT_

#include <RedshiftLibrary/spectrum/LSF.h>

namespace NSEpic
{

/**
 * \ingroup Redshift
 */
class CLSFConstantGaussian : public CLSF
{

public:

    CLSFConstantGaussian(const Float64 sigma=0.0);
    ~CLSFConstantGaussian();

    Float64             GetSigma() const;
    void                SetSigma(const Float64 sigma);
    bool                IsValid() const;

private:

    Float64             m_sigma=0.0;
};

}

#endif
