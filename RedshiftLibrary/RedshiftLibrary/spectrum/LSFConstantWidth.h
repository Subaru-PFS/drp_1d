#ifndef _REDSHIFT_SPECTRUM_LSFCONSTANTWIDTH_
#define _REDSHIFT_SPECTRUM_LSFCONSTANTWIDTH_

#include <RedshiftLibrary/spectrum/LSF.h>

namespace NSEpic
{

/**
 * \ingroup Redshift
 */
class CLSFConstantGaussianWidth : public CLSF
{

public:

    CLSFConstantGaussianWidth(const Float64 width=0.0);
    ~CLSFConstantGaussianWidth();

    virtual Float64             GetWidth(Float64 lambda=-1.0)const override;
    virtual void                SetWidth(const Float64 width) override;
    virtual bool                IsValid() const override;

private:

    Float64             m_width = 0.0;
};

}

#endif
