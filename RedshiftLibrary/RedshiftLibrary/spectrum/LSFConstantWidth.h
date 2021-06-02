#ifndef _REDSHIFT_SPECTRUM_LSFCONSTANTWIDTH_
#define _REDSHIFT_SPECTRUM_LSFCONSTANTWIDTH_

#include "RedshiftLibrary/spectrum/LSF.h"

namespace NSEpic
{

/**
 * \ingroup Redshift
 */
class CLSFGaussianConstantWidth : public CLSF
{

public:
    CLSFGaussianConstantWidth(const Float64 width = 0.);

    Float64               GetWidth(Float64 lambda=-1.0) const override;
    bool                  IsValid() const override;

    //static std::shared_ptr<CLSF> make_LSF();

private:
    const Float64             m_width = 0.0;
};
/*inline
std::shared_ptr<CLSF> CLSFGaussianConstantWidth::make_LSF()
{
    return std::make_shared<CLSFGaussianConstantWidth>();
}*/

}

#endif
