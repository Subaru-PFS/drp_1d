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

    Float64               GetWidth(Float64 lambda) const override;
    bool                  IsValid() const override;

    static std::shared_ptr<CLSF> make_LSF(const std::shared_ptr<const TLSFArguments>& args);

private:
    const Float64             m_width = 0.0;
};
inline
std::shared_ptr<CLSF> CLSFGaussianConstantWidth::make_LSF(const std::shared_ptr<const TLSFArguments>& args)
{
    const std::shared_ptr<const TLSFGaussianConstantWidthArgs>& args_ = std::dynamic_pointer_cast<const TLSFGaussianConstantWidthArgs>(args);
    return std::make_shared<CLSFGaussianConstantWidth>(args_->width);
}

}

#endif
