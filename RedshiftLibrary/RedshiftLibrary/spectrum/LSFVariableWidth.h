#ifndef _REDSHIFT_SPECTRUM_LSFVARIABLEWIDTH_
#define _REDSHIFT_SPECTRUM_LSFVARIABLEWIDTH_

#include "RedshiftLibrary/spectrum/LSF.h"
#include "RedshiftLibrary/spectrum/spectralaxis.h"
namespace NSEpic
{

/**
 * \ingroup Redshift
 */
class CLSFGaussianVariableWidth : public CLSF
{
    public:
        CLSFGaussianVariableWidth(const std::shared_ptr<const TLSFGaussianVarWidthArgs>& args);
        Float64               GetWidth(Float64 lambda) const override;

        bool                  IsValid() const override;

        static std::shared_ptr<CLSF> make_LSF(const std::shared_ptr<const TLSFArguments>& args);

    private:
        TFloat64List             m_width;
        CSpectrumSpectralAxis    m_spcAxis;
};

inline
std::shared_ptr<CLSF> CLSFGaussianVariableWidth::make_LSF(const std::shared_ptr<const TLSFArguments>& args)
{
    const std::shared_ptr<const TLSFGaussianVarWidthArgs>& args_ = std::dynamic_pointer_cast<const TLSFGaussianVarWidthArgs>(args);
    return std::make_shared<CLSFGaussianVariableWidth>(args_);
}


}

#endif
