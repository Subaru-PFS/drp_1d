#ifndef _REDSHIFT_SPECTRUM_LSFBYEXTRINSIC_FIXEDRESOLUTION_COMPONENTS_
#define _REDSHIFT_SPECTRUM_LSFBYEXTRINSIC_FIXEDRESOLUTION_COMPONENTS_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/spectrum/LSF.h"
#include "RedshiftLibrary/ray/lineprofile.h"
namespace NSEpic
{
  /**
   * \ingroup Redshift
   */
  class CLSFGaussianConstantResolution : public CLSF
  {
    public:
        CLSFGaussianConstantResolution(const Float64 resolution);

        Float64             GetWidth(Float64 lambda) const override;
        bool                IsValid() const override;

        static std::shared_ptr<CLSF>   make_LSF(const std::shared_ptr<const TLSFArguments>& args);

    private:
        const Float64 m_instrumentResolutionEmpiricalFactor = 230.0/325.0/2.35;
        //CLineProfile_ptr m_profile{std::make_shared<CLineProfileSYM>()}; // default to sym
        const Float64 m_Resolution;
  };
inline
std::shared_ptr<CLSF> CLSFGaussianConstantResolution::make_LSF(const std::shared_ptr<const TLSFArguments>& args)
{
  const std::shared_ptr<const TLSFGaussianConstantResolutionArgs>& args_ = std::dynamic_pointer_cast<const TLSFGaussianConstantResolutionArgs>(args);
  return std::make_shared<CLSFGaussianConstantResolution>(args_->resolution);
}

}

#endif