#ifndef _REDSHIFT_SPECTRUM_LSFBYEXTRINSIC_NISPVSSPSF2017_COMPONENTS_
#define _REDSHIFT_SPECTRUM_LSFBYEXTRINSIC_NISPVSSPSF2017_COMPONENTS_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/spectrum/LSF.h"
#include "RedshiftLibrary/ray/lineprofile.h"
namespace NSEpic
{
  /**
   * \ingroup Redshift
   */
  class CLSFGaussianNISPVSSPSF201707 : public CLSF
  {
    public:
        CLSFGaussianNISPVSSPSF201707(Float64 sourcesize=0.1);
        Float64             GetWidth(Float64 lambda) const override;
        bool                IsValid() const override;

        static std::shared_ptr<CLSF>   make_LSF(const std::shared_ptr<const TLSFArguments>& args);
    private:
        const Float64 m_SourceSizeDispersion;// = 0.1; // source size in the dispersion direction
                                        // (sigma arcsec)
  };

inline
std::shared_ptr<CLSF> CLSFGaussianNISPVSSPSF201707::make_LSF(const std::shared_ptr<const TLSFArguments>& args)
{
  const std::shared_ptr<const TLSFGaussianNISPVSSPSF201707Args>& args_ = std::dynamic_pointer_cast<const TLSFGaussianNISPVSSPSF201707Args>(args);
  return std::make_shared<CLSFGaussianNISPVSSPSF201707>(args_->sourcesize);
}

}

#endif