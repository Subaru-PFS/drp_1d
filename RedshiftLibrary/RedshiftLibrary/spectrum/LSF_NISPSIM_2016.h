#ifndef _REDSHIFT_SPECTRUM_LSFBYEXTRINSIC_NISPSIM2016_COMPONENTS_
#define _REDSHIFT_SPECTRUM_LSFBYEXTRINSIC_NISPSIM2016_COMPONENTS_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/spectrum/LSF.h"
#include "RedshiftLibrary/ray/lineprofile.h"
#include "RedshiftLibrary/processflow/parameterstore.h"
namespace NSEpic
{
  /**
   * \ingroup Redshift
   */
  class CLSFGaussianNISPSIM2016 : public CLSF
  {
    public:
        CLSFGaussianNISPSIM2016();

        Float64             GetWidth(Float64 lambda) const override;
        bool                IsValid() const override;
        void                SetSourcesizeDispersion(Float64 sigma);

        static std::shared_ptr<CLSF> make_LSF(const std::shared_ptr<const TLSFArguments>& args);
  };
inline
std::shared_ptr<CLSF>   CLSFGaussianNISPSIM2016::make_LSF(const std::shared_ptr<const TLSFArguments>& args)
{
     return std::make_shared<CLSFGaussianNISPSIM2016>();
}

}
#endif