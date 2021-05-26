#ifndef _REDSHIFT_SPECTRUM_LSFBYEXTRINSIC_NISPSIM2016_COMPONENTS_
#define _REDSHIFT_SPECTRUM_LSFBYEXTRINSIC_NISPSIM2016_COMPONENTS_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/spectrum/LSF.h"
#include "RedshiftLibrary/ray/lineprofile.h"

namespace NSEpic
{
  /**
   * \ingroup Redshift
   */
  class CLSFGaussianNISPSIM2016 : public CLSF
  {
    public:
        CLSFGaussianNISPSIM2016();

        virtual Float64             GetWidth(Float64 lambda=-1.0) const override;
        virtual bool                IsValid() const override;
        void                        SetSourcesizeDispersion(Float64 sigma);
        //static std::shared_ptr<CLSF> make_LSF();
  };
/*inline
std::shared_ptr<CLSF>   CLSFGaussianNISPSIM2016::make_LSF()
{
     return std::make_shared<CLSFGaussianNISPSIM2016>();
}*/

}
#endif