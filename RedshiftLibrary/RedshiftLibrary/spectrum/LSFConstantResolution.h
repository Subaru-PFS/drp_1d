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

        virtual Float64             GetWidth(Float64 lambda=-1.0) const override;
        virtual bool                IsValid() const override;

        //static std::shared_ptr<CLSF>   make_LSF(const Float64 resolution);

    private:
        const Float64 m_instrumentResolutionEmpiricalFactor = 230.0/325.0/2.35;
        //CLineProfile_ptr m_profile{std::make_shared<CLineProfileSYM>()}; // default to sym
        const Float64 m_Resolution;
  };
/*inline
std::shared_ptr<CLSF> CLSFGaussianConstantResolution::make_LSF(const Float64 resolution)
{
    return std::make_shared<CLSFGaussianConstantResolution>(resolution);
}*/

}

#endif