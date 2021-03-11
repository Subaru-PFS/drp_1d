#ifndef _REDSHIFT_SPECTRUM_LSFBYEXTRINSIC_COMPONENTS_
#define _REDSHIFT_SPECTRUM_LSFBYEXTRINSIC_COMPONENTS_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/spectrum/LSF.h"
#include "RedshiftLibrary/ray/lineprofile.h"

namespace NSEpic
{
  /**
   * \ingroup Redshift
   */
  class CLSFbyExtrinsicComponents : public CLSF
  {
    enum TLSFType {
      FIXED,
      NISPSIM2016,
      NISPVSSPSF201707,
      FROMSPECTRUMDATA
    };
    public:
        CLSFbyExtrinsicComponents(const std::string LSFType, 
                                  const Float64 resolution, 
                                  const Float64 nominalWidth=13.);
        ~CLSFbyExtrinsicComponents();
        virtual Float64             GetWidth(Float64 lambda=-1.0) const override;
        virtual void                SetWidth(const Float64 width) override;
        virtual bool                IsValid() const override;
        void                        SetSourcesizeDispersion(Float64 sigma);
        //Float64                     GetLineProfile(Float64 x, Float64 x0);
    private:
        Float64 m_width = 0.;
        TLSFType m_type;
        CLineProfile_ptr m_profile{std::make_shared<CLineProfileSYM>()}; // default to sym
        Float64 m_NominalWidth = 13.;
        Float64 m_Resolution;
        Float64 m_instrumentResolutionEmpiricalFactor;
        Float64 m_SourceSizeDispersion; // source size in the dispersion direction
                                        // (sigma arcsec)
  };
}

#endif
