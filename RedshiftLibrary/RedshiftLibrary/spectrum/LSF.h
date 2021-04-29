#ifndef _REDSHIFT_SPECTRUM_LSF_
#define _REDSHIFT_SPECTRUM_LSF_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/ray/lineprofile.h"
#include "RedshiftLibrary/ray/lineprofileSYM.h"
namespace NSEpic
{
class CLineProfile;

  typedef struct{
    Float64 resolution; 
    Float64 width;
    Float64 sourcesize;
  }TLSFArguments;
/**
 * \ingroup Redshift
 */
class CLSF
{

public:
    enum TLSFType {
      GaussianConstantWidth,
      GaussianConstantResolution,
      GaussianNISPSIM2016,
      GaussianNISPVSSPSF201707
    };
    CLSF(TLSFType name);
    CLSF(TLSFType name, CLineProfile_ptr profile);
    virtual ~CLSF() = default;

    CLSF(const CLSF & other) = default; 
    CLSF(CLSF && other) = default; 
    CLSF& operator=(const CLSF& other) = default;  
    CLSF& operator=(CLSF&& other) = default; 
    //GetWidth requires observed wavelength, not restframe
    virtual Float64             GetWidth(Float64 lambda=-1.0) const=0;
    virtual bool                IsValid() const = 0;
    Float64                     GetLineProfile(Float64 lambda, Float64 lambda0 = 0.) const;
    Float64                     GetLineProfile (Float64 lambda, Float64 lambda0, Float64 sigma0) const;
    void                        SetSourcesizeDispersion(Float64 sigma) const{};//empty default implementation
    std::shared_ptr<const CLineProfile>   GetProfile() const;                      
    const TLSFType  m_name;
    
protected:
    CLineProfile_ptr m_profile;
};

}
#endif
