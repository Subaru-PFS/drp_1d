#ifndef _REDSHIFT_SPECTRUM_LSF_
#define _REDSHIFT_SPECTRUM_LSF_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/ray/lineprofile.h"
#include "RedshiftLibrary/ray/lineprofileSYM.h"
#include "RedshiftLibrary/processflow/parameterstore.h"
namespace NSEpic
{
  class CLineProfile;

  struct TLSFArguments
  {
    //std::string type;
    virtual ~TLSFArguments(){};
    TLSFArguments()=default;
    TLSFArguments(const TLSFArguments & other) = default; 
    TLSFArguments(TLSFArguments && other) = default; 
    TLSFArguments& operator=(const TLSFArguments& other) = default;  
    TLSFArguments& operator=(TLSFArguments&& other) = default; 
 };

 struct TLSFGaussianVarWidthArgs : virtual TLSFArguments
 {
   //std::string type = "GaussianVariableWidth";
   TFloat64List lambdas; 
   TFloat64List width;
   TLSFGaussianVarWidthArgs(TFloat64List _lambdas, TFloat64List _width):lambdas(_lambdas), width(_width){}
 };

struct TLSFGaussianConstantWidthArgs : virtual TLSFArguments
 {
   Float64 width;
   TLSFGaussianConstantWidthArgs(const std::shared_ptr<const CParameterStore>& parameterStore)
   {
      width = parameterStore->GetScoped<Float64>("LSF.width");//13.
   }
 };

struct TLSFGaussianConstantResolutionArgs : virtual TLSFArguments
 {
   Float64 resolution;
   TLSFGaussianConstantResolutionArgs(const std::shared_ptr<const CParameterStore>& parameterStore)
   {
      resolution = parameterStore->GetScoped<Float64>("LSF.resolution");// 2350.0
   }
 };

struct TLSFGaussianNISPVSSPSF201707Args : virtual TLSFArguments
 {
   Float64 sourcesize;
   TLSFGaussianNISPVSSPSF201707Args(const std::shared_ptr<const CParameterStore>& parameterStore)
   {
      sourcesize = parameterStore->GetScoped<Float64>("LSF.sourcesize");//0.1
   }
 };
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
      GaussianNISPVSSPSF201707,
      GaussianVariableWidth
    };
    CLSF(TLSFType name);
    CLSF(TLSFType name, CLineProfile_ptr profile);
    virtual ~CLSF() = default;

    CLSF(const CLSF & other) = default; 
    CLSF(CLSF && other) = default; 
    CLSF& operator=(const CLSF& other) = default;  
    CLSF& operator=(CLSF&& other) = default; 
    //GetWidth requires observed wavelength, not restframe
    virtual Float64             GetWidth(Float64 lambda) const=0;
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
