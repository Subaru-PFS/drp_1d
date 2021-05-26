#ifndef _REDSHIFT_SPECTRUM_LSF_
#define _REDSHIFT_SPECTRUM_LSF_

#define __stdcall

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/ray/lineprofile.h>
#include <RedshiftLibrary/ray/lineprofileSYM.h>
namespace NSEpic
{
class CLineProfile;

  typedef struct{
    Float64 resolution; 
    Float64 width;
    Float64 sourcesize;
  }TLSFArguments;
    //typedef struct LSFArguments TLSFArguments;
    //TLSFArguments LSFDefaultArgs = { 2350.0, 0., 0.1};//default width to 0, i.e., invalid lsf
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

    //CLSF(const CLSF & other) = default; 
    CLSF(CLSF && other) = default; 
    //CLSF& operator=(const CLSF& other) = default;  
    CLSF& operator=(CLSF&& other) = default; 

    virtual Float64             GetWidth(Float64 lambda=-1.0) const=0;
    virtual bool                IsValid() const = 0;
    Float64                     GetLineProfile(Float64 lambda, Float64 lambda0 = 0.);
    void                        SetSourcesizeDispersion(Float64 sigma) const{};//empty default implementation
    //define a factory method (static by default) to return an instance of the subclass whenever needed
    static std::shared_ptr<CLSF>  make_LSF(const std::string lsfType, const TLSFArguments& args = { 2350.0, 0., 0.1});
                                          
    const TLSFType  m_name;
  /*  
    template<typename ...T>
    typedef std::shared_ptr<CLSF> ( *_CreateLSFFn)(T val, ...); //typedef a function pointer 
    static _CreateLSFFn CreateLSFFn; //define it as static to be called with no object instanciation*/
protected:
    CLineProfile_ptr m_profile;


};

}
#endif
