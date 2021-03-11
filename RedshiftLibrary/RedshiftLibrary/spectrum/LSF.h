#ifndef _REDSHIFT_SPECTRUM_LSF_
#define _REDSHIFT_SPECTRUM_LSF_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/ray/lineprofile.h>
#include <RedshiftLibrary/ray/lineprofileSYM.h>
namespace NSEpic
{
class CLineProfile;
/**
 * \ingroup Redshift
 */
class CLSF
{

public:

    virtual ~CLSF(){};

    virtual Float64             GetWidth(Float64 lambda=-1.0) const=0;
    virtual void                SetWidth(const Float64 width)=0;
    virtual bool                IsValid() const=0;
    Float64                     GetLineProfile(Float64 lambda, Float64 lambda0 = 0.);
    void                        SetSourcesizeDispersion(Float64 sigma) const{};//empty default implementation
protected:
    CLineProfile_ptr m_profile{std::make_shared<CLineProfileSYM>()}; // default to sym default to sym

};

inline
Float64 CLSF::GetLineProfile (Float64 lambda, Float64 lambda0){
    return m_profile->GetLineProfile(lambda, lambda0, GetWidth());
}
}
#endif
