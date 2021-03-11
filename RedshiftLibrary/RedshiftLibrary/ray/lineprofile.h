#ifndef _REDSHIFT_LINE_PROFILE_
#define _REDSHIFT_LINE_PROFILE_

#include <string>
#include "RedshiftLibrary/common/datatypes.h"
namespace NSEpic
{
    /**
     * \ingroup Redshift
     * Abstract class for different line profiles
     */
    class CLineProfile
    {
        /*
        enum TProfile
        {
        NONE,
        SYM,
        SYMXL,
        LOR,
        ASYM,
        ASYM2,
        ASYMFIT,
        ASYMFIXED,
        EXTINCT
        };
        */
        public:
            CLineProfile(const Float64 nsigmasupport=8.0);
            CLineProfile(const Float64 nsigmasupport=8.0, const std::string name="NONE");
            virtual Float64 GetLineProfile(Float64 x, Float64 x0, Float64 sigma)=0;
            virtual Float64 GetLineFlux( Float64 A, Float64 sigma)=0;
            virtual Float64 GetLineProfileDerivZ(Float64 x, Float64 lambda0, Float64 redshift, Float64 sigma)=0;
            virtual Float64 GetLineProfileDerivSigma(Float64 x, Float64 x0, Float64 sigma)=0;
            virtual Float64 GetNSigmaSupport();
            virtual TFloat64List GetLineProfileVector()=0;//TODO: equivalent to ::computeKernel
            const std::string& GetName();

            CLineProfile(const CLineProfile & other) = default; 
            CLineProfile(CLineProfile && other) = default; 
            CLineProfile& operator=(const CLineProfile& other) = default;  
            CLineProfile& operator=(CLineProfile&& other) = default; 

            virtual ~CLineProfile(){};//to make sure derived objects are correctly deleted from a pointer to the base class 
        protected:
            Float64 m_nsigmasupport;
            std::string m_name = "NONE";//hack to avoid using dynamic casting

    };
    typedef std::shared_ptr<CLineProfile> CLineProfile_ptr;
    typedef std::vector<CLineProfile_ptr> TProfileList;

    inline
    const std::string& CLineProfile::GetName(){
        return m_name;
    }

    //no need to define a constructor here
    inline 
    CLineProfile::CLineProfile(const Float64 nsigmasupport): 
    m_nsigmasupport(nsigmasupport)
    {}

    inline 
    CLineProfile::CLineProfile(const Float64 nsigmasupport, const std::string name): 
    m_nsigmasupport(nsigmasupport),
    m_name(name)
    {}

    inline
    Float64 CLineProfile::GetNSigmaSupport(){
        return m_nsigmasupport;
    }

}
#endif