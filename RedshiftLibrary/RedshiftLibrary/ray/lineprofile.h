#ifndef _REDSHIFT_LINE_PROFILE_
#define _REDSHIFT_LINE_PROFILE_

#include <string>
#include "RedshiftLibrary/common/datatypes.h"
#include <cmath>
namespace NSEpic
{
    /**
    * struct that holds ASYMFIXED profile parameters
    */
    typedef struct {
        Float64 sigma, alpha, delta;
    } TAsymParams;
    
    enum TProfile
    {
        NONE,
        SYM,
        LOR,
        ASYM,
        ASYMFIT,
        //ASYMFIXED,//doesnt exist anymore since merged with ASYM
    };
    /**
     * \ingroup Redshift
     * Abstract class for different line profiles
     */
    class CLineProfile
    {
        
        public:
            CLineProfile(const Float64 nsigmasupport=8.0);
            CLineProfile(const Float64 nsigmasupport=8.0, const TProfile=NONE);
            virtual Float64 GetLineProfile(Float64 x, Float64 x0, Float64 sigma)=0;
            virtual Float64 GetLineFlux( Float64 A, Float64 sigma)=0;
            virtual Float64 GetLineProfileDerivZ(Float64 x, Float64 lambda0, Float64 redshift, Float64 sigma)=0;
            virtual Float64 GetLineProfileDerivSigma(Float64 x, Float64 x0, Float64 sigma)=0;
            virtual Float64 GetNSigmaSupport() const;

            const TProfile& GetName();
            virtual const TAsymParams GetAsymParams(){return {NAN,NAN,NAN};};
            virtual Float64 GetAsymDelta();
            virtual Bool isAsymFit();
            virtual Bool isAsymFixed();
            virtual void SetAsymParams(TAsymParams params){};
            virtual void resetAsymFitParams();
            CLineProfile(const CLineProfile & other) = default; 
            CLineProfile(CLineProfile && other) = default; 
            CLineProfile& operator=(const CLineProfile& other) = default;  
            CLineProfile& operator=(CLineProfile&& other) = default; 

            virtual ~CLineProfile(){};//to make sure derived objects are correctly deleted from a pointer to the base class 
        protected:
            Float64 m_nsigmasupport;
            const TProfile m_name;//hack to avoid using dynamic casting

    };
    typedef std::shared_ptr<CLineProfile> CLineProfile_ptr;
    typedef std::vector<CLineProfile_ptr> TProfileList;

    inline
    const TProfile& CLineProfile::GetName(){
        return m_name;
    }

    //no need to define a constructor here
    inline
    CLineProfile::CLineProfile(const Float64 nsigmasupport): 
    m_nsigmasupport(nsigmasupport),
    m_name(NONE)
    {}

    inline
    CLineProfile::CLineProfile(const Float64 nsigmasupport, const TProfile name): 
    m_nsigmasupport(nsigmasupport),
    m_name(name)
    {}

    inline
    Float64 CLineProfile::GetNSigmaSupport()const{
        return m_nsigmasupport;
    }
    inline
    Float64 CLineProfile::GetAsymDelta(){
        return 0.;//default. Mainly used for asylfit/fixed
    }
    inline
    Bool CLineProfile::isAsymFit(){
        return 0; //default to no
    }
    inline
    Bool CLineProfile::isAsymFixed(){
        return 0; //default to no
    }
    inline
    void CLineProfile::resetAsymFitParams(){
        return;
    }
}
#endif