#ifndef _REDSHIFT_LINE_PROFILE_ASYM_
#define _REDSHIFT_LINE_PROFILE_ASYM_
#include <string>
#include <math.h>
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/ray/lineprofile.h"
namespace NSEpic
{
    /**
     * \ingroup Redshift
     */
    class CLineProfileASYM: public CLineProfile
    {
        public:
            CLineProfileASYM(const Float64 nsigmasupport = 8.0, 
                            const TAsymParams params = {1., 4.5, 0.}, 
                            const std::string centeringMethod = "none");
            CLineProfileASYM(const TProfile pltype,
                            const Float64 nsigmasupport = 8.0, 
                            const TAsymParams params = {2., 2.5, 0.}, 
                            const std::string centeringMethod = "mean");//mainly called by asymfit
 
            Float64 GetLineProfile(Float64 x, Float64 x0, const Float64 sigma) override;
            Float64 GetLineFlux(Float64 A, const Float64 sigma) override;
            Float64 GetLineProfileDerivZ(Float64 x, Float64 x0, Float64 redshift, const Float64 sigma) override;
            Float64 GetLineProfileDerivSigma(Float64 x, Float64 x0, const Float64 sigma) override;
            Float64 GetNSigmaSupport() override;

            Float64 GetAsymDelta();
            const TAsymParams  GetAsymParams() override;
            virtual Bool    isAsymFixed() override;
            virtual Bool    isAsymFit()   override;

            virtual ~CLineProfileASYM() = default;
            CLineProfileASYM(const CLineProfileASYM & other) = default; 
            CLineProfileASYM(CLineProfileASYM && other) = default; 
            CLineProfileASYM& operator=(const CLineProfileASYM& other) = default;  
            CLineProfileASYM& operator=(CLineProfileASYM&& other) = default; 
        protected: 
            Bool isValid();
            Float64 m_asym_sigma_coeff = 1.0;//vs 2. for asymFit/Fixed
            Float64 m_asym_alpha = 4.5;
            Float64 m_asym_delta = 0.;
            std::string m_centeringMethod = "none";
            Float64 m_constSigma = 1;//vs 2.5 for AsymFit and AsymFixed
        private:
            Float64 GetXSurc(Float64 xc, Float64& sigma, Float64& xsurc);
    };
}
#endif