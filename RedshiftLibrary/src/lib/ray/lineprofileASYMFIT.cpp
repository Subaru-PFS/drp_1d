#include "RedshiftLibrary/ray/lineprofileASYMFIT.h"
#include "RedshiftLibrary/log/log.h"
using namespace NSEpic;
using namespace std;

CLineProfileASYMFIT::CLineProfileASYMFIT(const Float64 nsigmasupport, TAsymParams params, const std::string centeringMethod):
CLineProfileASYM(ASYMFIT, nsigmasupport, params, centeringMethod)
{

}
 
void CLineProfileASYMFIT::SetAsymParams(TAsymParams params)
{
    if(std::isnan(params.sigma) || std::isnan(params.alpha) || std::isnan(params.delta)){
        Log.LogWarning("CLineProfileASYMFIT::setAsymParams AsymFit params are NaN");
    } 
    m_asym_sigma_coeff = params.sigma;
    m_asym_alpha = params.alpha;
    m_asym_delta = params.delta;
}

void CLineProfileASYMFIT::resetAsymFitParams()
{
    m_asym_sigma_coeff = 2.;
    m_asym_alpha = 0.;
    m_asym_delta = 0.;
}
 
Bool CLineProfileASYMFIT::isAsymFit()
{
    return 1;
}
 
Bool CLineProfileASYMFIT::isAsymFixed()
{
    return 0;
}