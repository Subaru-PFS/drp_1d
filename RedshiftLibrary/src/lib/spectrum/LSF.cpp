#include "RedshiftLibrary/spectrum/LSF.h"
#include "RedshiftLibrary/spectrum/LSF_NISPSIM_2016.h"
#include "RedshiftLibrary/spectrum/LSF_NISPVSSPSF_201707.h"
#include "RedshiftLibrary/spectrum/LSFConstantResolution.h"
#include "RedshiftLibrary/spectrum/LSFConstantWidth.h"
#include "RedshiftLibrary/log/log.h"

using namespace NSEpic;
using namespace std;

CLSF::CLSF(TLSFType name):m_name(name){} //To be deleted once child classes are adapted
CLSF::CLSF(TLSFType name, CLineProfile_ptr profile):
m_name(name),
m_profile(profile)
{}
Float64 CLSF::GetLineProfile (Float64 lambda, Float64 lambda0){
    return m_profile->GetLineProfile(lambda, lambda0, GetWidth());
}


std::shared_ptr<CLSF>  CLSF::make_LSF(const std::string lsfType, const TLSFArguments& args)
{
    if( lsfType == "GaussianConstantWidth"){
        return std::make_shared<CLSFGaussianConstantWidth>(args.width);
    } else if( lsfType == "GaussianConstantResolution"){
        return std::make_shared<CLSFGaussianConstantResolution>(args.resolution);
    }else if( lsfType == "GaussianNISPSIM2016"){
        return std::make_shared<CLSFGaussianNISPSIM2016>();
    } else if( lsfType == "GaussianNISPVSSPSF201707"){
        return std::make_shared<CLSFGaussianNISPVSSPSF201707>(args.sourcesize);
    } else {
        Log.LogError("Unknown lsfType %s", lsfType.c_str());
        throw std::runtime_error("Unknown lsfType");
    }
}
