#include "RedshiftLibrary/spectrum/LSF.h"
#include "RedshiftLibrary/spectrum/LSF_NISPSIM_2016.h"
#include "RedshiftLibrary/spectrum/LSF_NISPVSSPSF_201707.h"
#include "RedshiftLibrary/spectrum/LSFConstantResolution.h"
#include "RedshiftLibrary/spectrum/LSFConstantWidth.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/common/range.h"
using namespace NSEpic;
using namespace std;

CLSF::CLSF(TLSFType name):m_name(name){} //To be deleted once child classes are adapted
    CLSF::CLSF(TLSFType name, CLineProfile_ptr profile):
    m_name(name),
    m_profile(profile)
{

}

Float64 CLSF::GetLineProfile (Float64 lambda, Float64 lambda0){
    return m_profile->GetLineProfile(lambda, lambda0, GetWidth());
}
