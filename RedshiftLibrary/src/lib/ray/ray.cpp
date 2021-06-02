#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/ray/ray.h>

using namespace NSEpic;
using namespace std;
#include <fstream>

CRay::CRay():m_Offset(0.){}

CRay::CRay(const string& name,
           Float64 pos, UInt32 type,
           std::shared_ptr<CLineProfile> profile,
           UInt32 force,
           Float64 amp,
           Float64 width,
           Float64 cut ,
           Float64 posErr,
           Float64 sigmaErr,
           Float64 ampErr,
           const std::string& groupName,
           Float64 nominalAmp,
           const string &velGroupName,
	       Int32 id):
m_Name(name),
m_Pos(pos),
m_Type(type),
m_Force(force),
m_Amp(amp),
m_Width(width),
m_Cut(cut),
m_Profile(profile),
m_PosFitErr(posErr),
m_SigmaFitErr(sigmaErr),
m_AmpFitErr(ampErr),
m_GroupName(groupName),
m_NominalAmplitude(nominalAmp),
m_VelGroupName(velGroupName),
m_id(id),
m_Offset(0.),
m_OffsetFit(false)
{
}

bool CRay::operator < (const CRay& str) const
{
    if(m_Pos == str.m_Pos){
        return (m_Amp < str.m_Amp);
    }else{
        return (m_Pos < str.m_Pos);
    }
}

bool CRay::operator != (const CRay& str) const
{
    if(m_Pos == str.m_Pos){
        return (m_Amp != str.m_Amp);
    }else{
        return (m_Pos != str.m_Pos);
    }
}

void CRay::SetAsymParams(TAsymParams asymParams)
{
    if(!m_Profile)
        throw runtime_error("CRay::SetAsymParams: lineprofile is not initialized");
    m_Profile->SetAsymParams(asymParams);
}
void CRay::resetAsymFitParams()
{
    if(!m_Profile)
        throw runtime_error("CRay::resetAsymParams: lineprofile is not initialized");
    m_Profile->resetAsymFitParams();
}
const TAsymParams CRay::GetAsymParams()
{
    if(!m_Profile)
        throw std::runtime_error("CRay::GetAsymParams: lineprofile is not initialized");
    return m_Profile->GetAsymParams();
}

Bool CRay::GetIsStrong() const
{
    return m_Force == nForce_Strong;
}

Bool CRay::GetIsEmission() const
{
    return m_Type == nType_Emission;
}

Int32 CRay::GetType() const
{
    return m_Type;
}

std::shared_ptr<CLineProfile> CRay::GetProfile() const
{
    if(!m_Profile)
        throw runtime_error("Current Ray does not have a set profile ");
    return m_Profile;
     
}

bool CRay::SetProfile(const std::shared_ptr<CLineProfile>& profile)
{
    m_Profile = profile;
    return true;
}

Int32 CRay::GetForce() const
{
    return m_Force;
}

Float64 CRay::GetPosition() const
{
    return m_Pos;
}

Float64 CRay::GetOffset() const
{
    return m_Offset;
}

bool CRay::SetOffset(Float64 val)
{
    if(m_OffsetFit)
    {
        m_Offset = val;
        return true;
    }else{
        return false;
    }
}

bool CRay::GetOffsetFitEnabled() const
{
    return m_OffsetFit;
}

bool CRay::EnableOffsetFit(bool val)
{
    m_OffsetFit = val;
    return true;
}

Float64 CRay::GetAmplitude() const
{
    return m_Amp;
}

Float64 CRay::GetWidth() const
{
    return m_Width;
}

Float64 CRay::GetCut() const
{
    return m_Cut;
}

Float64 CRay::GetPosFitError() const
{
    return m_PosFitErr;
}

Float64 CRay::GetSigmaFitError() const
{
    return m_SigmaFitErr;
}

Float64 CRay::GetAmpFitError() const
{
    return m_AmpFitErr;
}

const std::string& CRay::GetName() const
{
    return m_Name;
}

const std::string& CRay::GetGroupName() const
{
    return m_GroupName;
}

const Float64 CRay::GetNominalAmplitude() const
{
    return m_NominalAmplitude;
}

const std::string& CRay::GetVelGroupName() const
{
    return m_VelGroupName;
}

Int32 CRay::GetID() const
{
    return m_id;
}

/**
 * This function converts lambda vacuum to lambda air
 *
 * see. Morton 2000 ApJS, 130, 403-436 (https://iopscience.iop.org/article/10.1086/317349)
 * 
 * Accomodate the equation to PFS by setting s to cte = 0.8,  i.e., average factor for PFS wavelength.
 * By doing so, we zero out the variability in wavelength caused by the vacuum-to-air conversion.
 * This allows the conversion on the input line-catalog and then the redshift 
 * (if not, we would have to apply the conversion after redshifting) 
 */
void CRay::ConvertVacuumToAir()
{
    Float64 s = 0.8; //(1e4)/m_Pos;
    Float64 coeff = 1 + 8.34254*1e-5 + (2.406147*1e-2)/(130-s*s) + (1.5998*1e-4)/(38.9-s*s) ;

    m_Pos = m_Pos/coeff;

    return;
}

void CRay::Save(  std::ostream& stream ) const
{
    stream << GetName() << "\t" << GetPosition() << "\t";
    if( GetIsStrong() )
        stream << "Strong"  << "\t";
    else
        stream << "Weak"  << "\t";
    stream << GetCut() << "\t" << GetWidth() << "\t" << GetAmplitude() << "\t" << GetPosFitError() << "\t" << GetSigmaFitError() << "\t" << GetAmpFitError() << "\t";
}

void CRay::SaveDescription(  std::ostream& stream ) const
{
    stream << "#";
    stream << "Name" << "\t" << "Position" << "\t";
    stream << "Force"  << "\t";
    stream << "Cut" << "\t" << "Width" << "\t" << "Amp" << "\t" << "PosFitErr" << "\t" << "SigmaFitErr" << "\t" << "AmpFitErr" << "\t";
}


