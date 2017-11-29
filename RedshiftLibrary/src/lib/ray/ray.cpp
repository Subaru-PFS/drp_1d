#include <RedshiftLibrary/ray/ray.h>

using namespace NSEpic;
using namespace std;
#include <fstream>

CRay::CRay()
{
    m_Offset = 0.0;
}

CRay::CRay(const string& name, Float64 pos, UInt32 type, std::string profile, UInt32 force, Float64 amp, Float64 width, Float64 cut , Float64 posErr, Float64 sigmaErr, Float64 ampErr, const std::string& groupName, Float64 nominalAmp)
{
    m_Name = name;
    m_Pos = pos;
    m_Type = type;
    m_Force = force;

    m_Profile = profile;

    m_Amp = amp;
    m_Width = width;
    m_Cut = cut;

    m_PosFitErr = posErr;
    m_SigmaFitErr = sigmaErr;
    m_AmpFitErr = ampErr;


    m_GroupName = groupName;
    m_NominalAmplitude = nominalAmp;

    m_Offset = 0.0;
    m_OffsetFit = false;
}

CRay::~CRay()
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

std::string CRay::GetProfile() const
{
    return m_Profile;
}

bool CRay::SetProfile(std::string profile)
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



/**
 * This function converts lambda vacuum to lambda air
 *
 * see. Morton 1991, ApJS, 77, 119.
 */
Void CRay::ConvertVacuumToAir()
{
    Float64 s = (1e-4)/m_Pos;
    Float64 coeff = 1 + 8.34254*1e-5 + (2.406147*1e-2)/(130-s*s) + (1.5998*1e-4)/(38.9-s*s) ;

    m_Pos = m_Pos/coeff;

    return;
}

Void CRay::Save(  std::ostream& stream ) const
{
    stream << GetName() << "\t" << GetPosition() << "\t";
    if( GetIsStrong() )
        stream << "Strong"  << "\t";
    else
        stream << "Weak"  << "\t";
    stream << GetCut() << "\t" << GetWidth() << "\t" << GetAmplitude() << "\t" << GetPosFitError() << "\t" << GetSigmaFitError() << "\t" << GetAmpFitError() << "\t";
}

Void CRay::SaveDescription(  std::ostream& stream ) const
{
    stream << "#";
    stream << "Name" << "\t" << "Position" << "\t";
    stream << "Force"  << "\t";
    stream << "Cut" << "\t" << "Width" << "\t" << "Amp" << "\t" << "PosFitErr" << "\t" << "SigmaFitErr" << "\t" << "AmpFitErr" << "\t";
}


