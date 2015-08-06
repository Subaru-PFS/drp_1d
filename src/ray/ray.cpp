#include <epic/redshift/ray/ray.h>

using namespace NSEpic;
using namespace std;
#include <fstream>

CRay::CRay()
{

}

CRay::CRay( const string& name, Float64 pos, UInt32 type, UInt32 force, Float64 amp, Float64 width, Float64 cut )
{
    m_Name = name;
    m_Pos = pos;
    m_Type = type;
    m_Force = force;

    m_Amp = amp;
    m_Width = width;
    m_Cut = cut;
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

Int32 CRay::GetForce() const
{
    return m_Force;
}

Float64 CRay::GetPosition() const
{
    return m_Pos;
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

const std::string& CRay::GetName() const
{
    return m_Name;
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
    stream << GetCut() << "\t" << GetWidth() << "\t";
}


