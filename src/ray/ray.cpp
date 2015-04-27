#include <epic/redshift/ray/ray.h>

using namespace NSEpic;
using namespace std;
#include <fstream>

CRay::CRay()
{

}

CRay::CRay( const string& name, Float64 pos, UInt32 type, UInt32 force )
{
    m_Name = name;
    m_Pos = pos;
    m_Type = type;
    m_Force = force;
}

CRay::~CRay()
{

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

const std::string& CRay::GetName() const
{
    return m_Name;
}

Void CRay::Save( std::ostream& stream ) const
{
    stream << GetName() << "\t" << GetPosition() << "\t";
    if( GetIsStrong() )
        stream << "Strong";
    else
        stream << "Weak";
}


