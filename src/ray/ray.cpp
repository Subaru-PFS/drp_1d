#include <epic/redshift/ray/ray.h>

using namespace __NS__;
using namespace std;

CRay::CRay()
{

}

CRay::CRay( const string& name, Float64 pos, UInt32 type )
{
    m_Name = name;
    m_Pos = pos;
    m_Type = type;
}

CRay::~CRay()
{

}

Bool CRay::GetIsStrong() const
{
    return m_Type & nType_Strong;
}

Bool CRay::GetIsEmission() const
{
    return m_Type & nType_Emission;
}

Float64 CRay::GetPosition() const
{
    return m_Pos;
}

const std::string& CRay::GetName() const
{
    return m_Name;
}
