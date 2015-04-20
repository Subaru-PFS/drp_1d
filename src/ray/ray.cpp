#include <epic/redshift/ray/ray.h>

using namespace NSEpic;
using namespace std;
#include <fstream>

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
    return m_Type == nForce_Strong;
}

Bool CRay::GetIsEmission() const
{
    return m_Type == nType_Emission;
}

Float64 CRay::GetPosition() const
{
    return m_Pos;
}

const std::string& CRay::GetName() const
{
    return m_Name;
}

const std::string CRay::GetDescription() const
{
    std::string strForce = "Weak";
    if(GetIsStrong()){
       strForce = "Strong";
    }
    char tmpChar[256];
    sprintf(tmpChar, "%f\t%s", GetPosition(),  strForce.c_str() );

    std::string tmpStr = tmpChar;

    return tmpStr;
}

