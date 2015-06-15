#include <epic/redshift/spectrum/template/template.h>

#include <epic/redshift/common/mask.h>

using namespace NSEpic;
using namespace std;

CTemplate::CTemplate() :
    m_Category( nCategory_None )
{

}

CTemplate::CTemplate( const char* name, ECategory category ) :
    m_Category( category ),
    m_Name( name )
{

}

CTemplate::~CTemplate()
{

}

const char* CTemplate::GetCategoryName( ECategory cat )
{
    static const char* emission = "emission";
    static const char* galaxy = "galaxy";
    static const char* star = "star";
    static const char* qso = "qso";
    static const char* none = "none";

    if( cat == nCategory_Emission )
        return emission;
    else if( cat == nCategory_Galaxy )
        return galaxy;
    else if( cat == nCategory_Star )
        return star;
    else if( cat == nCategory_Qso )
        return qso;

    return none;
}

const std::string& CTemplate::GetName() const
{
    return m_Name;
}

CTemplate::ECategory CTemplate::GetCategory() const
{
    return m_Category;
}

