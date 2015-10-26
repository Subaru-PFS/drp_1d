#include <epic/redshift/spectrum/template/template.h>

#include <epic/redshift/common/mask.h>

using namespace NSEpic;
using namespace std;

CTemplate::CTemplate( )
{

}

CTemplate::CTemplate( const std::string& name, const std::string& category ) :
    m_Category( category ),
    m_Name( name )
{

}

CTemplate::~CTemplate()
{

}

const std::string& CTemplate::GetName() const
{
    return m_Name;
}

const std::string& CTemplate::GetCategory() const
{
    return m_Category;
}

