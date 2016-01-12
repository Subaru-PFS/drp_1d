#include <epic/redshift/spectrum/template/template.h>

#include <epic/redshift/common/mask.h>

using namespace NSEpic;
using namespace std;

/**
 * Constructor, empty.
 */
CTemplate::CTemplate( )
{

}

/**
 * Constructor, assigns values to members.
 */
CTemplate::CTemplate( const std::string& name, const std::string& category ) :
    m_Category( category ),
    m_Name( name )
{

}

/**
 * Destructor, empty.
 */
CTemplate::~CTemplate()
{

}

/**
 * Returns the value stored in m_Name.
 */
const std::string& CTemplate::GetName() const
{
    return m_Name;
}

/**
 * Returns the value stored in m_Category.
 */
const std::string& CTemplate::GetCategory() const
{
    return m_Category;
}

