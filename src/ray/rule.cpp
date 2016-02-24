#include <epic/redshift/ray/rule.h>

using namespace NSEpic;
using namespace std;

/**
 * \brief Sets itself as disabled per default.
 **/
CRule::CRule ( )
{
  Enabled = false;
}
    
/**
 * Checks the elements and corrects if necessary.
 **/
void CRule::Apply( std::vector<boost::shared_ptr<CLineModelElement> > LinemodelElements )
{
  if ( ! Enabled )
    {
      return;
    }
  if ( Check( LinemodelElements ) )
    {
      return;
    }
  Correct( LinemodelElements );
}
