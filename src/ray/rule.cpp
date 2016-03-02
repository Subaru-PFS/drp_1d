#include <epic/core/common/datatypes.h>
#include <epic/redshift/ray/rule.h>

#include <cstdarg>

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
void CRule::Apply( CLineModelElementList& LineModelElementList )
{
  if ( ! Enabled )
    {
      return;
    }
  if ( Check( LineModelElementList ) )
    {
      return;
    }
  Correct( LineModelElementList );
}
