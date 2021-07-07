#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/ray/rule.h"

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

CRule::~CRule ( )
{

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

/**
 * Returns the Logs string.
 **/
std::string CRule::GetLogs( )
{
    std::string _logs = (std::string)Logs;
    Logs.clear();
    return _logs;
}
