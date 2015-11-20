#include <epic/redshift/continuum/standard.h>

using namespace NSEpic;
using namespace std;

/**
 * Empty constructor.
 */
CContinuumStandard::CContinuumStandard()
{

}

/**
 * Empty destructor.
 */
CContinuumStandard::~CContinuumStandard()
{

}

/**
 * Does nothing, but returns true.
 */
Bool CContinuumStandard::Remove( CSpectrum& s )
{
    return true;
}
