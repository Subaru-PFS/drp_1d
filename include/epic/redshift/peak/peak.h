#ifndef _REDSHIFT_PEAK_PEAK_
#define _REDSHIFT_PEAK_PEAK_

#include <epic/core/common/managedobject.h>
#include <epic/core/common/datatypes.h>

namespace __NS__
{

class CPeak
{

public:

    CPeak( );
    CPeak( Int32 startIndex, Int32 stopIndex );
    CPeak( const CPeak& other );

    Void Set( Int32 startIndex, Int32 stopIndex );

private:

    Int32 m_StartIndex;
    Int32 m_StopIndex;

};


typedef std::vector<CPeak> TPeakList;

}

#endif
