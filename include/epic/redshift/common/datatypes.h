#ifndef _CORE_REDSHIFT_DATATYPES_
#define _CORE_REDSHIFT_DATATYPES_

#include <epic/core/common/datatypes.h>

#include <vector>

namespace NSEpic
{

typedef UInt8                   Mask;
typedef Float64                 Redshift;
typedef Float64                 Sample;
typedef std::vector<Mask>     TMaskList;
typedef std::vector<Redshift>   TRedshiftList;
typedef std::vector<Sample>     TAxisSampleList;

}

#endif
