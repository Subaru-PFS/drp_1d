#ifndef _CORE_COMMON_DATATYPES_
#define _CORE_COMMON_DATATYPES_

#include <vector>
#include <string>

namespace NSEpic
{

/**
 * \ingroup Core
 * @{
 */

#ifndef NULL
#define NULL (0)
#endif

typedef int                     Int32;
typedef short                   Int16;
typedef signed char             Int8;
typedef long long               Int64;
typedef unsigned long long      UInt64;

typedef unsigned int            UInt32;
typedef unsigned short          UInt16;
typedef unsigned char           UInt8;
typedef	float                   Float32;
typedef	double                  Float64;
typedef char                    Char;
typedef unsigned char           Byte;
typedef void                    Void;
typedef unsigned int            Bool;
typedef const char*             String;

typedef std::vector<Float64>    TFloat64List;
typedef std::vector<Int64>      TInt64List;
typedef std::vector<Bool>       TBoolList;
typedef std::vector<UInt8>      TUInt8List;
typedef std::vector<std::string>      TStringList;

struct SPoint
{
    SPoint( )
    {
        X = 0.0;
        Y = 0.0;
    }

    SPoint( Float64 x, Float64 y)
    {
        X = x;
        Y = y;
    }
    Float64 X;
    Float64 Y;
};

typedef std::vector<SPoint>   TPointList;

/**@}*/

}

#endif
