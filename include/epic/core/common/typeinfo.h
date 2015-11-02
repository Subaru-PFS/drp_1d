#ifndef _CORE_COMMON_TYPEINFO_
#define _CORE_COMMON_TYPEINFO_

#include <epic/core/common/datatypes.h>

namespace NSEpic
{

class CSerializable;

/**
 * \ingroup Core
 * Provide Run time type information
 */
class CTypeInfo
{

public:

    typedef CSerializable* (*InstantiationCallback)();

    CTypeInfo( String typeName, InstantiationCallback cb );
    ~CTypeInfo();

    CSerializable* Allocate() const;

    Bool operator==( const CTypeInfo & id ) const;
    Bool operator!=( const CTypeInfo & id ) const;

    Bool IsA( String fullTypeName ) const;
    Bool IsA( const CTypeInfo & info ) const;


    String                  Name;
    InstantiationCallback   Callback;

private:

    friend class CTypeInfoRegistry;

    UInt16                  TypeID;
    CTypeInfo*              NextTypeInfo;

    static UInt16           TypeIDCount;
};

}

#endif
