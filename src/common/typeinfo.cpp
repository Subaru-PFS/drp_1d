#include <epic/core/common/typeinfo.h>


#include <cstring>
#include <typeinfo>
#include <iostream>

using namespace NSEpic;

UInt16 CTypeInfo::TypeIDCount = 0;

CTypeInfo::CTypeInfo( String typeName, InstantiationCallback cb ) :
    Name( typeName ),
    Callback( cb ),
    NextTypeInfo( NULL ),
    TypeID( TypeIDCount++ )
{
}

CTypeInfo::~CTypeInfo()
{

}

CSerializable* CTypeInfo::Allocate() const
{
    if( Callback )
        return Callback();

    return NULL;
}

Bool CTypeInfo::operator==( const CTypeInfo & id ) const
{
    return IsA( id );
}

Bool CTypeInfo::operator!=( const CTypeInfo & id ) const
{
    return ! IsA( id );
}

Bool CTypeInfo::IsA( String fullTypeName ) const
{
    return ( strcmp( fullTypeName, (String)Name ) == 0 );
}

Bool CTypeInfo::IsA( const CTypeInfo & info ) const
{
    return TypeID == info.TypeID;
}
