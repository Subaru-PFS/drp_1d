#ifndef _CORE_COMMON_MANAGEDOBJECT_
#define _CORE_COMMON_MANAGEDOBJECT_

#include <epic/core/common/datatypes.h>


#define DEFINE_MANAGED_OBJECT( T )


#define IMPLEMENT_MANAGED_OBJECT( T )

#define IMPLEMENT_MANAGED_OBJECT_NOT_INSTANCIABLE( T )


namespace NSEpic
{

class CTypeInfo;

/**
 * \ingroup Core
 * Base class for any managed object.
 *
 * Managed objects can be used in conjonction with smart point CRef and also provide Type information CTypeInfo
 */
class CManagedObject
{

public:
    
    CManagedObject();
    virtual ~CManagedObject();
    
    UInt32      GetRefCount() const; 
    
    virtual UInt32  AddPointer() ;
    virtual UInt32  ReleasePointer();
    
private:
    
    UInt32    m_PointerCounter;
};
    

}

#endif
