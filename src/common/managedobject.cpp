#include <epic/core/common/managedobject.h>

using namespace NSEpic;

CManagedObject::CManagedObject()
{	
    m_PointerCounter = 0;
}

CManagedObject::~CManagedObject()
{
}

UInt32 CManagedObject::GetRefCount()const 
{ 
    return m_PointerCounter; 
}
 
UInt32 CManagedObject::AddPointer() 
{
    m_PointerCounter++;
    return m_PointerCounter;
}

UInt32 CManagedObject::ReleasePointer()
{
    m_PointerCounter--;
    return m_PointerCounter;
}

