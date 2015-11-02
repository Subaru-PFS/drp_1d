#ifndef _CORE_COMMON_CONSTREF_
#define _CORE_COMMON_CONSTREF_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/managedobject.h>

namespace NSEpic
{

/**
 * \ingroup Core
 * Smart pointer to use in conjunction with CManagedObject
 *
 */
template < typename DataType >
class CConstRef
{

public:
    
    CConstRef();
    CConstRef( const DataType* n );
    CConstRef( const CConstRef<DataType>& n );
    ~CConstRef();

    CConstRef<DataType>&      operator = ( const CConstRef<DataType>& n );
    CConstRef<DataType>&      operator = ( const DataType* n );
    const DataType&           operator * () const;
    const DataType*           operator-> () const ;
    operator const DataType* () const ;
    
    const DataType*           GetPointer();
    const DataType*           GetPointer() const;

    Bool                        IsValid() const;

private:

    Void              ReleasePointer();
    Void              AcquirePointer( const Void* n );
    Void              SetPointer( const CManagedObject* p );

    const DataType*	m_Pointer;

};

#include <epic/core/common/constref.hpp>


}

#endif
