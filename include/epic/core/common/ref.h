#ifndef _CORE_COMMON_REF_
#define _CORE_COMMON_REF_

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
class CRef
{

public:
    
    CRef();
    CRef( DataType* n );
    CRef( const CRef<DataType>& n );
    ~CRef();
    
    CRef<DataType>&     operator = ( DataType* n );
    CRef<DataType>&     operator = ( const CRef<DataType>& n );
    DataType&           operator * () const;
    DataType*           operator-> () const ;
    operator DataType* () const ;

    DataType*           GetPointer() const;

    Bool                IsValid() const;

private:

    Void              ReleasePointer();
    Void              AcquirePointer( const Void* n );
    Void              SetPointer( CManagedObject* p );

    DataType*	m_Pointer;

};

#include <epic/core/common/ref.hpp>


}

#endif
