 

template < class DataType >
CRef<DataType>::CRef()
{
    m_Pointer = NULL;
}

template < class DataType >
CRef<DataType>::CRef( DataType* n )
{
    m_Pointer = NULL;
    AcquirePointer( n );
}

template < class DataType >
CRef<DataType>::CRef( const CRef<DataType>& n )
{
    m_Pointer = NULL;
    AcquirePointer( n.m_Pointer );
}

template < class DataType >
CRef<DataType>::~CRef()
{
    ReleasePointer();
    m_Pointer = NULL;
}

template < class DataType >
DataType* CRef<DataType>::GetPointer() const
{
    return m_Pointer;
}

template < class DataType >
void CRef<DataType>::ReleasePointer()
{
    if( m_Pointer == NULL ) 
        return;
    
    CManagedObject* m = (CManagedObject*) m_Pointer;
    
    m->ReleasePointer();
    if( m->GetRefCount() == 0 )   
    {
       delete m;
    }

    m_Pointer = NULL;
}

template < class DataType >
void CRef<DataType>::AcquirePointer( const void* n )
{
    if( m_Pointer != NULL ) 
        ReleasePointer(); 

    if( n == NULL ) 
    {
        m_Pointer = NULL;
        return;
    }

    m_Pointer = (DataType*) n;

    CManagedObject* m = (CManagedObject*) m_Pointer;

    m->AddPointer(); 
}

template <class DataType>
Bool CRef<DataType>::IsValid() const
{
    return m_Pointer != NULL;
}

template <class DataType>
DataType& CRef<DataType>::operator* () const
{
    //DebugAssertMsg( m_Pointer!=NULL, "CPointer::operator*\nNull pointer access");
    return *(m_Pointer);
}

template < class DataType >
DataType* CRef<DataType>::operator-> () const 
{ 
    //DebugAssertMsg( m_Pointer!=NULL, "CPointer::operator->\nNull pointer access");
    return (m_Pointer); 
}

template < class DataType >
CRef<DataType>& CRef<DataType>::operator = ( DataType* n )
{
    if( n == (m_Pointer)) 
        return *this; 
    
    ReleasePointer();
    AcquirePointer( n ); 
    return *this;
}

template < class DataType >
CRef<DataType>& CRef<DataType>::operator = ( const CRef<DataType>& n )
{
    if( n.m_Pointer == (m_Pointer))
        return *this;

    ReleasePointer();
    AcquirePointer( n.m_Pointer );
    return *this;
}

template < class DataType >
CRef<DataType>::operator DataType* () const
{
    return m_Pointer;
}
