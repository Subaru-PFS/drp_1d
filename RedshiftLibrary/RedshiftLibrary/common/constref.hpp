 

template < class DataType >
CConstRef<DataType>::CConstRef()
{
    m_Pointer = NULL;
}

template < class DataType >
CConstRef<DataType>::CConstRef( const DataType* n )
{
    m_Pointer = NULL;
    AcquirePointer( n );
}

template < class DataType >
CConstRef<DataType>::CConstRef( const CConstRef<DataType>& n )
{
    m_Pointer = NULL;
    AcquirePointer( n.m_Pointer );
}

template < class DataType >
CConstRef<DataType>& CConstRef<DataType>::operator = ( const CConstRef<DataType>& n )
{
    m_Pointer = NULL;
    AcquirePointer( n.m_Pointer );
}

template < class DataType >
CConstRef<DataType>::~CConstRef()
{
    ReleasePointer();
    m_Pointer = NULL;
}

template < class DataType >
const DataType* CConstRef<DataType>::GetPointer()
{
    return m_Pointer;
}

template < class DataType >
const DataType* CConstRef<DataType>::GetPointer() const
{
    return m_Pointer;
}

template < class DataType >
void CConstRef<DataType>::ReleasePointer()
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
void CConstRef<DataType>::AcquirePointer( const void* n )
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
Bool CConstRef<DataType>::IsValid() const
{
    return m_Pointer != NULL;
}

template <class DataType>
const DataType& CConstRef<DataType>::operator* () const
{
    //DebugAssertMsg( m_Pointer!=NULL, "CPointer::operator*\nNull pointer access");
    return *(m_Pointer);
}

template < class DataType >
const DataType* CConstRef<DataType>::operator-> () const
{ 
    //DebugAssertMsg( m_Pointer!=NULL, "CPointer::operator->\nNull pointer access");
    return (m_Pointer); 
}

template < class DataType >
CConstRef<DataType>& CConstRef<DataType>::operator = ( const DataType* n )
{
    if( n == (m_Pointer)) 
        return *this; 
    
    ReleasePointer();
    AcquirePointer( n ); 
    return *this;
}

template < class DataType >
CConstRef<DataType>::operator const DataType* () const
{
    return m_Pointer;
}
