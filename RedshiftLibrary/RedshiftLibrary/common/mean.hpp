template < class T >
CMean<T>::CMean()
{
}

template < class T >
CMean<T>::~CMean()
{
}

template < class T >
T CMean<T>::Find( const T* x, Int32 n ) 
{
    T sum = 0;

    for ( Int32 i=0; i<n ;i++ )
    {
        sum += x[i];
    }

    return (T) (sum/n);
}

