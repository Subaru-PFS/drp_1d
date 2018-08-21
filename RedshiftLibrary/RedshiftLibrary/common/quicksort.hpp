
template < class T >
CQuickSort<T>::CQuickSort()
{
}

template < class T >
CQuickSort<T>::~CQuickSort()
{
}

template < class T >
void CQuickSort<T>::Sort( T* value, Int32 n ) const
{
    Sort(value, 0, n );
}

template < class T >
void CQuickSort<T>::Sort( T* value, Int32 beg, Int32 end ) const
{
    if (end > beg + 1)
    {
        Float64 piv = value[beg];
        Int32 l = beg + 1, r = end;
        while (l < r)
        {
            if (value[l] <= piv)
                l++;
            else
                swap( value[l], value[--r] );
        }
        swap( value[--l], value[beg] );

        Sort(value, beg, l);
        Sort(value, r, end);
    }
}

template < class T >
void CQuickSort<T>::SortIndexes( const T* value, Int32* index, Int32 n ) const
{
    vector<T> tmp( n );

    for( Int32 i=0; i<n; i++ )
    {
        index[i] = i;
        tmp[i] = value[i];
    }

    SortIndexes(tmp.data(), index, 0, n );
}

template < class T >
void CQuickSort<T>::SortIndexes( T* value, Int32* index, Int32 beg, Int32 end ) const
{
    if (end > beg + 1)
    {
        Float64 piv = value[beg];
        Int32 l = beg + 1, r = end;
        while (l < r) 
        {
            if (value[l] <= piv)
                l++;
            else
            {
                --r;
                swap( value[l], value[r] );
                swap( index[l], index[r] );
            }
        }

        --l;
        swap( value[l], value[beg] );
        swap( index[l], index[beg] );

        SortIndexes( value, index, beg, l );
        SortIndexes( value, index, r, end );
    }
}
