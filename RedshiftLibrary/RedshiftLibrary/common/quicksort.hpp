// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================

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
