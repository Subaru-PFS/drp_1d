

#define SORT_PAIR(a, b)                                                        \
    {                                                                          \
        if ((a) > (b))                                                         \
            swap((a), (b));                                                    \
    }

template <class T> CMedian<T>::CMedian() {}

template <class T> CMedian<T>::~CMedian() {}

template <class T> T CMedian<T>::BeersFind(const T *x, Int32 n)
{
    T medVal;
    Int32 i, n2;

    vector<T> xwork(n);

    for (i = 0; i < n; i++)
        xwork[i] = x[i];

    CQuickSort<T> quickSort;
    quickSort.Sort(xwork.data(), n);

    n2 = (Int32)(n * 0.5);
    if (2 * n2 == n)
        medVal = (0.5 * (xwork[n2 - 1] + xwork[n2]));
    else
        medVal = xwork[n2];

    return medVal;
}

template <class T> T CMedian<T>::Find(const T *a, Int32 n)
{
    T medianVal;

    switch (n)
    {
    case 0:
    case 1:
        return a[0];
        break;

    case 3:
        medianVal = Opt3Find(a);
        break;

    case 5:
        medianVal = Opt5Find(a);
        break;

    case 7:
        medianVal = Opt7Find(a);
        break;

    case 9:
        medianVal = Opt9Find(a);
        break;

    default:
        medianVal = ((n > MEDIAN_FAST_OR_BEERS_THRESHOLD) ? FastFind(a, n)
                                                          : BeersFind(a, n));
        break;
    }

    return medianVal;
}

/**
   Implementation of the fast median search based upon Niklaus Wirth works.

   See paper:

   Fast median search: an ANSI C implementation
   Nicolas Devillard - ndevilla AT free DOT fr
   July 1998
*/
template <class T> T CMedian<T>::FastFind(const T *a, Int32 n)
{
    // Define custom K depending if n is odd or even
    Int32 k = ((n & 1) ? (n / 2.0) : ((n / 2.0) - 1));

    Int32 i, j, l, m;
    T x;

    vector<T> a_copy(n);

    for (Int32 i = 0; i < n; i++)
        a_copy[i] = a[i];

    l = 0;
    m = n - 1;
    while (l < m)
    {
        x = a_copy[k];
        i = l;
        j = m;
        do
        {
            while (a_copy[i] < x)
            {
                i++;
            }
            while (x < a_copy[j])
            {
                j--;
            }
            if (i <= j)
            {
                swap(a_copy[i], a_copy[j]);
                i++;
                j--;
            }
        } while (i <= j);

        if (j < k)
        {
            l = i;
        }
        if (k < i)
        {
            m = j;
        }
    }

    return (a_copy[k]);
}

template <class T> T CMedian<T>::Opt3Find(const T *ip)
{
    T p[3];
    for (Int32 i = 0; i < 3; i++)
    {
        p[i] = ip[i];
    }

    SORT_PAIR(p[0], p[1]);
    SORT_PAIR(p[1], p[2]);
    SORT_PAIR(p[0], p[1]);
    return (p[1]);
}

template <class T> T CMedian<T>::Opt5Find(const T *ip)
{
    T p[5];
    for (Int32 i = 0; i < 5; i++)
    {
        p[i] = ip[i];
    }
    SORT_PAIR(p[0], p[1]);
    SORT_PAIR(p[3], p[4]);
    SORT_PAIR(p[0], p[3]);
    SORT_PAIR(p[1], p[4]);
    SORT_PAIR(p[1], p[2]);
    SORT_PAIR(p[2], p[3]);
    SORT_PAIR(p[1], p[2]);
    return (p[2]);
}

template <class T> T CMedian<T>::Opt7Find(const T *ip)
{
    T p[7];
    for (Int32 i = 0; i < 7; i++)
    {
        p[i] = ip[i];
    }

    SORT_PAIR(p[0], p[5]);
    SORT_PAIR(p[0], p[3]);
    SORT_PAIR(p[1], p[6]);
    SORT_PAIR(p[2], p[4]);
    SORT_PAIR(p[0], p[1]);
    SORT_PAIR(p[3], p[5]);
    SORT_PAIR(p[2], p[6]);
    SORT_PAIR(p[2], p[3]);
    SORT_PAIR(p[3], p[6]);
    SORT_PAIR(p[4], p[5]);
    SORT_PAIR(p[1], p[4]);
    SORT_PAIR(p[1], p[3]);
    SORT_PAIR(p[3], p[4]);
    return (p[3]);
}

template <class T> T CMedian<T>::Opt9Find(const T *ip)
{
    T p[9];
    for (Int32 i = 0; i < 9; i++)
    {
        p[i] = ip[i];
    }

    SORT_PAIR(p[1], p[2]);
    SORT_PAIR(p[4], p[5]);
    SORT_PAIR(p[7], p[8]);
    SORT_PAIR(p[0], p[1]);
    SORT_PAIR(p[3], p[4]);
    SORT_PAIR(p[6], p[7]);
    SORT_PAIR(p[1], p[2]);
    SORT_PAIR(p[4], p[5]);
    SORT_PAIR(p[7], p[8]);
    SORT_PAIR(p[0], p[3]);
    SORT_PAIR(p[5], p[8]);
    SORT_PAIR(p[4], p[7]);
    SORT_PAIR(p[3], p[6]);
    SORT_PAIR(p[1], p[4]);
    SORT_PAIR(p[2], p[5]);
    SORT_PAIR(p[4], p[7]);
    SORT_PAIR(p[4], p[2]);
    SORT_PAIR(p[6], p[4]);
    SORT_PAIR(p[4], p[2]);
    return (p[4]);
}
