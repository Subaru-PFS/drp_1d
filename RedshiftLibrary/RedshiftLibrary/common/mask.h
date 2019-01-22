#ifndef _REDSHIFT_COMMON_WEIGHTS_
#define _REDSHIFT_COMMON_WEIGHTS_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>

namespace NSEpic
{

class CSpectrum;

/**
 * \ingroup Redshift
 * Wavelength interval objects.
 */
class CMask
{

public:

    CMask();
    explicit CMask( UInt32 weightsCount );
    ~CMask();

    const Mask*     GetMasks() const;
    CMask&          operator &= ( const CMask& other );
    UInt32          GetMasksCount() const;
    Mask            operator[]( const UInt32 i ) const;
    Mask&           operator[]( const UInt32 i );
    Float64         CompouteOverlapRate( const CMask& other ) const;
    Float64         IntersectAndComputeOverlapRate( const CMask& other ) const;

    Bool            IntersectWith( const CMask& other );
    UInt32          GetMaskedSampleCount() const;
    UInt32          GetUnMaskedSampleCount() const;
    void            SetSize( UInt32 s );

private:

    TMaskList    m_Mask;

};

inline
Mask CMask::operator[]( const UInt32 i ) const
{
    return m_Mask[i];
}

inline
Mask& CMask::operator[]( const UInt32 i )
{
    return m_Mask[i];
}

inline
UInt32 CMask::GetMasksCount() const
{
    return m_Mask.size();
}

inline
const Mask* CMask::GetMasks() const
{
    return m_Mask.data();
}

inline
UInt32 CMask::GetMaskedSampleCount() const
{
    return m_Mask.size()-GetUnMaskedSampleCount();
}

inline
void CMask::SetSize( UInt32 s )
{
    m_Mask.resize( s );
}

inline
UInt32 CMask::GetUnMaskedSampleCount() const
{
    UInt32 n = 0;
    for( UInt32 i=0;i<m_Mask.size();i++ )
    {
        n += m_Mask[i];
    }
    return n;
}

}

#endif
