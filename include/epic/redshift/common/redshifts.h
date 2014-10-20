#ifndef _REDSHIFT_COMMON_REDSHIFTS_
#define _REDSHIFT_COMMON_REDSHIFTS_

#include <epic/redshift/common/datatypes.h>
#include <epic/core/common/range.h>

namespace NSEpic
{

class CSpectrum;

class CRedshifts
{

public:

    CRedshifts( );
    CRedshifts( UInt32 redshiftsCount );
    CRedshifts( const TFloat64Range& range, Float64 delta );
    CRedshifts( Redshift startRS, Redshift endRS, Float64 delta );
    CRedshifts( Float64* redshifts, UInt32 count );
    ~CRedshifts();

    const TFloat64Range&    GetRange() const;
    UInt32                  GetRedshiftsCount() const;
    const Redshift*         GetRedshifts() const;
    Redshift                operator[]( const UInt32 i ) const;
    Redshift&               operator[]( const UInt32 i );

    Bool                    SpreadOver( Float64 startRS, Float64 endRS, Float64 delta );
    Bool                    SpreadOver( const TFloat64Range& range, Float64 delta );

private:

    TRedshiftList       m_Redshifts;
    TFloat64Range       m_Range;

};

inline
Redshift CRedshifts::operator[]( const UInt32 i ) const
{
    return m_Redshifts[i];
}

inline
Redshift& CRedshifts::operator[]( const UInt32 i )
{
    return m_Redshifts[i];
}

inline
const TFloat64Range& CRedshifts::GetRange() const
{
    return m_Range;
}

inline
UInt32 CRedshifts::GetRedshiftsCount() const
{
    return m_Redshifts.size();
}

inline
const Redshift* CRedshifts::GetRedshifts() const
{
    return m_Redshifts.data();
}

}

#endif
