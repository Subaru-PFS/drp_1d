#include <epic/redshift/common/redshifts.h>

#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/spectrum/axis.h>

#include <math.h>
#include <float.h>

using namespace NSEpic;

CRedshifts::CRedshifts()
{
}

CRedshifts::CRedshifts( Float64* redshifts, UInt32 count )
{
    Float64 min = DBL_MIN;
    Float64 max = DBL_MAX;

    m_Redshifts.resize( count );

    for( Int32 i=0;i<count;i++ )
    {
        m_Redshifts[i] = redshifts[i];

        if( redshifts[i] < min )
            min = redshifts[i];
        else if( redshifts[i] > max )
            max = redshifts[i];

    }

    m_Range.SetBegin( min );
    m_Range.SetEnd( max );
}

CRedshifts::CRedshifts( Redshift startRS, Redshift endRS, Float64 delta )
{
    SpreadOver( startRS, endRS, delta );
}

CRedshifts::CRedshifts( const TFloat64Range& range, Float64 delta )
{
    SpreadOver( range.GetBegin(), range.GetEnd(), delta );
}

CRedshifts::CRedshifts( UInt32 weightsCount ) :
    m_Redshifts( weightsCount )
{

}

CRedshifts::~CRedshifts()
{
}

Bool CRedshifts::SpreadOver( Float64 startRS, Float64 endRS, Float64 delta )
{
    TFloat64Range range( startRS, endRS );
    return SpreadOver( range, delta );
}

Bool CRedshifts::SpreadOver( const TFloat64Range& range, Float64 delta )
{
    if( range.GetIsEmpty() || delta == 0.0 )
        return false;

    if( range.GetLength() < delta )
        return false;

    Int32 count = range.GetLength() / delta;

    m_Range = range;
    m_Redshifts.resize( count+1 );
    for( UInt32 i=0; i<m_Redshifts.size() ;i++ )
    {
        m_Redshifts[i] = range.GetBegin() + delta * i;
    }

    return true;
}
