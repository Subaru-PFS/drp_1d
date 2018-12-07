#include <RedshiftLibrary/common/mask.h>

#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/axis.h>

#include <math.h>

using namespace NSEpic;

/**
 *
 */
CMask::CMask()
{

}

/**
 *
 */
CMask::CMask( UInt32 weightsCount ) :
    m_Mask( weightsCount )
{

}

/**
 *
 */
CMask::~CMask()
{

}

/**
 *
 */
CMask& CMask::operator &= ( const CMask& other )
{
    if( GetMasksCount() != other.GetMasksCount() )
        return *this;

    for( UInt32 i = 0; i<GetMasksCount(); i++ )
    {
        m_Mask[i] = other[i] & m_Mask[i];
    }

    return *this;
}

/**
 *
 */
Bool CMask::IntersectWith( const CMask& other )
{
    if( GetMasksCount() != other.GetMasksCount() )
        return false;

    Mask* selfWeight = m_Mask.data();
    const Mask* otherWeight = other.GetMasks();

    for( UInt32 j=0; j<GetMasksCount(); j++ )
    {
        selfWeight[j] = selfWeight[j] & otherWeight[j];
    }

    return true;
}

/**
 *
 */
Float64 CMask::CompouteOverlapRate( const CMask& other ) const
{
    if( other.GetMasksCount() != GetMasksCount() )
        return -1.0;

    Float64 selfRate=0;
    Float64 otherRate=0;

    /* method1
    selfRate = GetUnMaskedSampleCount();
    otherRate = other.GetUnMaskedSampleCount();
    //*/

    //* method2
    const Mask* selfWeight = GetMasks();
    const Mask* otherWeight = other.GetMasks();

    for( UInt32 i=0; i<GetMasksCount(); i++)
    {
        //selfRate+=(Float64) selfWeight[i];
        //otherRate+=(Float64) otherWeight[i];
        selfRate+=(Int32) selfWeight[i];
        otherRate+=(Int32) otherWeight[i];
    }
    //*/

    if( selfRate == 0.0 )
        return 0;

    return (Float64)otherRate/(Float64)selfRate;
}

/**
 *
 */
Float64 CMask::IntersectAndComputeOverlapRate( const CMask& other ) const
{
    if( other.GetMasksCount() != GetMasksCount() )
        return -1.0;

    Int32 selfRate=0;
    Int32 otherRate=0;

    const Mask* selfWeight = GetMasks();
    const Mask* otherWeight = other.GetMasks();

    for( UInt32 i=0; i<GetMasksCount(); i++)
    {
        selfRate+=(Int32) selfWeight[i];
        otherRate+=(Int32) (otherWeight[i]&selfWeight[i]);
    }

    if( selfRate == 0.0 )
        return 0;

    return (Float64)otherRate/(Float64)selfRate;
}
