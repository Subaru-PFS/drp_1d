#include <epic/redshift/common/mask.h>

#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/spectrum/axis.h>

#include <math.h>

using namespace NSEpic;

CMask::CMask()
{
}

CMask::CMask( UInt32 weightsCount ) :
    m_Mask( weightsCount )
{

}

CMask::~CMask()
{
}

Bool CMask::IntersectWith( const CMask& other )
{
    if( GetMasksCount() != other.GetMasksCount() )
        return false;

    Mask* selfWeight = m_Mask.data();
    const Mask* otherWeight = other.GetMasks();

    for( Int32 j=0; j<GetMasksCount(); j++ )
    {
        selfWeight[j] = selfWeight[j] & otherWeight[j];
    }

    return true;
}

Float64 CMask::CompouteOverlapRate( const CMask& other ) const
{
    if( other.GetMasksCount() != GetMasksCount() )
        return -1.0;

    Float64 selfRate=0;
    Float64 otherRate=0;

    const Mask* selfWeight = GetMasks();
    const Mask* otherWeight = other.GetMasks();

    for( UInt32 i=0; i<GetMasksCount(); i++)
    {
        selfRate+=(Float64) selfWeight[i];
        otherRate+=(Float64) otherWeight[i];
    }

    if( selfRate == 0.0 )
        return 0;

    return otherRate/selfRate;
}
