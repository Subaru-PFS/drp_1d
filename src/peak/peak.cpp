#include <epic/redshift/peak/peak.h>

#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/common/median.h>

#include <math.h>

using namespace __NS__;

CPeak::CPeak( ) :
    m_StartIndex( -1 ),
    m_StopIndex( -1 )
{

}

CPeak::CPeak( const CPeak& other ) :
    m_StartIndex( other.m_StartIndex ),
    m_StopIndex( other.m_StopIndex )
{

}

CPeak::CPeak( Int32 startIndex, Int32 stopIndex ) :
    m_StartIndex( startIndex ),
    m_StopIndex( stopIndex )
{

}

Void  CPeak::Set( Int32 startIndex, Int32 stopIndex )
{
    m_StartIndex = startIndex;
    m_StopIndex = stopIndex;
}
