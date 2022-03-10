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
#include "RedshiftLibrary/common/mask.h"


#include <math.h>

using namespace NSEpic;


/**
 *
 */
CMask::CMask( Int32 weightsCount ) :
    m_Mask( weightsCount )
{

}

/**
 *
 */
CMask& CMask::operator &= ( const CMask& other )
{
    if( GetMasksCount() != other.GetMasksCount() )
        return *this;

    for( Int32 i = 0; i<GetMasksCount(); i++ )
    {
        m_Mask[i] = other[i] & m_Mask[i];
    }

    return *this;
}

/**
 *
 */
bool CMask::IntersectWith( const CMask& other )
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

    for( Int32 i=0; i<GetMasksCount(); i++)
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

    for( Int32 i=0; i<GetMasksCount(); i++)
    {
        selfRate+=(Int32) selfWeight[i];
        otherRate+=(Int32) (otherWeight[i]&selfWeight[i]);
    }

    if( selfRate == 0.0 )
        return 0;

    return (Float64)otherRate/(Float64)selfRate;
}
