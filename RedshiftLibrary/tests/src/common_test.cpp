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

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/median.h"
#include "RedshiftLibrary/common/mean.h"

#include <time.h>
#include <iostream>
#include <stdlib.h>
#include <boost/test/unit_test.hpp>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(Common)

BOOST_AUTO_TEST_CASE(Median3)
{
    CMedian<Float64> median;

    const Int32 size = 3;

    TFloat64List data(size);
    data[0] = 12;
    data[1] = 64;
    data[2] = 32;

    BOOST_CHECK( median.Find( data ) == 32 );
}

BOOST_AUTO_TEST_CASE(Median5)
{
    CMedian<Float64> median;

    const Int32 size = 5;

    TFloat64List data(size);
    data[0] = 12;
    data[1] = 64;
    data[2] = 32;
    data[3] = 74;
    data[4] = 69;

    BOOST_CHECK( median.Find( data ) == 64 );
}

BOOST_AUTO_TEST_CASE(Median7)
{
    CMedian<Float64> median;

    const Int32 size = 7;

    TFloat64List data(size);
    data[0] = 12;
    data[1] = 125;
    data[2] = 32;
    data[3] = 74;
    data[4] = 850;
    data[5] = 101;
    data[6] = 100;  

    BOOST_CHECK( median.Find( data ) == 100 );
}

BOOST_AUTO_TEST_CASE(Median9)
{
    CMedian<Float64> median;

    const Int32 size = 9;

    TFloat64List data(size);
    data[0] = 12;
    data[1] = 125;
    data[2] = 32;
    data[3] = 74;
    data[4] = 850;
    data[5] = 100;
    data[6] = 101;
    data[7] = 452;
    data[8] = 0;    

    BOOST_CHECK( median.Find( data ) == 100 );
}

BOOST_AUTO_TEST_CASE(MedianFast)
{
    CMedian<Float64> median;

    const Int32 size = ( MEDIAN_FAST_OR_BEERS_THRESHOLD * 2 ) + 1 ;

    TFloat64List data(size);
    srand( time(0) );
    for( Int32 i=0; i<MEDIAN_FAST_OR_BEERS_THRESHOLD; i++ )
    {
        data[i] = ( (Float64) rand() / (Float64) (RAND_MAX) ) * 100.0;
    }

    for( Int32 i=0; i<MEDIAN_FAST_OR_BEERS_THRESHOLD; i++ )
    {
        data[i+MEDIAN_FAST_OR_BEERS_THRESHOLD] = 200.0 + ( (Float64) rand() / (Float64) (RAND_MAX) ) * 100.0;
    }

    data[ size -1 ] = 150.0;
    
    BOOST_CHECK( median.Find( data ) == 150.0 );
}

BOOST_AUTO_TEST_CASE(MedianBeers)
{
    CMedian<Float64> median;

    const Int32 size = ( ( MEDIAN_FAST_OR_BEERS_THRESHOLD - 1 ) | 1 ) ;
    const Int32 halfSize = size / 2 ;

    BOOST_CHECK(  size <= MEDIAN_FAST_OR_BEERS_THRESHOLD );
    BOOST_CHECK( halfSize * 2 + 1 == size );

    TFloat64List data(size);
    srand( time(0) );
    for( Int32 i=0; i<halfSize; i++ )
    {
        data[i] = ( (Float64) rand() / (Float64) (RAND_MAX) ) * 100.0;
    }

    for( Int32 i=0; i<halfSize; i++ )
    {
        data[i+halfSize] = 200.0 + ( (Float64) rand() / (Float64) (RAND_MAX) ) * 100.0;
    }

    data[ size - 1 ] = 150.0;

    BOOST_CHECK( median.Find( data ) == 150.0 );
}

BOOST_AUTO_TEST_CASE(Mean)
{
    CMean<Float64> mean;

    const Int32 size = 128;

    TFloat64List data(size);

    for( Int32 i=0; i<size; i++ )
    {
        data[i] = 150.0;
    }

    BOOST_CHECK( mean.Find( data ) == 150.0 );
}



BOOST_AUTO_TEST_SUITE_END()
