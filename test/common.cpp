
#include <epic/core/common/datatypes.h>
#include <epic/redshift/common/median.h>
#include <epic/redshift/common/mean.h>

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

    Float64 data[size];
    data[0] = 12;
    data[1] = 64;
    data[2] = 32;

    BOOST_CHECK( median.Find( data, size ) == 32 );
}

BOOST_AUTO_TEST_CASE(Median5)
{
    CMedian<Float64> median;

    const Int32 size = 5;

    Float64 data[size];
    data[0] = 12;
    data[1] = 64;
    data[2] = 32;
    data[3] = 74;
    data[4] = 69;

    BOOST_CHECK( median.Find( data, size ) == 64 );
}

BOOST_AUTO_TEST_CASE(Median7)
{
    CMedian<Float64> median;

    const Int32 size = 7;

    Float64 data[size];
    data[0] = 12;
    data[1] = 125;
    data[2] = 32;
    data[3] = 74;
    data[4] = 850;
    data[5] = 101;
    data[6] = 100;  

    BOOST_CHECK( median.Find( data, size ) == 100 );
}

BOOST_AUTO_TEST_CASE(Median9)
{
    CMedian<Float64> median;

    const Int32 size = 9;

    Float64 data[size];
    data[0] = 12;
    data[1] = 125;
    data[2] = 32;
    data[3] = 74;
    data[4] = 850;
    data[5] = 100;
    data[6] = 101;
    data[7] = 452;
    data[8] = 0;    

    BOOST_CHECK( median.Find( data, size ) == 100 );
}

BOOST_AUTO_TEST_CASE(MedianBeers)
{
    CMedian<Float64> median;

    const Int32 size = ( MEDIAN_FAST_OR_BEERS_THRESHOLD * 2 ) + 1 ;

    Float64 data[size];
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

    BOOST_CHECK( median.Find( data, size ) == 150.0 );
}

BOOST_AUTO_TEST_CASE(MedianFast)
{
    CMedian<Float64> median;

    const Int32 size = ( ( MEDIAN_FAST_OR_BEERS_THRESHOLD - 1 ) | 1 ) ;
    const Int32 halfSize = size / 2 ;

    BOOST_CHECK(  size <= MEDIAN_FAST_OR_BEERS_THRESHOLD );
    BOOST_CHECK( halfSize * 2 + 1 == size );

    Float64 data[size];
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

    BOOST_CHECK( median.Find( data, size ) == 150.0 );
}

BOOST_AUTO_TEST_CASE(Mean)
{
    CMean<Float64> mean;

    const Int32 size = 128;

    Float64 data[size];

    for( Int32 i=0; i<size; i++ )
    {
        data[i] = 150.0;
    }

    BOOST_CHECK( mean.Find( data, size ) == 150.0 );
}



BOOST_AUTO_TEST_SUITE_END()
