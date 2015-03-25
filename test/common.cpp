#include "common.h"

#include <epic/core/common/datatypes.h>
#include <time.h>
#include <stdlib.h>

using namespace NSEpicTest;
using namespace NSEpic;

void CRedshiftCommonTestCase::setUp()
{
}

void CRedshiftCommonTestCase::tearDown()
{
}

#include <epic/redshift/common/median.h>

void CRedshiftCommonTestCase::Median3()
{
    CMedian<Float64> median;

    const Int32 size = 3;

    Float64 data[size];
    data[0] = 12;
    data[1] = 64;
    data[2] = 32;

    CPPUNIT_ASSERT( median.Find( data, size ) == 32 );
}

void CRedshiftCommonTestCase::Median5()
{
    CMedian<Float64> median;

    const Int32 size = 5;

    Float64 data[size];
    data[0] = 12;
    data[1] = 64;
    data[2] = 32;
    data[3] = 74;
    data[4] = 69;

    CPPUNIT_ASSERT( median.Find( data, size ) == 64 );
}

void CRedshiftCommonTestCase::Median7()
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

    CPPUNIT_ASSERT( median.Find( data, size ) == 100 );
}

void CRedshiftCommonTestCase::Median9()
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

    CPPUNIT_ASSERT( median.Find( data, size ) == 100 );
}

void CRedshiftCommonTestCase::MedianBeers()
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

    CPPUNIT_ASSERT( median.Find( data, size ) == 150.0 );
}

void CRedshiftCommonTestCase::MedianFast()
{
    CMedian<Float64> median;

    const Int32 size = ( ( MEDIAN_FAST_OR_BEERS_THRESHOLD - 1 ) | 1 ) ;
    const Int32 halfSize = size / 2 ;

    CPPUNIT_ASSERT(  size <= MEDIAN_FAST_OR_BEERS_THRESHOLD );
    CPPUNIT_ASSERT( halfSize * 2 + 1 == size );

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

    CPPUNIT_ASSERT( median.Find( data, size ) == 150.0 );
}

#include <epic/redshift/common/mean.h>

void CRedshiftCommonTestCase::Mean()
{
    CMean<Float64> mean;

    const Int32 size = 128;

    Float64 data[size];

    for( Int32 i=0; i<size; i++ )
    {
        data[i] = 150.0;
    }

    CPPUNIT_ASSERT( mean.Find( data, size ) == 150.0 );
}
