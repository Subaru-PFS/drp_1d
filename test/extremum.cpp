#include "extremum.h"

#include <epic/core/common/datatypes.h>
#include <epic/redshift/extremum/extremum.h>

using namespace __NS__;

void CRedshiftExtremumTestCase::setUp()
{
}

void CRedshiftExtremumTestCase::tearDown()
{
}

void CRedshiftExtremumTestCase::Find()
{
    CExtremum peaks;

    Float64 x[] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
    Float64 y[] = { 1.0, 5.0, 0.0, 3.0, 9.0, 0.1, 0.2, 0.5, 8.0, 1.0, 4.0 };

    // Simple test
    {
        TPointList maxPoint;

        peaks.SetMaxPeakCount( 2 );
        peaks.SetRefreshCount( 2 );
        peaks.Find( x, y, 11, maxPoint);

        CPPUNIT_ASSERT( maxPoint.size() == 2 );

        CPPUNIT_ASSERT( maxPoint[0].X == x[4] );
        CPPUNIT_ASSERT( maxPoint[0].Y == y[4] );

        CPPUNIT_ASSERT( maxPoint[1].X == x[8] );
        CPPUNIT_ASSERT( maxPoint[1].Y == y[8] );
    }

    // Test with full range specified
    {
        TPointList maxPoint;

        peaks.SetMaxPeakCount( 2 );
        peaks.SetRefreshCount( 2 );
        peaks.SetXRange( TFloat64Range( 0.0, 1.0 ) );
        peaks.Find( x, y, 11, maxPoint);

        CPPUNIT_ASSERT( maxPoint.size() == 2 );

        CPPUNIT_ASSERT( maxPoint[0].X == x[4] );
        CPPUNIT_ASSERT( maxPoint[0].Y == y[4] );

        CPPUNIT_ASSERT( maxPoint[1].X == x[8] );
        CPPUNIT_ASSERT( maxPoint[1].Y == y[8] );
    }
}

