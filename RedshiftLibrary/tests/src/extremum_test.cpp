#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/extremum/extremum.h>

#include <boost/test/unit_test.hpp>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(Extremum)

void print_point(TPointList points)
{
  for (SPoint p : points) {
    BOOST_TEST_MESSAGE("(" << p.X << ", " << p.Y << "), ");
  }
  BOOST_TEST_MESSAGE("=======");
}

void check_points(TPointList points, TPointList test_points) {
  BOOST_CHECK(points.size() == test_points.size());
  for (int i=0; i<points.size(); ++i) {
    BOOST_CHECK((points[i].X == test_points[i].X &&
		 points[i].Y == test_points[i].Y));
  }
}

BOOST_AUTO_TEST_CASE(Extremum1)
{
    CExtremum peaks1;
    CExtremum peaks2(true);
    CExtremum peaks_empty;
    TPointList maxPoint;
    TFloat64List x = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
    TFloat64List y = { 1.0, 5.0, 0.0, 3.0, 9.0, 0.1, 0.2, 0.5, 8.0, 1.0, 4.0 };
    TFloat64List y_2 = { 5.0, 0.0, 0.0, 3.0, 9.0, 0.1, 0.2, 0.5, 8.0, 1.0, 4.0 };
    TFloat64List x_empty = {};

    peaks2 = CExtremum(TFloat64Range(-10.0, 10.0), 5, true, 2);
    peaks2 = CExtremum(TFloat64Range(-10.0, 10.0), 5, false, 2);

    peaks1.SetMaxPeakCount(2);
    peaks1.SetRefreshCount(2);
    peaks1.SetXRange( TFloat64Range(-10.0, 10.0) );

    peaks1.Find( x, y, maxPoint);
    check_points(maxPoint, TPointList({ {0.4,9}, {0.8,8} }));

    peaks1.Find( x, y_2, maxPoint);
    check_points(maxPoint, TPointList({ {0.4,9}, {0.8,8} }));

    BOOST_CHECK( peaks1.Find( x_empty, y, maxPoint) == false);
    print_point(maxPoint);

    peaks_empty = CExtremum(TFloat64Range(1.0, 1.0), 5, false, 2);
    peaks_empty.Find( x, y, maxPoint);

    // peaks1.Find( x_empty, x_empty, maxPoint);
    // print_point(maxPoint);
//    BOOST_ERROR("FIXME: fails with empty list");

/*
    // Simple test
    {
        TPointList maxPoint;

        peaks.SetMaxPeakCount( 2 );
        peaks.SetRefreshCount( 2 );
        peaks.Find( x, y, 11, maxPoint);

        CPPUNIT_ASSERT( maxPoint.size() == 1 );

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
    */
}

BOOST_AUTO_TEST_SUITE_END()
