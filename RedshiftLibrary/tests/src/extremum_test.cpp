#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/extremum/extremum.h>

#include <boost/test/unit_test.hpp>
#include <iostream>  
#include <numeric>  
#include <cmath>
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
  std::cout << " returned peaks: "<<points.size()<<"\n";
  for(Int32 i=0; i<points.size(); i++){
    std::cout << " returned points: " << i << " " <<points[i].X << " " << points[i].Y << " " << test_points[i].X << " " << test_points[i].Y<<"\n";
  }
  
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

    Float64 radius = 0.005;
    peaks2 = CExtremum(TFloat64Range(-10.0, 10.0), 5, radius, true);
    peaks2 = CExtremum(TFloat64Range(-10.0, 10.0), 5, radius, false);

    peaks1.SetMaxPeakCount(2);
    peaks1.SetXRange( TFloat64Range(-10.0, 10.0) );

    peaks1.Find( x, y, maxPoint);
    check_points(maxPoint, TPointList({ {0.4,9}, {0.8,8} }));
    maxPoint.clear();
    peaks1.Find( x, y_2, maxPoint);
    check_points(maxPoint, TPointList({ {0.4,9}, {0.8,8} }));
  
    maxPoint.clear();
    BOOST_CHECK( peaks1.Find( x_empty, y, maxPoint) == false);
    print_point(maxPoint);

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

BOOST_AUTO_TEST_CASE(Extremum_PlankDetection)
{
    CExtremum peaks1;
    TPointList maxPoint;
    TFloat64List x = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
    TFloat64List y_plank = { 5.0, 0.0, 9.0, 9.0, 9.0, 0.1, 0.2, 0.5, 8.0, 1.0, 4.0 };
    TFloat64List y_plank_chair = { 5.0, 9.5, 9.0, 9.0, 9.0, 0.1, 0.2, 0.5, 8.0, 1.0, 4.0 };
    TFloat64List y_plank_stair = { 5.0, 0.0, 9.0, 9.0, 9.0, 9.5, 0.2, 0.5, 8.0, 1.0, 4.0 };
    TFloat64List y_plank_low = { 5.0, 0.0, 0.0, 0.0, 0.0, 9.5, 0.2, 0.5, 8.0, 1.0, 4.0 };

    Float64 radius = 0.005;
    peaks1.SetMaxPeakCount(2);
    peaks1.SetXRange( TFloat64Range(-10.0, 10.0) );

    peaks1.Find( x, y_plank, maxPoint);
    check_points(maxPoint, TPointList({ {0.3,9}, {0.8,8} }));

    maxPoint.clear();
    peaks1.Find( x, y_plank_chair, maxPoint);
    check_points(maxPoint, TPointList({ {0.1,9.5}, {0.8,8} }));

    maxPoint.clear();
    peaks1.Find( x, y_plank_stair, maxPoint);
    check_points(maxPoint, TPointList({ {0.5,9.5}, {0.8,8} }));

    maxPoint.clear();
    peaks1.Find( x, y_plank_low, maxPoint);
    check_points(maxPoint, TPointList({ {0.5,9.5}, {0.8,8} }));
}

BOOST_AUTO_TEST_CASE(Extremum_cut_isolated)
{
    CExtremum peaks1;
    TPointList maxPoint;
    TFloat64List x = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
    TFloat64List y = { 1.0, 5.0, 0.0, 3.0, 1.0, 0.1, 6.0, 0.5, 8.0, 1.0, 4.0 };

    Float64 radius = 0.005;
    peaks1.SetMaxPeakCount(5);
    peaks1.SetXRange( TFloat64Range(-10.0, 10.0) );
    peaks1.SetMeritCut(3);

    TInt32List sortedIndexes(x.size());
    iota(sortedIndexes.begin(), sortedIndexes.end(), 0);
    sort(sortedIndexes.begin(), sortedIndexes.end(),
       [&y](Float64 i1, Float64 i2) {return y[i1] > y[i2];});

    peaks1.SetSortedIndexes(sortedIndexes);

    Bool v = peaks1.Cut_Threshold(x, y, 2);
    for (Int32 i = 0; i < x.size(); i++) {       
        maxPoint.push_back(SPoint(x[i],  y[i]) );
    }
    check_points(maxPoint, TPointList({  {0.8,8}, {0.6,6}, {0.1,5} }));
}

BOOST_AUTO_TEST_CASE(Extremum_FilterOutNeighboringPeaks)
{
    
    TPointList maxPoint;
    TFloat64List x = { 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
    TFloat64List y = { 1.0, 5.0, 0.0, 3.0, 1.0, 0.1, 6.0, 0.5, 8.0, 1.0, 4.0 };

    Float64 radius = 0.2;
    CExtremum peaks1( TFloat64Range(-10.0, 10.0), 5, radius, false);

    //testing only the sliding window algo
    peaks1.FilterOutNeighboringPeaks(x, y, 2);
    for (Int32 i = 0; i < x.size(); i++) {
        
        maxPoint.push_back(SPoint(x[i],  y[i]) );
    }
    check_points(maxPoint, TPointList({{0.1,5.0}, {0.8,8}})); //here we lose a good peak (0.6, 6) since it is close to 0.8, within sliding windows 
}

BOOST_AUTO_TEST_CASE(Extremum_FilterOutNeighboringPeaks_consecutifClosePeaks)
{
    
    TPointList maxPoint;
    TFloat64List x = { 0.01, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1 };
    TFloat64List y = { 1.0, 5.0, 0.0, 3.0, 1.0, 0.1, 6.0, 0.5, 8.0, 1.0, 4.0 };

    Float64 radius = 0.05;
    CExtremum peaks1( TFloat64Range(-10.0, 10.0), 5, radius, false);

    //testing only the sliding window algo
    peaks1.FilterOutNeighboringPeaks(x, y, 2);
    for (Int32 i = 0; i < x.size(); i++) { 
        maxPoint.push_back(SPoint(x[i],  y[i]) );
    }
    check_points(maxPoint, TPointList({{0.08,8}})); //here we lose a good peak (0.6, 6) since it is close to 0.8, within sliding windows 
}


BOOST_AUTO_TEST_CASE(Extremum_FilterOutNeighboringPeaks_consecutifClosePeaks2)
{
    
    TPointList maxPoint;
    TFloat64List x = { 0.01, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1 };
    TFloat64List y = { 1.0, 3.0, 0.0, 4.0, 1.0, 0.1, 6.0, 0.5, 8.0, 1.0, 4.0 };

    Float64 radius = 0.005;
    CExtremum peaks1( TFloat64Range(-10.0, 10.0), 5, radius, false);

    //testing only the sliding window algo
    peaks1.FilterOutNeighboringPeaks(x, y, 2);
    for (Int32 i = 0; i < x.size(); i++) {
        maxPoint.push_back(SPoint(x[i],  y[i]) );
    }
    check_points(maxPoint, TPointList({{0.03,4.0}, {0.06,6.0}, {0.08,8}, {0.1,4.0}})); //here we lose a good peak (0.6, 6) since it is close to 0.8, within sliding windows 
}

BOOST_AUTO_TEST_CASE(Extremum_FilterOutNeighboringPeaks_increasingPDF)
{
    
    TPointList maxPoint;
    TFloat64List x = { 0.01, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1 };
    TFloat64List y = { 1.0, 3.0, 4.0, 5.0, 6.0, 7.1, 8.0, 9.5, 10.0, 11.0, 14.0 };

    Float64 radius = 0.005;
    CExtremum peaks1( TFloat64Range(-10.0, 10.0), 5, radius, false);

    //testing only the sliding window algo
    peaks1.FilterOutNeighboringPeaks(x, y, 2);
    for (Int32 i = 0; i < x.size(); i++) {
        maxPoint.push_back(SPoint(x[i],  y[i]) );
    }
    check_points(maxPoint, TPointList({{0.03,5.0}, {0.06,8.0}, {0.1,14.0}})); //here we lose a good peak (0.6, 6) since it is close to 0.8, within sliding windows 

    //bigger radius
    radius = 0.1;
    //CExtremum peaks1( TFloat64Range(-10.0, 10.0), 5, radius, false);
    maxPoint = {};
    //testing only the sliding window algo
    peaks1.FilterOutNeighboringPeaks(x, y, 2);
    for (Int32 i = 0; i < x.size(); i++) {
        maxPoint.push_back(SPoint(x[i],  y[i]) );
    }
    check_points(maxPoint, TPointList({{0.03,5.0}, {0.06,8.0}, {0.1,14.0}})); //here we lose a good peak (0.6, 6) since it is close to 0.8, within sliding windows 

}
//testing case where nb of peaks = keepmin
BOOST_AUTO_TEST_CASE(Extremum_FilterOutNeighboringPeaks2)
{
    
    TPointList maxPoint;
    TFloat64List x = { 0.01, 0.1};
    TFloat64List y = { 1.0, 5.0};

    Float64 radius = 0.2;
    CExtremum peaks1( TFloat64Range(-10.0, 10.0), 5, radius, false);

    //testing only the sliding window algo
    Int32  keepmin = 2;
    peaks1.FilterOutNeighboringPeaks(x, y, keepmin);
    for (Int32 i = 0; i < x.size(); i++) {
        
        maxPoint.push_back(SPoint(x[i],  y[i]) );
    }
    check_points(maxPoint, TPointList({{0.01,1.0}, {0.1,5.0}})); //here we lose a good peak (0.6, 6) since it is close to 0.8, within sliding windows 
}

BOOST_AUTO_TEST_CASE(Extremum_Truncate)
{
    
    TPointList maxPoint;
    TFloat64List x = { 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
    TFloat64List y = { 1.0, 5.0, 0.0, 3.0, 1.0, 0.1, 6.0, 0.5, 8.0, 1.0, 4.0 };

    Int32 maxCount = 2;
    CExtremum peaks1( false); //TFloat64Range(-10.0, 10.0), 5, radius, false);
    peaks1.Truncate( x, y, maxCount, maxPoint);
    
    check_points(maxPoint, TPointList({{0.8,8}, {0.6,6}}));
}
BOOST_AUTO_TEST_SUITE_END()
