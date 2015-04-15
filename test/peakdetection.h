#ifndef _TEST_REDSHIFT_PEAKDETECTION_
#define _TEST_REDSHIFT_PEAKDETECTION_

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
  

namespace NSEpicTest
{

class CRedshiftPeakDetectionTestCase : public CppUnit::TestCase
{
    CPPUNIT_TEST_SUITE(CRedshiftPeakDetectionTestCase);
    CPPUNIT_TEST(Compute);
    CPPUNIT_TEST_SUITE_END();
 
public:

    void setUp();
    void tearDown();
 
private:

    void Compute();

};

CPPUNIT_TEST_SUITE_REGISTRATION( CRedshiftPeakDetectionTestCase );
}
#endif
