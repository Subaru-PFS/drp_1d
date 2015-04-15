#ifndef _TEST_REDSHIFT_SMOOTH_
#define _TEST_REDSHIFT_SMOOTH_

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
  
namespace NSEpicTest
{

class CRedshiftSmoothTestCase : public CppUnit::TestCase
{
    CPPUNIT_TEST_SUITE(CRedshiftSmoothTestCase);
    CPPUNIT_TEST(Mean);
    CPPUNIT_TEST_SUITE_END();
 
public:

    void setUp();
    void tearDown();
 
private:

    void Mean();

};

CPPUNIT_TEST_SUITE_REGISTRATION( CRedshiftSmoothTestCase );

}

#endif
