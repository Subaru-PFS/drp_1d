#ifndef _TEST_REDSHIFT_RAY_
#define _TEST_REDSHIFT_RAY_

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
  
namespace NSEpicTest
{

class CRedshiftRayTestCase : public CppUnit::TestCase
{
    CPPUNIT_TEST_SUITE(CRedshiftRayTestCase);
    CPPUNIT_TEST(LoadCatalog);
    CPPUNIT_TEST(MatchingTest1);
    CPPUNIT_TEST_SUITE_END();
 
public:

    void setUp();
    void tearDown();
 
private:

    void LoadCatalog();
    void MatchingTest1();

};

CPPUNIT_TEST_SUITE_REGISTRATION( CRedshiftRayTestCase );

}

#endif
