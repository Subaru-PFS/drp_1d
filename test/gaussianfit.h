#ifndef _TEST_REDSHIFT_GAUSSIANFIT_
#define _TEST_REDSHIFT_GAUSSIANFIT_

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
  
namespace NSEpicTest
{

class CRedshiftGaussianFitTestCase : public CppUnit::TestCase
{
    CPPUNIT_TEST_SUITE(CRedshiftGaussianFitTestCase);
    CPPUNIT_TEST(TestFit1);
    CPPUNIT_TEST(TestFit2);
    CPPUNIT_TEST_SUITE_END();
 
public:

    void setUp();
    void tearDown();
 
private:

    void TestFit1();
    void TestFit2();

};

CPPUNIT_TEST_SUITE_REGISTRATION( CRedshiftGaussianFitTestCase );

}

#endif
