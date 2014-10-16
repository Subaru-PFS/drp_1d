#ifndef _TEST_REDSHIFT_COMMON_
#define _TEST_REDSHIFT_COMMON_

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
  
class CRedshiftCommonTestCase : public CppUnit::TestCase
{
    CPPUNIT_TEST_SUITE(CRedshiftCommonTestCase);

    CPPUNIT_TEST(Median3);
    CPPUNIT_TEST(Median5);
    CPPUNIT_TEST(Median7);
    CPPUNIT_TEST(Median9);
    CPPUNIT_TEST(MedianBeers);
    CPPUNIT_TEST(MedianFast);

    CPPUNIT_TEST(Mean);

    CPPUNIT_TEST_SUITE_END();
 
public:

    void setUp();
    void tearDown();
 
private:

    void Median3();
    void Median5();
    void Median7();
    void Median9();
    void MedianBeers();
    void MedianFast();
    void Mean();

};

CPPUNIT_TEST_SUITE_REGISTRATION( CRedshiftCommonTestCase );
 
#endif
