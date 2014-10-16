#ifndef _TEST_REDSHIFT_OPERATOR_
#define _TEST_REDSHIFT_OPERATOR_

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

class CRedshiftOperatorTestCase : public CppUnit::TestCase
{
    CPPUNIT_TEST_SUITE(CRedshiftOperatorTestCase);
    CPPUNIT_TEST(Correlation1);
    CPPUNIT_TEST(Correlation2);
    CPPUNIT_TEST(Correlation3);
    CPPUNIT_TEST(Correlation4);
    CPPUNIT_TEST_SUITE_END();

public:

    void setUp();
    void tearDown();

private:

    void Correlation1();
    void Correlation2();
    void Correlation3();
    void Correlation4();
};

CPPUNIT_TEST_SUITE_REGISTRATION( CRedshiftOperatorTestCase );

#endif
