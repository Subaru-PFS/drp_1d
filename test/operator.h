#ifndef _TEST_REDSHIFT_OPERATOR_
#define _TEST_REDSHIFT_OPERATOR_

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

namespace NSEpicTest
{

class CRedshiftOperatorTestCase : public CppUnit::TestCase
{
    CPPUNIT_TEST_SUITE(CRedshiftOperatorTestCase);
    CPPUNIT_TEST(CorrelationAtZEqualZero);
    CPPUNIT_TEST(CorrelationAtGivenZ);
    CPPUNIT_TEST_SUITE_END();

public:

    void setUp();
    void tearDown();

private:

    void CorrelationAtZEqualZero();
    void CorrelationAtGivenZ();
};

CPPUNIT_TEST_SUITE_REGISTRATION( CRedshiftOperatorTestCase );

}

#endif
