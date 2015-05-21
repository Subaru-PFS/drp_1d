#ifndef _TEST_REDSHIFT_OPERATOR_
#define _TEST_REDSHIFT_OPERATOR_

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

namespace NSEpicTest
{

class CRedshiftOperatorTestCase : public CppUnit::TestCase
{
    CPPUNIT_TEST_SUITE(CRedshiftOperatorTestCase);
    //CPPUNIT_TEST(CorrelationAtZEqualZero);
    //CPPUNIT_TEST(CorrelationAtGivenZ);
    CPPUNIT_TEST(CorrelationMatchWithEZ);
    CPPUNIT_TEST_SUITE_END();

public:

    void setUp();
    void tearDown();

private:

    void CorrelationMatchWithEZ();
    void CorrelationAtZEqualZero();
    void CorrelationAtGivenZ();

    void CorrelationMatchWithEZ( const char* spectraPath, const char* noisePath, const char* tplPath, const char* resultPath );

};

CPPUNIT_TEST_SUITE_REGISTRATION( CRedshiftOperatorTestCase );

}

#endif
