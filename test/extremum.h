#ifndef _TEST_REDSHIFT_PEAK_
#define _TEST_REDSHIFT_PEAK_

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

class CRedshiftExtremumTestCase : public CppUnit::TestCase
{
    CPPUNIT_TEST_SUITE(CRedshiftExtremumTestCase);
    CPPUNIT_TEST(Find);
    CPPUNIT_TEST_SUITE_END();

public:

    void setUp();
    void tearDown();

private:

    void Find();

};

CPPUNIT_TEST_SUITE_REGISTRATION( CRedshiftExtremumTestCase );

#endif
