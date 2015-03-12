#ifndef _TEST_REDSHIFT_PROCESSFLOW_
#define _TEST_REDSHIFT_PROCESSFLOW_

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

class CRedshiftProcessFlowTestCase : public CppUnit::TestCase
{
    CPPUNIT_TEST_SUITE(CRedshiftProcessFlowTestCase);
    CPPUNIT_TEST(ProcessShifted);
    CPPUNIT_TEST(ProcessShiftedDecimated);
    CPPUNIT_TEST_SUITE_END();

public:

    void setUp();
    void tearDown();

private:

    void ProcessShifted();
    void ProcessShiftedDecimated();

};

CPPUNIT_TEST_SUITE_REGISTRATION( CRedshiftProcessFlowTestCase );

#endif
