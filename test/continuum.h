#ifndef _TEST_REDSHIFT_CONTINUUM_
#define _TEST_REDSHIFT_CONTINUUM_

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
  
namespace NSEpicTest
{

class CRedshiftContinuumTestCase : public CppUnit::TestCase
{
    CPPUNIT_TEST_SUITE(CRedshiftContinuumTestCase);
    CPPUNIT_TEST(Compute);
    CPPUNIT_TEST_SUITE_END();

public:

    void setUp();
    void tearDown();

private:

    void Compute();

};

CPPUNIT_TEST_SUITE_REGISTRATION( CRedshiftContinuumTestCase );
}

#endif
