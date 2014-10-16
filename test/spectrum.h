#ifndef _TEST_REDSHIFT_SPECTRUM_
#define _TEST_REDSHIFT_SPECTRUM_

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
  
class CRedshiftSpectrumTestCase : public CppUnit::TestCase
{
    CPPUNIT_TEST_SUITE(CRedshiftSpectrumTestCase);
    CPPUNIT_TEST(Load);
    CPPUNIT_TEST_SUITE_END();
 
public:

    void setUp();
    void tearDown();
 
private:

    void Load();

};

CPPUNIT_TEST_SUITE_REGISTRATION( CRedshiftSpectrumTestCase );
 
#endif
