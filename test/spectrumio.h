#ifndef _TEST_REDSHIFT_SPECTRUM_
#define _TEST_REDSHIFT_SPECTRUM_

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
  
namespace NSEpicTest
{

class CCRedshiftSpectrumioTestCase : public CppUnit::TestCase
{
    CPPUNIT_TEST_SUITE(CCRedshiftSpectrumioTestCase);
    CPPUNIT_TEST( VVDSReadValidFile );
    CPPUNIT_TEST( VVDSReadInvalidFile );
    CPPUNIT_TEST_SUITE_END();
 
public:

    void setUp();
    void tearDown();
 
private:

    void VVDSReadValidFile();
    void VVDSReadInvalidFile();

};

CPPUNIT_TEST_SUITE_REGISTRATION( CCRedshiftSpectrumioTestCase );

}

#endif
