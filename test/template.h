#ifndef _TEST_REDSHIFT_MERIT_
#define _TEST_REDSHIFT_MERIT_

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

namespace NSEpicTest
{

class CRedshiftTemplateTestCase : public CppUnit::TestCase
{
    CPPUNIT_TEST_SUITE(CRedshiftTemplateTestCase);
    CPPUNIT_TEST( LoadCatalog );
    CPPUNIT_TEST_SUITE_END();

public:

    void setUp();
    void tearDown();

private:

    void LoadCatalog();

};

CPPUNIT_TEST_SUITE_REGISTRATION( CRedshiftTemplateTestCase );

}

#endif
