#ifndef RAYDETECTION
#define RAYDETECTION


#include <epic/core/common/datatypes.h>
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

namespace NSEpicTest
{

class CRedshiftRayDetectionTestCase : public CppUnit::TestCase
{
    CPPUNIT_TEST_SUITE(CRedshiftRayDetectionTestCase);
    CPPUNIT_TEST(EzValidationTest);
    CPPUNIT_TEST_SUITE_END();

public:

    void setUp();
    void tearDown();

private:

    NSEpic::TFloat64List LoadDetectedRayPositions( const char* filePath );
    void EzValidationTest();

};

CPPUNIT_TEST_SUITE_REGISTRATION( CRedshiftRayDetectionTestCase );

}


#endif // RAYDETECTION

