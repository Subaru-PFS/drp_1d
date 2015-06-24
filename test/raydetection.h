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
    CPPUNIT_TEST(SyntheticValidationTest);
    //CPPUNIT_TEST(EzValidationTest); //deactivated, 20150624, due to irregular sampling compatibility implementation (differs from EZ)
    CPPUNIT_TEST_SUITE_END();

public:

    void setUp();
    void tearDown();

private:

    NSEpic::TFloat64List LoadDetectedRayPositions( const char* filePath );
    void EzValidationTest();
    void SyntheticValidationTest();

};

CPPUNIT_TEST_SUITE_REGISTRATION( CRedshiftRayDetectionTestCase );

}


#endif // RAYDETECTION

