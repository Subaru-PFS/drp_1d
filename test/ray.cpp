#include "ray.h"

#include <epic/core/common/datatypes.h>
#include <epic/redshift/ray/catalog.h>

using namespace NSEpic;

void CRedshiftRayTestCase::setUp()
{
}

void CRedshiftRayTestCase::tearDown()
{
}

void CRedshiftRayTestCase::LoadCatalog()
{
    CRayCatalog catalog;
    Bool returnValue;
    
    returnValue = catalog.Load( "../test/redshift/data/raycatalog_OK1.txt" );
    CPPUNIT_ASSERT_MESSAGE( "Failled to load or parse raycatalog_OK1.txt", returnValue == true );

    returnValue = catalog.Load( "../test/redshift/data/raycatalog_NOK1.txt" );
    CPPUNIT_ASSERT_MESSAGE( "Succeeded to parse invalid raycatalog_NOK1.txt", returnValue == false );

    returnValue = catalog.Load( "../test/redshift/data/raycatalog_NOK1.txt" );
    CPPUNIT_ASSERT_MESSAGE( "Succeeded to parse invalid raycatalog_NOK1.txt", returnValue == false );

}

