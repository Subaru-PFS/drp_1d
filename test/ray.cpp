#include "ray.h"

#include <epic/core/common/datatypes.h>
#include <epic/core/common/ref.h>
#include <epic/redshift/ray/catalog.h>
#include <epic/redshift/ray/matching.h>
#include <epic/redshift/operator/raymatchingresult.h>

#include <math.h>

using namespace NSEpic;

using namespace NSEpicTest;

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

    returnValue = catalog.Load( "../test/data/RayTestCase/raycatalog_OK1.txt" );
    CPPUNIT_ASSERT_MESSAGE( "Failled to load or parse raycatalog_OK1.txt", returnValue == true );

    returnValue = catalog.Load( "../test/data/RayTestCase/raycatalog_NOK1.txt" );
    CPPUNIT_ASSERT_MESSAGE( "Succeeded to parse invalid raycatalog_NOK1.txt", returnValue == false );

    returnValue = catalog.Load( "../test/data/RayTestCase/raycatalog_NOK1.txt" );
    CPPUNIT_ASSERT_MESSAGE( "Succeeded to parse invalid raycatalog_NOK1.txt", returnValue == false );

}

void CRedshiftRayTestCase::MatchingTest1()
// load a simple EL catalog and test the match with a redshifted version of itself
{
    CRayCatalog restFrameCatalog;
    Bool returnValue;

    returnValue = restFrameCatalog.Load( "../test/data/RayTestCase/raycatalog_testMatch1.txt" );
    CPPUNIT_ASSERT_MESSAGE( "Failed to load or parse raycatalog_testMatch1.txt", returnValue == true );

    CRayCatalog detectedCatalog;
    Float64 shiftLambda = 1.5;
    CRayCatalog::TRayVector cataloglist = restFrameCatalog.GetList();
    CRayCatalog::TRayVector ::iterator it;
    for( it = cataloglist.begin(); it != cataloglist.end(); ++it )
    {
        detectedCatalog.Add( CRay( (*it).GetName(), (*it).GetPosition()*shiftLambda, 0 ) );
    }

    CRayMatching rayMatching;
    TFloat64Range redshiftrange( 0.0, 5.0);
    CRef<CRayMatchingResult> result = rayMatching.Compute(detectedCatalog, restFrameCatalog, redshiftrange, 2, 0.002 );
    CPPUNIT_ASSERT_MESSAGE( "Failed to match ray catalogs for MatchingTest1.txt", result != NULL );

    Float64 res = result->GetMeanRedshiftSolutionByIndex(0);
    CPPUNIT_ASSERT_MESSAGE( "Failed to find redshift accurately for MatchingTest1", fabs(res-(shiftLambda-1)) < 0.0001 );

}

