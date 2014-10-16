#include "processflow.h"

#include <epic/core/common/datatypes.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/processflow/processflow.h>
#include <epic/redshift/processflow/context.h>
#include <epic/core/common/ref.h>

using namespace __NS__;

void CRedshiftProcessFlowTestCase::setUp()
{
}

void CRedshiftProcessFlowTestCase::tearDown()
{
}

void CRedshiftProcessFlowTestCase::Process()
{
    CProcessFlowContext ctx;
    CProcessFlow processFlow;

    CProcessFlowContext::SParam params;
    params.lambdaRange = TFloat64Range( 3800.0, 12500.0 );
    params.templateCategoryList.clear();
    params.templateCategoryList.push_back( CTemplate::nCategory_Galaxy );
    params.templateCategoryList.push_back( CTemplate::nCategory_Qso );

    Bool retVal = ctx.Init( "../test/redshift/data/spectrum1_z_1.2299.fits", NULL, "../test/redshift/data/templatecatalog/", "../test/redshift/data/raycatalog.txt", params );
    CPPUNIT_ASSERT( retVal == true );

    retVal = processFlow.Process( ctx );
    CPPUNIT_ASSERT( retVal == true );

    Float64 redshift;
    Float64 merit;
    std::string tplName;

    ctx.GetBestCorrelationResult( redshift, merit, tplName );
    CPPUNIT_ASSERT( tplName == "elliptical.txt" );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.2299, redshift, 0.0001 );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.1955, merit, 0.0001 );


}
