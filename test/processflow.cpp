#include "processflow.h"

#include <epic/core/common/datatypes.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/operator/blindsolveresult.h>
#include <epic/redshift/processflow/processflow.h>
#include <epic/redshift/processflow/context.h>
#include <epic/core/common/ref.h>

using namespace NSEpic;

using namespace NSEpicTest;

void CRedshiftProcessFlowTestCase::setUp()
{
}

void CRedshiftProcessFlowTestCase::tearDown()
{
}

void CRedshiftProcessFlowTestCase::ProcessShifted1 ()
{
    CProcessFlowContext ctx;
    CProcessFlow processFlow;

    CProcessFlowContext::SParam params;
    params.lambdaRange = TFloat64Range( 3800.0, 12500.0 );
    params.redshiftStep = 0.00001;
    params.smoothWidth = 0;
    params.templateCategoryList.clear();
    params.templateCategoryList.push_back( CTemplate::nCategory_Galaxy );

    Bool retVal = ctx.Init( "../test/data/ProcessFlowTestCase/lbgabs_1K_2z3_20J22.5__EZ_fits-W-F_0.fits", NULL, "../test/data/ProcessFlowTestCase/template_shifted1/", NULL, params );
    CPPUNIT_ASSERT( retVal == true );

    retVal = processFlow.Process( ctx );
    CPPUNIT_ASSERT( retVal == true );

    Float64 redshift;
    Float64 merit;
    std::string tplName;

    const CBlindSolveResult* blindSolveResult = (CBlindSolveResult*)ctx.GetGlobalResult( "blindsolve" );
    blindSolveResult->GetBestCorrelationResult( ctx, redshift, merit, tplName );

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.02952, redshift, 0.00001 );

}

void CRedshiftProcessFlowTestCase::ProcessShifted2()
{
    CProcessFlowContext ctx;
    CProcessFlow processFlow;

    CProcessFlowContext::SParam params;
    params.lambdaRange = TFloat64Range( 3800.0, 12500.0 );
    params.redshiftStep = 0.00001;
    params.smoothWidth = 0;
    params.templateCategoryList.clear();
    params.templateCategoryList.push_back( CTemplate::nCategory_Galaxy );

    Bool retVal = ctx.Init( "../test/data/ProcessFlowTestCase/lbgabs_1K_2z3_20J22.5__EZ_fits-W-F_206.fits", NULL, "../test/data/ProcessFlowTestCase/template_shifted2/", NULL, params );
    CPPUNIT_ASSERT( retVal == true );

    retVal = processFlow.Process( ctx );
    CPPUNIT_ASSERT( retVal == true );

    Float64 redshift;
    Float64 merit;
    std::string tplName;

    const CBlindSolveResult* blindSolveResult = (CBlindSolveResult*)ctx.GetGlobalResult( "blindsolve" );
    blindSolveResult->GetBestCorrelationResult( ctx, redshift, merit, tplName );

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.77570, redshift, 0.00001 );

}

void CRedshiftProcessFlowTestCase::ProcessShiftedDecimated()
{
    CProcessFlowContext ctx;
    CProcessFlow processFlow;

    CProcessFlowContext::SParam params;
    params.lambdaRange = TFloat64Range( 3800.0, 12500.0 );
    params.redshiftStep = 0.00001;
    params.smoothWidth = 0;
    params.templateCategoryList.clear();
    params.templateCategoryList.push_back( CTemplate::nCategory_Galaxy );

    Bool retVal = ctx.Init( "../test/data/ProcessFlowTestCase/lbgabs_1K_2z3_20J22.5__EZ_fits-W-F_0.fits", NULL, "../test/data/ProcessFlowTestCase/template_shifted_decimated/", NULL, params );
    CPPUNIT_ASSERT( retVal == true );

    retVal = processFlow.Process( ctx );
    CPPUNIT_ASSERT( retVal == true );

    Float64 redshift;
    Float64 merit;
    std::string tplName;

    const CBlindSolveResult* blindSolveResult = (CBlindSolveResult*)ctx.GetGlobalResult( "blindsolve" );
    blindSolveResult->GetBestCorrelationResult( ctx, redshift, merit, tplName );

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.02952, redshift, 0.00001 );

}
