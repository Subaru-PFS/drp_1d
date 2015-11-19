#include "processflow.h"

#include <epic/core/common/datatypes.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/method/blindsolveresult.h>
#include <epic/redshift/processflow/processflow.h>
#include <epic/redshift/processflow/parameterstore.h>
#include <epic/redshift/processflow/context.h>

#include <boost/property_tree/ptree.hpp>

using namespace NSEpic;

using namespace NSEpicTest;

void CRedshiftProcessFlowTestCase::setUp()
{

}

void CRedshiftProcessFlowTestCase::tearDown()
{

}

void CRedshiftProcessFlowTestCase::ProcessShifted1()
{
    CProcessFlowContext ctx;
    CProcessFlow processFlow;

    std::shared_ptr<CParameterStore> params = std::shared_ptr<CParameterStore>( new CParameterStore() );
    params->Set( "lambdaRange", TFloat64Range( 3800.0, 12500.0 ) );
    params->Set( "redshiftRange", TFloat64Range( 0.0, 5.0 ) );
    params->Set( "redshiftStep", 0.0001);
    params->Set( "smoothWidth", (Int64)0 );
    params->Set( "templateCategoryList", TStringList { "galaxy" } );
    params->Set( "method", "blindsolve");


    Bool retVal = ctx.Init( "../test/data/ProcessFlowTestCase/lbgabs_1K_2z3_20J22.5__EZ_fits-W-F_0.fits", NULL, "../test/data/ProcessFlowTestCase/template_shifted1/", NULL, params );
    CPPUNIT_ASSERT( retVal == true );

    retVal = processFlow.Process( ctx );
    CPPUNIT_ASSERT( retVal == true );

    Float64 redshift;
    Float64 merit;
    std::string tplName;

    auto blindSolveResult = std::dynamic_pointer_cast<const CBlindSolveResult>( ctx.GetDataStore().GetGlobalResult( "blindsolve" ).lock() );
    blindSolveResult->GetBestFitResult( ctx.GetDataStore(), redshift, merit, tplName );

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.0295, redshift, 0.0001 );
}

void CRedshiftProcessFlowTestCase::ProcessShifted2()
{
    CProcessFlowContext ctx;
    CProcessFlow processFlow;

    std::shared_ptr<CParameterStore> params = std::shared_ptr<CParameterStore>( new CParameterStore() );
    params->Set( "lambdaRange", TFloat64Range( 3800.0, 12500.0 ) );
    params->Set( "redshiftRange", TFloat64Range( 0.0, 5.0 ) );
    params->Set( "redshiftStep", 0.0001);
    params->Set( "smoothWidth", (Int64)0 );
    params->Set( "templateCategoryList",  TStringList { "galaxy" } );
    params->Set( "method", "blindsolve");

    Bool retVal = ctx.Init( "../test/data/ProcessFlowTestCase/lbgabs_1K_2z3_20J22.5__EZ_fits-W-F_206.fits", NULL, "../test/data/ProcessFlowTestCase/template_shifted2/", NULL, params );
    CPPUNIT_ASSERT( retVal == true );

    retVal = processFlow.Process( ctx );
    CPPUNIT_ASSERT( retVal == true );

    Float64 redshift=0.0;
    Float64 merit=0.0;
    std::string tplName="";

    auto blindSolveResult = std::dynamic_pointer_cast< const CBlindSolveResult>( ctx.GetDataStore().GetGlobalResult( "blindsolve" ).lock() );
    blindSolveResult->GetBestFitResult( ctx.GetDataStore(), redshift, merit, tplName );

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.7757, redshift, 0.0001 );

}

void CRedshiftProcessFlowTestCase::ProcessShiftedDecimated()
{
    CProcessFlowContext ctx;
    CProcessFlow processFlow;

    std::shared_ptr<CParameterStore> params = std::shared_ptr<CParameterStore>( new CParameterStore() );
    params->Set( "lambdaRange", TFloat64Range( 3800.0, 12500.0 ) );
    params->Set( "redshiftStep", 0.0001);
    params->Set( "smoothWidth", (Int64)0 );
    params->Set( "templateCategoryList", TStringList { "galaxy" } );
    params->Set( "method", "blindsolve");

    Bool retVal = ctx.Init( "../test/data/ProcessFlowTestCase/lbgabs_1K_2z3_20J22.5__EZ_fits-W-F_0.fits", NULL, "../test/data/ProcessFlowTestCase/template_shifted_decimated/", NULL, params );
    CPPUNIT_ASSERT( retVal == true );

    retVal = processFlow.Process( ctx );
    CPPUNIT_ASSERT( retVal == true );

    Float64 redshift;
    Float64 merit;
    std::string tplName;

    auto blindSolveResult = std::dynamic_pointer_cast<const CBlindSolveResult>( ctx.GetDataStore().GetGlobalResult( "blindsolve" ).lock() );
    blindSolveResult->GetBestFitResult( ctx.GetDataStore(), redshift, merit, tplName );

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.02952, redshift, 0.00001 );

}
