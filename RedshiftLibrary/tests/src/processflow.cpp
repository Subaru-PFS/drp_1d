#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/method/blindsolveresult.h>
#include <RedshiftLibrary/processflow/processflow.h>
#include <RedshiftLibrary/processflow/parameterstore.h>
#include <RedshiftLibrary/processflow/context.h>

#include <boost/test/unit_test.hpp>

#include <boost/property_tree/ptree.hpp>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(ProcessFlow)

BOOST_AUTO_TEST_CASE( ProcessShifted1 )
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

    std::string procID = "lbgabs_1K_2z3_20J22.5__EZ_fits-W-F_0";
    Bool retVal = ctx.Init( "../RedshiftLibrary/tests/src/data/ProcessFlowTestCase/lbgabs_1K_2z3_20J22.5__EZ_fits-W-F_0.fits", NULL, procID, "../RedshiftLibrary/tests/src/data/ProcessFlowTestCase/template_shifted1/", NULL, params, NULL );
    BOOST_CHECK( retVal == true );

    retVal = processFlow.Process( ctx );
    BOOST_CHECK( retVal == true );

    Float64 redshift;
    Float64 merit;
    std::string tplName;

    auto blindSolveResult = std::dynamic_pointer_cast<const CBlindSolveResult>( ctx.GetDataStore().GetGlobalResult( "blindsolve" ).lock() );
    blindSolveResult->GetBestFitResult( ctx.GetDataStore(), redshift, merit, tplName );

    BOOST_CHECK_CLOSE_FRACTION( 2.0295, redshift, 0.0001 );
}

BOOST_AUTO_TEST_CASE( ProcessShifted2 )
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

    std::string procID = "processing_id_unused";
    Bool retVal = ctx.Init( "../RedshiftLibrary/tests/src/data/ProcessFlowTestCase/lbgabs_1K_2z3_20J22.5__EZ_fits-W-F_206.fits", NULL, procID, "../RedshiftLibrary/tests/src/data/ProcessFlowTestCase/template_shifted2/", NULL, params, NULL );
    BOOST_CHECK( retVal == true );

    retVal = processFlow.Process( ctx );
    BOOST_CHECK( retVal == true );

    Float64 redshift=0.0;
    Float64 merit=0.0;
    std::string tplName="";

    auto blindSolveResult = std::dynamic_pointer_cast< const CBlindSolveResult>( ctx.GetDataStore().GetGlobalResult( "blindsolve" ).lock() );
    blindSolveResult->GetBestFitResult( ctx.GetDataStore(), redshift, merit, tplName );

    BOOST_CHECK_CLOSE_FRACTION( 2.7757, redshift, 0.0001 );
}

BOOST_AUTO_TEST_CASE( ProcessShiftedDecimated )
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


    std::string procID = "processing_id_unused";
    Bool retVal = ctx.Init( "../RedshiftLibrary/tests/src/data/ProcessFlowTestCase/lbgabs_1K_2z3_20J22.5__EZ_fits-W-F_0.fits", NULL, procID, "../RedshiftLibrary/tests/src/data/ProcessFlowTestCase/template_shifted_decimated/", NULL, params, NULL );
    BOOST_CHECK( retVal == true );

    retVal = processFlow.Process( ctx );
    BOOST_CHECK( retVal == true );

    Float64 redshift;
    Float64 merit;
    std::string tplName;

    auto blindSolveResult = std::dynamic_pointer_cast<const CBlindSolveResult>( ctx.GetDataStore().GetGlobalResult( "blindsolve" ).lock() );
    blindSolveResult->GetBestFitResult( ctx.GetDataStore(), redshift, merit, tplName );

    BOOST_CHECK_CLOSE_FRACTION( 2.02952, redshift, 0.00001 );
}

BOOST_AUTO_TEST_SUITE_END()
