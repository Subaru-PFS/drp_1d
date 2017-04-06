#include <epic/core/common/datatypes.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/method/blindsolveresult.h>
#include <epic/redshift/processflow/processflow.h>
#include <epic/redshift/processflow/parameterstore.h>
#include <epic/redshift/processflow/context.h>

#include <epic/redshift/method/linemodelsolve.h>
#include <epic/redshift/method/linemodelsolveresult.h>
#include <epic/redshift/linemodel/modelfittingresult.h>

#include <boost/test/unit_test.hpp>

#include <boost/property_tree/ptree.hpp>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(LinemodelRules)

Float64 getLinemodelDoubletRatio(std::string spc, std::string noise, bool enableRatioRule){

    CProcessFlowContext ctx;
    CProcessFlow processFlow;


    TFloat64Range redshiftRange = TFloat64Range( 0.0, 0.0 );
    TFloat64Range spcLambdaRange = TFloat64Range( 3800.0, 12000.0 );

    std::shared_ptr<CParameterStore> params = std::shared_ptr<CParameterStore>( new CParameterStore() );
    params->Set( "lambdaRange", spcLambdaRange);
    params->Set( "redshiftRange",  redshiftRange);
    params->Set( "redshiftStep", 0.1);
    params->Set( "smoothWidth", (Int64)0 );
    params->Set( "templateCategoryList", TStringList { "galaxy" } );
    params->Set( "method", "linemodel");

    std::string procID = "processing_id_unused";
    Bool retVal = ctx.Init( spc.c_str(), noise.c_str(), procID, NULL, "../test/data/LinemodelRulesTestCase/raycatalog_test_elratiorules.txt",params );
    BOOST_CHECK( retVal == true );

    //these tplcatalog related variables are unused here.
    CTemplateCatalog tplCatalog;
    Bool rValue = tplCatalog.Load( "../test/data/templatecatalog/" );
    TStringList tplCategories;


    if(enableRatioRule){
        ctx.GetDataStore().SetScopedParam("linemodelsolve.linemodel.rules", "oiiratio");
    }else{
        ctx.GetDataStore().SetScopedParam("linemodelsolve.linemodel.rules", "no");
    }

    retVal = processFlow.Process( ctx );
    BOOST_CHECK( retVal == true );

    // Create redshift initial list by spanning redshift acdross the given range, with the given delta
    Float64 redshiftStep = 0.001;
    TFloat64List redshifts = redshiftRange.SpreadOver( redshiftStep );

    CLineModelSolve Solve;
    std::shared_ptr<const CLineModelSolveResult> solveResult = Solve.Compute(ctx.GetDataStore(),
                                                                             ctx.GetSpectrum(),
                                                                             ctx.GetSpectrumWithoutContinuum(),
                                                                             tplCatalog,
                                                                             tplCategories,
                                                                             ctx.GetRayCatalog(),
                                                                             spcLambdaRange,
                                                                             redshifts);


    std::string scope = "linemodelsolve.linemodel_fit_extrema_0";
    auto results = ctx.GetDataStore().GetGlobalResult(scope.c_str());

    Float64 a1 = -1.0;
    Float64 a2 = -1.0;
    if(!results.expired()){
        std::shared_ptr<const CModelFittingResult> result = std::dynamic_pointer_cast<const CModelFittingResult>( results.lock() );

        a1 = result->GetLineModelSolution().Amplitudes[0];
        a2 = result->GetLineModelSolution().Amplitudes[1];
    }


    Float64 ratio = a2/a1;
    return ratio;
}

BOOST_AUTO_TEST_CASE( OIIRatioRange1 )
{
    std::string spc, noise;
    Float64 ratio;

    Float64 OIIRatioRangeLimit = 2.5; //this threshold has to be set according the regulament.cpp value

    //TEST 1 : initially A1 = 3*A2, oiiratiorange rule should modify both amplitudes to get a ratio = OIIRatioRangeLimit
    spc = "../test/data/LinemodelRulesTestCase/simu_rules_ratiorange_1.fits";
    noise = "../test/data/LinemodelRulesTestCase/simu_rules_ratiorange_1_noise.fits";
    ratio = getLinemodelDoubletRatio(spc, noise, 0); //rule disabled, get ratio
    BOOST_CHECK_MESSAGE( ratio > OIIRatioRangeLimit, "RatioRange: 1st spectrum control test failed = leads to a ratio not > OIIRatioRangeLimit" );
    ratio = getLinemodelDoubletRatio(spc, noise, 1); //rule enabled, get ratio
    BOOST_CHECK_CLOSE_FRACTION( OIIRatioRangeLimit, ratio, 0.1);

    //TEST 2 : same as Test 1, but initially A2 = 3*A1
    spc = "../test/data/LinemodelRulesTestCase/simu_rules_ratiorange_2.fits";
    noise = "../test/data/LinemodelRulesTestCase/simu_rules_ratiorange_2_noise.fits";
    ratio = getLinemodelDoubletRatio(spc, noise, 1);
    BOOST_CHECK_CLOSE_FRACTION( 1./OIIRatioRangeLimit, ratio, 0.1 );

    //TEST 3 : A1 = 1.5*A2, so that the ratio should be unchanged
    spc = "../test/data/LinemodelRulesTestCase/simu_rules_ratiorange_3.fits";
    noise = "../test/data/LinemodelRulesTestCase/simu_rules_ratiorange_3_noise.fits";
    Float64 ratioIN = getLinemodelDoubletRatio(spc, noise, 0);
    Float64 ratioMOD = getLinemodelDoubletRatio(spc, noise, 1);
    BOOST_CHECK_MESSAGE( ratioIN == ratioMOD, "RatioRange: 3rd spectrum test failed = rule should not have been applied in this case" );
}


std::vector<Float64> getLinemodelFittedAmplitudes(std::string spc, std::string noise, std::string ctlgPath, bool enableSuperStrongRule){

    CProcessFlowContext ctx;
    CProcessFlow processFlow;


    TFloat64Range redshiftRange = TFloat64Range( 0.0, 0.0 );
    TFloat64Range spcLambdaRange = TFloat64Range( 2000.0, 12000.0 );

    std::shared_ptr<CParameterStore> params = std::shared_ptr<CParameterStore>( new CParameterStore() );
    params->Set( "lambdaRange", spcLambdaRange);
    params->Set( "redshiftRange",  redshiftRange);
    params->Set( "redshiftStep", 0.01);
    params->Set( "smoothWidth", (Int64)0 );
    params->Set( "templateCategoryList", TStringList { "galaxy" } );
    params->Set( "method", "linemodel");


    std::string procID = "processing_id_unused";
    Bool retVal = ctx.Init( spc.c_str(), noise.c_str(), procID, NULL, ctlgPath.c_str(),params );
    BOOST_CHECK( retVal == true );


    if(enableSuperStrongRule){
        ctx.GetDataStore().SetScopedParam("linemodelsolve.linemodel.rules", "superstrong");
    }else{
        ctx.GetDataStore().SetScopedParam("linemodelsolve.linemodel.rules", "no");
    }

    retVal = processFlow.Process( ctx );
    BOOST_CHECK( retVal == true );


    //these tplcatalog related variables are unused here.
    CTemplateCatalog tplCatalog;
    Bool rValue = tplCatalog.Load( "../test/data/templatecatalog/" );
    TStringList tplCategories;

    // Create redshift initial list by spanning redshift acdross the given range, with the given delta
    Float64 redshiftStep = 0.01;
    TFloat64List redshifts = redshiftRange.SpreadOver( redshiftStep );

    CLineModelSolve Solve;
    std::shared_ptr<const CLineModelSolveResult> solveResult = Solve.Compute(ctx.GetDataStore(),
                                                                             ctx.GetSpectrum(),
                                                                             ctx.GetSpectrumWithoutContinuum(),
                                                                             tplCatalog,
                                                                             tplCategories,
                                                                             ctx.GetRayCatalog(),
                                                                             spcLambdaRange,
                                                                             redshifts);


    std::string scope = "linemodelsolve.linemodel_fit_extrema_0";
    auto results = ctx.GetDataStore().GetGlobalResult(scope.c_str());

    std::vector<Float64> amps;
    if(!results.expired()){
        std::shared_ptr<const CModelFittingResult> result = std::dynamic_pointer_cast<const CModelFittingResult>( results.lock() );

        amps = result->GetLineModelSolution().Amplitudes;
    }
    return amps;
}


BOOST_AUTO_TEST_CASE( OIIIMultilineSuperstrongRule )
{
    std::string spc, noise, ctlg;
    Float64 oiii_ratio;
    std::vector<Float64> amplis;

    //TEST 1 : initially oii = 2.0, oiiia = 30.0 and oiiib = 10.0
    spc = "../test/data/LinemodelRulesTestCase/simu_rules_multiline_superstong_1.fits";
    noise = "../test/data/LinemodelRulesTestCase/simu_rules_multiline_superstong_1_noise.fits";
    ctlg = "../test/data/LinemodelRulesTestCase/raycatalog_test_elmultilinesuperstong.txt";
    amplis = getLinemodelFittedAmplitudes(spc, noise, ctlg, 0); //rule disabled, get amplitudes
    oiii_ratio = amplis[0]/amplis[1];
    BOOST_CHECK_MESSAGE( oiii_ratio < 3.01 && oiii_ratio > 2.99, "Multiline-Superstrong: 1st test (no rule) failed = leads to a oiii nominal ratio not conserved" );

    amplis = getLinemodelFittedAmplitudes(spc, noise, ctlg, 1); //rule enabled, get amplitudes
    oiii_ratio = amplis[0]/amplis[1];
    BOOST_CHECK_MESSAGE( oiii_ratio < 3.01 && oiii_ratio > 2.99, "Multiline-Superstrong: 2nd test (superstrong rule) failed = leads to a oiii nominal ratio not conserved" );
    //ratio = getLinemodelDoubletRatio(spc, noise, 1); //rule enabled, get ratio
    //BOOST_CHECK_CLOSE_FRACTION( 2.0, ratio, 0.1);

}

BOOST_AUTO_TEST_SUITE_END()
