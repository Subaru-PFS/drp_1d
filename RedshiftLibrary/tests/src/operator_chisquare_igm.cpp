#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/operator/chisquare2.h>
#include <RedshiftLibrary/operator/chisquareloglambda.h>
#include <RedshiftLibrary/operator/chisquareresult.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/io/fitsreader.h>
#include <RedshiftLibrary/spectrum/io/genericreader.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/continuum/irregularsamplingmedian.h>
#include <RedshiftLibrary/noise/fromfile.h>

#include <fstream>
#include <boost/math/special_functions.hpp>

#include <boost/test/unit_test.hpp>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(Operator_Chisquare_igm)


void UtilChisquareTestFit( const char* spectraPath,
                           const char* noisePath,
                           const char* tplPath,
                           bool disableMask,
                           const Float64 redshift,
                           const Float64 targetFittedAmplitude,
                           const Float64 targetFittedMeiksinIdx,
                           std::string chisquareOperator)
{
    Bool retVal;
    CSpectrum s;
    CTemplate t;

    BOOST_TEST_MESSAGE( "\n  Operator = " << chisquareOperator );

    Float64 z = redshift;
    BOOST_TEST_MESSAGE( "  Redshift = " << redshift );

    // Load spectrum and templates
    CSpectrumIOGenericReader reader;
    retVal = reader.Read( spectraPath, s );
    BOOST_CHECK( retVal );

    if( noisePath )
    {
        CNoiseFromFile noise;
        noise.SetNoiseFilePath( noisePath );
        noise.AddNoise( s );
    }

    retVal = reader.Read( tplPath, t );
    BOOST_CHECK( retVal );

    Float64 redshiftDelta = 0.0001;
    TFloat64List redshifts = TFloat64Range( z, z ).SpreadOver( redshiftDelta );

    TFloat64Range lambdaRange = TFloat64Range( 300*(1+z), 4000*(1+z) );

    //building the mask
    Int32 sampleCount = s.GetFluxAxis().GetSamplesCount();
    std::vector<CMask> additional_spcMasks;
    CMask spcMask = Mask();
    spcMask.SetSize(sampleCount);

    std::string calibrationPath = "../RedshiftLibrary/tests/src/data/Operator_Chisquare_igmTestCase/calibration";

    std::shared_ptr<CChisquareResult> r;

    if(chisquareOperator=="chisquare2")
    {
        COperatorChiSquare2 chi(calibrationPath);
        r = std::dynamic_pointer_cast<CChisquareResult>( chi.Compute( s, t, TFloat64Range( 200, 20000 ), redshifts, 1.0, additional_spcMasks, "precomputedfinegrid", 1 ) );
        BOOST_CHECK( r != NULL );
        BOOST_CHECK( r->Status[0] == COperatorChiSquare2::nStatus_OK );
    }
    else
    {
        bool enableRebinLog = true;
        COperatorChiSquareLogLambda chi(calibrationPath);
        chi.enableSpcLogRebin(enableRebinLog);
        r = std::dynamic_pointer_cast<CChisquareResult>( chi.Compute( s, t, lambdaRange, redshifts, 1.0, additional_spcMasks, "precomputedfinegrid", 1 ) );
        BOOST_CHECK( r != NULL );
        BOOST_CHECK( r->Status[0] == COperatorChiSquareLogLambda::nStatus_OK );
    }

    Float64 fit_amplitude = r->FitAmplitude[0];
    BOOST_CHECK_CLOSE_FRACTION( targetFittedAmplitude, fit_amplitude, 0.1 );

    Float64 fit_lstSquare = r->ChiSquare[0];
    BOOST_CHECK( fit_lstSquare < 2.0 );
    BOOST_TEST_MESSAGE( "  Lst-Square = " << fit_lstSquare );

    Float64 fit_meiksinIdx = r->FitMeiksinIdx[0];
    BOOST_TEST_MESSAGE( "  RefIdx = " << targetFittedMeiksinIdx << ", Meiksin Idx = " << fit_meiksinIdx );
    BOOST_CHECK( targetFittedMeiksinIdx == fit_meiksinIdx );

}

BOOST_AUTO_TEST_CASE(ChisquareTestCstFlux)
{
    std::string chisquareOperator="chisquare2";

    for(UInt32 k=0; k<2; k++)
    {
        if(k==0)
        {
            chisquareOperator="chisquare2";
        }else{
            chisquareOperator="chisquarelog";
        }

        Int32 targetFitMeiksinId = 5;
        Float64 z = 0.0;

        //*
        //z=0 test
        targetFitMeiksinId = 5;
        z = 0.0;
        UtilChisquareTestFit( "../RedshiftLibrary/tests/src/data/Operator_Chisquare_igmTestCase/fits_chisquare_igm_cstflux/spc_synth_z0_constant1p0_meiksin-z2c5.fits",
                              "../RedshiftLibrary/tests/src/data/Operator_Chisquare_igmTestCase/fits_chisquare_igm_cstflux/spc_synth_z0_constant0p5_ErrF.fits",
                              "../RedshiftLibrary/tests/src/data/Operator_Chisquare_igmTestCase/fits_chisquare_igm_cstflux/templates/galaxy/template_constantflux.dat", //constant template f=1.0 for all lambda
                              true,
                              z,
                              1.0,
                              targetFitMeiksinId,
                              chisquareOperator);
        //*/

        //z=2.75 test
        targetFitMeiksinId = 5;
        z = 2.75;
        UtilChisquareTestFit( "../RedshiftLibrary/tests/src/data/Operator_Chisquare_igmTestCase/fits_chisquare_igm_cstflux/spc_synth_z2p75_constant1p0_meiksin-z3c5.fits",
                              "../RedshiftLibrary/tests/src/data/Operator_Chisquare_igmTestCase/fits_chisquare_igm_cstflux/spc_synth_z2p75_constant0p5_ErrF.fits",
                              "../RedshiftLibrary/tests/src/data/Operator_Chisquare_igmTestCase/fits_chisquare_igm_cstflux/templates/galaxy/template_constantflux.dat", //constant template f=1.0 for all lambda
                              true,
                              z,
                              1.0,
                              targetFitMeiksinId,
                              chisquareOperator);

        //z=4.75 test
        targetFitMeiksinId = 1;
        z = 4.75;
        UtilChisquareTestFit( "../RedshiftLibrary/tests/src/data/Operator_Chisquare_igmTestCase/fits_chisquare_igm_cstflux/spc_synth_z4p75_constant1p0_meiksin-z5c1.fits",
                              "../RedshiftLibrary/tests/src/data/Operator_Chisquare_igmTestCase/fits_chisquare_igm_cstflux/spc_synth_z4p75_constant0p5_ErrF.fits",
                              "../RedshiftLibrary/tests/src/data/Operator_Chisquare_igmTestCase/fits_chisquare_igm_cstflux/templates/galaxy/template_constantflux.dat", //constant template f=1.0 for all lambda
                              true,
                              z,
                              1.0,
                              targetFitMeiksinId,
                              chisquareOperator);
    }
}


BOOST_AUTO_TEST_SUITE_END()
