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
#include "test-config.h"

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(Operator_Chisquare_ism)


void UtilChisquareTestFit( const char* spectraPath,
                           const char* noisePath,
                           const char* tplPath,
                           bool disableMask,
                           const Float64 redshift,
                           const Float64 targetFittedAmplitude,
                           const Float64 targetFitIsmCalzettiCoeff,
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

    std::string calibrationPath = DATA_ROOT_DIR "Operator_Chisquare_ismTestCase/calibration";

    std::shared_ptr<CChisquareResult> r;
    Int32 enable_IGM = 0;
    Int32 enable_ISM = 1;
    if(chisquareOperator=="chisquare2")
    {
        COperatorChiSquare2 chi(calibrationPath);
        r = std::dynamic_pointer_cast<CChisquareResult>( chi.Compute( s, t, TFloat64Range( 200, 20000 ), redshifts, 1.0, additional_spcMasks, "precomputedfinegrid", enable_IGM, enable_ISM ) );
        BOOST_CHECK( r != NULL );
        BOOST_CHECK( r->Status[0] == COperatorChiSquare2::nStatus_OK );
    }
    else
    {
        bool enableRebinLog = true;
        COperatorChiSquareLogLambda chi(calibrationPath);
        chi.enableSpcLogRebin(enableRebinLog);
        r = std::dynamic_pointer_cast<CChisquareResult>( chi.Compute( s, t, lambdaRange, redshifts, 1.0, additional_spcMasks, "precomputedfinegrid", enable_IGM, enable_ISM  ) );
        BOOST_CHECK( r != NULL );
        BOOST_CHECK( r->Status[0] == COperatorChiSquareLogLambda::nStatus_OK );
    }

    Float64 fit_amplitude = r->FitAmplitude[0];
    BOOST_CHECK_CLOSE_FRACTION( targetFittedAmplitude, fit_amplitude, 0.1 );

    Float64 fit_lstSquare = r->ChiSquare[0];
    BOOST_CHECK( fit_lstSquare < 2.0 );
    BOOST_TEST_MESSAGE( "  Lst-Square = " << fit_lstSquare );

    Float64 fit_dustCoeff = r->FitDustCoeff[0];
    BOOST_TEST_MESSAGE( "  RefIdx = " << targetFitIsmCalzettiCoeff << ", Calzetti coeff = " << fit_dustCoeff );
    BOOST_CHECK( targetFitIsmCalzettiCoeff == fit_dustCoeff );

    r.reset();
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

        Float64 targetFitIsmCalzetti = 0.0;
        Float64 z = 0.0;

        //z=0 test 0
        targetFitIsmCalzetti = 0.0;
        z = 0.0;
        UtilChisquareTestFit( DATA_ROOT_DIR "Operator_Chisquare_ismTestCase/fits_chisquare_ism_cstflux/spc_synth_z0_a1_constantflux_gmu6500gw15ga025_ismCalzetti0_200A-10kA_TF.fits",
                              DATA_ROOT_DIR "Operator_Chisquare_ismTestCase/fits_chisquare_ism_cstflux/spc_synth_z0_a1_constantflux_gmu6500gw15ga025_200A-10kA_ErrF.fits",
                              DATA_ROOT_DIR "Operator_Chisquare_ismTestCase/fits_chisquare_ism_cstflux/templates/galaxy/template_constantflux_gmu6500gw15ga025.dat",
                              true,
                              z,
                              1.0,
                              targetFitIsmCalzetti,
                              chisquareOperator);

        //z=0 test 1
        targetFitIsmCalzetti = 0.2;
        z = 0.0;
        UtilChisquareTestFit( DATA_ROOT_DIR "Operator_Chisquare_ismTestCase/fits_chisquare_ism_cstflux/spc_synth_z0_a1_constantflux_gmu6500gw15ga025_ismCalzetti02_200A-10kA_TF.fits",
                              DATA_ROOT_DIR "Operator_Chisquare_ismTestCase/fits_chisquare_ism_cstflux/spc_synth_z0_a1_constantflux_gmu6500gw15ga025_200A-10kA_ErrF.fits",
                              DATA_ROOT_DIR "Operator_Chisquare_ismTestCase/fits_chisquare_ism_cstflux/templates/galaxy/template_constantflux_gmu6500gw15ga025.dat",
                              true,
                              z,
                              1.0,
                              targetFitIsmCalzetti,
                              chisquareOperator);


        //z=0 test 2
        targetFitIsmCalzetti = 0.8;
        z = 0.0;
        UtilChisquareTestFit( DATA_ROOT_DIR "Operator_Chisquare_ismTestCase/fits_chisquare_ism_cstflux/spc_synth_z0_a1_constantflux_gmu6500gw15ga025_ismCalzetti08_200A-10kA_TF.fits",
                              DATA_ROOT_DIR "Operator_Chisquare_ismTestCase/fits_chisquare_ism_cstflux/spc_synth_z0_a1_constantflux_gmu6500gw15ga025_200A-10kA_ErrF.fits",
                              DATA_ROOT_DIR "Operator_Chisquare_ismTestCase/fits_chisquare_ism_cstflux/templates/galaxy/template_constantflux_gmu6500gw15ga025.dat",
                              true,
                              z,
                              1.0,
                              targetFitIsmCalzetti,
                              chisquareOperator);

    }
}


BOOST_AUTO_TEST_SUITE_END()
