#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/operator/chisquare2.h>
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

BOOST_AUTO_TEST_SUITE(Operator_Chisquare_igm)

void UtilChisquareTestFit( const char* spectraPath, const char* noisePath, const char* tplPath, bool disableMask, const Float64 redshift, const Float64 targetFittedAmplitude, const Float64 targetFittedMeiksinIdx )
{
    Bool retVal;
    CSpectrum s;
    CTemplate t;

    Float64 z = redshift;
    BOOST_TEST_MESSAGE( "  Redshift = " << redshift );

    // Load spectrum and templates
    CSpectrumIOGenericReader reader;
    retVal = reader.Read( spectraPath, std::shared_ptr<CSpectrum>(&s) );
    BOOST_CHECK( retVal );

    if( noisePath )
    {
        CNoiseFromFile noise;
        noise.SetNoiseFilePath( noisePath );
        noise.AddNoise( s );
    }

    retVal = reader.Read( tplPath, std::shared_ptr<CTemplate>(&t) );
    BOOST_CHECK( retVal );

    Float64 redshiftDelta = 0.0001;
    TFloat64List redshifts = TFloat64Range( z, z ).SpreadOver( redshiftDelta );

    //building the mask
    Int32 sampleCount = s.GetFluxAxis().GetSamplesCount();
    std::vector<CMask> additional_spcMasks;
    CMask spcMask = Mask();
    spcMask.SetSize(sampleCount);

    std::string calibrationPath = DATA_ROOT_DIR "Operator_Chisquare_igmTestCase/calibration";

    COperatorChiSquare2 chi(calibrationPath);
    auto r = std::dynamic_pointer_cast<CChisquareResult>( chi.Compute( s, t, TFloat64Range( 200, 20000 ), redshifts, 1.0, additional_spcMasks, "precomputedfinegrid", 1 ) );
    BOOST_CHECK( r != NULL );
    BOOST_CHECK( r->Status[0] == COperatorChiSquare2::nStatus_OK );

    Float64 fit_amplitude = r->FitAmplitude[0];
    BOOST_CHECK_CLOSE_FRACTION( targetFittedAmplitude, fit_amplitude, 0.1 );

    Float64 fit_lstSquare = r->ChiSquare[0];
    BOOST_CHECK( fit_lstSquare < 0.01 );
    BOOST_TEST_MESSAGE( "  Lst-Square = " << fit_lstSquare );

    Float64 fit_meiksinIdx = r->FitMeiksinIdx[0];
    BOOST_TEST_MESSAGE( "  RefIdx = " << targetFittedMeiksinIdx << ", Meiksin Idx = " << fit_meiksinIdx );
    BOOST_CHECK( targetFittedMeiksinIdx == fit_meiksinIdx );



}

BOOST_AUTO_TEST_CASE(ChisquareTestCstFlux)
{
    Int32 targetFitMeiksinId = 5;
    Float64 z = 0.0;

    //*
    //z=0 test
    targetFitMeiksinId = 5;
    z = 0.0;
    UtilChisquareTestFit( DATA_ROOT_DIR "Operator_Chisquare_igmTestCase/fits_chisquare_igm_cstflux/spc_synth_z0_constant1p0_meiksin-z2c5.fits",
                            DATA_ROOT_DIR "Operator_Chisquare_igmTestCase/fits_chisquare_igm_cstflux/spc_synth_z0_constant0p5_ErrF.fits",
                            DATA_ROOT_DIR "Operator_Chisquare_igmTestCase/fits_chisquare_igm_cstflux/templates/galaxy/template_constantflux.dat", //constant template f=1.0 for all lambda
                            true,
                            z,
                            1.0,
                            targetFitMeiksinId);
    //*/

    //z=2.75 test
    targetFitMeiksinId = 5;
    z = 2.75;
    UtilChisquareTestFit( DATA_ROOT_DIR "Operator_Chisquare_igmTestCase/fits_chisquare_igm_cstflux/spc_synth_z2p75_constant1p0_meiksin-z3c5.fits",
                            DATA_ROOT_DIR "Operator_Chisquare_igmTestCase/fits_chisquare_igm_cstflux/spc_synth_z2p75_constant0p5_ErrF.fits",
                            DATA_ROOT_DIR "Operator_Chisquare_igmTestCase/fits_chisquare_igm_cstflux/templates/galaxy/template_constantflux.dat", //constant template f=1.0 for all lambda
                            true,
                            z,
                            1.0,
                            targetFitMeiksinId);

    //z=4.75 test
    targetFitMeiksinId = 1;
    z = 4.75;
    UtilChisquareTestFit( DATA_ROOT_DIR "Operator_Chisquare_igmTestCase/fits_chisquare_igm_cstflux/spc_synth_z4p75_constant1p0_meiksin-z5c1.fits",
                            DATA_ROOT_DIR "Operator_Chisquare_igmTestCase/fits_chisquare_igm_cstflux/spc_synth_z4p75_constant0p5_ErrF.fits",
                            DATA_ROOT_DIR "Operator_Chisquare_igmTestCase/fits_chisquare_igm_cstflux/templates/galaxy/template_constantflux.dat", //constant template f=1.0 for all lambda
                            true,
                            z,
                            1.0,
                            targetFitMeiksinId);

}


BOOST_AUTO_TEST_SUITE_END()
