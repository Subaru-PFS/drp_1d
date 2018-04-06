#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/operator/correlation.h>
#include <RedshiftLibrary/operator/correlationresult.h>
#include <RedshiftLibrary/operator/chisquare.h>
#include <RedshiftLibrary/operator/chisquareresult.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/io/fitsreader.h>
#include <RedshiftLibrary/spectrum/io/genericreader.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/continuum/median.h>
#include <RedshiftLibrary/noise/fromfile.h>

#include <fstream>
#include <boost/math/special_functions.hpp>

#include <boost/test/unit_test.hpp>
#include "test-config.h"

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(Operator)

/**
 * Correlate two spectrum over a given Z range: [0 - 3]
 * and assert that correlation is maximized at Z = 0.0
 */
BOOST_AUTO_TEST_CASE(CorrelationAtZEqualZero)
{
    Bool retVal;
    std::shared_ptr<CSpectrum> spectrum = std::shared_ptr<CSpectrum>( new CSpectrum() );
    std::shared_ptr<CTemplate> _template = std::shared_ptr<CTemplate>( new CTemplate() );

    CSpectrumIOFitsReader reader;

    retVal = reader.Read( DATA_ROOT_DIR "OperatorTestCase/spectrum1_z_1.2299.fits",
			  spectrum );
    BOOST_CHECK( retVal );
    retVal = reader.Read( DATA_ROOT_DIR "OperatorTestCase/spectrum1_z_1.2299.fits",
			  _template );
    BOOST_CHECK( retVal );

    spectrum->ConvertToLogScale();
    _template->ConvertToLogScale();

    TFloat64Range lambdaRange( spectrum->GetLambdaRange().GetBegin(),
			       spectrum->GetLambdaRange().GetEnd() );

    COperatorCorrelation correlation;
    TFloat64List redshifts;

    Float64 redshiftDelta = 0.0001;
    redshifts = TFloat64Range( 0.0, 3.0 ).SpreadOver( redshiftDelta );
    // prepare the unused masks
    std::vector<CMask> maskList;
    auto r = std::dynamic_pointer_cast<CCorrelationResult>( correlation.Compute( *spectrum, *_template, lambdaRange, redshifts, 0.99, maskList ) );
    BOOST_CHECK( r != NULL );


    CExtremum extremum;
    TPointList extremumList;
    extremum.Find( r->Redshifts, r->Correlation, extremumList );

    BOOST_CHECK_CLOSE_FRACTION( 0.0, extremumList[0].X, 0.000001 );


}

/**
 * Shift back a spectrum to it's rest pose (knowing it's z)
 * cross correlate between the shifted version and the unshifted one over a Z range of [0-3]
 * and check that the correlation factor is maximized at the expected Z
 */
BOOST_AUTO_TEST_CASE(CorrelationAtGivenZ)
{
    Bool retVal;
    std::shared_ptr<CSpectrum> spectrum = std::shared_ptr<CSpectrum>( new CSpectrum() );
    std::shared_ptr<CTemplate> _template = std::shared_ptr<CTemplate>( new CTemplate() );

    Float64 z = 1.2299;

    CSpectrumIOFitsReader reader;

    retVal = reader.Read( DATA_ROOT_DIR "OperatorTestCase/spectrum1_z_1.2299.fits",
			  spectrum );
    BOOST_CHECK( retVal );
    retVal = reader.Read( DATA_ROOT_DIR "OperatorTestCase/spectrum1_z_1.2299.fits",
			  _template );
    BOOST_CHECK( retVal );

    // Shift template back to rest pose
    CSpectrumSpectralAxis& tplSpectralAxis = _template->GetSpectralAxis();
    tplSpectralAxis.ShiftByWaveLength( 1.0 + z,  CSpectrumSpectralAxis::nShiftBackward );

    spectrum->ConvertToLogScale();
    _template->ConvertToLogScale();

    TFloat64Range lambdaRange( spectrum->GetLambdaRange().GetBegin(),
			       spectrum->GetLambdaRange().GetEnd() );

    Float64 redshiftDelta = 0.0001;
    TFloat64List redshifts = TFloat64Range( 0.0, 3.0 ).SpreadOver( redshiftDelta );

    //CRedshifts redshifts( &z, 1 );

    // prepare the unused masks
    std::vector<CMask> maskList;
    COperatorCorrelation correlation;
    auto r = std::dynamic_pointer_cast<CCorrelationResult>( correlation.Compute( *spectrum, *_template, lambdaRange, redshifts, 0.7, maskList ) );
    BOOST_CHECK( r != NULL );

    const TFloat64List& results = r->Correlation;
    const COperatorCorrelation::TStatusList& status = r->Status;

    BOOST_CHECK( results.size() == status.size() );

    CExtremum extremum;
    TPointList extremumList;
    extremum.Find( r->Redshifts, r->Correlation, extremumList );

    BOOST_CHECK_CLOSE_FRACTION( z, extremumList[0].X, redshiftDelta*2 );


}

void UtilCorrelationMatchWithEZ( const char* spectraPath, const char* noisePath, const char* tplPath, const char* resultPath )
{

    Bool retVal;
    std::shared_ptr<CSpectrum> spectrum = std::shared_ptr<CSpectrum>( new CSpectrum() );
    std::shared_ptr<CTemplate> _template = std::shared_ptr<CTemplate>( new CTemplate() );

    Float64 z = 1.2299;

    // Load spectrum and templates
    CSpectrumIOGenericReader reader;
    retVal = reader.Read( spectraPath, spectrum );
    BOOST_CHECK( retVal );

    if( noisePath )
    {
        CNoiseFromFile noise;
        noise.SetNoiseFilePath( noisePath );
        noise.AddNoise( *spectrum );
    }



    retVal = reader.Read( tplPath, _template );
    BOOST_CHECK( retVal );

    {
        CContinuumMedian continuum;
        spectrum->RemoveContinuum( continuum );
    }
    {
        CContinuumMedian continuum;
        _template->RemoveContinuum(continuum);
    }

    spectrum->ConvertToLogScale();
    _template->ConvertToLogScale();


    Float64 redshiftDelta = 0.0001;
    TFloat64List redshifts = TFloat64Range( 0.0, 2.0 ).SpreadOver( redshiftDelta );

    // prepare the unused masks
    std::vector<CMask> maskList;
    COperatorCorrelation correlation;
    auto r = std::dynamic_pointer_cast<CCorrelationResult>( correlation.Compute( *spectrum, *_template, TFloat64Range( 5600, 7000 ), redshifts, 1.0, maskList ) );
    BOOST_CHECK( r != NULL );

    CCorrelationResult referenceResult;

    std::ifstream input( resultPath );
    BOOST_CHECK( input.is_open() );

    referenceResult.Load( input );

    for( Int32 i=0; i<referenceResult.Correlation.size(); i++ )
    {
        if( boost::math::isnan( referenceResult.Correlation[i] ) )
        {
            BOOST_CHECK( boost::math::isnan( r->Correlation[i] ) );
        }
        else
        {
            BOOST_CHECK_CLOSE_FRACTION( referenceResult.Correlation[i], r->Correlation[i], 0.00001 );
        }


        BOOST_CHECK_CLOSE_FRACTION( referenceResult.Redshifts[i], r->Redshifts[i], 0.00001 );
        BOOST_CHECK_CLOSE_FRACTION( referenceResult.Overlap[i], r->Overlap[i], 0.00001 );
    }

}

BOOST_AUTO_TEST_CASE(CorrelationMatchWithEZ)
{

    UtilCorrelationMatchWithEZ( DATA_ROOT_DIR "OperatorTestCase/fromVVDSDeep/spectra/sc_020086397_F02P016_vmM1_red_31_1_atm_clean.fits",
                            NULL,
                            DATA_ROOT_DIR "OperatorTestCase/fromVVDSDeep/template/galaxy/zcosmos_red.txt",
                            DATA_ROOT_DIR "OperatorTestCase/fromVVDSDeep/results_nonoise/sc_020086397_F02P016_vmM1_red_31_1_atm_clean.csv" );
/*
    UtilCorrelationMatchWithEZ( DATA_ROOT_DIR "OperatorTestCase/fromVVDSDeep/spectra/sc_020100776_F02P017_vmM1_red_129_1_atm_clean.fits",
                            NULL,
                            DATA_ROOT_DIR "OperatorTestCase/fromVVDSDeep/template/galaxy/zcosmos_red.txt",
                            DATA_ROOT_DIR "OperatorTestCase/fromVVDSDeep/results_nonoise/sc_020100776_F02P017_vmM1_red_129_1_atm_clean.csv" );


    UtilCorrelationMatchWithEZ( DATA_ROOT_DIR "OperatorTestCase/fromVVDSDeep/spectra/sc_020088501_F02P017_vmM1_red_82_1_atm_clean.fits",
                            DATA_ROOT_DIR "OperatorTestCase/fromVVDSDeep/spectra/sc_020088501_F02P017_vmM1_red_82_1_noise.fits",
                            DATA_ROOT_DIR "OperatorTestCase/fromVVDSDeep/template/galaxy/zcosmos_red.txt",
                            DATA_ROOT_DIR "OperatorTestCase/fromVVDSDeep/results_withnoise/sc_020088501_F02P017_vmM1_red_82_1_atm_clean.csv" );


    UtilCorrelationMatchWithEZ( DATA_ROOT_DIR "OperatorTestCase/fromVVDSDeep/spectra/sc_020123432_F02P019_vmM1_red_72_1_atm_clean.fits",
                            DATA_ROOT_DIR "OperatorTestCase/fromVVDSDeep/spectra/sc_020123432_F02P019_vmM1_red_72_1_noise.fits",
                            DATA_ROOT_DIR "OperatorTestCase/fromVVDSDeep/template/galaxy/zcosmos_red.txt",
                            DATA_ROOT_DIR "OperatorTestCase/fromVVDSDeep/results_withnoise/sc_020123432_F02P019_vmM1_red_72_1_atm_clean.csv" );
                            */

}


void UtilChisquareMatchWithEZ( const char* spectraPath, const char* noisePath, const char* tplPath, const char* resultPath )
{
    Bool retVal;
    std::shared_ptr<CSpectrum> spectrum = std::shared_ptr<CSpectrum>( new CSpectrum() );
    std::shared_ptr<CTemplate> _template = std::shared_ptr<CTemplate>( new CTemplate() );

    Float64 z = 1.2299;

    // Load spectrum and templates
    CSpectrumIOGenericReader reader;
    retVal = reader.Read( spectraPath, spectrum);
    BOOST_CHECK( retVal );

    if( noisePath )
    {
        CNoiseFromFile noise;
        noise.SetNoiseFilePath( noisePath );
        noise.AddNoise( *spectrum );
    }

    retVal = reader.Read( tplPath, _template );
    BOOST_CHECK( retVal );

    Float64 redshiftDelta = 0.0001;
    TFloat64List redshifts = TFloat64Range( 0.0, 5.0 ).SpreadOver( redshiftDelta );
    Float64 overlapThres = 1.0;

    // prepare the unused masks
    std::vector<CMask> maskList;
    COperatorChiSquare chi;
    auto r = std::dynamic_pointer_cast<CChisquareResult>( chi.Compute( *spectrum, *_template, TFloat64Range( 0., 10000. ), redshifts, overlapThres, maskList ) );
    BOOST_CHECK( r != NULL );

    CChisquareResult referenceResult;

    std::ifstream input( resultPath );
    BOOST_CHECK( input.is_open() );

    referenceResult.Load( input );

    Bool redshiftsOK = true;
    Bool overlapsOK = true;
    Int32 nOverlaptested = 0;
    Int32 nRedshiftstested = 0;
    Int32 nChisquaretested = 0;

    for( Int32 i=0; i<referenceResult.ChiSquare.size(); i++ )
    {
        if(std::abs(referenceResult.Redshifts[i] - r->Redshifts[i])>redshiftDelta)
        {
            redshiftsOK = false;
            break;
        }else
        {
            nOverlaptested++;
        }

        if(std::abs(referenceResult.Overlap[i] - r->Overlap[i])>0.5)
        {
            overlapsOK = false;
            break;
        }else
        {
            nRedshiftstested++;
        }

        if(r->Overlap[i]>=overlapThres && referenceResult.Overlap[i]>=overlapThres)
        {
            if( boost::math::isnan( referenceResult.ChiSquare[i] ) )
            {
                BOOST_CHECK( boost::math::isnan( r->ChiSquare[i] ) );
            }
            else
            {
                BOOST_CHECK_CLOSE_FRACTION( referenceResult.ChiSquare[i], r->ChiSquare[i], 0.00001 );
                nChisquaretested++;
            }
        }

    }
    BOOST_CHECK_MESSAGE( redshiftsOK, "redshift values shift comparison with reference, FAILED" );
    BOOST_CHECK_MESSAGE( overlapsOK, "overlap values shift comparison with reference, FAILED" );
    BOOST_CHECK_MESSAGE( nOverlaptested>50000, "PB, less than 50000 overlap values tested" );
    BOOST_CHECK_MESSAGE( nRedshiftstested>50000, "PB, less than 50000 redshift values tested" );
    BOOST_CHECK_MESSAGE( nChisquaretested>20000, "PB, less than 30000 merit values tested" );
}

BOOST_AUTO_TEST_CASE(ChisquareMatchWithEZ)
{
    UtilChisquareMatchWithEZ( DATA_ROOT_DIR "OperatorTestCase/fromVVDSDeep/spectra/sc_020088501_F02P017_vmM1_red_82_1_atm_clean.fits",
                            DATA_ROOT_DIR "OperatorTestCase/fromVVDSDeep/spectra/sc_020088501_F02P017_vmM1_red_82_1_noise.fits",
                            DATA_ROOT_DIR "OperatorTestCase/fromVVDSDeep/template/galaxy/zcosmos_red.txt",
                            DATA_ROOT_DIR "OperatorTestCase/fromVVDSDeep/chisquare_results_withnoise/sc_020088501_F02P017_vmM1_red_82_1_atm_clean.csv" );


    UtilChisquareMatchWithEZ( DATA_ROOT_DIR "OperatorTestCase/fromVVDSDeep/spectra/sc_020123432_F02P019_vmM1_red_72_1_atm_clean.fits",
                            DATA_ROOT_DIR "OperatorTestCase/fromVVDSDeep/spectra/sc_020123432_F02P019_vmM1_red_72_1_noise.fits",
                            DATA_ROOT_DIR "OperatorTestCase/fromVVDSDeep/template/galaxy/zcosmos_red.txt",
                            DATA_ROOT_DIR "OperatorTestCase/fromVVDSDeep/chisquare_results_withnoise/sc_020123432_F02P019_vmM1_red_72_1_atm_clean.csv" );
}


BOOST_AUTO_TEST_SUITE_END()
