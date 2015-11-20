#include "operator.h"

#include <epic/core/common/datatypes.h>
#include <epic/redshift/operator/correlation.h>
#include <epic/redshift/operator/correlationresult.h>
#include <epic/redshift/operator/chisquare.h>
#include <epic/redshift/operator/chisquareresult.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/spectrum/io/fitsreader.h>
#include <epic/redshift/spectrum/io/genericreader.h>
#include <epic/redshift/extremum/extremum.h>
#include <epic/redshift/continuum/irregularsamplingmedian.h>
#include <epic/redshift/noise/fromfile.h>

#include <boost/math/special_functions.hpp>

using namespace NSEpic;
using namespace NSEpicTest;

void CRedshiftOperatorTestCase::setUp()
{
}

void CRedshiftOperatorTestCase::tearDown()
{
}

/**
 * Correlate two spectrum over a given Z range: [0 - 3]
 * and assert that correlation is maximized at Z = 0.0
 */
void CRedshiftOperatorTestCase::CorrelationAtZEqualZero()
{
    Bool retVal;
    CSpectrum s;
    CTemplate t;

    CSpectrumIOFitsReader reader;

    retVal = reader.Read( "../test/data/OperatorTestCase/spectrum1_z_1.2299.fits", s );
    CPPUNIT_ASSERT( retVal );
    retVal = reader.Read( "../test/data/OperatorTestCase/spectrum1_z_1.2299.fits", t );
    CPPUNIT_ASSERT( retVal );

    s.ConvertToLogScale();
    t.ConvertToLogScale();

    TFloat64Range lambdaRange( s.GetLambdaRange().GetBegin(), s.GetLambdaRange().GetEnd() );

    COperatorCorrelation correlation;
    TFloat64List redshifts;

    Float64 redshiftDelta = 0.0001;
    redshifts = TFloat64Range( 0.0, 3.0 ).SpreadOver( redshiftDelta );
    auto r = std::dynamic_pointer_cast<CCorrelationResult>( correlation.Compute( s, t, lambdaRange, redshifts, 0.99 ) );
    CPPUNIT_ASSERT( r != NULL );


    CExtremum extremum;
    TPointList extremumList;
    extremum.Find( r->Redshifts, r->Correlation, extremumList );

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, extremumList[0].X, 0.000001 );


}

/**
 * Shift back a spectrum to it's rest pose (knowing it's z)
 * cross correlate between the shifted version and the unshifted one over a Z range of [0-3]
 * and check that the correlation factor is maximized at the expected Z
 */
void CRedshiftOperatorTestCase::CorrelationAtGivenZ()
{
    Bool retVal;
    CSpectrum s;
    CTemplate t;

    Float64 z = 1.2299;

    CSpectrumIOFitsReader reader;

    retVal = reader.Read( "../test/data/OperatorTestCase/spectrum1_z_1.2299.fits", s );
    CPPUNIT_ASSERT( retVal );
    retVal = reader.Read( "../test/data/OperatorTestCase/spectrum1_z_1.2299.fits", t );
    CPPUNIT_ASSERT( retVal );

    // Shift template back to rest pose
    CSpectrumSpectralAxis& tplSpectralAxis = t.GetSpectralAxis();
    tplSpectralAxis.ShiftByWaveLength( 1.0 + z,  CSpectrumSpectralAxis::nShiftBackward );

    s.ConvertToLogScale();
    t.ConvertToLogScale();

    TFloat64Range lambdaRange( s.GetLambdaRange().GetBegin(), s.GetLambdaRange().GetEnd() );


    Float64 redshiftDelta = 0.0001;
    TFloat64List redshifts = TFloat64Range( 0.0, 3.0 ).SpreadOver( redshiftDelta );

    //CRedshifts redshifts( &z, 1 );

    COperatorCorrelation correlation;
    auto r = std::dynamic_pointer_cast<CCorrelationResult>( correlation.Compute( s, t, lambdaRange, redshifts, 0.7 ) );
    CPPUNIT_ASSERT( r != NULL );

    const TFloat64List& results = r->Correlation;
    const COperatorCorrelation::TStatusList& status = r->Status;

    CPPUNIT_ASSERT( results.size() == status.size() );

    CExtremum extremum;
    TPointList extremumList;
    extremum.Find( r->Redshifts, r->Correlation, extremumList );

    CPPUNIT_ASSERT_DOUBLES_EQUAL( z, extremumList[0].X, redshiftDelta*2 );


}

void CRedshiftOperatorTestCase::CorrelationMatchWithEZ()
{

    CorrelationMatchWithEZ( "../test/data/OperatorTestCase/fromVVDSDeep/spectra/sc_020086397_F02P016_vmM1_red_31_1_atm_clean.fits",
                            NULL,
                            "../test/data/OperatorTestCase/fromVVDSDeep/template/galaxy/zcosmos_red.txt",
                            "../test/data/OperatorTestCase/fromVVDSDeep/results_nonoise/sc_020086397_F02P016_vmM1_red_31_1_atm_clean.csv" );
/*
    CorrelationMatchWithEZ( "../test/data/OperatorTestCase/fromVVDSDeep/spectra/sc_020100776_F02P017_vmM1_red_129_1_atm_clean.fits",
                            NULL,
                            "../test/data/OperatorTestCase/fromVVDSDeep/template/galaxy/zcosmos_red.txt",
                            "../test/data/OperatorTestCase/fromVVDSDeep/results_nonoise/sc_020100776_F02P017_vmM1_red_129_1_atm_clean.csv" );


    CorrelationMatchWithEZ( "../test/data/OperatorTestCase/fromVVDSDeep/spectra/sc_020088501_F02P017_vmM1_red_82_1_atm_clean.fits",
                            "../test/data/OperatorTestCase/fromVVDSDeep/spectra/sc_020088501_F02P017_vmM1_red_82_1_noise.fits",
                            "../test/data/OperatorTestCase/fromVVDSDeep/template/galaxy/zcosmos_red.txt",
                            "../test/data/OperatorTestCase/fromVVDSDeep/results_withnoise/sc_020088501_F02P017_vmM1_red_82_1_atm_clean.csv" );


    CorrelationMatchWithEZ( "../test/data/OperatorTestCase/fromVVDSDeep/spectra/sc_020123432_F02P019_vmM1_red_72_1_atm_clean.fits",
                            "../test/data/OperatorTestCase/fromVVDSDeep/spectra/sc_020123432_F02P019_vmM1_red_72_1_noise.fits",
                            "../test/data/OperatorTestCase/fromVVDSDeep/template/galaxy/zcosmos_red.txt",
                            "../test/data/OperatorTestCase/fromVVDSDeep/results_withnoise/sc_020123432_F02P019_vmM1_red_72_1_atm_clean.csv" );
                            */

}

void CRedshiftOperatorTestCase::CorrelationMatchWithEZ( const char* spectraPath, const char* noisePath, const char* tplPath, const char* resultPath )
{

    Bool retVal;
    CSpectrum s;
    CTemplate t;

    Float64 z = 1.2299;

    // Load spectrum and templates
    CSpectrumIOGenericReader reader;
    retVal = reader.Read( spectraPath, s );
    CPPUNIT_ASSERT( retVal );

    if( noisePath )
    {
        CNoiseFromFile noise;
        noise.SetNoiseFilePath( noisePath );
        noise.AddNoise( s );
    }



    retVal = reader.Read( tplPath, t );
    CPPUNIT_ASSERT( retVal );

    {
        CContinuumIrregularSamplingMedian continuum;
        s.RemoveContinuum( continuum );
    }
    {
        CContinuumIrregularSamplingMedian continuum;
        t.RemoveContinuum(continuum);
    }

    s.ConvertToLogScale();
    t.ConvertToLogScale();


    Float64 redshiftDelta = 0.0001;
    TFloat64List redshifts = TFloat64Range( 0.0, 2.0 ).SpreadOver( redshiftDelta );

    COperatorCorrelation correlation;
    auto r = std::dynamic_pointer_cast<CCorrelationResult>( correlation.Compute( s, t, TFloat64Range( 5600, 7000 ), redshifts, 1.0 ) );
    CPPUNIT_ASSERT( r != NULL );

    CCorrelationResult referenceResult;

    std::ifstream input( resultPath );
    CPPUNIT_ASSERT( input.is_open() );

    referenceResult.Load( input );

    for( Int32 i=0; i<referenceResult.Correlation.size(); i++ )
    {
        if( boost::math::isnan( referenceResult.Correlation[i] ) )
        {
            CPPUNIT_ASSERT( boost::math::isnan( r->Correlation[i] ) );
        }
        else
        {
            CPPUNIT_ASSERT_DOUBLES_EQUAL( referenceResult.Correlation[i], r->Correlation[i], 0.00001 );
        }


        CPPUNIT_ASSERT_DOUBLES_EQUAL( referenceResult.Redshifts[i], r->Redshifts[i], 0.00001 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL( referenceResult.Overlap[i], r->Overlap[i], 0.00001 );
    }

}

void CRedshiftOperatorTestCase::ChisquareMatchWithEZ()
{
    ChisquareMatchWithEZ( "../test/data/OperatorTestCase/fromVVDSDeep/spectra/sc_020088501_F02P017_vmM1_red_82_1_atm_clean.fits",
                            "../test/data/OperatorTestCase/fromVVDSDeep/spectra/sc_020088501_F02P017_vmM1_red_82_1_noise.fits",
                            "../test/data/OperatorTestCase/fromVVDSDeep/template/galaxy/zcosmos_red.txt",
                            "../test/data/OperatorTestCase/fromVVDSDeep/chisquare_results_withnoise/sc_020088501_F02P017_vmM1_red_82_1_atm_clean.csv" );


    ChisquareMatchWithEZ( "../test/data/OperatorTestCase/fromVVDSDeep/spectra/sc_020123432_F02P019_vmM1_red_72_1_atm_clean.fits",
                            "../test/data/OperatorTestCase/fromVVDSDeep/spectra/sc_020123432_F02P019_vmM1_red_72_1_noise.fits",
                            "../test/data/OperatorTestCase/fromVVDSDeep/template/galaxy/zcosmos_red.txt",
                            "../test/data/OperatorTestCase/fromVVDSDeep/chisquare_results_withnoise/sc_020123432_F02P019_vmM1_red_72_1_atm_clean.csv" );
}


void CRedshiftOperatorTestCase::ChisquareMatchWithEZ( const char* spectraPath, const char* noisePath, const char* tplPath, const char* resultPath )
{
    Bool retVal;
    CSpectrum s;
    CTemplate t;

    Float64 z = 1.2299;

    // Load spectrum and templates
    CSpectrumIOGenericReader reader;
    retVal = reader.Read( spectraPath, s );
    CPPUNIT_ASSERT( retVal );

    if( noisePath )
    {
        CNoiseFromFile noise;
        noise.SetNoiseFilePath( noisePath );
        noise.AddNoise( s );
    }

    retVal = reader.Read( tplPath, t );
    CPPUNIT_ASSERT( retVal );

    Float64 redshiftDelta = 0.0001;
    TFloat64List redshifts = TFloat64Range( 0.0, 2.0 ).SpreadOver( redshiftDelta );

    COperatorChiSquare chi;
    auto r = std::dynamic_pointer_cast<CChisquareResult>( chi.Compute( s, t, TFloat64Range( 5600, 7000 ), redshifts, 1.0 ) );
    CPPUNIT_ASSERT( r != NULL );

    CChisquareResult referenceResult;

    std::ifstream input( resultPath );
    CPPUNIT_ASSERT( input.is_open() );

    referenceResult.Load( input );

    for( Int32 i=0; i<referenceResult.ChiSquare.size(); i++ )
    {
        if( boost::math::isnan( referenceResult.ChiSquare[i] ) )
        {
            CPPUNIT_ASSERT( boost::math::isnan( r->ChiSquare[i] ) );
        }
        else
        {
            CPPUNIT_ASSERT_DOUBLES_EQUAL( referenceResult.ChiSquare[i], r->ChiSquare[i], 0.00001 );
        }


        CPPUNIT_ASSERT_DOUBLES_EQUAL( referenceResult.Redshifts[i], r->Redshifts[i], 0.00001 );
        CPPUNIT_ASSERT_DOUBLES_EQUAL( referenceResult.Overlap[i], r->Overlap[i], 0.00001 );
    }

}
