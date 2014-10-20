#include "operator.h"

#include <epic/core/common/datatypes.h>
#include <epic/redshift/operator/correlation.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/spectrum/io/fitsreader.h>
#include <epic/redshift/extremum/extremum.h>
#include <epic/redshift/continuum/median.h>

using namespace NSEpic;

void CRedshiftOperatorTestCase::setUp()
{
}

void CRedshiftOperatorTestCase::tearDown()
{
}

void CRedshiftOperatorTestCase::Correlation4()
{
    Bool retVal;
    CSpectrum s;
    CTemplate t;

    CSpectrumIOFitsReader reader;

    retVal = reader.Read( "../test/redshift/data/spectrum1_z_1.2299.fits", s );
    CPPUNIT_ASSERT( retVal );
    retVal = reader.Read( "../test/redshift/data/spectrum1_z_1.2299.fits", t );
    CPPUNIT_ASSERT( retVal );

    s.ConvertToLogScale();
    t.ConvertToLogScale();

    TFloat64Range lambdaRange( s.GetLambdaRange().GetBegin(), s.GetLambdaRange().GetEnd() );

    COperatorCorrelation correlation;
    CRedshifts redshifts;

    Float64 redshiftDelta = 0.0001;
    redshifts.SpreadOver( 0.0, 3.0, redshiftDelta );
    Bool r = correlation.Compute( s, t, lambdaRange, redshifts, 1.0 );
    CPPUNIT_ASSERT( r == true );

    const TFloat64List& results = correlation.GetResults();
    const COperatorCorrelation::TStatusList& status = correlation.GetStatus();

    CPPUNIT_ASSERT( results.size() == status.size() );

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, results[0], redshiftDelta*2 );

    for( Int32 i = 1; i<results.size(); i++ )
    {
        CPPUNIT_ASSERT( status[i] != COperatorCorrelation::nStatus_OK );
    }

}

void CRedshiftOperatorTestCase::Correlation1()
{
    Bool retVal;
    CSpectrum s;
    CTemplate t;

    Float64 z = 1.2299;

    CSpectrumIOFitsReader reader;

    // Those spectrum are already at Z = 1.2299
    retVal = reader.Read( "../test/redshift/data/spectrum1_z_1.2299.fits", s );
    CPPUNIT_ASSERT( retVal );
    retVal = reader.Read( "../test/redshift/data/spectrum1_z_1.2299.fits", t );
    CPPUNIT_ASSERT( retVal );

    // Shift template back to rest pose
    CSpectrumAxis& tplSpectralAxis = t.GetSpectralAxis();
    for( Int32 i=0;i<tplSpectralAxis.GetSamplesCount();i++ )
    {
        tplSpectralAxis[i] = tplSpectralAxis[i] / (1.0 + z );
    }
    s.RemoveContinuum<CContinuumMedian>();
    s.ConvertToLogScale();

    t.RemoveContinuum<CContinuumMedian>();
    t.ConvertToLogScale();

    TFloat64Range lambdaRange( s.GetLambdaRange().GetBegin(), s.GetLambdaRange().GetEnd() );
    CRedshifts redshifts( &z, 1 );

    // Compute correlation for Z = 1.2299
    // So template is moved back to its original position, and correlation result should be maximized ( i.e: close to 1.0 )
    COperatorCorrelation correlation;
    Bool r = correlation.Compute( s, t, lambdaRange, redshifts, 0.99 );
    CPPUNIT_ASSERT( r == true );

    const TFloat64List& results = correlation.GetResults();
    const COperatorCorrelation::TStatusList& status = correlation.GetStatus();
    CPPUNIT_ASSERT( results.size() == status.size() );

    // test that correlation result is maximized ( i.e: close to 1.0 )
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, results[0], 0.0002 );


}

void CRedshiftOperatorTestCase::Correlation3()
{
    Bool retVal;
    CSpectrum s;
    CTemplate t;

    Float64 z = 1.2299;

    CSpectrumIOFitsReader reader;

    retVal = reader.Read( "../test/redshift/data/spectrum1_z_1.2299.fits", s );
    CPPUNIT_ASSERT( retVal );
    retVal = reader.Read( "../test/redshift/data/spectrum1_z_1.2299.fits", t );
    CPPUNIT_ASSERT( retVal );

    CSpectrumSpectralAxis& tplSpectralAxis = t.GetSpectralAxis();
    for( Int32 i=0;i<tplSpectralAxis.GetSamplesCount();i++ )
    {
        tplSpectralAxis[i] = tplSpectralAxis[i] / (1.0 + z );
    }

    s.RemoveContinuum<CContinuumMedian>();
    s.ConvertToLogScale();

    t.RemoveContinuum<CContinuumMedian>();
    t.ConvertToLogScale();

    TFloat64Range lambdaRange( s.GetLambdaRange().GetBegin(), s.GetLambdaRange().GetEnd() );


    Float64 redshiftDelta = 0.0001;
    CRedshifts redshifts( 0.0, 3.0, redshiftDelta );

    //CRedshifts redshifts( &z, 1 );

    COperatorCorrelation correlation;
    Bool r = correlation.Compute( s, t, lambdaRange, redshifts, 0.7 );
    CPPUNIT_ASSERT( r == true );

    const TFloat64List& results = correlation.GetResults();
    const COperatorCorrelation::TStatusList& status = correlation.GetStatus();

    CPPUNIT_ASSERT( results.size() == status.size() );

    CExtremum extremum;
    TPointList extremumList;
    extremum.Find( redshifts.GetRedshifts(), correlation.GetResults().data(), redshifts.GetRedshiftsCount(), extremumList );

    CPPUNIT_ASSERT_DOUBLES_EQUAL( z, extremumList[0].X, redshiftDelta*2 );


}

void CRedshiftOperatorTestCase::Correlation2()
{
    Bool retVal;
    CSpectrum s;
    CTemplate t;

    Float64 z = 1.11219;

    CSpectrumIOFitsReader reader;

    retVal = reader.Read( "../test/redshift/data/spectrum2_z_1.11219.fits", s );
    CPPUNIT_ASSERT( retVal );

    retVal = reader.Read( "../test/redshift/data/spectrum2_z_1.11219.fits", t );
    CPPUNIT_ASSERT( retVal );

    CSpectrumAxis& tplSpectralAxis = t.GetSpectralAxis();
    for( Int32 i=0;i<tplSpectralAxis.GetSamplesCount();i++ )
    {
        tplSpectralAxis[i] = tplSpectralAxis[i] / (1.0 + z );
    }

    t.RemoveContinuum<CContinuumMedian>();
    s.RemoveContinuum<CContinuumMedian>();

    t.ConvertToLogScale();
    s.ConvertToLogScale();

    TFloat64Range lambdaRange( s.GetLambdaRange().GetBegin(), s.GetLambdaRange().GetEnd() );


    Float64 redshiftDelta = 0.0001;
    CRedshifts redshifts( 0.0, 3.0, redshiftDelta );

    //CRedshifts redshifts( &z, 1 );

    COperatorCorrelation correlation;
    Bool r = correlation.Compute( s, t, lambdaRange, redshifts, 0.7 );
    CPPUNIT_ASSERT( r == true );

    const TFloat64List& results = correlation.GetResults();
    const COperatorCorrelation::TStatusList& status = correlation.GetStatus();

    CPPUNIT_ASSERT( results.size() == status.size() );

    CExtremum extremum;
    TPointList extremumList;
    extremum.Find( redshifts.GetRedshifts(), correlation.GetResults().data(), redshifts.GetRedshiftsCount(), extremumList );

    CPPUNIT_ASSERT_DOUBLES_EQUAL( z, extremumList[0].X, redshiftDelta*2 );


}
