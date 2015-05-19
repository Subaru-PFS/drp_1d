#include "operator.h"

#include <epic/core/common/datatypes.h>
#include <epic/core/common/ref.h>
#include <epic/redshift/operator/correlation.h>
#include <epic/redshift/operator/correlationresult.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/spectrum/io/fitsreader.h>
#include <epic/redshift/extremum/extremum.h>
#include <epic/redshift/continuum/median.h>

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
    CRef<CCorrelationResult> r = (CCorrelationResult*) correlation.Compute( s, t, lambdaRange, redshifts, 0.99 );
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
    CRef<CCorrelationResult> r = (CCorrelationResult*) correlation.Compute( s, t, lambdaRange, redshifts, 0.7 );
    CPPUNIT_ASSERT( r != NULL );

    const TFloat64List& results = r->Correlation;
    const COperatorCorrelation::TStatusList& status = r->Status;

    CPPUNIT_ASSERT( results.size() == status.size() );

    CExtremum extremum;
    TPointList extremumList;
    extremum.Find( r->Redshifts, r->Correlation, extremumList );

    CPPUNIT_ASSERT_DOUBLES_EQUAL( z, extremumList[0].X, redshiftDelta*2 );


}
