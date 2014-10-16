#include "spectrumio.h"

#include <epic/core/common/datatypes.h>
#include <epic/redshift/spectrum/spectrum.h>

using namespace __NS__;
using namespace std;

void CCRedshiftSpectrumioTestCase::setUp()
{
}

void CCRedshiftSpectrumioTestCase::tearDown()
{
}

#include <epic/redshift/spectrum/io/fitsreader.h>

void CCRedshiftSpectrumioTestCase::VVDSReadValidFile()
{
    CSpectrumIOFitsReader reader;

    CSpectrum s;

    Bool retVal = reader.Read( "../test/redshift/data/spectrum1_z_1.2299.fits", s );
    CPPUNIT_ASSERT( retVal == true );

    CPPUNIT_ASSERT( s.GetSampleCount() == 11391 );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 3800.0, s.GetLambdaRange().GetBegin(), 0.01 );

}

void CCRedshiftSpectrumioTestCase::VVDSReadInvalidFile()
{
    CSpectrumIOFitsReader reader;

    CSpectrum s;

    Bool rValue = reader.Read( "../test/redshift/data/invalidspectrum1.fits", s );

    CPPUNIT_ASSERT( rValue == false );
    CPPUNIT_ASSERT( s.GetSampleCount() == 0 );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, s.GetLambdaRange().GetBegin(), 0.01 );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, s.GetResolution(), 0.01  );

}


