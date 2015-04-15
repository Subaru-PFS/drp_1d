#include "smooth.h"

#include <epic/core/common/datatypes.h>

using namespace NSEpic;

using namespace NSEpicTest;

void CRedshiftSmoothTestCase::setUp()
{
}

void CRedshiftSmoothTestCase::tearDown()
{
}

#include <epic/redshift/spectrum/fluxaxis.h>
#include <epic/redshift/spectrum/spectrum.h>

void CRedshiftSmoothTestCase::Mean()
{
    CSpectrum spectrum;

    spectrum.GetFluxAxis().SetSize( 3 );
    Float64* data = spectrum.GetFluxAxis().GetSamples();
    data[0] = 10;
    data[1] = 5;
    data[2] = 10;

    spectrum.GetFluxAxis().ApplyMeanSmooth( 1 );

    CPPUNIT_ASSERT_DOUBLES_EQUAL( 7.5, spectrum.GetFluxAxis()[0], 0.0001 );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 8.3333, spectrum.GetFluxAxis()[1], 0.0001 );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 7.5, spectrum.GetFluxAxis()[2], 0.0001 );
}

