#include "gaussianfit.h"

#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>
#include <epic/redshift/gaussianfit/gaussianfit.h>
#include <epic/redshift/spectrum/io/fitsreader.h>
#include <epic/redshift/spectrum/spectrum.h>

using namespace NSEpic;

using namespace NSEpicTest;

void CRedshiftGaussianFitTestCase::setUp()
{
}

void CRedshiftGaussianFitTestCase::tearDown()
{
}

void CRedshiftGaussianFitTestCase::TestFit1()
{
    CSpectrumIOFitsReader reader;
    CSpectrum spectrum;

    Bool retVal = reader.Read( "../test/data/GaussianfitTestCase/gaussian.fits", spectrum );
    CPPUNIT_ASSERT_MESSAGE( "Failed to load test data", retVal == true );


    CGaussianFit fitter;
    CGaussianFit::EStatus status = fitter.Compute( spectrum, TInt32Range( 3000, 5000 ) );
    //CPPUNIT_ASSERT_MESSAGE( "Failed to fit gaussian", status & CGaussianFit::nStatus_Success );

    Float64 gaussAmp;
    Float64 gaussPos;
    Float64 gaussWidth;
    fitter.GetResults( gaussAmp, gaussPos, gaussWidth );

    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE( "Invalid gaussian amplitude", 1.0, gaussAmp, 0.01 );
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE( "Invalid gaussian position", 4000.0, gaussPos, 2 );
    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE( "Invalid gaussian width", 565.0/1.4142, gaussWidth, 1.0);
}


void CRedshiftGaussianFitTestCase::TestFit2()
{
    CSpectrumIOFitsReader reader;
    CSpectrum spectrum;

    Bool retVal = reader.Read( "../test/data/GaussianfitTestCase/flat.fits", spectrum );
    CPPUNIT_ASSERT_MESSAGE( "Failed to load test data", retVal == true );


    CGaussianFit fitter;
    CGaussianFit::EStatus status = fitter.Compute( spectrum, TInt32Range( 3000, 5000 ) );
    CPPUNIT_ASSERT_MESSAGE( "Invalid return code", status == CGaussianFit::nStatus_IterationHasNotConverged );

    Float64 gaussAmp;
    Float64 gaussPos;
    Float64 gaussWidth;
    fitter.GetResults( gaussAmp, gaussPos, gaussWidth );

    CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE( "Invalid gaussian amplitude", 0.0, gaussAmp, 0.01 );


}


