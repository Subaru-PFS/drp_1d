#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>
#include <epic/redshift/gaussianfit/gaussianfit.h>
#include <epic/redshift/spectrum/io/fitsreader.h>
#include <epic/redshift/spectrum/spectrum.h>

#include <boost/test/unit_test.hpp>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(GaussianFit)

BOOST_AUTO_TEST_CASE(TestFit1)
{
    CSpectrumIOFitsReader reader;
    CSpectrum spectrum;

    Bool retVal = reader.Read( "../test/data/GaussianfitTestCase/gaussian.fits", spectrum );
    BOOST_CHECK( retVal == true );


    CGaussianFit fitter;
    CGaussianFit::EStatus status = fitter.Compute( spectrum, TInt32Range( 3000, 5000 ) );
    //CPPUNIT_ASSERT_MESSAGE( "Failed to fit gaussian", status & CGaussianFit::nStatus_Success );

    Float64 gaussAmp;
    Float64 gaussPos;
    Float64 gaussWidth;
    fitter.GetResults( gaussAmp, gaussPos, gaussWidth );


    Float64 gaussAmpErr;
    Float64 gaussPosErr;
    Float64 gaussWidthErr;
    fitter.GetResultsError( gaussAmpErr, gaussPosErr, gaussWidthErr );

    BOOST_CHECK_CLOSE_FRACTION( 1.0, gaussAmp, 0.01 );
    BOOST_CHECK_CLOSE_FRACTION( 4000.0, gaussPos, 2 );
    BOOST_CHECK_CLOSE_FRACTION( 565.0/1.4142, gaussWidth, 1.0);
}

BOOST_AUTO_TEST_CASE(TestFit1_1em17)
{
    CSpectrumIOFitsReader reader;
    CSpectrum spectrum;

    Bool retVal = reader.Read( "../test/data/GaussianfitTestCase/gaussian_1em17.fits", spectrum );
    BOOST_CHECK( retVal == true );


    CGaussianFit fitter;
    CGaussianFit::EStatus status = fitter.Compute( spectrum, TInt32Range( 3000, 5000 ) );
    //CPPUNIT_ASSERT_MESSAGE( "Failed to fit gaussian", status & CGaussianFit::nStatus_Success );

    Float64 gaussAmp;
    Float64 gaussPos;
    Float64 gaussWidth;
    fitter.GetResults( gaussAmp, gaussPos, gaussWidth );


    Float64 gaussAmpErr;
    Float64 gaussPosErr;
    Float64 gaussWidthErr;
    fitter.GetResultsError( gaussAmpErr, gaussPosErr, gaussWidthErr );

    BOOST_CHECK_CLOSE_FRACTION( 1e-17, gaussAmp, 0.01 );
    BOOST_CHECK_CLOSE_FRACTION( 4000.0, gaussPos, 2 );
    BOOST_CHECK_CLOSE_FRACTION( 565.0/1.4142, gaussWidth, 1.0);
}

BOOST_AUTO_TEST_CASE(TestFit1_1em18_noisy)
{
    CSpectrumIOFitsReader reader;
    CSpectrum spectrum;

    Bool retVal = reader.Read( "../test/data/GaussianfitTestCase/gaussian_1em18_noisy.fits", spectrum );
    BOOST_CHECK( retVal == true );


    CGaussianFit fitter;
    CGaussianFit::EStatus status = fitter.Compute( spectrum, TInt32Range( 3000, 5000 ) );
    //CPPUNIT_ASSERT_MESSAGE( "Failed to fit gaussian", status & CGaussianFit::nStatus_Success );

    Float64 gaussAmp;
    Float64 gaussPos;
    Float64 gaussWidth;
    fitter.GetResults( gaussAmp, gaussPos, gaussWidth );


    Float64 gaussAmpErr;
    Float64 gaussPosErr;
    Float64 gaussWidthErr;
    fitter.GetResultsError( gaussAmpErr, gaussPosErr, gaussWidthErr );

    BOOST_CHECK_CLOSE_FRACTION( 1e-18, gaussAmp, 0.2 );
    BOOST_CHECK_CLOSE_FRACTION( 4000.0, gaussPos, 2 );
    BOOST_CHECK_CLOSE_FRACTION( 565.0/1.4142, gaussWidth, 1.0);
}


BOOST_AUTO_TEST_CASE(TestFit2)
{
    CSpectrumIOFitsReader reader;
    CSpectrum spectrum;

    Bool retVal = reader.Read( "../test/data/GaussianfitTestCase/flat.fits", spectrum );
    BOOST_CHECK( retVal == true );


    CGaussianFit fitter;
    CGaussianFit::EStatus status = fitter.Compute( spectrum, TInt32Range( 3000, 5000 ) );
    //BOOST_CHECK( status == CGaussianFit::nStatus_IterationHasNotConverged );

    Float64 gaussAmp;
    Float64 gaussPos;
    Float64 gaussWidth;
    fitter.GetResults( gaussAmp, gaussPos, gaussWidth );

    BOOST_CHECK_CLOSE_FRACTION(   0.0, gaussAmp, 0.01 );

}

BOOST_AUTO_TEST_SUITE_END()


