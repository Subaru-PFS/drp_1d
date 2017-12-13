#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/io/fitsreader.h>

#include <boost/test/unit_test.hpp>

using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(SpectrumIO)


BOOST_AUTO_TEST_CASE(VVDSReadValidFile)
{
    CSpectrumIOFitsReader reader;

    CSpectrum s;

    Bool retVal = reader.Read( "./RedshiftLibrary/tests/src/data/SpectrumioTestCase/spectrum1_z_1.2299.fits", s );
    BOOST_CHECK( retVal == true );

    BOOST_CHECK( s.GetSampleCount() == 11391 );
    BOOST_CHECK_CLOSE_FRACTION( 3800.0, s.GetLambdaRange().GetBegin(), 0.01 );

}

BOOST_AUTO_TEST_CASE(VVDSReadInvalidFile)
{
    CSpectrumIOFitsReader reader;

    CSpectrum s;

    Bool rValue = reader.Read( "./RedshiftLibrary/tests/src/data/SpectrumioTestCase/invalidspectrum1.fits", s );

    BOOST_CHECK( rValue == false );
    BOOST_CHECK( s.GetSampleCount() == 0 );
    BOOST_CHECK_CLOSE_FRACTION( 0.0, s.GetLambdaRange().GetBegin(), 0.01 );
    BOOST_CHECK_CLOSE_FRACTION( 0.0, s.GetResolution(), 0.01  );

}


BOOST_AUTO_TEST_SUITE_END()

