#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/io/fitsreader.h>

#include <boost/test/unit_test.hpp>
#include "test-config.h"

using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(SpectrumIO)


BOOST_AUTO_TEST_CASE(VVDSReadValidFile)
{
    CSpectrumIOFitsReader reader;

    CSpectrum spectrum;

    BOOST_CHECK_NO_THROW(reader.Read( DATA_ROOT_DIR "SpectrumioTestCase/spectrum1_z_1.2299.fits",
				      spectrum ));
    BOOST_CHECK( spectrum.GetSampleCount() == 11391 );
    BOOST_CHECK_CLOSE_FRACTION( 3800.0, spectrum.GetLambdaRange().GetBegin(), 0.01 );

}

BOOST_AUTO_TEST_CASE(VVDSReadInvalidFile)
{
    CSpectrumIOFitsReader reader;

    CSpectrum spectrum;

    BOOST_CHECK_THROW(reader.Read( DATA_ROOT_DIR "SpectrumioTestCase/invalidspectrum1.fits",
				   spectrum ), std::string);
    //BOOST_CHECK( spectrum.GetSampleCount() == 0 );
    //BOOST_CHECK_CLOSE_FRACTION( 0.0, spectrum.GetLambdaRange().GetBegin(), 0.01 );
    //BOOST_CHECK_CLOSE_FRACTION( 0.0, spectrum.GetResolution(), 0.01  );

}


BOOST_AUTO_TEST_SUITE_END()

