#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/io/fitsreader.h>

#include <boost/filesystem.hpp>

#include <boost/test/unit_test.hpp>
#include "test-tools.h"
#include "test-config.h"

using namespace NSEpic;
using namespace CPFTest;
using namespace std;

BOOST_AUTO_TEST_SUITE(SpectrumIO)

BOOST_AUTO_TEST_CASE(VVDSReadValidFile)
{
    CSpectrumIOFitsReader reader;

    CSpectrum spectrum;

    boost::filesystem::path tempfile = generate_spectrum_fits(123, 3800, 12600);

    BOOST_CHECK_NO_THROW(reader.Read( tempfile.c_str(), spectrum ));
    BOOST_CHECK( spectrum.GetSampleCount() == 123 );
    BOOST_CHECK_CLOSE_FRACTION( 3800.0, spectrum.GetLambdaRange().GetBegin(), 0.01 );

    boost::filesystem::remove(tempfile.native());

}

BOOST_AUTO_TEST_CASE(VVDSReadInvalidFile)
{
    CSpectrumIOFitsReader reader;

    CSpectrum spectrum;

    BOOST_CHECK_THROW(reader.Read( DATA_ROOT_DIR "SpectrumioTestCase/invalidspectrum1.fits",
                                   spectrum ), std::runtime_error);
    //BOOST_CHECK( spectrum.GetSampleCount() == 0 );
    //BOOST_CHECK_CLOSE_FRACTION( 0.0, spectrum.GetLambdaRange().GetBegin(), 0.01 );
    //BOOST_CHECK_CLOSE_FRACTION( 0.0, spectrum.GetResolution(), 0.01  );

}


BOOST_AUTO_TEST_SUITE_END()

