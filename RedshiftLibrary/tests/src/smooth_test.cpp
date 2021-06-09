#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/fluxaxis.h>
#include <RedshiftLibrary/spectrum/spectrum.h>

#include <boost/test/unit_test.hpp>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(Smooth)


BOOST_AUTO_TEST_CASE(Mean)
{
    CSpectrumFluxAxis flux(3);
    CSpectrum spectrum;

    Float64* data = flux.GetSamples();
    data[0] = 10;
    data[1] = 5;
    data[2] = 10;

    flux.ApplyMeanSmooth( 1 );

    BOOST_CHECK_CLOSE_FRACTION( 7.5, flux[0], 0.0001 );
    BOOST_CHECK_CLOSE_FRACTION( 8.3333, flux[1], 0.0001 );
    BOOST_CHECK_CLOSE_FRACTION( 7.5, flux[2], 0.0001 );
}


BOOST_AUTO_TEST_SUITE_END()
