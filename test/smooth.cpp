#include <epic/core/common/datatypes.h>
#include <epic/redshift/spectrum/fluxaxis.h>
#include <epic/redshift/spectrum/spectrum.h>

#include <boost/test/unit_test.hpp>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(Smooth)


BOOST_AUTO_TEST_CASE(Mean)
{
    CSpectrum spectrum;

    spectrum.GetFluxAxis().SetSize( 3 );
    Float64* data = spectrum.GetFluxAxis().GetSamples();
    data[0] = 10;
    data[1] = 5;
    data[2] = 10;

    spectrum.GetFluxAxis().ApplyMeanSmooth( 1 );

    BOOST_CHECK_CLOSE_FRACTION( 7.5, spectrum.GetFluxAxis()[0], 0.0001 );
    BOOST_CHECK_CLOSE_FRACTION( 8.3333, spectrum.GetFluxAxis()[1], 0.0001 );
    BOOST_CHECK_CLOSE_FRACTION( 7.5, spectrum.GetFluxAxis()[2], 0.0001 );
}


BOOST_AUTO_TEST_SUITE_END()
