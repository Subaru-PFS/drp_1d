#include <boost/test/unit_test.hpp>
#include <RedshiftLibrary/noise/flat.h>
#include <RedshiftLibrary/spectrum/axis.h>
#include <RedshiftLibrary/spectrum/spectrum.h>

using namespace NSEpic;
using namespace std;



BOOST_AUTO_TEST_SUITE(CNoiseFlat_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(AddNoise_test)
{
  CNoiseFlat ONoiseFlat= CNoiseFlat();
  CSpectrum OSpectrum = CSpectrum();
  CSpectrumFluxAxis& fluxAxis = OSpectrum.GetFluxAxis();
  fluxAxis.SetSize(3);
  fluxAxis[0]= 1.0;
  fluxAxis[1]= 1.5;
  fluxAxis[2]= 2.5;
  BOOST_CHECK(ONoiseFlat.AddNoise(OSpectrum) == true);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
