#include <boost/test/unit_test.hpp>
#include <RedshiftLibrary/noise/flat.h>
#include <RedshiftLibrary/noise/fromfile.h>
#include <RedshiftLibrary/spectrum/axis.h>
#include <RedshiftLibrary/spectrum/io/genericreader.h>
#include <RedshiftLibrary/spectrum/io/fitsreader.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include "test-config.h"

using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(CNoiseFromFile_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(AddNoise_test)
{
  CNoiseFromFile noiseFromFile= CNoiseFromFile();
  CSpectrum OSpectrum = CSpectrum();
  CSpectrumFluxAxis& fluxAxis = OSpectrum.GetFluxAxis();
  std::shared_ptr<CSpectrumIOGenericReader> reader = std::shared_ptr<CSpectrumIOGenericReader>( new CSpectrumIOGenericReader() );

  fluxAxis.SetSize(3);
  fluxAxis[0]= 1.0;
  fluxAxis[1]= 1.5;
  fluxAxis[2]= 2.5;
  noiseFromFile.AddNoise(OSpectrum);

  BOOST_CHECK_THROW(noiseFromFile.SetNoiseFilePath("/this/file/should/not/exist", reader), string);
  BOOST_CHECK_NO_THROW(noiseFromFile.SetNoiseFilePath(DATA_ROOT_DIR "SpectrumioTestCase/spectrum1_z_1.2299.fits", reader));

  /*
  CSpectrumIOFitsReader reader;

  Bool retVal = reader.Read( "01", OSpectrum );
  BOOST_CHECK( retVal == true );
  BOOST_CHECK(noiseFromFile.SetNoiseFilePath("02") == true);
  */
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END ()
