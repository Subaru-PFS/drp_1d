#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>


#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "test-config.h"
#include <math.h>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(Template)

BOOST_AUTO_TEST_CASE(Constructor)
{
  Float64 array[] = {0.,2.,3.,6.};
  CSpectrumSpectralAxis spectralAxis(array, 4, false) ;
  CSpectrumFluxAxis fluxAxis(array, 4);

  BOOST_CHECK_NO_THROW(CTemplate tmpl("name", "category"));
  BOOST_CHECK_NO_THROW(CTemplate tmpl2("name", "category", spectralAxis, fluxAxis));
}


BOOST_AUTO_TEST_CASE(Save)
{
  Float64 array[] = {0.,2.,3.,6.};
  Float64 flux[] = {0.1,0.2,0.3,0.4};
  CSpectrumSpectralAxis spectralAxis(array, 4, false) ;
  CSpectrumSpectralAxis spectralAxisLog(array, 4, true) ;
  CSpectrumFluxAxis fluxAxis(flux, 4);
  CTemplate tmpl("name", "category", spectralAxis, fluxAxis);
  CTemplate tmplLog("name", "category", spectralAxisLog, fluxAxis);

  CSpectrum spectrum;

  // linear
  boost::filesystem::path tempfile = boost::filesystem::unique_path("tst_%%%%%%%%%%.txt");
  const char* filename = tempfile.c_str();
  BOOST_CHECK_NO_THROW(tmpl.Save(filename));
  spectrum.LoadSpectrum(filename, NULL);
  BOOST_CHECK_CLOSE(spectrum.GetSpectralAxis()[0], 0.0, 1e-8);
  BOOST_CHECK_CLOSE(spectrum.GetSpectralAxis()[3], 6.0, 1e-8);
  BOOST_CHECK_CLOSE(spectrum.GetFluxAxis()[0], 0.1, 1e-8);
  BOOST_CHECK_CLOSE(spectrum.GetFluxAxis()[3], 0.4, 1e-8);
  boost::filesystem::remove(filename);

  // logscale
  tempfile = boost::filesystem::unique_path("%%%%%%%%%%.txt");
  const char* filenameLog = tempfile.native().c_str();
  BOOST_CHECK_NO_THROW(tmplLog.Save(filenameLog));
  spectrum.LoadSpectrum(filenameLog, NULL);
  BOOST_CHECK_CLOSE(spectrum.GetSpectralAxis()[0], exp(0.0), 1e-8);
  BOOST_CHECK_CLOSE(spectrum.GetSpectralAxis()[3], exp(6.0), 1e-8);
  BOOST_CHECK_CLOSE(spectrum.GetFluxAxis()[0], 0.1, 1e-8);
  BOOST_CHECK_CLOSE(spectrum.GetFluxAxis()[3], 0.4, 1e-8);
  boost::filesystem::remove(filename);

  // bad filename
  BOOST_CHECK(tmpl.Save(NULL) == false);

}

BOOST_AUTO_TEST_SUITE_END()

