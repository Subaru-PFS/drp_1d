#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"


#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "test-config.h"
#include <math.h>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(Template)

BOOST_AUTO_TEST_CASE(Constructor)
{
  TFloat64List array = {0.,2.,3.,6.};
  CSpectrumSpectralAxis spectralAxis(array,false) ;
  CSpectrumFluxAxis fluxAxis(array);

  BOOST_CHECK_NO_THROW(CTemplate tmpl("name", "category"));
  BOOST_CHECK_NO_THROW(CTemplate tmpl2("name", "category", spectralAxis, fluxAxis));
}


BOOST_AUTO_TEST_CASE(Save)
{
  TFloat64List array = {0.,2.,3.,6.};
  TFloat64List flux = {0.1,0.2,0.3,0.4};
  CSpectrumSpectralAxis spectralAxis(array, false) ;
  CSpectrumSpectralAxis spectralAxisLog(array, true) ;
  CSpectrumFluxAxis fluxAxis(flux);
  CTemplate tmpl("name", "category", spectralAxis, fluxAxis);
  CTemplate tmplLog("name", "category", spectralAxisLog, fluxAxis);

  // linear
  boost::filesystem::path tempfile = boost::filesystem::unique_path("tst_%%%%%%%%%%.txt");
  const char* filename = tempfile.c_str();
  BOOST_CHECK_NO_THROW(tmpl.Save(filename)); 
  boost::filesystem::remove(filename);

  // logscale
  tempfile = boost::filesystem::unique_path("%%%%%%%%%%.txt");
  const char* filenameLog = tempfile.native().c_str();
  BOOST_CHECK_NO_THROW(tmplLog.Save(filenameLog));
  boost::filesystem::remove(filename);

  // bad filename
  BOOST_CHECK(tmpl.Save(NULL) == false);

}

BOOST_AUTO_TEST_SUITE_END()

