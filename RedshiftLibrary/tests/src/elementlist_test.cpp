#include <RedshiftLibrary/linemodel/elementlist.h>
#include <RedshiftLibrary/noise/flat.h>
#include <RedshiftLibrary/noise/fromfile.h>
#include <RedshiftLibrary/spectrum/io/genericreader.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/log/consolehandler.h>
#include <RedshiftLibrary/tests/test-tools.h>

#include <time.h>
#include <iostream>
#include <stdlib.h>
#include <boost/test/unit_test.hpp>
#include "test-config.h"

namespace bfs = boost::filesystem;
using namespace NSEpic;
using namespace std;
using namespace CPFTest;

BOOST_AUTO_TEST_SUITE(test_elementlist)

BOOST_AUTO_TEST_CASE(Constructor)
{
  bfs::path noisePath;
  bfs::path linecatalogPath;
  bfs::path calibrationPath;
  string unused_calibrationPath = "";
  Int32 lineTypeFilter = CRay::nType_Emission;
  Int32 forceFilter = CRay::nForce_Strong;
  string opt_lineWidthType = "velocitydriven";
  Float64 opt_nsigmasupport = 8.;
  Float64 opt_resolution = 2350; //unused with velocity driven linewidth
  Float64 initVelocity = 50.0;
  Float64 opt_velocityEmission = initVelocity;
  Float64 opt_velocityAbsorption = initVelocity;
  string opt_rules = "no";
  string opt_rigidity = "rules";

  CSpectrum spectrum;
  CNoiseFromFile noise;
  CTemplateCatalog tplCatalog;
  TStringList tplCategories;
  CRayCatalog lineCatalog;
  CSpectrumIOGenericReader reader;
  TFloat64Range range(12500,18500);
  CLineModelSolution solution;
  CContinuumModelSolution c_solution;
  int iterations = 1;

  generate_spectrum(spectrum, 1000, 3500, 12500);
  noisePath = generate_noise_fits(1000, 3500, 12500);
  noise.SetNoiseFilePath(noisePath.c_str(), reader);
  noise.AddNoise(spectrum);

  CSpectrum spectrumContinuum = spectrum;
  CSpectrumFluxAxis& continuumFluxAxis = spectrumContinuum.GetFluxAxis();
  for(UInt32 i=0; i<continuumFluxAxis.GetSamplesCount(); i++){
    continuumFluxAxis[i] = 0.0;
  }

  generate_template_catalog(tplCatalog, 100, 3500., 12500.);

  calibrationPath = generate_calibration_dir();

  linecatalogPath = generate_linecatalog_file(FULL);
  lineCatalog.Load( linecatalogPath.c_str() );
  CRayCatalog::TRayVector lineList = lineCatalog.GetFilteredList(lineTypeFilter, forceFilter);

  // no continuum
  CLineModelElementList model_nocontinuum(spectrum, spectrumContinuum,
                                          tplCatalog, tplCatalog, tplCategories,
					  calibrationPath.c_str(), lineList,
					  "lmfit", "nocontinuum",
                                          opt_lineWidthType, opt_nsigmasupport, opt_resolution, opt_velocityEmission,
					  opt_velocityAbsorption, opt_rules, opt_rigidity);

  // continuum from spectrum
  CLineModelElementList model_fromspectrum(spectrum, spectrumContinuum,
                                           tplCatalog, tplCatalog, tplCategories,
					   calibrationPath.c_str(), lineList,
					   "lmfit", "fromspectrum",
                                           opt_lineWidthType, opt_nsigmasupport, opt_resolution, opt_velocityEmission,
					   opt_velocityAbsorption, opt_rules, opt_rigidity);

  model_fromspectrum.fit(0.5, range, solution, c_solution, iterations, false);

  // tplfit
  BOOST_TEST_MESSAGE("TODO : tplfit doesn't work. Bad Meiksin generation ?");
  CLineModelElementList model_tplfit(spectrum, spectrumContinuum,
                                     tplCatalog, tplCatalog, tplCategories,
   				     calibrationPath.c_str(), lineList,
   				     "lmfit", "tplfit",
                                     opt_lineWidthType, opt_nsigmasupport, opt_resolution, opt_velocityEmission,
   				     opt_velocityAbsorption, opt_rules, opt_rigidity);

  /*
  model_tplfit.fit(0.5, range, solution, iterations, false);
  */

  // GetSpectrumModelContinuum
  CSpectrum continuum = model_tplfit.GetSpectrumModelContinuum();

  // setPassMode
  model_tplfit.setPassMode(1);
  //BOOST_CHECK( model_tplfit.m_forceDisableLyaFitting == true );
  model_tplfit.setPassMode(2);
  //BOOST_CHECK( model_tplfit.m_forceDisableLyaFitting == false );

  // GetModelSpectrum
  CSpectrum modelspectrum = model_tplfit.GetModelSpectrum();

  // GetObservedSpectrumWithLinesRemoved
  CSpectrum emission = model_tplfit.GetObservedSpectrumWithLinesRemoved(CRay::nType_Emission);
  CSpectrum absorption = model_tplfit.GetObservedSpectrumWithLinesRemoved(CRay::nType_Absorption);

  bfs::remove(noisePath);
  bfs::remove(linecatalogPath);
  bfs::remove_all(calibrationPath);
}

BOOST_AUTO_TEST_SUITE_END()
