#include <RedshiftLibrary/linemodel/elementlist.h>
#include <RedshiftLibrary/noise/flat.h>
#include <RedshiftLibrary/noise/fromfile.h>
#include <RedshiftLibrary/spectrum/io/genericreader.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/log/consolehandler.h>

#include <time.h>
#include <iostream>
#include <stdlib.h>
#include <boost/test/unit_test.hpp>
#include "test-config.h"

using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(test_elementlist)

BOOST_AUTO_TEST_CASE(Constructor)
{
  CLog log;
  std::shared_ptr<CLogConsoleHandler> logConsoleHandler = std::shared_ptr<CLogConsoleHandler>( new CLogConsoleHandler( CLog::GetInstance() ) );
  logConsoleHandler->SetLevelMask ( CLog::nLevel_Debug );

  string spectrumPath = DATA_ROOT_DIR "ElementListTestCase/spc_15_56957008_SIR_F.fits";
  string noisePath    = DATA_ROOT_DIR "ElementListTestCase/spc_15_56957008_SIR_ErrF.fits";
  string linecatalogPath = DATA_ROOT_DIR "ElementListTestCase/linecatalogamazedvacuum_C1_noHepsilon.txt";
  //string calibrationPath="";
  string calibrationPath = DATA_ROOT_DIR "ElementListTestCase/calibration";
  Int32 lineTypeFilter = CRay::nType_Emission;
  Int32 forceFilter = CRay::nForce_Strong;
  string opt_lineWidthType = "velocitydriven";
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
  int iterations = 1;

  reader.Read( spectrumPath.c_str(), spectrum);
  noise.SetNoiseFilePath(noisePath.c_str(), reader);
  noise.AddNoise(spectrum);

  CSpectrum spectrumContinuum = spectrum;
  CSpectrumFluxAxis& continuumFluxAxis = spectrumContinuum.GetFluxAxis();
  for(UInt32 i=0; i<continuumFluxAxis.GetSamplesCount(); i++){
    continuumFluxAxis[i] = 0.0;
  }

  tplCatalog.Load( DATA_ROOT_DIR "ElementListTestCase/calibration/templates/BC03_sdss_tremonti21/" );

  lineCatalog.Load( linecatalogPath.c_str() );
  CRayCatalog::TRayVector lineList = lineCatalog.GetFilteredList(lineTypeFilter, forceFilter);
  /*
  // no continuum
  CLineModelElementList model_nocontinuum(spectrum, spectrumContinuum,
					  tplCatalog, tplCategories,
					  unused_calibrationPath, lineList,
					  "lmfit", "nocontinuum",
					  opt_lineWidthType, opt_resolution, opt_velocityEmission,
					  opt_velocityAbsorption, opt_rules, opt_rigidity);

  // continuum from spectrum
  CLineModelElementList model_fromspectrum(spectrum, spectrumContinuum,
					   tplCatalog, tplCategories,
					   unused_calibrationPath, lineList,
					   "lmfit", "fromspectrum",
					   opt_lineWidthType, opt_resolution, opt_velocityEmission,
					   opt_velocityAbsorption, opt_rules, opt_rigidity);

  model_fromspectrum.fit(0.5, range, solution, iterations, false);

  */

  // tplfit
  CLineModelElementList model_tplfit(spectrum, spectrumContinuum,
   				     tplCatalog, tplCategories,
   				     calibrationPath, lineList,
   				     "lmfit", "tplfit",
   				     opt_lineWidthType, opt_resolution, opt_velocityEmission,
   				     opt_velocityAbsorption, opt_rules, opt_rigidity);

  model_tplfit.fit(0.5, range, solution, iterations, false);

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

}

BOOST_AUTO_TEST_SUITE_END()
