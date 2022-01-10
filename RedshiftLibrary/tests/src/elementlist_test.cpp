// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#include "RedshiftLibrary/linemodel/linemodelfitting.h"
#include "RedshiftLibrary/noise/flat.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/log/consolehandler.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/tests/test-tools.h"

#include <time.h>
#include <iostream>
#include <stdlib.h>
#include <boost/test/unit_test.hpp>
#include "test-config.h"
#include "RedshiftLibrary/spectrum/LSFFactory.h"
#include "RedshiftLibrary/spectrum/LSF.h"
#include "RedshiftLibrary/spectrum/LSF_NISPSIM_2016.h"
#include "RedshiftLibrary/spectrum/LSF_NISPVSSPSF_201707.h"
#include "RedshiftLibrary/spectrum/LSFConstantResolution.h"
#include "RedshiftLibrary/spectrum/LSFConstantWidth.h"

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
  Float64 initVelocity = 50.0;
  Float64 opt_velocityEmission = initVelocity;
  Float64 opt_velocityAbsorption = initVelocity;
  string opt_rules = "no";
  string opt_rigidity = "rules";

  CSpectrum spectrum;

  CTemplateCatalog tplCatalog;
  TStringList tplCategories;
  CRayCatalog lineCatalog;
  TFloat64Range range(12500,18500);
  CLineModelSolution solution;
  CContinuumModelSolution c_solution;
  int iterations = 1;

  generate_spectrum(spectrum, 1000, 3500, 12500);
  /*  noisePath = generate_noise_fits(1000, 3500, 12500);
  noise.SetNoiseFilePath(noisePath.c_str(), reader);
  noise.AddNoise(spectrum);
  */
  
  std::string lsfType="GaussianConstantWidth";
  Float64 width = 13.;
  TScopeStack scopeStack;
  std::shared_ptr<CParameterStore> store = std::make_shared<CParameterStore>(scopeStack);
  store->Set( "LSF.width", width );
  std::shared_ptr<TLSFArguments> args = std::make_shared<TLSFGaussianConstantWidthArgs>(store);
  std::shared_ptr<CLSF> lsf = LSFFactory.Create(lsfType, args);
  spectrum.SetLSF(lsf);

  generate_template_catalog(tplCatalog, 100, 3500., 12500.);

  calibrationPath = generate_calibration_dir();

  linecatalogPath = generate_linecatalog_file(FULL);
  lineCatalog.Load( linecatalogPath.c_str() );
  CRayCatalog::TRayVector lineList = lineCatalog.GetFilteredList(lineTypeFilter, forceFilter);

  // no continuum
  CLineModelFitting model_nocontinuum(    spectrum,
                                          range,
                                          tplCatalog, tplCategories,
					                                calibrationPath.c_str(), lineList,
					                                "lmfit", "nocontinuum",-INFINITY,
                                          opt_lineWidthType, opt_nsigmasupport, opt_velocityEmission,
					                                opt_velocityAbsorption, opt_rules, opt_rigidity);

  // continuum from spectrum
  CLineModelFitting model_fromspectrum(    spectrum,
                                           range,
                                           tplCatalog, tplCategories,
					                                 calibrationPath.c_str(), lineList,
					                                 "lmfit", "fromspectrum", -INFINITY,
                                           opt_lineWidthType, opt_nsigmasupport, opt_velocityEmission,
					                                 opt_velocityAbsorption, opt_rules, opt_rigidity);

  model_fromspectrum.fit(0.5, range, solution, c_solution, iterations, false);

  // tplfit
  BOOST_TEST_MESSAGE("TODO : tplfit doesn't work. Bad Meiksin generation ?");
  CLineModelFitting model_tplfit(    spectrum,
                                     range,
                                     tplCatalog, tplCategories,
   				                           calibrationPath.c_str(), lineList,
   				                           "lmfit", "tplfit", -5.0,
                                     opt_lineWidthType, opt_nsigmasupport, opt_velocityEmission,
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
