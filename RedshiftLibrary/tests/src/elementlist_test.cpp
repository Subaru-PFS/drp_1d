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
#include <iostream>

#include <boost/test/unit_test.hpp>

#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/linemodel/linemodelfitting.h"
#include "RedshiftLibrary/log/consolehandler.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/spectrum/LSF.h"
#include "RedshiftLibrary/spectrum/LSFConstantResolution.h"
#include "RedshiftLibrary/spectrum/LSFConstantWidth.h"
#include "RedshiftLibrary/spectrum/LSFFactory.h"
#include "RedshiftLibrary/spectrum/LSF_NISPSIM_2016.h"
#include "RedshiftLibrary/spectrum/LSF_NISPVSSPSF_201707.h"
#include "test-config.h"

namespace bfs = boost::filesystem;
using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(test_elementlist)

BOOST_AUTO_TEST_CASE(Constructor) {
  bfs::path noisePath;
  bfs::path linecatalogPath;

  auto lineTypeFilter = CLine::EType::nType_Emission;
  auto forceFilter = CLine::EForce::nForce_Strong;
  string opt_lineWidthType = "velocityDriven";
  Float64 opt_nsigmasupport = 8.;
  Float64 initVelocity = 50.0;
  Float64 opt_velocityEmission = initVelocity;
  Float64 opt_velocityAbsorption = initVelocity;
  string opt_rules = "no";
  string opt_rigidity = "rules";

  CSpectrum spectrum;

  // CTemplateCatalog tplCatalog;
  TStringList tplCategories;
  CLineCatalog lineCatalog;
  TFloat64Range range(12500, 18500);
  CLineModelSolution solution;
  CTplModelSolution c_solution;
  int iterations = 1;

  // generate_spectrum(spectrum, 1000, 3500, 12500);
  /*  noisePath = generate_noise_fits(1000, 3500, 12500);
  noise.SetNoiseFilePath(noisePath.c_str(), reader);
  noise.AddNoise(spectrum);
  */

  std::string lsfType = "gaussianConstantWidth";
  Float64 width = 13.;
std:
  shared_ptr<CScopeStack> scopeStack = std::make_shared<CScopeStack>();
  std::shared_ptr<CParameterStore> store =
      std::make_shared<CParameterStore>(scopeStack);
  store->Set("lsf.width", width);
  std::shared_ptr<TLSFArguments> args =
      std::make_shared<TLSFGaussianConstantWidthArgs>(store);
  std::shared_ptr<CLSF> lsf = LSFFactory.Create(lsfType, args);
  spectrum.SetLSF(lsf);

  // generate_template_catalog(tplCatalog, 100, 3500., 12500.);
  TAsymParams asymP;
  /* //TODO restore this test , for the moment if linecatalog is not empty, test
    will fail with follogin message: fatal error: in
    "test_elementlist/Constructor": NSEpic::GlobalException: Could not find
    template with name fromspectrum

    lineCatalog.AddLineFromParams("Halpha",6562.8,"E","S","SYM",asymP,"",1.,"E1",INFINITY,false,0,"Halpha_,6562.8_E");
    lineCatalog.AddLineFromParams("Hbeta",4861.3,"E","S","SYM",asymP,"",1.,"E1",INFINITY,false,1,"Hbeta_4861.3_E");
    lineCatalog.AddLineFromParams("Hgamma",4340.4,"E","W","SYM",asymP,"",1.,"E1",INFINITY,false,2,"Hgamma_4340.4_E");
    lineCatalog.AddLineFromParams("Hdelta",4101.7,"E","W","SYM",asymP,"",1.,"E1",INFINITY,false,3,"Hdelta_4101.7_E");
  */
  // BUILT FROM REGEX replace :
  //"(\d+\.\d+)\\t([a-z0-9_\[\]\(\)\-\.\/]+)\\t([AE])\\t([WS])\\t([A-Z]+)\\t([A-Z_\-1]+)\\t(-?\d\.?\d*)\\t*-1\\n"
  // lineCatalog.AddLineFromParams("$2",$1,"$3","$4","$5",asymP,"$6",$7,"$3
  // 1",INFINITY,false,);
  /*
  lineCatalog.AddLineFromParams("P5A",12821.59,"A","W","SYM",asymP,undefStr,1,"A1",0,false,0);
  lineCatalog.AddLineFromParams("P6A",10941.09,"A","W","SYM",asymP,undefStr,1,"A1",0,false,1);
  lineCatalog.AddLineFromParams("P7A",10052.13,"A","W","SYM",asymP,undefStr,1,"A1",0,false,2);
  lineCatalog.AddLineFromParams("P8A",9548.59,"A","W","SYM",asymP,undefStr,1,"A1",0,false,3);
  lineCatalog.AddLineFromParams("P9A",9231.55,"A","W","SYM",asymP,undefStr,1,"A1",0,false,4);
  lineCatalog.AddLineFromParams("P10A",9017.38,"A","W","SYM",asymP,undefStr,1,"A1",0,false,5);
  lineCatalog.AddLineFromParams("P11A",8865.22,"A","W","SYM",asymP,undefStr,1,"A1",0,false,6);
  lineCatalog.AddLineFromParams("CaII_t3A",8664.50,"A","S","SYM",asymP,undefStr,1,"A1",0,false,7);
  lineCatalog.AddLineFromParams("CaII_t2A",8544.42,"A","S","SYM",asymP,undefStr,1,"A1",0,false,8);
  lineCatalog.AddLineFromParams("CaII_t1A",8500.35,"A","S","SYM",asymP,undefStr,1,"A1",0,false,9);
  lineCatalog.AddLineFromParams("HalphaA",6564.61,"A","S","SYM",asymP,undefStr,1,"A1",0,false,10);
  lineCatalog.AddLineFromParams("NaD",5895.6,"A","W","SYM",asymP,undefStr,1,"A1",0,false,11);
  lineCatalog.AddLineFromParams("MgI5175",5176.71,"A","S","SYM",asymP,undefStr,1,"A1",0,false,12);
  lineCatalog.AddLineFromParams("HbetaA",4862.72,"A","W","SYM",asymP,undefStr,1,"A1",0,false,13);
  lineCatalog.AddLineFromParams("HgammaA",4341.58,"A","W","SYM",asymP,undefStr,1,"A1",0,false,14);
  lineCatalog.AddLineFromParams("GBand",4304.57,"A","W","SYM",asymP,undefStr,1,"A1",0,false,15);
  lineCatalog.AddLineFromParams("HdeltaA",4102.81,"A","W","SYM",asymP,undefStr,1,"A1",0,false,16);
  lineCatalog.AddLineFromParams("HepsilonA",3971.15,"A","W","SYM",asymP,undefStr,1,"A1",0,false,17);
  lineCatalog.AddLineFromParams("CaII_H",3969.55,"A","S","SYM",asymP,"A_Ca",22,"A1",0,false,18);
  lineCatalog.AddLineFromParams("CaII_K",3934.73,"A","S","SYM",asymP,"A_Ca",23,"A1",0,false,19);
  lineCatalog.AddLineFromParams("H8A",3890.11,"A","W","SYM",asymP,undefStr,1,"A1",0,false,20);
  lineCatalog.AddLineFromParams("H9A",3836.43,"A","W","SYM",asymP,undefStr,1,"A1",0,false,21);
  lineCatalog.AddLineFromParams("H10A",3798.93,"A","W","SYM",asymP,undefStr,1,"A1",0,false,22);
  lineCatalog.AddLineFromParams("H11A",3771.65,"A","W","SYM",asymP,undefStr,1,"A1",0,false,23);
  lineCatalog.AddLineFromParams("FeI",3582.19,"A","W","SYM",asymP,undefStr,1,"A1",0,false,24);
  lineCatalog.AddLineFromParams("MgI2852",2853.73,"A","S","SYM",asymP,undefStr,1,"A1",0,false,25);
  lineCatalog.AddLineFromParams("MgII2803",2804.29,"A","S","SYM",asymP,"A_MgII",0.3054,"A1",0,false,26);
  lineCatalog.AddLineFromParams("MgII2796",2797.11,"A","S","SYM",asymP,"A_MgII",0.6123,"A1",0,false,27);
  lineCatalog.AddLineFromParams("FeII2600",2600.87,"A","S","SYM",asymP,undefStr,1,"A1",0,false,28);
  lineCatalog.AddLineFromParams("FeII2586",2587.35,"A","W","SYM",asymP,undefStr,1,"A1",0,false,29);
  lineCatalog.AddLineFromParams("[SII]6731",6732.68,"E","W","SYM",asymP,undefStr,1,"E1",0,false,30);
  lineCatalog.AddLineFromParams("[SII]6716",6718.29,"E","W","SYM",asymP,undefStr,1,"E1",0,false,31);
  lineCatalog.AddLineFromParams("[NII](doublet-1)",6585.27,"E","W","SYM",asymP,"E_NII",2.95,"E1",0,false,32);
  lineCatalog.AddLineFromParams("Halpha",6564.61,"E","S","SYM",asymP,undefStr,1,"E1",0,false,33);
  lineCatalog.AddLineFromParams("[NII](doublet-1/2.95)",6549.86,"E","W","SYM",asymP,"E_NII",1,"E1",0,false,34);
  lineCatalog.AddLineFromParams("[OIII](doublet-1)",5008.24,"E","S","SYM",asymP,"E_OIII",3,"E1",0,false,35);
  lineCatalog.AddLineFromParams("[OIII](doublet-1/3)",4960.29,"E","W","SYM",asymP,"E_OIII",1,"E1",0,false,36);
  lineCatalog.AddLineFromParams("Hbeta",4862.72,"E","S","SYM",asymP,undefStr,1,"E1",0,false,37);
  lineCatalog.AddLineFromParams("Hgamma",4341.58,"E","W","SYM",asymP,undefStr,1,"E1",0,false,38);
  lineCatalog.AddLineFromParams("Hdelta",4102.81,"E","W","SYM",asymP,undefStr,1,"E1",0,false,39);
  lineCatalog.AddLineFromParams("Hepsilon",3971.15,"E","W","SYM",asymP,undefStr,1,"E1",0,false,40);
  lineCatalog.AddLineFromParams("H8",3890.11,"E","W","SYM",asymP,undefStr,1,"E1",0,false,41);
  lineCatalog.AddLineFromParams("NeIII",3869.05,"E","W","SYM",asymP,undefStr,1,"E1",0,false,42);
  lineCatalog.AddLineFromParams("H9",3836.43,"E","W","SYM",asymP,undefStr,1,"E1",0,false,43);
  lineCatalog.AddLineFromParams("H10",3798.93,"E","W","SYM",asymP,undefStr,1,"E1",0,false,44);
  lineCatalog.AddLineFromParams("H11",3771.65,"E","W","SYM",asymP,undefStr,1,"E1",0,false,45);
  lineCatalog.AddLineFromParams("[OII]3729",3729.88,"E","S","SYM",asymP,undefStr,1,"E1",0,false,46);
  lineCatalog.AddLineFromParams("[OII]3726",3727.09,"E","S","SYM",asymP,undefStr,1,"E1",0,false,47);
  lineCatalog.AddLineFromParams("[NeVa]",3426.73,"E","W","SYM",asymP,undefStr,1,"E1",0,false,48);
  lineCatalog.AddLineFromParams("[NeVb]",3346.81,"E","W","SYM",asymP,undefStr,1,"E1",0,false,49);
  lineCatalog.AddLineFromParams("MgII",2799.12,"E","W","SYM",asymP,undefStr,1,"E1",0,false,50);
  lineCatalog.AddLineFromParams("Fe2632",2632.82,"E","W","SYM",asymP,undefStr,1,"E1",0,false,51);;
  */
  /*
 lineCatalog.AddLineFromParams("CaH",3698.5,"A","S","SYM",asymP,"",1.,"A1",0,false,0);
 lineCatalog.AddLineFromParams("MgI",5175,"A","S","SYM",asymP,"",1.,"A1",0,false,1);
 lineCatalog.AddLineFromParams("[SII]",10320,"E","W","SYM",asymP,"",1.,"E1",0,false,2);
 lineCatalog.AddLineFromParams("Halpha",6562.8,"E","S","SYM",asymP,"",1.,"E1",0,false,3);
 lineCatalog.AddLineFromParams("Hdelta",4101.7,"E","W","SYM",asymP,"",1.,"E1",0,false,4);
  */

  /*
  model_tplfit.fit(0.5, range, solution, iterations, false);
  */

  // GetSpectrumModelContinuum

  // setPassMode

  // BOOST_CHECK( model_tplfit.m_forceDisableLyaFitting == true );

  // BOOST_CHECK( model_tplfit.m_forceDisableLyaFitting == false );

  // GetModelSpectrum

  // GetObservedSpectrumWithLinesRemoved

  bfs::remove(noisePath);
  bfs::remove(linecatalogPath);
}

BOOST_AUTO_TEST_SUITE_END()
