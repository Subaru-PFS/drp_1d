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
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/linemodel/modelfittingresult.h"
#include "RedshiftLibrary/processflow/resultstore.h"
#include "RedshiftLibrary/processflow/parameterstore.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(ModelFittingResult)

BOOST_AUTO_TEST_CASE(Constructor)
{
  CLineModelSolution lineModelSolution;
  std::shared_ptr<CLineProfile> profilesym{std::make_shared<CLineProfileSYM>()};
  CRay ray1 = CRay("Abs",5500, 1, profilesym, 2, 10.2, 10.3, 10.4 ,10.5 , 10.6, 10.7, "group", 10.8);
  CRay ray2 = CRay("Abs",5500, 2, profilesym, 2, 10.2, 10.3, 10.4 ,10.5 , 10.6, 10.7, "group", 10.8);
  CRay ray3 = CRay("Abs",5500, 2, profilesym, 1, 10.2, 10.3, 10.4 ,10.5 , 10.6, 10.7, "group", 10.8);
  CRay ray4 = CRay("Em",5520, 2, profilesym, 2, 10.2, 10.3, 20.4 ,10.5 , 10.6, 10.7, "group", 10.8);
  CRayCatalog::TRayVector _restRayList = {ray1, ray2, ray3, ray4};

  TScopeStack scopeStack;
  COperatorResultStore resultStore = COperatorResultStore(scopeStack);
  CParameterStore parameStore = CParameterStore(scopeStack);
  

  boost::filesystem::path temp = boost::filesystem::unique_path();

  ofstream stream(temp.c_str());

  CModelFittingResult result_empty = CModelFittingResult();

  CLineModelResult linemodelResult = CLineModelResult();
  linemodelResult.Redshifts.push_back(0.6);
  linemodelResult.Redshifts.push_back(0.8);

  lineModelSolution.ElementId.push_back(0);
  lineModelSolution.ElementId.push_back(1);
  lineModelSolution.ElementId.push_back(2);
  lineModelSolution.ElementId.push_back(1);
  lineModelSolution.Amplitudes.push_back(2.1);
  lineModelSolution.Amplitudes.push_back(2.1);
  lineModelSolution.Amplitudes.push_back(2.1);
  lineModelSolution.Amplitudes.push_back(1.2);
  lineModelSolution.LambdaObs.push_back(2.1);
  lineModelSolution.LambdaObs.push_back(2.1);
  lineModelSolution.LambdaObs.push_back(2.1);
  lineModelSolution.LambdaObs.push_back(1.2);
  lineModelSolution.Rays.push_back(ray1);
  lineModelSolution.Rays.push_back(ray2);
  lineModelSolution.Rays.push_back(ray3);
  lineModelSolution.Rays.push_back(ray4);
  lineModelSolution.Errors.push_back(0.5);
  lineModelSolution.Errors.push_back(0.5);
  lineModelSolution.Errors.push_back(0.5);
  lineModelSolution.Errors.push_back(0.3);
  lineModelSolution.FittingError.push_back(1.5);
  lineModelSolution.FittingError.push_back(1.5);
  lineModelSolution.FittingError.push_back(1.5);
  lineModelSolution.FittingError.push_back(2.5);
  lineModelSolution.fittingGroupInfo.push_back("foo");
  lineModelSolution.fittingGroupInfo.push_back("bar");
  lineModelSolution.fittingGroupInfo.push_back("baz");
  lineModelSolution.fittingGroupInfo.push_back("quux");
  lineModelSolution.Velocity.push_back(0.5);
  lineModelSolution.Velocity.push_back(0.5);
  lineModelSolution.Velocity.push_back(0.5);
  lineModelSolution.Velocity.push_back(0.3);
  lineModelSolution.Offset.push_back(0.001);
  lineModelSolution.Offset.push_back(0.001);
  lineModelSolution.Offset.push_back(0.002);
  lineModelSolution.Offset.push_back(0.003);
  lineModelSolution.Sigmas.push_back(0.01);
  lineModelSolution.Sigmas.push_back(0.01);
  lineModelSolution.Sigmas.push_back(0.02);
  lineModelSolution.Sigmas.push_back(0.03);
  lineModelSolution.Fluxs.push_back(0.1);
  lineModelSolution.Fluxs.push_back(0.1);
  lineModelSolution.Fluxs.push_back(0.2);
  lineModelSolution.Fluxs.push_back(0.3);
  lineModelSolution.FluxErrors.push_back(0.001);
  lineModelSolution.FluxErrors.push_back(0.001);
  lineModelSolution.FluxErrors.push_back(0.002);
  lineModelSolution.FluxErrors.push_back(0.003);
  lineModelSolution.FluxDirectIntegration.push_back(0.0001);
  lineModelSolution.FluxDirectIntegration.push_back(0.0001);
  lineModelSolution.FluxDirectIntegration.push_back(0.0002);
  lineModelSolution.FluxDirectIntegration.push_back(0.0003);
  lineModelSolution.CenterContinuumFlux.push_back(0.1);
  lineModelSolution.CenterContinuumFlux.push_back(0.1);
  lineModelSolution.CenterContinuumFlux.push_back(0.2);
  lineModelSolution.CenterContinuumFlux.push_back(0.3);
  lineModelSolution.ContinuumError.push_back(0.1);
  lineModelSolution.ContinuumError.push_back(0.1);
  lineModelSolution.ContinuumError.push_back(0.2);
  lineModelSolution.ContinuumError.push_back(0.3);

  linemodelResult.LineModelSolutions.push_back(lineModelSolution);

  CModelFittingResult result = CModelFittingResult(lineModelSolution, 0.5, 1.2,
						   _restRayList, -1.0, -1.0 );

  CLineModelSolution solution = result.GetLineModelSolution();

  //BOOST_CHECK_CLOSE( result_loaded.Redshift, 0.5, 1e-18 );
  //BOOST_CHECK_CLOSE( result_loaded.VelocityEmission, -1.0, 1e-18 );
  //BOOST_CHECK_CLOSE( result_loaded.VelocityAbsorption, -1.0, 1e-18 );
  //BOOST_CHECK_CLOSE( result_loaded.Merit, 1.2, 1e-18 );
  BOOST_CHECK( solution.Amplitudes.size() == 4 );

  boost::filesystem::remove(temp);
}

BOOST_AUTO_TEST_SUITE_END()
