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
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/line/catalog.h"
#include "RedshiftLibrary/line/lineRatioCatalog.h"
#include "RedshiftLibrary/operator/linematching.h"
#include "RedshiftLibrary/operator/linematchingresult.h"
#include "test-config.h"

#include <boost/lexical_cast.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/tokenizer.hpp>

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

using namespace NSEpic;

using namespace std;
using namespace boost;

BOOST_AUTO_TEST_SUITE(Line)

//

BOOST_AUTO_TEST_CASE(LoadLineRatioCatalog) {
  CLineCatalog catalog;
  TAsymParams asymP;
  catalog.AddLineFromParams("Halpha", 6562.8, "E", "S", "SYM", asymP, "", 1.,
                            "E1", INFINITY, false, 0, "Halpha_,6562.8_E");
  catalog.AddLineFromParams("Hbeta", 4861.3, "E", "S", "SYM", asymP, "", 1.,
                            "E1", INFINITY, false, 1, "Hbeta_4861.3_E");
  catalog.AddLineFromParams("Hgamma", 4340.4, "E", "W", "SYM", asymP, "", 1.,
                            "E1", INFINITY, false, 2, "Hgamma_4340.4_E");
  catalog.AddLineFromParams("Hdelta", 4101.7, "E", "W", "SYM", asymP, "", 1.,
                            "E1", INFINITY, false, 3, "Hdelta_4101.7_E");

  // TODO this test should be moved to python
  //    BOOST_CHECK_NO_THROW(catalog.Load( DATA_ROOT_DIR
  //    "LineTestCase/linecatalog_OK1.txt" )); BOOST_CHECK_THROW(catalog.Load(
  //    DATA_ROOT_DIR "LineTestCase/linecatalog_NOK1.txt" ),
  //		      GlobalException);
}

// load a simple EL catalog and test the match with a redshifted version of
// itself
BOOST_AUTO_TEST_CASE(MatchingTest1) {
  CLineCatalog restFrameCatalog;
  TAsymParams asymP;
  restFrameCatalog.AddLineFromParams("Halpha", 6562.8, "E", "S", "SYM", asymP,
                                     "", 1., "E1", INFINITY, false, 0,
                                     "Halpha_6562.8_E");
  restFrameCatalog.AddLineFromParams("Hbeta", 4861.3, "E", "S", "SYM", asymP,
                                     "", 1., "E1", INFINITY, false, 1,
                                     "Hbeta_4861.3_E");
  restFrameCatalog.AddLineFromParams("Hgamma", 4340.4, "E", "W", "SYM", asymP,
                                     "", 1., "E1", INFINITY, false, 2,
                                     "Hgamma_4340.4_E");
  restFrameCatalog.AddLineFromParams("Hdelta", 4101.7, "E", "W", "SYM", asymP,
                                     "", 1., "E1", INFINITY, false, 3,
                                     "Hdelta_4101.7_E");

  CLineDetectedCatalog detectedCatalog;
  Float64 shiftLambda = 1.5;
  CLineMap const &cataloglist = restFrameCatalog.GetList();
  CLineMap::iterator it;
  CLineProfile_ptr profilesym = std::make_unique<CLineProfileSYM>();
  for (auto const &[id, line] : cataloglist) {
    detectedCatalog.Add(
        CLineDetected(line.GetName(), line.GetPosition() * shiftLambda,
                      CLine::EType::nType_Emission, profilesym->Clone(),
                      CLine::EForce::nForce_Strong));
  }

  CLineMatching lineMatching;
  TFloat64Range redshiftrange(0.0, 5.0);
  auto result = lineMatching.Compute(detectedCatalog, restFrameCatalog,
                                     redshiftrange, 2, 0.002);
  BOOST_CHECK(result != NULL);

  Float64 res = result->GetMeanRedshiftSolutionByIndex(0);
  BOOST_CHECK(fabs(res - (shiftLambda - 1)) < 0.0001);
}

BOOST_AUTO_TEST_CASE(BuilLineRatioCatalog) {
  CLineCatalog restFrameCatalog;
  TAsymParams asymP;
  restFrameCatalog.AddLineFromParams("Halpha", 6562.8, "E", "S", "SYM", asymP,
                                     "", 1., "E1", INFINITY, false, 0,
                                     "Halpha_6562.8_E");
  restFrameCatalog.AddLineFromParams("Hbeta", 4861.3, "E", "S", "SYM", asymP,
                                     "", 1., "E1", INFINITY, false, 1,
                                     "Hbeta_4861.3_E");
  restFrameCatalog.AddLineFromParams("Hgamma", 4340.4, "E", "W", "SYM", asymP,
                                     "", 1., "E1", INFINITY, false, 2,
                                     "Hgamma_4340.4_E");
  restFrameCatalog.AddLineFromParams("Hdelta", 4101.7, "E", "W", "SYM", asymP,
                                     "", 1., "E1", INFINITY, false, 3,
                                     "Hdelta_4101.7_E");

  CLineRatioCatalog lrCatalog("H", restFrameCatalog);
  lrCatalog.addVelocity("velA", 200);
  lrCatalog.addVelocity("velE", 100);
  lrCatalog.setPrior(0.2);
  lrCatalog.setIsmIndex(2);
}

/*

BOOST_AUTO_TEST_CASE(MatchingTest2_EzValidationTest)
// load linedetection results from VVDS DEEP and compare results with EZ python
EZELMatch results
{
    //load restframe catalog
    CLineCatalog restFrameCatalog;

    BOOST_CHECK_NO_THROW(restFrameCatalog.Load( DATA_ROOT_DIR
"LineTestCase/LineMatchingVVDS/linecatalog.txt" ));

    //load detected lines results
    TFloat64List linePosList = UtilLoadDetectedLinePositions(DATA_ROOT_DIR
"LineTestCase/LineMatchingVVDS/sc_020100776_F02P017_vmM1_red_129_1_atm_clean.fits/detectedLineCatalog.csv");
    BOOST_CHECK( linePosList.size()>0 );
    CLineProfile_ptr profilesym{std::unique_ptr<CLineProfileSYM>(new
CLineProfileSYM()) }; CLineCatalog detectedCatalog; for( int i=0;
i<linePosList.size(); i++)
    {
        char buffer [64];
        sprintf(buffer,"loaded_%d",i);
        detectedCatalog.Add( CLine("", linePosList[i], 2, profilesym->Clone(), 2
) );
    }

    CLineMatching lineMatching;
    TFloat64Range redshiftrange( 0.0, 2.0);
    auto result = lineMatching.Compute(detectedCatalog, restFrameCatalog,
redshiftrange, 1, 0.002 ); BOOST_CHECK( result != NULL );

    //Load LineMatching reference results
    TFloat64List zListRef = UtilLoadLineMatchingResults(DATA_ROOT_DIR
"LineTestCase/LineMatchingVVDS/sc_020100776_F02P017_vmM1_red_129_1_atm_clean.fits/lineMatching.csv");
    BOOST_CHECK( zListRef.size()>0 );

    // Check number of matching results
    NSEpic::CLineMatchingResult::TSolutionSetList sol =
result->GetSolutionsListOverNumber(0); BOOST_CHECK( sol.size() ==
zListRef.size() );

    // Check that all the redshifts values are present
    for( int i=0; i<zListRef.size(); i++)
    {
        Float64 zref= zListRef[i];
        bool found = 0;
        for( int j=0; j<sol.size(); j++)
        {
            Float64 zsol= result->GetMeanRedshiftSolution(sol[i]);
            if(abs(zsol - zref)<1e-10){
                found = true;
                break;
            }
        }
        BOOST_CHECK( found );
    }
}
 */
BOOST_AUTO_TEST_SUITE_END()
