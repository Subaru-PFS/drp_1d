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
#include "RedshiftLibrary/ray/catalog.h"
#include "RedshiftLibrary/ray/lineRatioCatalog.h"
#include "RedshiftLibrary/operator/raymatching.h"
#include "RedshiftLibrary/operator/raymatchingresult.h"

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iostream>

#include <math.h>
#include <boost/test/unit_test.hpp>
#include "test-config.h"

using namespace NSEpic;

using namespace std;
using namespace boost;

BOOST_AUTO_TEST_SUITE(Ray)

//


BOOST_AUTO_TEST_CASE(LoadLineRatioCatalog)
{
    CRayCatalog catalog;
    TAsymParams asymP;
    catalog.AddRayFromParams("Halpha",6562.8,"E","S","SYM",asymP,"",1.,"E1",INFINITY,false,0);
    catalog.AddRayFromParams("Hbeta",4861.3,"E","S","SYM",asymP,"",1.,"E1",INFINITY,false,1);
    catalog.AddRayFromParams("Hgamma",4340.4,"E","W","SYM",asymP,"",1.,"E1",INFINITY,false,2);
    catalog.AddRayFromParams("Hdelta",4101.7,"E","W","SYM",asymP,"",1.,"E1",INFINITY,false,3);

    // TODO this test should be moved to python
    //    BOOST_CHECK_NO_THROW(catalog.Load( DATA_ROOT_DIR "RayTestCase/raycatalog_OK1.txt" ));
    //    BOOST_CHECK_THROW(catalog.Load( DATA_ROOT_DIR "RayTestCase/raycatalog_NOK1.txt" ),
    //		      GlobalException);

}

// load a simple EL catalog and test the match with a redshifted version of itself
BOOST_AUTO_TEST_CASE(MatchingTest1)
{
    CRayCatalog restFrameCatalog;
    TAsymParams asymP;
    restFrameCatalog.AddRayFromParams("Halpha",6562.8,"E","S","SYM",asymP,"",1.,"E1",INFINITY,false,0);
    restFrameCatalog.AddRayFromParams("Hbeta",4861.3,"E","S","SYM",asymP,"",1.,"E1",INFINITY,false,1);
    restFrameCatalog.AddRayFromParams("Hgamma",4340.4,"E","W","SYM",asymP,"",1.,"E1",INFINITY,false,2);
    restFrameCatalog.AddRayFromParams("Hdelta",4101.7,"E","W","SYM",asymP,"",1.,"E1",INFINITY,false,3);

    
    CRayCatalog detectedCatalog;
    Float64 shiftLambda = 1.5;
    CRayCatalog::TRayVector cataloglist = restFrameCatalog.GetList();
    CRayCatalog::TRayVector ::iterator it;
    CLineProfile_ptr profilesym{std::unique_ptr<CLineProfileSYM>(new CLineProfileSYM()) };
    for( it = cataloglist.begin(); it != cataloglist.end(); ++it )
    {
        detectedCatalog.Add( CRay( (*it).GetName(), (*it).GetPosition()*shiftLambda, 2, profilesym->Clone(), 2 ) );
    }

    CRayMatching rayMatching;
    TFloat64Range redshiftrange( 0.0, 5.0);
    auto result = rayMatching.Compute(detectedCatalog, restFrameCatalog, redshiftrange, 2, 0.002 );
    BOOST_CHECK( result != NULL );

    Float64 res = result->GetMeanRedshiftSolutionByIndex(0);
    BOOST_CHECK( fabs(res-(shiftLambda-1)) < 0.0001 );

}

BOOST_AUTO_TEST_CASE(BuilLineRatioCatalog)
{
    CRayCatalog restFrameCatalog;
    TAsymParams asymP;
    restFrameCatalog.AddRayFromParams("Halpha",6562.8,"E","S","SYM",asymP,"",1.,"E1",INFINITY,false,0);
    restFrameCatalog.AddRayFromParams("Hbeta",4861.3,"E","S","SYM",asymP,"",1.,"E1",INFINITY,false,1);
    restFrameCatalog.AddRayFromParams("Hgamma",4340.4,"E","W","SYM",asymP,"",1.,"E1",INFINITY,false,2);
    restFrameCatalog.AddRayFromParams("Hdelta",4101.7,"E","W","SYM",asymP,"",1.,"E1",INFINITY,false,3);

    CLineRatioCatalog lrCatalog("H",restFrameCatalog);
    lrCatalog.addVelocity("velA",200);
    lrCatalog.addVelocity("velE",100);
    lrCatalog.setPrior(0.2);
    lrCatalog.setIsmIndex(2);

}


/*

BOOST_AUTO_TEST_CASE(MatchingTest2_EzValidationTest)
// load raydetection results from VVDS DEEP and compare results with EZ python EZELMatch results
{
    //load restframe catalog
    CRayCatalog restFrameCatalog;

    BOOST_CHECK_NO_THROW(restFrameCatalog.Load( DATA_ROOT_DIR "RayTestCase/RayMatchingVVDS/raycatalog.txt" ));

    //load detected lines results
    TFloat64List rayPosList = UtilLoadDetectedRayPositions(DATA_ROOT_DIR "RayTestCase/RayMatchingVVDS/sc_020100776_F02P017_vmM1_red_129_1_atm_clean.fits/detectedRayCatalog.csv");
    BOOST_CHECK( rayPosList.size()>0 );
    CLineProfile_ptr profilesym{std::unique_ptr<CLineProfileSYM>(new CLineProfileSYM()) };
    CRayCatalog detectedCatalog;
    for( int i=0; i<rayPosList.size(); i++)
    {
        char buffer [64];
        sprintf(buffer,"loaded_%d",i);
        detectedCatalog.Add( CRay("", rayPosList[i], 2, profilesym->Clone(), 2 ) );
    }

    CRayMatching rayMatching;
    TFloat64Range redshiftrange( 0.0, 2.0);
    auto result = rayMatching.Compute(detectedCatalog, restFrameCatalog, redshiftrange, 1, 0.002 );
    BOOST_CHECK( result != NULL );

    //Load RayMatching reference results
    TFloat64List zListRef = UtilLoadRayMatchingResults(DATA_ROOT_DIR "RayTestCase/RayMatchingVVDS/sc_020100776_F02P017_vmM1_red_129_1_atm_clean.fits/rayMatching.csv");
    BOOST_CHECK( zListRef.size()>0 );

    // Check number of matching results
    NSEpic::CRayMatchingResult::TSolutionSetList sol = result->GetSolutionsListOverNumber(0);
    BOOST_CHECK( sol.size() == zListRef.size() );

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
