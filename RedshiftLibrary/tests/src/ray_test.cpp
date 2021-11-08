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
NSEpic::TFloat64List UtilLoadDetectedRayPositions( const char* filePath ){
    TFloat64List posList;

    ifstream file;
    file.open( filePath, ifstream::in );
    if( file.rdstate() & ios_base::failbit )
        return posList;

    string line;

    // Read file line by line
    while( getline( file, line ) )
    {
        char_separator<char> sep(" \t");

        // Tokenize each line
        typedef tokenizer< char_separator<char> > ttokenizer;
        ttokenizer tok( line, sep );

        // Check if it's not a comment
        ttokenizer::iterator it = tok.begin();
        if( it != tok.end() && *it != "#" )
        {
            // Parse name
            string name;
            if( it != tok.end() )
            {
                name = *it;
            }
            else
            {
                return posList;
            }

            ++it;
            // Parse position
            double pos = 0.0;
            try
            {
                pos = lexical_cast<double>(*it);
            }
            catch (bad_lexical_cast&)
            {
                pos = 0.0;
                return posList;
            }
            posList.push_back(pos);

        }
    }
    file.close();

    return posList;
}


//
NSEpic::TFloat64List UtilLoadRayMatchingResults( const char* filePath ){
    TFloat64List zList;

    ifstream file;
    file.open( filePath, ifstream::in );
    if( file.rdstate() & ios_base::failbit )
        return zList;

    string line;

    // Read file line by line
    while( getline( file, line ) )
    {
        char_separator<char> sep("\t");

        // Tokenize each line
        typedef tokenizer< char_separator<char> > ttokenizer;
        ttokenizer tok( line, sep );

        // Check if it's not a comment
        double z = -1.0;
        ttokenizer::iterator it = tok.begin();
        std::string str = *it;
        std::string str1 = str.substr(0,1);
        int comment = strcmp(str1.c_str(), "#");
        if(comment != 0){
            if( it != tok.end() )
            {
                while(it != tok.end()){


                    try
                    {
                        z = lexical_cast<double>(*it);
                    }
                    catch (bad_lexical_cast&)
                    {;
                    }
                    ++it;

                }

            }

            zList.push_back(z);
        }
    }

    file.close();
    return zList;
}


BOOST_AUTO_TEST_CASE(LoadCatalog)
{
    CRayCatalog catalog;

    BOOST_CHECK_NO_THROW(catalog.Load( DATA_ROOT_DIR "RayTestCase/raycatalog_OK1.txt" ));
    BOOST_CHECK_THROW(catalog.Load( DATA_ROOT_DIR "RayTestCase/raycatalog_NOK1.txt" ),
		      GlobalException);

}

// load a simple EL catalog and test the match with a redshifted version of itself
BOOST_AUTO_TEST_CASE(MatchingTest1)
{
    CRayCatalog restFrameCatalog;

    BOOST_CHECK_NO_THROW(restFrameCatalog.Load( DATA_ROOT_DIR "RayTestCase/raycatalog_testMatch1.txt" ));

    CRayCatalog detectedCatalog;
    Float64 shiftLambda = 1.5;
    CRayCatalog::TRayVector cataloglist = restFrameCatalog.GetList();
    CRayCatalog::TRayVector ::iterator it;
    std::shared_ptr<CLineProfile> profilesym{std::make_shared<CLineProfileSYM>()};
    for( it = cataloglist.begin(); it != cataloglist.end(); ++it )
    {
        detectedCatalog.Add( CRay( (*it).GetName(), (*it).GetPosition()*shiftLambda, 2, profilesym, 2 ) );
    }

    CRayMatching rayMatching;
    TFloat64Range redshiftrange( 0.0, 5.0);
    auto result = rayMatching.Compute(detectedCatalog, restFrameCatalog, redshiftrange, 2, 0.002 );
    BOOST_CHECK( result != NULL );

    Float64 res = result->GetMeanRedshiftSolutionByIndex(0);
    BOOST_CHECK( fabs(res-(shiftLambda-1)) < 0.0001 );

}

BOOST_AUTO_TEST_CASE(MatchingTest2_EzValidationTest)
// load raydetection results from VVDS DEEP and compare results with EZ python EZELMatch results
{
    //load restframe catalog
    CRayCatalog restFrameCatalog;

    BOOST_CHECK_NO_THROW(restFrameCatalog.Load( DATA_ROOT_DIR "RayTestCase/RayMatchingVVDS/raycatalog.txt" ));

    //load detected lines results
    TFloat64List rayPosList = UtilLoadDetectedRayPositions(DATA_ROOT_DIR "RayTestCase/RayMatchingVVDS/sc_020100776_F02P017_vmM1_red_129_1_atm_clean.fits/detectedRayCatalog.csv");
    BOOST_CHECK( rayPosList.size()>0 );
    std::shared_ptr<CLineProfile> profilesym{std::make_shared<CLineProfileSYM>()};
    CRayCatalog detectedCatalog;
    for( int i=0; i<rayPosList.size(); i++)
    {
        char buffer [64];
        sprintf(buffer,"loaded_%d",i);
        detectedCatalog.Add( CRay("", rayPosList[i], 2, profilesym, 2 ) );
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

BOOST_AUTO_TEST_SUITE_END()
