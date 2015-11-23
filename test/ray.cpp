#include "ray.h"

#include <epic/core/common/datatypes.h>
#include <epic/redshift/ray/catalog.h>
#include <epic/redshift/operator/raymatching.h>
#include <epic/redshift/operator/raymatchingresult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iostream>

#include <math.h>

using namespace NSEpic;
using namespace NSEpicTest;

using namespace std;
using namespace boost;

void CRedshiftRayTestCase::setUp()
{
}

void CRedshiftRayTestCase::tearDown()
{
}

void CRedshiftRayTestCase::LoadCatalog()
{
    CRayCatalog catalog;
    Bool returnValue;

    returnValue = catalog.Load( "../test/data/RayTestCase/raycatalog_OK1.txt" );
    CPPUNIT_ASSERT_MESSAGE( "Failled to load or parse raycatalog_OK1.txt", returnValue == true );

    returnValue = catalog.Load( "../test/data/RayTestCase/raycatalog_NOK1.txt" );
    CPPUNIT_ASSERT_MESSAGE( "Succeeded to parse invalid raycatalog_NOK1.txt", returnValue == false );

    returnValue = catalog.Load( "../test/data/RayTestCase/raycatalog_NOK1.txt" );
    CPPUNIT_ASSERT_MESSAGE( "Succeeded to parse invalid raycatalog_NOK1.txt", returnValue == false );

}

void CRedshiftRayTestCase::MatchingTest1()
// load a simple EL catalog and test the match with a redshifted version of itself
{
    CRayCatalog restFrameCatalog;
    Bool returnValue;

    returnValue = restFrameCatalog.Load( "../test/data/RayTestCase/raycatalog_testMatch1.txt" );
    CPPUNIT_ASSERT_MESSAGE( "Failed to load or parse raycatalog_testMatch1.txt", returnValue == true );

    CRayCatalog detectedCatalog;
    Float64 shiftLambda = 1.5;
    CRayCatalog::TRayVector cataloglist = restFrameCatalog.GetList();
    CRayCatalog::TRayVector ::iterator it;
    for( it = cataloglist.begin(); it != cataloglist.end(); ++it )
    {
        detectedCatalog.Add( CRay( (*it).GetName(), (*it).GetPosition()*shiftLambda, 2, 2 ) );
    }

    CRayMatching rayMatching;
    TFloat64Range redshiftrange( 0.0, 5.0);
    auto result = rayMatching.Compute(detectedCatalog, restFrameCatalog, redshiftrange, 2, 0.002 );
    CPPUNIT_ASSERT_MESSAGE( "Failed to match line catalogs for MatchingTest1.txt", result != NULL );

    Float64 res = result->GetMeanRedshiftSolutionByIndex(0);
    CPPUNIT_ASSERT_MESSAGE( "Failed to find redshift accurately for MatchingTest1", fabs(res-(shiftLambda-1)) < 0.0001 );

}

void CRedshiftRayTestCase::MatchingTest2_EzValidationTest()
// load raydetection results from VVDS DEEP and compare results with EZ python EZELMatch results
{
    //load restframe catalog
    CRayCatalog restFrameCatalog;
    Bool returnValue;

    returnValue = restFrameCatalog.Load( "../test/data/RayTestCase/RayMatchingVVDS/raycatalog.txt" );
    CPPUNIT_ASSERT_MESSAGE( "Failed to load or parse raycatalog.txt", returnValue == true );

    //load detected lines results
    TFloat64List rayPosList = LoadDetectedRayPositions("../test/data/RayTestCase/RayMatchingVVDS/sc_020100776_F02P017_vmM1_red_129_1_atm_clean.fits/detectedRayCatalog.csv");
    CPPUNIT_ASSERT_MESSAGE( "Failed to load or parse input detected lines results", rayPosList.size()>0 );

    CRayCatalog detectedCatalog;
    for( int i=0; i<rayPosList.size(); i++)
    {
        char buffer [64];
        sprintf(buffer,"loaded_%d",i);
        detectedCatalog.Add( CRay("", rayPosList[i], 2, 2 ) );
    }

    CRayMatching rayMatching;
    TFloat64Range redshiftrange( 0.0, 2.0);
    auto result = rayMatching.Compute(detectedCatalog, restFrameCatalog, redshiftrange, 1, 0.002 );
    CPPUNIT_ASSERT_MESSAGE( "Failed to match line catalogs for MatchingTest1.txt", result != NULL );

    //Load RayMatching reference results
    TFloat64List zListRef = LoadRayMatchingResults("../test/data/RayTestCase/RayMatchingVVDS/sc_020100776_F02P017_vmM1_red_129_1_atm_clean.fits/rayMatching.csv");
    CPPUNIT_ASSERT_MESSAGE( "Failed to load or parse raycatalog.txt", zListRef.size()>0 );

    // Check number of matching results
    NSEpic::CRayMatchingResult::TSolutionSetList sol = result->GetSolutionsListOverNumber(0);
    CPPUNIT_ASSERT_MESSAGE( "Failed in comparison of the number of matching solutions", sol.size() == zListRef.size() );

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
        CPPUNIT_ASSERT_MESSAGE( "Not all Matching solutions redshifts values comparison successfull at the given tolerance !", found );
    }
}

//
NSEpic::TFloat64List CRedshiftRayTestCase::LoadDetectedRayPositions( const char* filePath ){
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
            catch (bad_lexical_cast)
            {
                pos = 0.0;
                return posList;
            }
            posList.push_back(pos);

        }
    }


    return posList;
}


//
NSEpic::TFloat64List CRedshiftRayTestCase::LoadRayMatchingResults( const char* filePath ){
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
                    catch (bad_lexical_cast)
                    {;
                    }
                    ++it;

                }

            }

            zList.push_back(z);
        }
    }


    return zList;
}
