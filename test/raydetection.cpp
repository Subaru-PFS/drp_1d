#include "raydetection.h"

#include <boost/filesystem.hpp>
#include <epic/core/common/datatypes.h>
#include <epic/core/common/ref.h>
#include <epic/redshift/ray/catalog.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iostream>

#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/spectrum/io/fitsreader.h>
#include <epic/redshift/noise/fromfile.h>

#include <epic/redshift/operator/peakdetection.h>
#include <epic/redshift/operator/peakdetectionresult.h>
#include <epic/redshift/operator/raydetectionresult.h>
#include <epic/redshift/operator/raydetection.h>


#include <math.h>

using namespace NSEpic;
using namespace std;
using namespace NSEpicTest;
namespace bfs = boost::filesystem;
using namespace boost;

void CRedshiftRayDetectionTestCase::setUp()
{
}

void CRedshiftRayDetectionTestCase::tearDown()
{
}

void CRedshiftRayDetectionTestCase::EzValidationTest()
// load spectra from the VVDS DEEP and compare results with EZ python EZELFind results
{
    UInt32 nSpectraToBeTested = 5;
    std::string spectraPath = "../test/data/RayDetectionTestCase/fromVVDSDeep/spectra/";
    std::string refresults_nonoise_Path = "../test/data/RayDetectionTestCase/fromVVDSDeep/results_nonoise/";
    std::string refresults_withnoise_Path = "../test/data/RayDetectionTestCase/fromVVDSDeep/results_withnoise/";


    for(int kspectrum=0; kspectrum<nSpectraToBeTested; kspectrum++){

        // define test
        std::string inputspectrumStr = "";
        std::string inputnoiseStr = "";
        std::string inputRayDetectionResultsFileStr = "detectedRayCatalog.csv";
        if(kspectrum == 0){
            // ...
            inputspectrumStr.append( "sc_020086397_F02P016_vmM1_red_31_1_atm_clean.fits" );
            inputnoiseStr.append("");
        }else if(kspectrum == 1){
            // spectrum for testing retest function
            inputspectrumStr.append( "sc_020100776_F02P017_vmM1_red_129_1_atm_clean.fits" );
            inputnoiseStr.append("");
        }else if(kspectrum == 2){
            // ...
            inputspectrumStr.append( "sc_020123432_F02P019_vmM1_red_72_1_atm_clean.fits" );
            inputnoiseStr.append("");
        }else if(kspectrum == 3){
            // same as previous with noise
            inputspectrumStr.append( "sc_020123432_F02P019_vmM1_red_72_1_atm_clean.fits" );
            inputnoiseStr.append("sc_020123432_F02P019_vmM1_red_72_1_noise.fits");
        }else if(kspectrum == 4){
            // ...
            inputspectrumStr.append( "sc_020088501_F02P017_vmM1_red_82_1_atm_clean.fits" );
            inputnoiseStr.append("sc_020088501_F02P017_vmM1_red_82_1_noise.fits");
        }

        bfs::path inputspectrum = bfs::path( spectraPath ).append( inputspectrumStr.c_str() );
        bfs::path inputnoise = "";
        bfs::path inputrefdetectionresults = "";
        if(inputnoiseStr.size()>0){
            inputrefdetectionresults = bfs::path( refresults_withnoise_Path ).append( inputspectrumStr.c_str() ).append( inputRayDetectionResultsFileStr.c_str() );
            inputnoise = bfs::path( spectraPath ).append( inputnoiseStr.c_str() );
        }else{
            inputrefdetectionresults = bfs::path( refresults_nonoise_Path ).append( inputspectrumStr.c_str() ).append( inputRayDetectionResultsFileStr.c_str() );;
        }

        // load spectrum
        CSpectrumIOFitsReader reader;
        CSpectrum s;

        Bool retVal = reader.Read( inputspectrum.c_str(), s );
        CPPUNIT_ASSERT_MESSAGE(  "load fits", retVal == true);

        // load noise
        if(inputnoiseStr.size()>0){
            CNoiseFromFile noise;
            CPPUNIT_ASSERT_MESSAGE(  "load noise",  noise.SetNoiseFilePath( inputnoise.c_str() ) );
            CPPUNIT_ASSERT_MESSAGE(  "add noise",  noise.AddNoise( s ) );
        }

        // detect possible peaks
        Float64 winsize = 250.0;
        Float64 minsize = 3.0;
        Float64 maxsize = 70.0;
        Float64 cut = 5.0;
        Float64 strongcut = 2.0;

        CPeakDetection peakDetection(winsize, cut);
        CConstRef<CPeakDetectionResult> peakDetectionResult = peakDetection.Compute( s, s.GetLambdaRange() );


        // detected rays
        CRayDetection rayDetection(cut, strongcut, winsize, minsize, maxsize);
        CConstRef<CRayDetectionResult> rayDetectionResult = rayDetection.Compute( s, s.GetLambdaRange(), peakDetectionResult->PeakList, peakDetectionResult->EnlargedPeakList);

        //load reference results
        TFloat64List rayPosList = LoadDetectedRayPositions(inputrefdetectionresults.c_str());

        // Check results
        CPPUNIT_ASSERT_MESSAGE(  "N Ray Detected", rayDetectionResult->RayCatalog.GetList().size() == rayPosList.size());
        for(int i=0; i<rayDetectionResult->RayCatalog.GetList().size(); i++){
            Float64 pos = rayDetectionResult->RayCatalog.GetList()[i].GetPosition();
            Float64 posRef = rayPosList[i];
            CPPUNIT_ASSERT_DOUBLES_EQUAL( pos, posRef, 1e-6);
        }

    }
    return;
}


//
NSEpic::TFloat64List CRedshiftRayDetectionTestCase::LoadDetectedRayPositions( const char* filePath ){
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

void CRedshiftRayDetectionTestCase::SyntheticValidationTest()
// load synthetic spectra and check if the ray are correctly detected
{
    std::string spectraPath = "../test/data/RayDetectionTestCase/raydetection_simu_3700A40FWHM_7000A100FWHM.fits";
    bfs::path inputspectrum = bfs::path( spectraPath );


    // load spectrum
    CSpectrumIOFitsReader reader;
    CSpectrum s;

    Bool retVal = reader.Read( inputspectrum.c_str(), s );
    CPPUNIT_ASSERT_MESSAGE(  "load fits", retVal == true);

    // detect possible peaks
    Float64 winsize = 350.0;
    Float64 minsize = 3.0;
    Float64 maxsize = 200.0;
    Float64 cut = 4.0;
    Float64 strongcut = 2.0;

    CPeakDetection peakDetection(winsize, cut);
    CConstRef<CPeakDetectionResult> peakDetectionResult = peakDetection.Compute( s, s.GetLambdaRange() );


    // detected rays
    CRayDetection rayDetection(cut, strongcut, winsize, minsize, maxsize);
    CConstRef<CRayDetectionResult> rayDetectionResult = rayDetection.Compute( s, s.GetLambdaRange(), peakDetectionResult->PeakList, peakDetectionResult->EnlargedPeakList);



    // Check results
    CPPUNIT_ASSERT_MESSAGE(  "1 Ray Detected", rayDetectionResult->RayCatalog.GetList().size() == 2);
    Float64 pos1 = rayDetectionResult->RayCatalog.GetList()[0].GetPosition();
    CPPUNIT_ASSERT_DOUBLES_EQUAL( pos1, 3702, 1e-1);
    Float64 pos2 = rayDetectionResult->RayCatalog.GetList()[1].GetPosition();
    CPPUNIT_ASSERT_DOUBLES_EQUAL( pos2, 7002, 1e-1);


    return;
}
