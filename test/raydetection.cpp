#include <boost/filesystem.hpp>
#include <epic/core/common/datatypes.h>
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
#include <boost/test/unit_test.hpp>

using namespace NSEpic;
using namespace std;
using namespace boost;
namespace bfs = boost::filesystem;

BOOST_AUTO_TEST_SUITE(RayDetection)


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


// load spectra from the VVDS DEEP and compare results with EZ python EZELFind results
BOOST_AUTO_TEST_CASE(EzValidationTest)
{
     //deactivated, 20150624, due to irregular sampling compatibility implementation (differs from EZ)
    return;

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

        bfs::path inputspectrum = bfs::path( spectraPath ) / inputspectrumStr;
        bfs::path inputnoise = "";
        bfs::path inputrefdetectionresults = "";
        if(inputnoiseStr.size()>0){
            inputrefdetectionresults = bfs::path( refresults_withnoise_Path ) / inputspectrumStr / inputRayDetectionResultsFileStr;
            inputnoise = bfs::path( spectraPath ) / inputnoiseStr;
        }else{
            inputrefdetectionresults = bfs::path( refresults_nonoise_Path ) / inputspectrumStr / inputRayDetectionResultsFileStr;
            //inputrefdetectionresults /= bfs::path( inputspectrumStr );
            //inputrefdetectionresults /= bfs::path( inputRayDetectionResultsFileStr );
            //.append( bfs::path( inputspectrumStr.c_str() ) ).append( inputRayDetectionResultsFileStr.c_str() );;
        }

        // load spectrum
        CSpectrumIOFitsReader reader;
        CSpectrum s;

        Bool retVal = reader.Read( inputspectrum.c_str(), s );
        BOOST_CHECK( retVal == true);

        // load noise
        if(inputnoiseStr.size()>0){
            CNoiseFromFile noise;
            BOOST_CHECK( noise.SetNoiseFilePath( inputnoise.c_str() ) );
            BOOST_CHECK( noise.AddNoise( s ) );
        }

        // detect possible peaks
        Float64 winsize = 250.0;
        Float64 minsize = 3.0;
        Float64 maxsize = 70.0;
        Float64 cut = 5.0;
        Float64 strongcut = 2.0;

        CPeakDetection peakDetection(winsize, cut);
        auto peakDetectionResult = peakDetection.Compute( s, s.GetLambdaRange() );


        // detected rays
        CLineDetection lineDetection(CRay::nType_Emission, cut, strongcut, winsize, minsize, maxsize);
        auto lineDetectionResult = lineDetection.Compute( s, s.GetLambdaRange(), peakDetectionResult->PeakList, peakDetectionResult->EnlargedPeakList);

        //load reference results
        TFloat64List rayPosList = UtilLoadDetectedRayPositions(inputrefdetectionresults.c_str());

        // Check results
        BOOST_CHECK( lineDetectionResult->RayCatalog.GetList().size() == rayPosList.size());
        for(int i=0; i<lineDetectionResult->RayCatalog.GetList().size(); i++){
            Float64 pos = lineDetectionResult->RayCatalog.GetList()[i].GetPosition();
            Float64 posRef = rayPosList[i];
            BOOST_CHECK_CLOSE_FRACTION( pos, posRef, 1e-6);
        }

    }
    return;
}


//
BOOST_AUTO_TEST_CASE(SyntheticValidationTest)
// load synthetic spectra and check if the lines are correctly detected
{
    std::string spectraPath = "../test/data/RayDetectionTestCase/raydetection_simu_7lines.fits";
//    0#1000 : reference line = detected
//    #2000 : negative line = NOT detected
//    #3000 : very thin line = NOT detected
//    #4000 : very large line = NOT detected
//    1#5000 : ok line = detected
//    #6000+6050 : double/deformed gaussian line = should be rejected, NOT detected
//    #7000: too weak, NOT detected
    bfs::path inputspectrum = bfs::path( spectraPath );


    // load spectrum
    CSpectrumIOFitsReader reader;
    CSpectrum s;

    Bool retVal = reader.Read( inputspectrum.c_str(), s );
    BOOST_CHECK( retVal == true);

    // detect possible peaks
    Float64 winsize = 250.0;
    Float64 minsize = 3.0;
    Float64 maxsize = 80.0;
    Float64 cut = 15.0;
    Float64 strongcut = 2.0;

    CPeakDetection peakDetection(winsize, cut, 1, 2.0 , 0.0);
    auto peakDetectionResult = peakDetection.Compute( s, s.GetLambdaRange() );


    // detected rays
    CLineDetection lineDetection(CRay::nType_Emission, cut, strongcut, winsize, minsize, maxsize);
    auto lineDetectionResult = lineDetection.Compute( s, s.GetLambdaRange(), peakDetectionResult->PeakList, peakDetectionResult->EnlargedPeakList);



    // Check results
    Float64 tol = s.GetResolution();
    BOOST_CHECK( lineDetectionResult->RayCatalog.GetList().size() == 2);
    Float64 pos1 = lineDetectionResult->RayCatalog.GetList()[0].GetPosition();
    BOOST_CHECK_CLOSE_FRACTION( pos1, 1000, tol);
    Float64 pos2 = lineDetectionResult->RayCatalog.GetList()[1].GetPosition();
    BOOST_CHECK_CLOSE_FRACTION( pos2, 5000, tol);


    return;
}

BOOST_AUTO_TEST_SUITE_END()
