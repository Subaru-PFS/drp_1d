#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/spectrum.h>

#include <RedshiftLibrary/continuum/irregularsamplingmedian.h>

#include <RedshiftLibrary/spectrum/io/fitsreader.h>

#include <RedshiftLibrary/noise/flat.h>
#include <RedshiftLibrary/noise/fromfile.h>
#include <RedshiftLibrary/linemodel/modelfittingresult.h>
#include <RedshiftLibrary/linemodel/elementlist.h>

#include <boost/test/unit_test.hpp>

#include <boost/property_tree/ptree.hpp>
#include <math.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>

#include "test-config.h"


using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(LinemodelMultilineProfile)

/**
 * @brief checkProfileValue
 * @param linecatalogPath
 * @param opt_velocityEmission
 * @param lambdas: a set of wavelengths where the flux will be tested with regard to a reference
 * @param refProfileFluxValues: n reference flux values corresponding to the wavelengths in lambdas
 * @param enableGradientCheck
 * @param refProfileDerivAmplitudeValues: n*m reference gradient values corresponding to each element in the catalog (first d_amps, then d_velocity)
 * @param enableGradientMeansquareCheck
 * @param refProfileMeansquareValue: reference meansquare residual float value
 * @param refProfileGradientMeansquareValues: m reference meansquare gradient values corresponding to each element in the catalog (first d_amps, then d_velocity)
*/
void checkProfileValue(std::string linecatalogPath,
                       Float64 opt_velocityEmission,
                       std::vector<Float64> lambdas,
                       std::vector<Float64> refProfileFluxValues,
                       bool enableGradientCheck,
                       std::vector<std::vector<Float64> >  refProfileGradientValues,
                       bool enableGradientMeansquareCheck,
                       Float64   refProfileMeansquareValue,
                       std::vector<Float64>   refProfileGradientMeansquareValues)
{
    //some input params
    std::string spectrumPath = DATA_ROOT_DIR "LinemodelProfileTestCase/signalnoise_4lines_sig400_6000A_10000A.fits";       //unused, only needed for linemodel initialization
    std::string noisePath    = DATA_ROOT_DIR "LinemodelProfileTestCase/noise_4lines_sig400_6000A_10000A.fits";     //unused, only needed for linemodel initialization
    Int32 lineTypeFilter = CRay::nType_Emission;
    Int32 forceFilter = CRay::nForce_Strong;
    std::string opt_fittingmethod = "ones"; //all the elements amplitudes set to 1.0
    Float64 z = 0.0;

    // load spectrum
    CSpectrumIOFitsReader reader;

    CSpectrum spectrum;

    BOOST_CHECK_NO_THROW(reader.Read( spectrumPath.c_str(), spectrum));
    CNoiseFromFile noise;
    BOOST_CHECK_NO_THROW(noise.SetNoiseFilePath( noisePath.c_str(), reader ));
    BOOST_CHECK_NO_THROW(noise.AddNoise( spectrum ));

    // get continuum
    //CContinuumIrregularSamplingMedian continuum;
    //CSpectrumFluxAxis fluxAxisWithoutContinuumCalc;
    //Int32 retValCont = continuum.RemoveContinuum( spectrum, fluxAxisWithoutContinuumCalc );
    CSpectrum spectrumContinuum = spectrum;
    CSpectrumFluxAxis& continuumFluxAxis = spectrum.GetFluxAxis();
    for(UInt32 i=0; i<continuumFluxAxis.GetSamplesCount(); i++){
        //continuumFluxAxis[i] -= fluxAxisWithoutContinuumCalc[i];
        continuumFluxAxis[i] = 0.0;
    }


    //get line catalog
    CRayCatalog lineCatalog;
    BOOST_CHECK_NO_THROW(lineCatalog.Load( linecatalogPath.c_str() ));

    CRayCatalog::TRayVector lineList = lineCatalog.GetFilteredList(lineTypeFilter, forceFilter);
    BOOST_CHECK( lineList.size()>0);


    std::string opt_continuumcomponent = "fromspectrum";
    std::string opt_lineWidthType = "velocitydriven";

    Float64 opt_resolution = 2350; //unused with velocity driven linewidth
    Float64 opt_velocityAbsorption = 100.0; //unused
    std::string opt_rules = "no";
    std::string opt_rigidity = "rules";
    std::string unused_calibrationPath="";


    //these tplcatalog related variables are unused here.
    CTemplateCatalog tplCatalog;
    BOOST_CHECK_NO_THROW(tplCatalog.Load( DATA_ROOT_DIR "templatecatalog/" ));
    TStringList tplCategories;

    CLineModelElementList model(spectrum, spectrumContinuum, tplCatalog, tplCategories, unused_calibrationPath, lineList, opt_fittingmethod, opt_continuumcomponent, opt_lineWidthType, opt_resolution, opt_velocityEmission, opt_velocityAbsorption, opt_rules, opt_rigidity);
    TFloat64Range lambdaRange = TFloat64Range( 100.0, 12000.0 );
    CLineModelSolution modelSolution;
    model.fit(z, lambdaRange, modelSolution);

    CSpectrumSpectralAxis spcSpectralAxis = spectrum.GetSpectralAxis();

    std::vector<UInt32> validEltsIdx;
    for (UInt32 iElt=0; iElt<model.m_Elements.size(); iElt++)
    {
        validEltsIdx.push_back(iElt);
    }

    model.refreshModelDerivVelUnderElements(validEltsIdx);

    //check the flux values
    for(Int32 k=0; k<lambdas.size(); k++)
    {
        Float64 iLambda = spcSpectralAxis.GetIndexAtWaveLength(lambdas[k]);
        Float64 profileValue = model.getModelFluxVal(iLambda);
        Float64 refFluxValue = refProfileFluxValues[k];
        //BOOST_CHECK_CLOSE_FRACTION( profileValue, refFluxValue, 1e-5);
        BOOST_CHECK_SMALL( std::abs(profileValue - refFluxValue), 1e-3);
    }

    if(enableGradientCheck)
    {
        //check the deriv-amplitude values
        for(Int32 iElt=0; iElt<validEltsIdx.size(); iElt++)
        {
            for(Int32 k=0; k<lambdas.size(); k++)
            {
                Float64 iLambda = spcSpectralAxis.GetIndexAtWaveLength(lambdas[k]);
                Float64 val = model.getModelFluxDerivEltVal(validEltsIdx[iElt], iLambda);
                Float64 refDerivAmpValue = refProfileGradientValues[k][iElt];
                BOOST_CHECK_SMALL( std::abs(val - refDerivAmpValue), 1e-5);
            }

        }

        //check the deriv-velocity values
        Int32 idxVelGradient = validEltsIdx.size();
        for(Int32 k=0; k<lambdas.size(); k++)
        {
            Float64 iLambda = spcSpectralAxis.GetIndexAtWaveLength(lambdas[k]);
            Float64 val = model.getModelFluxDerivVelVal(iLambda);
            Float64 refDerivVelValue = refProfileGradientValues[k][idxVelGradient];
            BOOST_CHECK_SMALL( std::abs(val - refDerivVelValue), 1e-4);
        }
    }


    if(enableGradientMeansquareCheck)
    {

        CSpectrumFluxAxis spcFluxAxis = spectrum.GetFluxAxis();
        //some allocations
        Int32 nsamples = lambdas.size();
        Int32 nddl = validEltsIdx.size()+1; //amps+velocity
        Float64 normFactor = 1.0;
        //input model variables
        Float64* margs = (Float64*) calloc( nddl, sizeof( Float64 ) );
        for (Int32 i = 0; i < nddl-1; i++)
        {
            margs[i] = 1.0;
        }
        margs[nddl-1] = opt_velocityEmission;
        //output variables
        Float64 f=0.0;
        Float64* g = (Float64*) calloc( nddl, sizeof( Float64 ) );
        //buffer for calculation of the gradient
        Float64* mmy = (Float64*) calloc( nsamples, sizeof( Float64 ) );
        //prepare fluxdata
        Float64* fluxdata = (Float64*) calloc( nsamples, sizeof( Float64 ) );
        std::vector<UInt32> xInds;
        const Float64* flux = spcFluxAxis.GetSamples();
        for (Int32 i = 0; i < nsamples; i++)
        {
            xInds.push_back(i); //the support is the full wavelength range for this test
            fluxdata[i] = flux[i]*normFactor;
        }


        model.estimateMeanSqFluxAndGradient(margs, normFactor, validEltsIdx, xInds, lineTypeFilter, fluxdata, mmy, f, g);
        //BOOST_CHECK_SMALL( abs(f - refProfileMeansquareValue), 1e-4);
        BOOST_CHECK_CLOSE_FRACTION( f, refProfileMeansquareValue, 1e-3);



        //check if esimated meansquare gradient residuals have the same values as the reference
        for (Int32 i = 0; i < refProfileGradientMeansquareValues.size(); i++)
        {
            Float64 val = g[i];
            Float64 refValue = refProfileGradientMeansquareValues[i];
            //BOOST_CHECK_SMALL( abs(val - refValue), 1e-4);
            BOOST_CHECK_CLOSE_FRACTION( val, refValue, 1e-3);
        }
    }
}

//*
//test the sym profile on its center wavelength
//-> flux value must be 1.0
//-> deriv amp value must be 1.0
BOOST_AUTO_TEST_CASE( Linemodel_multiline_profile_sym_value_center )
{
    std::string linecatalogPath = DATA_ROOT_DIR "LinemodelProfileTestCase/linecatalog_test_linemodel_profile_4EL.txt";

    std::vector<Float64> lambda;
    lambda.push_back(8600.00);
    lambda.push_back(8400.00);
    lambda.push_back(7600.00);
    lambda.push_back(7400.00);

    std::vector<Float64> refProfileValue(4, 1.0);
    std::vector<std::vector<Float64> >  gradResidualMatrix( lambda.size(), std::vector<Float64>(5, 0.0) );
    std::vector<Float64> gradResidualMeansquareMatrix(4, 1.0); //unused
    Float64 refProfileMeansquareValue=1.0; //unused

    Float64 opt_velocityEmission = 100.0;
    checkProfileValue(linecatalogPath, opt_velocityEmission, lambda, refProfileValue, false, gradResidualMatrix, false, refProfileMeansquareValue, gradResidualMeansquareMatrix);
}
//*/

//*
//test the sym profile on a the full 6000A to 10000A wavelength range
BOOST_AUTO_TEST_CASE( Linemodel_multiline_profile_sym_value_fullwavelengthrange )
{
    std::string linecatalogPath = DATA_ROOT_DIR "LinemodelProfileTestCase/linecatalog_test_linemodel_profile_4EL.txt";

    //init lambda table
    std::vector<Float64> lambda;
    for(Int32 k=6000; k<10000; k++)
    {
        lambda.push_back((Float64)k);
    }

    //init refProfileValue from fits file
    std::string signalRefFluxPath = DATA_ROOT_DIR "LinemodelProfileTestCase/signal_4lines_sig100_6000A_10000A_a1.0.fits";
    CSpectrumIOFitsReader reader;
    CSpectrum spectrum;

    BOOST_CHECK_NO_THROW(reader.Read( signalRefFluxPath.c_str(), spectrum ));
    std::vector<Float64> refProfileValue(lambda.size());
    CSpectrumFluxAxis& spcFluxAxis = spectrum.GetFluxAxis();
    for(UInt32 i=0; i<spcFluxAxis.GetSamplesCount(); i++){
        refProfileValue[i] = spcFluxAxis[i];
    }

    //load gradient reference values from CSV file
    bool enableGradientCheck = true;
    std::vector<std::vector<Float64> >  gradResidualMatrix( lambda.size(), std::vector<Float64>(5, 0.0) );
    std::string gradResidual_filePath = DATA_ROOT_DIR "LinemodelProfileTestCase/gradResidual_signalnoise4lines_sig400_6000A_10000A_m100_a1.0.txt";
    std::ifstream file;
    file.open( gradResidual_filePath, std::ifstream::in );
    bool fileOpenFailed = file.rdstate() & std::ios_base::failbit;
    BOOST_CHECK( !fileOpenFailed );
    std::string line;
    // Read file line by line
    Int32 kLine = 0;
    while( getline( file, line ) )
    {
        boost::char_separator<char> sep(" \t");
        // Tokenize each line
        typedef boost::tokenizer< boost::char_separator<char> > ttokenizer;
        ttokenizer tok( line, sep );
        ttokenizer::iterator it = tok.begin();
        if( it != tok.end() )
        {
            Int32 kCol = 0;
            // Parse line
            while( it != tok.end() )
            {
	      double val;
                val = boost::lexical_cast<double>(*it);
                gradResidualMatrix[kLine][kCol] = val;
                ++it;
                kCol++;
            }
            BOOST_CHECK( kCol==5 );
            kLine++;
        }
    }
    file.close();
    BOOST_CHECK( kLine==lambda.size() );

    //init meansquare gradient values (hardcoded)
    bool enableGradientMeansquareCheck=true;
    //load meansquare value
    Float64 refProfileMeansquareValue=0.0;
    std::string meansquare_filePath = DATA_ROOT_DIR "LinemodelProfileTestCase/residualMeanSquare_signalnoise4lines_sig400_6000A_10000A_m100_a1.0.txt";
    file.open( meansquare_filePath, std::ifstream::in );
    fileOpenFailed = file.rdstate() & std::ios_base::failbit;
    BOOST_CHECK( !fileOpenFailed );
    kLine = 0;
    while( getline( file, line ) )
    {
        boost::char_separator<char> sep(" \t");
        // Tokenize each line
        typedef boost::tokenizer< boost::char_separator<char> > ttokenizer;
        ttokenizer tok( line, sep );
        ttokenizer::iterator it = tok.begin();
        if( it != tok.end() )
        {
            // Parse line
            while( it != tok.end() )
            {
	      double val;
                val = boost::lexical_cast<double>(*it);
                refProfileMeansquareValue = val;
                ++it;
            }
            kLine++;
        }
    }
    file.close();

    //load gradientMeansquare value
    std::vector<Float64> refProfileGradientMeansquareValues(5, 0.0);
    std::string gradientMeansquare_filePath = DATA_ROOT_DIR "LinemodelProfileTestCase/gradResidualMeanSquare_signalnoise4lines_sig400_6000A_10000A_m100_a1.0.txt";
    file.open( gradientMeansquare_filePath, std::ifstream::in );
    fileOpenFailed = file.rdstate() & std::ios_base::failbit;
    BOOST_CHECK( !fileOpenFailed );
    kLine = 0;
    while( getline( file, line ) )
    {
        boost::char_separator<char> sep(" \t");
        // Tokenize each line
        typedef boost::tokenizer< boost::char_separator<char> > ttokenizer;
        ttokenizer tok( line, sep );
        ttokenizer::iterator it = tok.begin();
        if( it != tok.end() )
        {
            // Parse line
            while( it != tok.end() )
            {
	        double val;
                val = boost::lexical_cast<double>(*it);
                refProfileGradientMeansquareValues[kLine] = val;
                ++it;
            }
            kLine++;
        }
    }
    file.close();


    Float64 opt_velocityEmission = 100.0;
    checkProfileValue(linecatalogPath,
                      opt_velocityEmission,
                      lambda,
                      refProfileValue,
                      enableGradientCheck,
                      gradResidualMatrix,
                      enableGradientMeansquareCheck,
                      refProfileMeansquareValue,
                      refProfileGradientMeansquareValues);
}
//*/

BOOST_AUTO_TEST_SUITE_END()
