#include <epic/core/common/datatypes.h>
#include <epic/redshift/spectrum/spectrum.h>

#include <epic/redshift/continuum/irregularsamplingmedian.h>

#include <epic/redshift/spectrum/io/fitsreader.h>

#include <epic/redshift/noise/flat.h>
#include <epic/redshift/noise/fromfile.h>
#include <epic/redshift/linemodel/modelfittingresult.h>
#include <epic/redshift/linemodel/elementlist.h>

#include <boost/test/unit_test.hpp>

#include <boost/property_tree/ptree.hpp>
#include <math.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>


using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(LinemodelMultilineProfile)

/**
 * @brief checkProfileValue
 * @param linecatalogPath
 * @param opt_velocityEmission
 * @param lambdas: a set of wavelengths where the flux will be tested with regard to a reference
 * @param refProfileFluxValues: n reference flux values corresponding to the wavelengths in lambdas
 * @param enableGradientCheck
 * @param refProfileDerivAmplitudeValues: m reference gradient values corresponding to each element in the catalog (first d_amps, then d_velocity)
 */
void checkProfileValue(std::string linecatalogPath,
                       Float64 opt_velocityEmission,
                       std::vector<Float64> lambdas,
                       std::vector<Float64> refProfileFluxValues,
                       bool enableGradientCheck,
                       std::vector<std::vector<Float64> >  refProfileGradientValues)
{
    //some input params
    std::string spectrumPath = "../test/data/LinemodelProfileTestCase/signalnoise_4lines_sig400_6000A_10000A.fits";       //unused, only needed for linemodel initialization
    std::string noisePath    = "../test/data/LinemodelProfileTestCase/noise_4lines_sig400_6000A_10000A.fits";     //unused, only needed for linemodel initialization
    Int32 lineTypeFilter = CRay::nType_Emission;
    Int32 forceFilter = CRay::nForce_Strong;
    std::string opt_fittingmethod = "ones"; //all the elements amplitudes set to 1.0
    Float64 z = 0.0;

    // load spectrum
    CSpectrumIOFitsReader reader;
    CSpectrum spectrum;

    Bool retVal = reader.Read( spectrumPath.c_str(), spectrum);
    BOOST_CHECK( retVal == true);
    CNoiseFromFile noise;
    retVal = noise.SetNoiseFilePath( noisePath.c_str() );
    BOOST_CHECK( retVal == true);
    retVal = noise.AddNoise( spectrum ) ;
    BOOST_CHECK( retVal == true);


    // get continuum
    //CContinuumIrregularSamplingMedian continuum;
    //CSpectrumFluxAxis fluxAxisWithoutContinuumCalc;
    //Int32 retValCont = continuum.RemoveContinuum( spectrum, fluxAxisWithoutContinuumCalc );
    CSpectrum spectrumContinuum = spectrum;
    CSpectrumFluxAxis& continuumFluxAxis = spectrumContinuum.GetFluxAxis();
    for(UInt32 i=0; i<continuumFluxAxis.GetSamplesCount(); i++){
        //continuumFluxAxis[i] -= fluxAxisWithoutContinuumCalc[i];
        continuumFluxAxis[i] = 0.0;
    }


    //get line catalog
    CRayCatalog lineCatalog;
    Bool rValue = lineCatalog.Load( linecatalogPath.c_str() );
    BOOST_CHECK( rValue == true);
    CRayCatalog::TRayVector lineList = lineCatalog.GetFilteredList(lineTypeFilter, forceFilter);
    BOOST_CHECK( lineList.size()>0);


    std::string opt_continuumcomponent = "fromspectrum";
    std::string opt_lineWidthType = "velocitydriven";

    Float64 opt_resolution = 2350; //unused with velocity driven linewidth
    Float64 opt_velocityAbsorption = 100.0; //unused
    std::string opt_rules = "no";


    //these tplcatalog related variables are unused here.
    CTemplateCatalog tplCatalog;
    Bool retValue = tplCatalog.Load( "../test/data/templatecatalog/" );
    TStringList tplCategories;

    CLineModelElementList model(spectrum, spectrumContinuum, tplCatalog, tplCategories, lineList, opt_fittingmethod, opt_continuumcomponent, opt_lineWidthType, opt_resolution, opt_velocityEmission, opt_velocityAbsorption, opt_rules);
    TFloat64Range lambdaRange = TFloat64Range( 100.0, 12000.0 );
    CLineModelResult::SLineModelSolution modelSolution;
    Float64 merit = model.fit(z, lambdaRange, modelSolution);

    CSpectrumSpectralAxis spcSpectralAxis = spectrum.GetSpectralAxis();

    std::vector<Int32> validEltsIdx;
    for (Int32 iElt=0; iElt<model.m_Elements.size(); iElt++)
    {
        validEltsIdx.push_back(iElt);
    }

    model.refreshModelDerivSigmaUnderElements(validEltsIdx);

    //check the flux values
    for(Int32 k=0; k<lambdas.size(); k++)
    {
        Float64 iLambda = spcSpectralAxis.GetIndexAtWaveLength(lambdas[k]);
        Float64 profileValue = model.getModelFluxVal(iLambda);
        Float64 refFluxValue = refProfileFluxValues[k];
        //BOOST_CHECK_CLOSE_FRACTION( profileValue, refFluxValue, 1e-5);
        BOOST_CHECK_SMALL( abs(profileValue - refFluxValue), 1e-3);
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
                BOOST_CHECK_SMALL( abs(val - refDerivAmpValue), 1e-5);
            }

        }

        //check the deriv-velocity values
        Int32 idxVelGradient = validEltsIdx.size();
        for(Int32 k=0; k<lambdas.size(); k++)
        {
            Float64 iLambda = spcSpectralAxis.GetIndexAtWaveLength(lambdas[k]);
            Float64 val = model.getModelFluxDerivSigmaVal(iLambda);
            Float64 refDerivVelValue = refProfileGradientValues[k][idxVelGradient];
            BOOST_CHECK_SMALL( abs(val - refDerivVelValue), 1e-4);
        }
    }
}

//*
//test the sym profile on its center wavelength
//-> flux value must be 1.0
//-> deriv amp value must be 1.0
BOOST_AUTO_TEST_CASE( Linemodel_multiline_profile_sym_value_center )
{
    std::string linecatalogPath = "../test/data/LinemodelProfileTestCase/linecatalog_test_linemodel_profile_4EL.txt";

    std::vector<Float64> lambda;
    lambda.push_back(8600.00);
    lambda.push_back(8400.00);
    lambda.push_back(7600.00);
    lambda.push_back(7400.00);

    std::vector<Float64> refProfileValue(4, 1.0);
    std::vector<std::vector<Float64> >  gradResidualMatrix( lambda.size(), std::vector<Float64>(5, 0.0) );

    Float64 opt_velocityEmission = 100.0;
    checkProfileValue(linecatalogPath, opt_velocityEmission, lambda, refProfileValue, false, gradResidualMatrix);
}
//*/

//*
//test the sym profile on a the full 6000A to 10000A wavelength range
BOOST_AUTO_TEST_CASE( Linemodel_multiline_profile_sym_value_fullwavelengthrange )
{
    std::string linecatalogPath = "../test/data/LinemodelProfileTestCase/linecatalog_test_linemodel_profile_4EL.txt";

    //init lambda table
    std::vector<Float64> lambda;
    for(Int32 k=6000; k<10000; k++)
    {
        lambda.push_back((Float64)k);
    }

    //init refProfileValue from fits file
    std::string signalRefFluxPath = "../test/data/LinemodelProfileTestCase/signal_4lines_sig100_6000A_10000A_a1.0.fits";
    CSpectrumIOFitsReader reader;
    CSpectrum spectrum;
    Bool retVal = reader.Read( signalRefFluxPath.c_str(), spectrum);
    BOOST_CHECK( retVal == true);
    std::vector<Float64> refProfileValue(lambda.size());
    CSpectrumFluxAxis& spcFluxAxis = spectrum.GetFluxAxis();
    for(UInt32 i=0; i<spcFluxAxis.GetSamplesCount(); i++){
        refProfileValue[i] = spcFluxAxis[i];
    }

    //load gradient reference values from CSV file
    std::vector<std::vector<Float64> >  gradResidualMatrix( lambda.size(), std::vector<Float64>(5, 0.0) );
    std::string gradResidual_filePath = "../test/data/LinemodelProfileTestCase/gradResidual_signalnoise4lines_sig400_6000A_10000A_m100_a1.0.txt";
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
            double val = 0.0;
            while( it != tok.end() )
            {
                val = boost::lexical_cast<double>(*it);
                gradResidualMatrix[kLine][kCol] = val;
                ++it;
                kCol++;
            }
            BOOST_CHECK( kCol==5 );
            kLine++;
        }
    }
    BOOST_CHECK( kLine==lambda.size() );

    bool enableGradientCheck = true;
    Float64 opt_velocityEmission = 100.0;
    checkProfileValue(linecatalogPath, opt_velocityEmission, lambda, refProfileValue, enableGradientCheck, gradResidualMatrix);
}
//*/

BOOST_AUTO_TEST_SUITE_END()
