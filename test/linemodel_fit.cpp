#include <epic/core/common/datatypes.h>
#include <epic/redshift/spectrum/spectrum.h>


#include <epic/redshift/spectrum/io/fitsreader.h>

#include <epic/redshift/linemodel/modelfittingresult.h>
#include <epic/redshift/linemodel/elementlist.h>

#include <boost/test/unit_test.hpp>

#include <boost/property_tree/ptree.hpp>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(LinemodelFit)


BOOST_AUTO_TEST_CASE( LinemodelFitDumb )
{
    std::string spectrumPath = "../test/data/LinemodelRulesTestCase/simu_rules_ratiorange_1.fits";
    std::string noisePath = "../test/data/LinemodelRulesTestCase/simu_rules_ratiorange_1_noise.fits";
    std::string linecatalogPath = "../test/data/LinemodelRulesTestCase/raycatalog_test_elratiorules.txt";

    // load spectrum
    CSpectrumIOFitsReader reader;
    CSpectrum spectrum;

    Bool retVal = reader.Read( spectrumPath.c_str(), spectrum);
    BOOST_CHECK( retVal == true);

    // get continuum (zero for the test)
    CSpectrum spectrumContinuum = spectrum;
    CSpectrumFluxAxis& continuumFluxAxis = spectrumContinuum.GetFluxAxis();
    for(UInt32 i=0; i<continuumFluxAxis.GetSamplesCount(); i++){
        continuumFluxAxis[i] = 0.0;
    }

    //get line catalog
    CRayCatalog lineCatalog;
    Bool rValue = lineCatalog.Load( linecatalogPath.c_str() );
    BOOST_CHECK( rValue == true);
    Int32 typeFilter = typeFilter = CRay::nType_Emission; //CRay::nType_Absorption;
    Int32 forceFilter = CRay::nForce_Strong;
    CRayCatalog::TRayVector lineList = lineCatalog.GetFilteredList(typeFilter, forceFilter);

    std::string opt_fittingmethod = "lmfit";
    std::string opt_continuumcomponent = "fromspectrum";
    std::string opt_lineWidthType = "fixedvelocity";
    Float64 opt_resolution = 2350;
    Float64 opt_velocityEmission = 400;
    Float64 opt_velocityAbsorption = 300;
    std::string opt_rules = "all";


    CLineModelElementList model(spectrum, spectrumContinuum, lineList, opt_fittingmethod, opt_continuumcomponent, opt_lineWidthType, opt_resolution, opt_velocityEmission, opt_velocityAbsorption, opt_rules);
    TFloat64Range lambdaRange = TFloat64Range( 3800.0, 12000.0 );
    CLineModelResult::SLineModelSolution modelSolution;
    Float64 merit = model.fit(0.0, lambdaRange, modelSolution);

}


BOOST_AUTO_TEST_SUITE_END()
