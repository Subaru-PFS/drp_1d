#include <epic/core/common/datatypes.h>
#include <epic/redshift/spectrum/spectrum.h>

#include <epic/redshift/continuum/irregularsamplingmedian.h>

#include <epic/redshift/spectrum/io/fitsreader.h>

#include <epic/redshift/linemodel/modelfittingresult.h>
#include <epic/redshift/linemodel/elementlist.h>

#include <boost/test/unit_test.hpp>

#include <boost/property_tree/ptree.hpp>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(LinemodelFit)


BOOST_AUTO_TEST_CASE( LinemodelFit_lmfit_velocity_synthetic2lines )
{
    std::string spectrumPath = "../test/data/LinemodelFitTestCase/simu_fit_synth_1.fits";
    std::string noisePath    = "../test/data/LinemodelFitTestCase/simu_fit_synth_1_noise.fits";
    std::string linecatalogPath = "../test/data/LinemodelFitTestCase/linecatalog_test_linemodel_fit_synth.txt";

    // load spectrum
    CSpectrumIOFitsReader reader;
    CSpectrum spectrum;

    Bool retVal = reader.Read( spectrumPath.c_str(), spectrum);
    BOOST_CHECK( retVal == true);

    // get continuum
    CContinuumIrregularSamplingMedian continuum;
    CSpectrumFluxAxis fluxAxisWithoutContinuumCalc;
    Int32 retValCont = continuum.RemoveContinuum( spectrum, fluxAxisWithoutContinuumCalc );
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
    Int32 typeFilter = CRay::nType_Emission; //CRay::nType_Absorption;
    Int32 forceFilter = CRay::nForce_Strong;
    CRayCatalog::TRayVector lineList = lineCatalog.GetFilteredList(typeFilter, forceFilter);

    std::string opt_fittingmethod = "lmfit";
    std::string opt_continuumcomponent = "fromspectrum";
    std::string opt_lineWidthType = "fixedvelocity";
    Float64 opt_resolution = 2350;
    Float64 opt_velocityEmission = 100;
    Float64 opt_velocityAbsorption = 300;
    std::string opt_rules = "no";


    CLineModelElementList model(spectrum, spectrumContinuum, lineList, opt_fittingmethod, opt_continuumcomponent, opt_lineWidthType, opt_resolution, opt_velocityEmission, opt_velocityAbsorption, opt_rules);
    TFloat64Range lambdaRange = TFloat64Range( 100.0, 12000.0 );
    CLineModelResult::SLineModelSolution modelSolution;
    Float64 z = 0.0;
    Float64 merit = model.fit(z, lambdaRange, modelSolution);

    Float64 amplitude_line1 = modelSolution.Amplitudes[0];
    Float64 amplitude_line2 = modelSolution.Amplitudes[1];
    Float64 velocity_emission = model.GetVelocityEmission();

    BOOST_CHECK_CLOSE_FRACTION( amplitude_line1, 1.0, 0.1);
    BOOST_CHECK_CLOSE_FRACTION( amplitude_line2, 3.0, 0.1);
    BOOST_CHECK_CLOSE_FRACTION( velocity_emission, 380.0, 0.1);
}

BOOST_AUTO_TEST_CASE( LinemodelFit_lmfit_velocity_pfsbatch6_emission_sym_single )
{
    std::string spectrumPath = "../test/data/LinemodelFitTestCase/55016588000024vacLine_TF.fits";
    std::string noisePath = "../test/data/LinemodelFitTestCase/55016588000024vacLine_ErrF.fits";
    std::string linecatalogPath = "../test/data/LinemodelFitTestCase/linecatalog_b9_emission_sym_single.txt";

    // load spectrum
    CSpectrumIOFitsReader reader;
    CSpectrum spectrum;

    Bool retVal = reader.Read( spectrumPath.c_str(), spectrum);
    BOOST_CHECK( retVal == true);

    // get continuum (zero for the test)
    CContinuumIrregularSamplingMedian continuum;
    CSpectrumFluxAxis fluxAxisWithoutContinuumCalc;
    Int32 retValCont = continuum.RemoveContinuum( spectrum, fluxAxisWithoutContinuumCalc );
    CSpectrum spectrumContinuum = spectrum;
    CSpectrumFluxAxis& continuumFluxAxis = spectrumContinuum.GetFluxAxis();
    for(UInt32 i=0; i<continuumFluxAxis.GetSamplesCount(); i++){
        continuumFluxAxis[i] -= fluxAxisWithoutContinuumCalc[i];
    }


    //get line catalog
    CRayCatalog lineCatalog;
    Bool rValue = lineCatalog.Load( linecatalogPath.c_str() );
    BOOST_CHECK( rValue == true);
    Int32 typeFilter = -1;// = CRay::nType_Emission; //CRay::nType_Absorption;
    Int32 forceFilter = -1;
    CRayCatalog::TRayVector lineList = lineCatalog.GetFilteredList(typeFilter, forceFilter);

    std::string opt_fittingmethod = "lmfit";
    std::string opt_continuumcomponent = "fromspectrum";
    std::string opt_lineWidthType = "fixedvelocity";
    Float64 opt_resolution = 2350;
    Float64 opt_velocityEmission = 100;
    Float64 opt_velocityAbsorption = 300;
    std::string opt_rules = "all";


    CLineModelElementList model(spectrum, spectrumContinuum, lineList, opt_fittingmethod, opt_continuumcomponent, opt_lineWidthType, opt_resolution, opt_velocityEmission, opt_velocityAbsorption, opt_rules);
    TFloat64Range lambdaRange = TFloat64Range( 3800.0, 12000.0 );
    CLineModelResult::SLineModelSolution modelSolution;
    Float64 z = 0.954114;
    Float64 merit = model.fit(z, lambdaRange, modelSolution);

}


BOOST_AUTO_TEST_SUITE_END()
