#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/spectrum.h>

#include <RedshiftLibrary/continuum/irregularsamplingmedian.h>

#include <RedshiftLibrary/spectrum/io/fitsreader.h>

#include <RedshiftLibrary/noise/flat.h>
#include <RedshiftLibrary/noise/fromfile.h>
#include <RedshiftLibrary/linemodel/modelfittingresult.h>
#include <RedshiftLibrary/linemodel/elementlist.h>

#include <boost/test/unit_test.hpp>
#include "test-config.h"

#include <boost/property_tree/ptree.hpp>
#include <math.h>


using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(LinemodelFit)

void checkAmplitudeAndVelocityFit(std::string spectrumPath, std::string noisePath, std::string linecatalogPath, std::string opt_fittingmethod, Int32 lineTypeFilter, Int32 forceFilter, Float64 initVelocity, Float64 z, std::vector<Float64> ampsRef, Float64 fittedVelocityRef)
{
    // load spectrum
    CSpectrumIOFitsReader reader;
    CSpectrum spectrum;

    Bool retVal = reader.Read( spectrumPath.c_str(), std::shared_ptr<CSpectrum>(&spectrum));
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
    Float64 opt_velocityEmission = initVelocity;
    Float64 opt_velocityAbsorption = initVelocity;
    std::string opt_rules = "no";
    std::string opt_rigidity = "rules";
    std::string unused_calibrationPath="";



    //these tplcatalog related variables are unused here.
    CTemplateCatalog tplCatalog;
    Bool retValue = tplCatalog.Load( DATA_ROOT_DIR "templatecatalog/" );
    TStringList tplCategories;

    CLineModelElementList model(spectrum, spectrumContinuum, tplCatalog, tplCategories, unused_calibrationPath, lineList, opt_fittingmethod, opt_continuumcomponent, opt_lineWidthType, opt_resolution, opt_velocityEmission, opt_velocityAbsorption, opt_rules, opt_rigidity);
    TFloat64Range lambdaRange = TFloat64Range( 100.0, 12000.0 );
    CLineModelSolution modelSolution;
    Float64 merit = model.fit(z, lambdaRange, modelSolution);


    for(Int32 k=0; k<ampsRef.size(); k++)
    {
        Float64 amp_fitted = modelSolution.Amplitudes[k];
        BOOST_CHECK_CLOSE_FRACTION( amp_fitted, ampsRef[k], 0.05);
    }


    Float64 velocity;
    if(lineTypeFilter==CRay::nType_Emission)
    {
        velocity = model.GetVelocityEmission();
    }else
    {
        velocity = model.GetVelocityAbsorption();
    }
    BOOST_CHECK_SMALL( fabs(velocity - fittedVelocityRef), 15.0); //todo: check with more precision ?

}


//BOOST_AUTO_TEST_CASE( LinemodelFit_lmfit_velocity_pfsbatch6_emission_sym_single )
//{
//    std::string spectrumPath = DATA_ROOT_DIR "LinemodelFitTestCase/55016588000024vacLine_TF.fits";
//    std::string noisePath = DATA_ROOT_DIR "LinemodelFitTestCase/55016588000024vacLine_ErrF.fits";
//    std::string linecatalogPath = DATA_ROOT_DIR "LinemodelFitTestCase/linecatalog_b9_emission_sym_single.txt";

//    std::string opt_fittingmethod = "lmfit";
//    Int32 lineTypeFilter = CRay::nType_Emission;
//    Float64 initialVelocity = 100.0;
//    Int32 forceFilter = CRay::nForce_Strong;

//    std::vector<Float64> ampsRef;
//    ampsRef.push_back(-1);
//    ampsRef.push_back(1.25*1e-18);
//    ampsRef.push_back(2.79*1e-19);
//    ampsRef.push_back(6.9*1e-19);
//    ampsRef.push_back(3.8*1e-19);
//    ampsRef.push_back(3.8*1e-19);

//    Float64 emissionVelocityRef = 250.0;
//    Float64 z = 0.954114;
//    checkAmplitudeAndVelocityFit(spectrumPath, noisePath, linecatalogPath, opt_fittingmethod, lineTypeFilter, forceFilter, initialVelocity, z, ampsRef, emissionVelocityRef);

//}

//BOOST_AUTO_TEST_CASE( LinemodelFit_lmfit_velocity_pfsbatch6_absorption_sym_single )
//{
//    std::string spectrumPath = DATA_ROOT_DIR "LinemodelFitTestCase/10000663000008vacLine_TF.fits";
//    std::string noisePath = DATA_ROOT_DIR "LinemodelFitTestCase/10000663000008vacLine_ErrF.fits";
//    std::string linecatalogPath = DATA_ROOT_DIR "LinemodelFitTestCase/linecatalog_b9_absorption_sym_single.txt";

//    std::string opt_fittingmethod = "lmfit";
//    Int32 lineTypeFilter = CRay::nType_Absorption;
//    Float64 initialVelocity = 100.0;
//    Int32 forceFilter = -1;

//    std::vector<Float64> ampsRef;
//    ampsRef.push_back(1.97e-16);
//    ampsRef.push_back(2.38e-16);
//    ampsRef.push_back(1.8e-16);
//    ampsRef.push_back(1.97e-16);
//    ampsRef.push_back(1.86e-16);
//    ampsRef.push_back(1.51e-16);
//    ampsRef.push_back(1.43e-16);

//    Float64 absorptionVelocityRef = 775.0;
//    Float64 z = 0.0772;
//    checkAmplitudeAndVelocityFit(spectrumPath, noisePath, linecatalogPath, opt_fittingmethod, lineTypeFilter, forceFilter, initialVelocity, z, ampsRef, absorptionVelocityRef);

//}

//BOOST_AUTO_TEST_CASE( LinemodelFit_lmfit_velocity_pfshighz_emission_asym_single )
//{
//    std::string spectrumPath = DATA_ROOT_DIR "LinemodelFitTestCase/lya/spec_lines_oc_z5.70_ew50_fwhm10_nb25.0mag_TF.fits";
//    std::string noisePath = DATA_ROOT_DIR "LinemodelFitTestCase/lya/spec_lines_oc_z5.70_ew50_fwhm10_nb25.0mag_ErrF.fits";
//    std::string linecatalogPath = DATA_ROOT_DIR "LinemodelFitTestCase/linecatalog_b9_emission_asym_single.txt";

//    std::string opt_fittingmethod = "lmfit";
//    //std::string opt_fittingmethod = "individual";
//    Int32 lineTypeFilter = CRay::nType_Emission;
//    Float64 initialVelocity = 200.0;
//    Int32 forceFilter = -1;

//    std::vector<Float64> ampsRef;
//    ampsRef.push_back(1e-18);
//    Float64 emissionVelocityRef = 50;
//    //Float64 z = 5.702; //sym profile
//    Float64 z = 5.698; //asym profile
//    checkAmplitudeAndVelocityFit(spectrumPath, noisePath, linecatalogPath, opt_fittingmethod, lineTypeFilter, forceFilter, initialVelocity, z, ampsRef, emissionVelocityRef);

//}


BOOST_AUTO_TEST_CASE( LinemodelFit_lmfit_1EmissionLineF )
{
    std::string spectrumPath = DATA_ROOT_DIR "LinemodelFitTestCase/signalnoise_4lines_sig400_6000A_10000A.fits";
    std::string noisePath    = DATA_ROOT_DIR "LinemodelFitTestCase/noise_4lines_sig400_6000A_10000A.fits";
    std::string linecatalogPath = DATA_ROOT_DIR "LinemodelFitTestCase/linecatalog_test_linemodel_profile_1EL.txt";

    std::string opt_fittingmethod = "lmfit";
    Int32 lineTypeFilter = CRay::nType_Emission;
    Int32 forceFilter = CRay::nForce_Strong;

    std::vector<Float64> ampsRef;
    ampsRef.push_back(1.0);

    Float64 emissionVelocityRef = 400.0;
    Float64 z = 0.0;

    Float64 initialVelocity = 50.0;
    for(Int32 kV=0; kV<5; kV++)
    {
        BOOST_TEST_MESSAGE( "LinemodelFit_lmfit_1EmissionLineF: initial velocity: " << initialVelocity );
        checkAmplitudeAndVelocityFit(spectrumPath, noisePath, linecatalogPath, opt_fittingmethod, lineTypeFilter, forceFilter, initialVelocity, z, ampsRef, emissionVelocityRef);
        initialVelocity += 100.0;
    }

}

BOOST_AUTO_TEST_CASE( LinemodelFit_lmfit_4UnitaryEmissionLinesTF )
{
    std::string spectrumPath = DATA_ROOT_DIR "LinemodelFitTestCase/signal_4lines_sig100_6000A_10000A_a1.0.fits";
    std::string noisePath    = DATA_ROOT_DIR "LinemodelFitTestCase/noise_4lines_sig400_6000A_10000A.fits";
    std::string linecatalogPath = DATA_ROOT_DIR "LinemodelFitTestCase/linecatalog_test_linemodel_profile_4EL.txt";

    std::string opt_fittingmethod = "lmfit";
    Int32 lineTypeFilter = CRay::nType_Emission;
    Int32 forceFilter = CRay::nForce_Strong;

    std::vector<Float64> ampsRef;
    ampsRef.push_back(1.0);
    ampsRef.push_back(1.0);
    ampsRef.push_back(1.0);
    ampsRef.push_back(1.0);

    Float64 emissionVelocityRef = 100.0;
    Float64 z = 0.0;

    Float64 initialVelocity = 50.0;
    for(Int32 kV=0; kV<5; kV++)
    {
        BOOST_TEST_MESSAGE( "LinemodelFit_lmfit_4UnitaryEmissionLinesTF: initial velocity: " << initialVelocity );
        checkAmplitudeAndVelocityFit(spectrumPath, noisePath, linecatalogPath, opt_fittingmethod, lineTypeFilter, forceFilter, initialVelocity, z, ampsRef, emissionVelocityRef);
        initialVelocity += 100.0;
    }
}

BOOST_AUTO_TEST_CASE( LinemodelFit_lmfit_4UnitaryEmissionLinesTF_normalisation )
{
    std::string spectrumPath = DATA_ROOT_DIR "LinemodelFitTestCase/signal_4lines_sig100_6000A_10000A_a1.0_1em17.fits";
    std::string noisePath    = DATA_ROOT_DIR "LinemodelFitTestCase/noise_4lines_sig400_6000A_10000A_1em17.fits";
    std::string linecatalogPath = DATA_ROOT_DIR "LinemodelFitTestCase/linecatalog_test_linemodel_profile_4EL.txt";

    std::string opt_fittingmethod = "lmfit";
    Int32 lineTypeFilter = CRay::nType_Emission;
    Int32 forceFilter = CRay::nForce_Strong;

    std::vector<Float64> ampsRef;
    ampsRef.push_back(1.0e-17);
    ampsRef.push_back(1.0e-17);
    ampsRef.push_back(1.0e-17);
    ampsRef.push_back(1.0e-17);

    Float64 emissionVelocityRef = 100.0;
    Float64 z = 0.0;

    Float64 initialVelocity = 50.0;
    for(Int32 kV=0; kV<5; kV++)
    {
        BOOST_TEST_MESSAGE( "LinemodelFit_lmfit_4UnitaryEmissionLinesTF_normalisation: initial velocity: " << initialVelocity );
        checkAmplitudeAndVelocityFit(spectrumPath, noisePath, linecatalogPath, opt_fittingmethod, lineTypeFilter, forceFilter, initialVelocity, z, ampsRef, emissionVelocityRef);
        initialVelocity += 100.0;
    }
}

BOOST_AUTO_TEST_CASE( LinemodelFit_lmfit_4EmissionLinesF )
{
    std::string spectrumPath = DATA_ROOT_DIR "LinemodelFitTestCase/signalnoise_4lines_sig400_6000A_10000A.fits";
    std::string noisePath    = DATA_ROOT_DIR "LinemodelFitTestCase/noise_4lines_sig400_6000A_10000A.fits";
    std::string linecatalogPath = DATA_ROOT_DIR "LinemodelFitTestCase/linecatalog_test_linemodel_profile_4EL.txt";

    std::string opt_fittingmethod = "lmfit";
    Int32 lineTypeFilter = CRay::nType_Emission;
    Int32 forceFilter = CRay::nForce_Strong;

    std::vector<Float64> ampsRef;
    ampsRef.push_back(1.0);
    ampsRef.push_back(3.0);
    ampsRef.push_back(4.0);
    ampsRef.push_back(0.1);

    Float64 emissionVelocityRef = 400.0;
    Float64 z = 0.0;

    Float64 initialVelocity = 50.0;
    for(Int32 kV=0; kV<5; kV++)
    {
        BOOST_TEST_MESSAGE( "LinemodelFit_lmfit_4EmissionLinesF: initial velocity: " << initialVelocity );
        checkAmplitudeAndVelocityFit(spectrumPath, noisePath, linecatalogPath, opt_fittingmethod, lineTypeFilter, forceFilter, initialVelocity, z, ampsRef, emissionVelocityRef);
        initialVelocity += 100.0;
    }

}


//BOOST_AUTO_TEST_CASE( LinemodelFit_lbfgsfit_4UnitaryEmissionLinesTF )
//{
//    std::string spectrumPath = DATA_ROOT_DIR "LinemodelFitTestCase/signal_4lines_sig100_6000A_10000A_a1.0.fits";
//    std::string noisePath    = DATA_ROOT_DIR "LinemodelFitTestCase/noise_4lines_sig400_6000A_10000A.fits";
//    std::string linecatalogPath = DATA_ROOT_DIR "LinemodelFitTestCase/linecatalog_test_linemodel_profile_4EL.txt";

//    std::string opt_fittingmethod = "lbfgsfit";
//    Int32 lineTypeFilter = CRay::nType_Emission;
//    Int32 forceFilter = CRay::nForce_Strong;

//    std::vector<Float64> ampsRef;
//    ampsRef.push_back(1.0);
//    ampsRef.push_back(1.0);
//    ampsRef.push_back(1.0);
//    ampsRef.push_back(1.0);

//    Float64 emissionVelocityRef = 100.0;
//    Float64 z = 0.0;

//    Float64 initialVelocity = 50.0;
//    for(Int32 kV=0; kV<5; kV++)
//    {
//        BOOST_TEST_MESSAGE( "LinemodelFit_lbfgsfit_4UnitaryEmissionLinesTF: initial velocity: " << initialVelocity );
//        checkAmplitudeAndVelocityFit(spectrumPath, noisePath, linecatalogPath, opt_fittingmethod, lineTypeFilter, forceFilter, initialVelocity, z, ampsRef, emissionVelocityRef);
//        initialVelocity += 100.0;
//    }
//}

//BOOST_AUTO_TEST_CASE( LinemodelFit_lbfgsfit_1EmissionLineF )
//{
//    std::string spectrumPath = DATA_ROOT_DIR "LinemodelFitTestCase/signalnoise_4lines_sig400_6000A_10000A.fits";
//    std::string noisePath    = DATA_ROOT_DIR "LinemodelFitTestCase/noise_4lines_sig400_6000A_10000A.fits";
//    std::string linecatalogPath = DATA_ROOT_DIR "LinemodelFitTestCase/linecatalog_test_linemodel_profile_1EL.txt";

//    std::string opt_fittingmethod = "lbfgsfit";//"hybrid";//
//    Int32 lineTypeFilter = CRay::nType_Emission;
//    Int32 forceFilter = CRay::nForce_Strong;

//    std::vector<Float64> ampsRef;
//    ampsRef.push_back(1.0);

//    Float64 emissionVelocityRef = 400.0;
//    Float64 z = 0.0;

//    Float64 initialVelocity = 450.0;
//    for(Int32 kV=0; kV<1; kV++)
//    {
//        BOOST_TEST_MESSAGE( "LinemodelFit_lbfgsfit_1EmissionLineF: initial velocity: " << initialVelocity );
//        checkAmplitudeAndVelocityFit(spectrumPath, noisePath, linecatalogPath, opt_fittingmethod, lineTypeFilter, forceFilter, initialVelocity, z, ampsRef, emissionVelocityRef);
//        initialVelocity += 400.0;
//    }

//}


//BOOST_AUTO_TEST_CASE( LinemodelFit_lbfgsfit_4EmissionLinesF )
//{
//    std::string spectrumPath = DATA_ROOT_DIR "LinemodelFitTestCase/signalnoise_4lines_sig400_6000A_10000A.fits";
//    std::string noisePath    = DATA_ROOT_DIR "LinemodelFitTestCase/noise_4lines_sig400_6000A_10000A.fits";
//    std::string linecatalogPath = DATA_ROOT_DIR "LinemodelFitTestCase/linecatalog_test_linemodel_profile_4EL.txt";

//    std::string opt_fittingmethod = "lbfgsfit";
//    Int32 lineTypeFilter = CRay::nType_Emission;
//    Int32 forceFilter = CRay::nForce_Strong;

//    std::vector<Float64> ampsRef;
//    ampsRef.push_back(1.0);
//    ampsRef.push_back(3.0);
//    ampsRef.push_back(4.0);
//    ampsRef.push_back(0.1);

//    Float64 emissionVelocityRef = 400.0;
//    Float64 z = 0.0;

//    Float64 initialVelocity = 100.0;
//    for(Int32 kV=0; kV<1; kV++)
//    {
//        BOOST_TEST_MESSAGE( "LinemodelFit_lbfgsfit_4EmissionLinesF: initial velocity: " << initialVelocity );
//        checkAmplitudeAndVelocityFit(spectrumPath, noisePath, linecatalogPath, opt_fittingmethod, lineTypeFilter, forceFilter, initialVelocity, z, ampsRef, emissionVelocityRef);
//        initialVelocity += 100.0;
//    }

//}

BOOST_AUTO_TEST_SUITE_END()
