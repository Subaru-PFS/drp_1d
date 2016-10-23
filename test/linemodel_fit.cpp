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


using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(LinemodelFit)

void checkAmplitudeAndVelocityFit(std::string spectrumPath, std::string noisePath, std::string linecatalogPath, std::string opt_fittingmethod, Int32 lineTypeFilter, Int32 forceFilter, Float64 initVelocity, Float64 z, std::vector<Float64> ampsRef, Float64 fittedVelocityRef)
{
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
    CContinuumIrregularSamplingMedian continuum;
    CSpectrumFluxAxis fluxAxisWithoutContinuumCalc;
    Int32 retValCont = continuum.RemoveContinuum( spectrum, fluxAxisWithoutContinuumCalc );
    CSpectrum spectrumContinuum = spectrum;
    CSpectrumFluxAxis& continuumFluxAxis = spectrumContinuum.GetFluxAxis();
    for(UInt32 i=0; i<continuumFluxAxis.GetSamplesCount(); i++){
        continuumFluxAxis[i] -= fluxAxisWithoutContinuumCalc[i];
        //continuumFluxAxis[i] = 0.0;
    }


    //get line catalog
    CRayCatalog lineCatalog;
    Bool rValue = lineCatalog.Load( linecatalogPath.c_str() );
    BOOST_CHECK( rValue == true);
    CRayCatalog::TRayVector lineList = lineCatalog.GetFilteredList(lineTypeFilter, forceFilter);
    BOOST_CHECK( lineList.size()>0);


    std::string opt_continuumcomponent = "fromspectrum";
    std::string opt_lineWidthType = "combined";
    Float64 opt_resolution = 2350;
    Float64 opt_velocityEmission = initVelocity;
    Float64 opt_velocityAbsorption = initVelocity;
    std::string opt_rules = "no";


    //these tplcatalog related variables are unused here.
    CTemplateCatalog tplCatalog;
    Bool retValue = tplCatalog.Load( "../test/data/templatecatalog/" );
    TStringList tplCategories;

    CLineModelElementList model(spectrum, spectrumContinuum, tplCatalog, tplCategories, lineList, opt_fittingmethod, opt_continuumcomponent, opt_lineWidthType, opt_resolution, opt_velocityEmission, opt_velocityAbsorption, opt_rules);
    TFloat64Range lambdaRange = TFloat64Range( 100.0, 12000.0 );
    CLineModelResult::SLineModelSolution modelSolution;
    Float64 merit = model.fit(z, lambdaRange, modelSolution);


    for(Int32 k=0; k<ampsRef.size(); k++)
    {
        Float64 amp_fitted = modelSolution.Amplitudes[k];
        BOOST_CHECK_CLOSE_FRACTION( amp_fitted, ampsRef[k], 0.15);
    }


    Float64 velocity;
    if(lineTypeFilter==CRay::nType_Emission)
    {
        velocity = model.GetVelocityEmission();
    }else
    {
        velocity = model.GetVelocityAbsorption();
    }
    BOOST_CHECK_SMALL( fabs(velocity - fittedVelocityRef), 50.0); //todo: check with more precision ?

}

BOOST_AUTO_TEST_CASE( LinemodelFit_lmfit_velocity_synthetic2lines_no_normalization )
{
    std::string spectrumPath = "../test/data/LinemodelFitTestCase/simu_fit_synth_3.fits";
    std::string noisePath    = "../test/data/LinemodelFitTestCase/simu_fit_synth_3_noise.fits";
    std::string linecatalogPath = "../test/data/LinemodelFitTestCase/linecatalog_test_linemodel_fit_synth_b.txt";

    std::string opt_fittingmethod = "lmfit";
    Int32 lineTypeFilter = CRay::nType_Emission;
    Float64 initialVelocity = 100.0;
    Int32 forceFilter = CRay::nForce_Strong;

    std::vector<Float64> ampsRef;
    ampsRef.push_back(1.0);
    ampsRef.push_back(3.0);
    ampsRef.push_back(4.0);
    ampsRef.push_back(1.0);

    Float64 emissionVelocityRef = 377.0;
    Float64 z = 0.0;
    checkAmplitudeAndVelocityFit(spectrumPath, noisePath, linecatalogPath, opt_fittingmethod, lineTypeFilter, forceFilter, initialVelocity, z, ampsRef, emissionVelocityRef);
}


BOOST_AUTO_TEST_CASE( LinemodelFit_lmfit_velocity_synthetic2lines_normalization )
{
    std::string spectrumPath = "../test/data/LinemodelFitTestCase/simu_fit_synth_3_weak.fits";
    std::string noisePath    = "../test/data/LinemodelFitTestCase/simu_fit_synth_3_weak_noise.fits";
    std::string linecatalogPath = "../test/data/LinemodelFitTestCase/linecatalog_test_linemodel_fit_synth_b.txt";

    std::string opt_fittingmethod = "lmfit";
    Int32 lineTypeFilter = CRay::nType_Emission;
    Float64 initialVelocity = 100.0;
    Int32 forceFilter = CRay::nForce_Strong;

    std::vector<Float64> ampsRef;
    ampsRef.push_back(1.0*1e-17);
    ampsRef.push_back(3.0*1e-17);
    ampsRef.push_back(4.0*1e-17);
    ampsRef.push_back(1.0*1e-17);

    Float64 emissionVelocityRef = 377.0;
    Float64 z = 0.0;
    checkAmplitudeAndVelocityFit(spectrumPath, noisePath, linecatalogPath, opt_fittingmethod, lineTypeFilter, forceFilter, initialVelocity, z, ampsRef, emissionVelocityRef);
}

BOOST_AUTO_TEST_CASE( LinemodelFit_lmfit_velocity_pfsbatch6_emission_sym_single )
{
    std::string spectrumPath = "../test/data/LinemodelFitTestCase/55016588000024vacLine_TF.fits";
    std::string noisePath = "../test/data/LinemodelFitTestCase/55016588000024vacLine_ErrF.fits";
    std::string linecatalogPath = "../test/data/LinemodelFitTestCase/linecatalog_b9_emission_sym_single.txt";

    std::string opt_fittingmethod = "lmfit";
    Int32 lineTypeFilter = CRay::nType_Emission;
    Float64 initialVelocity = 100.0;
    Int32 forceFilter = CRay::nForce_Strong;

    std::vector<Float64> ampsRef;
    ampsRef.push_back(-1);
    ampsRef.push_back(1.25*1e-18);
    ampsRef.push_back(2.79*1e-19);
    ampsRef.push_back(6.9*1e-19);
    ampsRef.push_back(3.8*1e-19);
    ampsRef.push_back(3.8*1e-19);

    Float64 emissionVelocityRef = 250.0;
    Float64 z = 0.954114;
    checkAmplitudeAndVelocityFit(spectrumPath, noisePath, linecatalogPath, opt_fittingmethod, lineTypeFilter, forceFilter, initialVelocity, z, ampsRef, emissionVelocityRef);

}

BOOST_AUTO_TEST_CASE( LinemodelFit_lmfit_velocity_pfsbatch6_absorption_sym_single )
{
    std::string spectrumPath = "../test/data/LinemodelFitTestCase/10000663000008vacLine_TF.fits";
    std::string noisePath = "../test/data/LinemodelFitTestCase/10000663000008vacLine_ErrF.fits";
    std::string linecatalogPath = "../test/data/LinemodelFitTestCase/linecatalog_b9_absorption_sym_single.txt";

    std::string opt_fittingmethod = "lmfit";
    Int32 lineTypeFilter = CRay::nType_Absorption;
    Float64 initialVelocity = 100.0;
    Int32 forceFilter = -1;

    std::vector<Float64> ampsRef;
    ampsRef.push_back(1.97e-16);
    ampsRef.push_back(2.38e-16);
    ampsRef.push_back(1.8e-16);
    ampsRef.push_back(1.97e-16);
    ampsRef.push_back(1.86e-16);
    ampsRef.push_back(1.51e-16);
    ampsRef.push_back(1.43e-16);

    Float64 absorptionVelocityRef = 775.0;
    Float64 z = 0.0772;
    checkAmplitudeAndVelocityFit(spectrumPath, noisePath, linecatalogPath, opt_fittingmethod, lineTypeFilter, forceFilter, initialVelocity, z, ampsRef, absorptionVelocityRef);

}

//BOOST_AUTO_TEST_CASE( LinemodelFit_lmfit_velocity_pfshighz_emission_asym_single )
//{
//    std::string spectrumPath = "../test/data/LinemodelFitTestCase/lya/spec_lines_oc_z5.70_ew50_fwhm10_nb25.0mag_TF.fits";
//    std::string noisePath = "../test/data/LinemodelFitTestCase/lya/spec_lines_oc_z5.70_ew50_fwhm10_nb25.0mag_ErrF.fits";
//    std::string linecatalogPath = "../test/data/LinemodelFitTestCase/linecatalog_b9_emission_asym_single.txt";

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


BOOST_AUTO_TEST_CASE( LinemodelFit_lbfgsfit_1EmissionLine )
{
    std::string spectrumPath = "../test/data/LinemodelFitTestCase/simu_fit_synth_3.fits";
    std::string noisePath    = "../test/data/LinemodelFitTestCase/simu_fit_synth_3_noise.fits";
    std::string linecatalogPath = "../test/data/LinemodelFitTestCase/linecatalog_test_linemodel_fit_synth_1line.txt";

    std::string opt_fittingmethod = "lbfgsfit";//"hybrid";//
    Int32 lineTypeFilter = CRay::nType_Emission;
    Float64 initialVelocity = 250.0;
    Int32 forceFilter = CRay::nForce_Strong;

    std::vector<Float64> ampsRef;
    ampsRef.push_back(1.0);

    Float64 emissionVelocityRef = 380;//377.0;
    Float64 z = 0.0;
    checkAmplitudeAndVelocityFit(spectrumPath, noisePath, linecatalogPath, opt_fittingmethod, lineTypeFilter, forceFilter, initialVelocity, z, ampsRef, emissionVelocityRef);
}

//BOOST_AUTO_TEST_CASE( LinemodelFit_lbfgsfit_2EmissionLines )
//{
//    std::string spectrumPath = "../test/data/LinemodelFitTestCase/simu_fit_synth_3.fits";
//    std::string noisePath    = "../test/data/LinemodelFitTestCase/simu_fit_synth_3_noise.fits";
//    std::string linecatalogPath = "../test/data/LinemodelFitTestCase/linecatalog_test_linemodel_fit_synth_2lines.txt";

//    std::string opt_fittingmethod = "lbfgsfit";//"hybrid";//
//    Int32 lineTypeFilter = CRay::nType_Emission;
//    Float64 initialVelocity = 250.0;
//    Int32 forceFilter = CRay::nForce_Strong;

//    std::vector<Float64> ampsRef;
//    ampsRef.push_back(3.0);
//    ampsRef.push_back(1.0);

//    Float64 emissionVelocityRef = 377.0;
//    Float64 z = 0.0;
//    checkAmplitudeAndVelocityFit(spectrumPath, noisePath, linecatalogPath, opt_fittingmethod, lineTypeFilter, forceFilter, initialVelocity, z, ampsRef, emissionVelocityRef);
//}

//BOOST_AUTO_TEST_CASE( LinemodelFit_lbfgsfit_4EmissionLines )
//{
//    std::string spectrumPath = "../test/data/LinemodelFitTestCase/simu_fit_synth_3.fits";
//    std::string noisePath    = "../test/data/LinemodelFitTestCase/simu_fit_synth_3_noise.fits";
//    std::string linecatalogPath = "../test/data/LinemodelFitTestCase/linecatalog_test_linemodel_fit_synth_b.txt";

//    std::string opt_fittingmethod = "lbfgsfit";//"hybrid";//
//    Int32 lineTypeFilter = CRay::nType_Emission;
//    Float64 initialVelocity = 250.0;
//    Int32 forceFilter = CRay::nForce_Strong;

//    std::vector<Float64> ampsRef;
//    ampsRef.push_back(1.0);
//    ampsRef.push_back(3.0);
//    ampsRef.push_back(4.0);
//    ampsRef.push_back(1.0);

//    Float64 emissionVelocityRef = 377.0;
//    Float64 z = 0.0;
//    checkAmplitudeAndVelocityFit(spectrumPath, noisePath, linecatalogPath, opt_fittingmethod, lineTypeFilter, forceFilter, initialVelocity, z, ampsRef, emissionVelocityRef);
//}

BOOST_AUTO_TEST_SUITE_END()
