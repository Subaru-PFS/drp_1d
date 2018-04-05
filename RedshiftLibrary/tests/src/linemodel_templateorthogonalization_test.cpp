#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/spectrum.h>

#include <RedshiftLibrary/continuum/irregularsamplingmedian.h>

#include <RedshiftLibrary/spectrum/io/fitsreader.h>

#include <RedshiftLibrary/noise/flat.h>
#include <RedshiftLibrary/noise/fromfile.h>
#include <RedshiftLibrary/linemodel/modelfittingresult.h>
#include <RedshiftLibrary/linemodel/elementlist.h>
#include <RedshiftLibrary/linemodel/templatesortho.h>


#include <boost/test/unit_test.hpp>
#include "test-config.h"

#include <boost/property_tree/ptree.hpp>
#include <math.h>


using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(Linemodel_tplorthogonalization)

/***
 *
**/
Float64 processOrtho(std::string spectrumPath, std::string noisePath, std::string linecatalogPath,
                     std::string opt_fittingmethod,
                     Int32 lineTypeFilter,
                     Int32 forceFilter,
                     Float64 initVelocity,
                     Float64 z,
                     bool enableOrtho)
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


    // get continuum from Median in case of opt_continuumcomponent==fromspectrum
    CSpectrum spectrumContinuum = spectrum;
    CSpectrumFluxAxis& continuumFluxAxis = spectrumContinuum.GetFluxAxis();
    for(UInt32 i=0; i<continuumFluxAxis.GetSamplesCount(); i++){
        continuumFluxAxis[i] = 0.0; //put zero as continuum in case of "tplfit" continuum for linemodel
    }


    //get line catalog
    CRayCatalog lineCatalog;
    Bool rValueLoadLineCatalog = lineCatalog.Load( linecatalogPath.c_str() );
    BOOST_CHECK( rValueLoadLineCatalog == true);
    CRayCatalog::TRayVector lineList = lineCatalog.GetFilteredList(lineTypeFilter, forceFilter);
    BOOST_CHECK( lineList.size()>0);


    std::string opt_continuumcomponent = "tplfit";
    std::string opt_lineWidthType = "velocitydriven";
    Float64 opt_resolution = 2350; //unused with velocity driven linewidth
    Float64 opt_velocityEmission = initVelocity;
    Float64 opt_velocityAbsorption = initVelocity;
    std::string opt_rules = "no";
    std::string opt_rigidity = "rules";
    std::string opt_calibrationPath= DATA_ROOT_DIR "Linemodel_tplorthogalization/calibration";


    //these tplcatalog related variables are unused here.
    CTemplateCatalog tplCatalog;
    std::string templatesPath= DATA_ROOT_DIR "Linemodel_tplorthogalization/templates/";
    BOOST_TEST_MESSAGE( "Loading templates from " << templatesPath );
    Bool retValue = tplCatalog.Load(templatesPath.c_str());
    BOOST_CHECK( retValue == true);
    TStringList tplCategories = TStringList { "galaxy" };

    CTemplateCatalog finalTplCatalog;
    if(enableOrtho)
    {
        //prepare continuum templates catalog
        CTemplatesOrthogonalization tplOrtho(tplCatalog,
                                             tplCategories,
                                             opt_calibrationPath,
                                             lineList,
                                             opt_fittingmethod,
                                             opt_continuumcomponent,
                                             opt_lineWidthType,
                                             opt_resolution,
                                             opt_velocityEmission,
                                             opt_velocityAbsorption,
                                             opt_rules,
                                             opt_rigidity);
        finalTplCatalog = tplOrtho.getOrthogonalTplCatalog();
    }else{
        finalTplCatalog = tplCatalog;
    }


    CLineModelElementList model(spectrum,
                                spectrumContinuum,
                                finalTplCatalog,
                                tplCategories,
                                opt_calibrationPath,
                                lineList,
                                opt_fittingmethod,
                                opt_continuumcomponent,
                                opt_lineWidthType,
                                opt_resolution,
                                opt_velocityEmission,
                                opt_velocityAbsorption,
                                opt_rules,
                                opt_rigidity);

    TFloat64Range lambdaRange = TFloat64Range( 1000.0, 20000.0 );
    CLineModelSolution modelSolution;
    Float64 merit = model.fit(z, lambdaRange, modelSolution);
    BOOST_TEST_MESSAGE( "enableOrtho=" << enableOrtho << ", Merit = " << merit );

    Float64 cAmp = model.getFitContinuum_tplAmplitude();
    BOOST_TEST_MESSAGE( "enableOrtho=" << enableOrtho << ",cAmp=" << cAmp );


    /*//debug:
    CSpectrum spcModel = model.GetModelSpectrum();
    FILE* f = fopen( "tplortho_model_dbg.txt", "w+" );
    for( Int32 t=0;t<spcModel.GetSampleCount();t++)
    {
        fprintf( f, "%f %e\n", spcModel.GetSpectralAxis()[t], spcModel.GetFluxAxis()[t]);
    }
    fclose( f );
    //*/

    CLineModelSolution solution = model.GetModelSolution();
    BOOST_TEST_MESSAGE( "linemodel solution amplitude=" << solution.Amplitudes[0] );

    return cAmp;
}

/***
 * test the fullmodel - templates orthogonalization
 *
 *
***/
BOOST_AUTO_TEST_CASE( Linemodel_tplorthogonalization )
{
    std::string spectrumPath_z0 = DATA_ROOT_DIR "Linemodel_tplorthogalization/simu_fm_tplortho_synth_5k8k_ha.fits";
    std::string spectrumPath_z0p9 = DATA_ROOT_DIR "Linemodel_tplorthogalization/simu_fm_tplortho_synth_5k8k_ha_z0.9.fits";
    std::string noisePath = DATA_ROOT_DIR "Linemodel_tplorthogalization/simu_fm_tplortho_synth_ha_5k8k_noise.fits";
    std::string linecatalogPath = DATA_ROOT_DIR "Linemodel_tplorthogalization/linecatalog_b.txt";

    std::string opt_fittingmethod = "individual";
    Int32 lineTypeFilter = CRay::nType_Emission;
    Float64 initialVelocity = 1150.0;
    Int32 forceFilter = CRay::nForce_Strong;

    Float64 z = 0.0;

    Bool enableOrtho = false;
    Float64 continuumAmp = -1.0;

    //test the fullmodel without ortho: amp should be too high due to the presence of high amplitude line
    enableOrtho = false;
    continuumAmp = processOrtho(spectrumPath_z0, noisePath, linecatalogPath, opt_fittingmethod, lineTypeFilter, forceFilter, initialVelocity, z, enableOrtho);
    BOOST_CHECK( continuumAmp<1e-21 == false);

    /* //deactivate for now, until the ortho is fully implemented
    //test the fullmodel with ortho at z=0
    enableOrtho = true;
    continuumAmp = processOrtho(spectrumPath_z0, noisePath, linecatalogPath, opt_fittingmethod, lineTypeFilter, forceFilter, initialVelocity, z, enableOrtho);
    BOOST_CHECK( continuumAmp<1e-21 == true);

    //test the fullmodel with ortho at z=0.9
    z = 0.9;
    enableOrtho = true;
    continuumAmp = processOrtho(spectrumPath_z0p9, noisePath, linecatalogPath, opt_fittingmethod, lineTypeFilter, forceFilter, initialVelocity, z, enableOrtho);
    BOOST_CHECK( continuumAmp<1e-21 == true);
    //*/
}

BOOST_AUTO_TEST_SUITE_END()
