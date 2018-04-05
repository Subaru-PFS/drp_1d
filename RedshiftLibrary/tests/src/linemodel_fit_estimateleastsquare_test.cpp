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

BOOST_AUTO_TEST_SUITE(LinemodelFit_EstimateLeastSquare)

/***
 * Check that the least-square-fast provides the same value as the least-square method
 *
**/
void checkLeastSquareFast(std::string spectrumPath, std::string noisePath, std::string linecatalogPath, std::string opt_fittingmethod, std::string opt_continuumcomponent, Int32 lineTypeFilter, Int32 forceFilter, Float64 initVelocity, Float64 z)
{
    // load spectrum
    CSpectrumIOFitsReader reader;
    std::shared_ptr<CSpectrum> spectrum = std::shared_ptr<CSpectrum>( new CSpectrum() );
    Bool retVal = reader.Read( spectrumPath.c_str(), spectrum);
    BOOST_CHECK( retVal == true);
    CNoiseFromFile noise;
    retVal = noise.SetNoiseFilePath( noisePath.c_str() );
    BOOST_CHECK( retVal == true);
    retVal = noise.AddNoise( *spectrum ) ;
    BOOST_CHECK( retVal == true);


    // get continuum from Median in case of opt_continuumcomponent==fromspectrum
    CContinuumIrregularSamplingMedian continuum;
    CSpectrumFluxAxis fluxAxisWithoutContinuumCalc;
    Int32 retValContinuumEstimation = continuum.RemoveContinuum( *spectrum, fluxAxisWithoutContinuumCalc );
    BOOST_CHECK( retValContinuumEstimation);
    CSpectrumFluxAxis& continuumFluxAxis = spectrum->GetFluxAxis();
    for(UInt32 i=0; i<continuumFluxAxis.GetSamplesCount(); i++){
        if(opt_continuumcomponent == "fromspectrum")
        {
            continuumFluxAxis[i] -= fluxAxisWithoutContinuumCalc[i];
        }else{
            continuumFluxAxis[i] = 0.0; //put zero as continuum in case of "tplfit" continuum for linemodel
        }
    }


    //get line catalog
    CRayCatalog lineCatalog;
    Bool rValueLoadLineCatalog = lineCatalog.Load( linecatalogPath.c_str() );
    BOOST_CHECK( rValueLoadLineCatalog == true);
    CRayCatalog::TRayVector lineList = lineCatalog.GetFilteredList(lineTypeFilter, forceFilter);
    BOOST_CHECK( lineList.size()>0);



    std::string opt_lineWidthType = "velocitydriven";
    Float64 opt_resolution = 2350; //unused with velocity driven linewidth
    Float64 opt_velocityEmission = initVelocity;
    Float64 opt_velocityAbsorption = initVelocity;
    std::string opt_rules = "no";
    std::string opt_rigidity = "tplshape";
    std::string opt_calibrationPath= DATA_ROOT_DIR "LinemodelFitEstimateLeastSquareTestCase/calibration";



    //these tplcatalog related variables are unused here.
    CTemplateCatalog tplCatalog;
    Bool retValue = tplCatalog.Load( DATA_ROOT_DIR "LinemodelFitEstimateLeastSquareTestCase/ContinuumTemplates_simulm2016Extended_dustfree201702_1/" );
    BOOST_CHECK( retValue == true);
    TStringList tplCategories = TStringList { "galaxy" };

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
    CTemplateCatalog orthoTplCatalog = tplOrtho.getOrthogonalTplCatalog();


    CLineModelElementList model(*spectrum, *spectrum, orthoTplCatalog, tplCategories, opt_calibrationPath, lineList, opt_fittingmethod, opt_continuumcomponent, opt_lineWidthType, opt_resolution, opt_velocityEmission, opt_velocityAbsorption, opt_rules, opt_rigidity);

    bool tplratioInitRet = model.initTplratioCatalogs();
    BOOST_CHECK_MESSAGE( tplratioInitRet, "Unable to intialize tpl-ratio catalog");

    TFloat64Range lambdaRange = TFloat64Range( 3900.0, 12500.0 );
    CLineModelSolution modelSolution;
    Float64 contreest_iterations = 0;
    Bool enableLogging=true;
    Float64 merit = model.fit(z, lambdaRange, modelSolution, contreest_iterations, enableLogging);
    BOOST_TEST_MESSAGE( "Merit = " << merit );

    /*//debug:
    model.refreshModel();
    CSpectrum spcModel = model.GetModelSpectrum();
    FILE* f = fopen( "estimateleastsquarefast_model_dbg.txt", "w+" );
    for( Int32 t=0;t<spcModel.GetSampleCount();t++)
    {
        fprintf( f, "%f %e\n", spcModel.GetSpectralAxis()[t], spcModel.GetFluxAxis()[t]);
    }
    fclose( f );
    CLineModelSolution solution = model.GetModelSolution();
    //*/


    Float64 lstSqFast = model.getLeastSquareMeritFast();
    BOOST_TEST_MESSAGE( "rigidity=" << opt_rigidity << ", continuum=" << opt_continuumcomponent << ", Lst-sq-fast = " << lstSqFast );

    Float64 lstSq = model.getLeastSquareMerit(lambdaRange);
    BOOST_TEST_MESSAGE( "rigidity=" << opt_rigidity << ", continuum=" << opt_continuumcomponent << ", Lst-sq = " << lstSq );

    Float64 tol = 0.1; //for now, tol is very weak, TODO improve algos and lower this tolerance. Limiting method is linemodel-tplfit=fullmodel
    BOOST_CHECK_CLOSE_FRACTION( lstSq, lstSqFast, tol);

}

/***
 * test the linemodel tpl-shape for the estimation of the least-square value.
 * Warning as of 2017-03-27: using a hardcoded tplshape catalog: linecatalogs_tplshape_ExtendedTemplatesMarch2016_B13D_mod2
 *
 *
***/
BOOST_AUTO_TEST_CASE( LinemodelFit_EstimateLstSq_tplshape_pfsbatch6 )
{
    std::string spectrumPath = DATA_ROOT_DIR "LinemodelFitEstimateLeastSquareTestCase/10000663000008vacLine_TF.fits";
    std::string noisePath = DATA_ROOT_DIR "LinemodelFitEstimateLeastSquareTestCase/10000663000008vacLine_ErrF.fits";
    std::string linecatalogPath = DATA_ROOT_DIR "LinemodelFitEstimateLeastSquareTestCase/linecatalogamazedvacuum_B13D.txt";
    std::string linecatalogPath_NoLinesInLbdaRange = DATA_ROOT_DIR "LinemodelFitEstimateLeastSquareTestCase/linecatalogamazedvacuum_NoLinesInWavelengthRange.txt"; //only 1 line outside the wavelength range
    std::string linecatalogPath_Ha = DATA_ROOT_DIR "LinemodelFitEstimateLeastSquareTestCase/linecatalogamazedvacuum_Ha.txt"; //only 1 line in the wavelength range


    std::string opt_fittingmethod = "individual";
    Int32 lineTypeFilter = CRay::nType_Emission;
    Float64 initialVelocity = 100.0;
    Int32 forceFilter = CRay::nForce_Strong;

    Float64 z = 0.077164;
    std::string opt_continuumcomponent = "fromspectrum";

    //test the linemodel tplshape, without lines fitted (catalog has no lines in the lbda range)
    opt_continuumcomponent = "fromspectrum";
    checkLeastSquareFast(spectrumPath, noisePath, linecatalogPath_NoLinesInLbdaRange, opt_fittingmethod, opt_continuumcomponent, lineTypeFilter, forceFilter, initialVelocity, z);

    //test the linemodel tplshape, with only 1 line in the catalog
    opt_continuumcomponent = "fromspectrum";
    checkLeastSquareFast(spectrumPath, noisePath, linecatalogPath_Ha, opt_fittingmethod, opt_continuumcomponent, lineTypeFilter, forceFilter, initialVelocity, z);

    //test the linemodel tplshape
    opt_continuumcomponent = "fromspectrum";
    checkLeastSquareFast(spectrumPath, noisePath, linecatalogPath, opt_fittingmethod, opt_continuumcomponent, lineTypeFilter, forceFilter, initialVelocity, z);

    //test the fullmodel tplshape, without lines fitted (catalog has no lines in the lbda range)
    opt_continuumcomponent = "tplfit";
    checkLeastSquareFast(spectrumPath, noisePath, linecatalogPath_NoLinesInLbdaRange, opt_fittingmethod, opt_continuumcomponent, lineTypeFilter, forceFilter, initialVelocity, z);

    //test the fullmodel tplshape
    opt_continuumcomponent = "tplfit";
    checkLeastSquareFast(spectrumPath, noisePath, linecatalogPath, opt_fittingmethod, opt_continuumcomponent, lineTypeFilter, forceFilter, initialVelocity, z);

}

BOOST_AUTO_TEST_SUITE_END()
