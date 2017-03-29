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

BOOST_AUTO_TEST_SUITE(LinemodelFit_EstimateLeastSquare)

/***
 * Check that the least-square-fast provides the same value as the least-square method
 *
**/
void checkLeastSquareFast(std::string spectrumPath, std::string noisePath, std::string linecatalogPath, std::string opt_fittingmethod, std::string opt_continuumcomponent, Int32 lineTypeFilter, Int32 forceFilter, Float64 initVelocity, Float64 z)
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


    // get continuum from Median in case of opt_continuumcomponent==fromspectrum
    CContinuumIrregularSamplingMedian continuum;
    CSpectrumFluxAxis fluxAxisWithoutContinuumCalc;
    Int32 retValContinuumEstimation = continuum.RemoveContinuum( spectrum, fluxAxisWithoutContinuumCalc );
    BOOST_CHECK( retValContinuumEstimation);
    CSpectrum spectrumContinuum = spectrum;
    CSpectrumFluxAxis& continuumFluxAxis = spectrumContinuum.GetFluxAxis();
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
    std::string opt_calibrationPath= "../test/data/LinemodelFitEstimateLeastSquareTestCase/calibration";



    //these tplcatalog related variables are unused here.
    CTemplateCatalog tplCatalog;
    Bool retValue = tplCatalog.Load( "../test/data/LinemodelFitEstimateLeastSquareTestCase/ContinuumTemplates_simulm2016Extended_dustfree201702_1/" );
    BOOST_CHECK( retValue == true);
    TStringList tplCategories;

    CLineModelElementList model(spectrum, spectrumContinuum, tplCatalog, tplCategories, opt_calibrationPath, lineList, opt_fittingmethod, opt_continuumcomponent, opt_lineWidthType, opt_resolution, opt_velocityEmission, opt_velocityAbsorption, opt_rules, opt_rigidity);
    TFloat64Range lambdaRange = TFloat64Range( 3800.0, 12600.0 );
    CLineModelResult::SLineModelSolution modelSolution;
    Float64 merit = model.fit(z, lambdaRange, modelSolution);
    //BOOST_TEST_MESSAGE( "Merit = " << merit );

    Float64 lstSqFast = model.getLeastSquareMeritFast();
    BOOST_TEST_MESSAGE( "rigidity=" << opt_rigidity << ", continuum=" << opt_continuumcomponent << ", Lst-sq-fast = " << lstSqFast );

    Float64 lstSq = model.getLeastSquareMerit(lambdaRange);
    BOOST_TEST_MESSAGE( "rigidity=" << opt_rigidity << ", continuum=" << opt_continuumcomponent << ", Lst-sq = " << lstSq );

    BOOST_CHECK_SMALL( fabs(lstSq - lstSqFast), 1e-5); //todo: check with more precision ?

}

/***
 * test the linemodel tpl-shape for the estimation of the least-square value.
 * Warning as of 2017-03-27: using a hardcoded tplshape catalog: linecatalogs_tplshape_ExtendedTemplatesMarch2016_B13D_mod2
 *
 *
***/
BOOST_AUTO_TEST_CASE( LinemodelFit_EstimateLstSq_tplshape_pfsbatch6 )
{
    std::string spectrumPath = "../test/data/LinemodelFitEstimateLeastSquareTestCase/10000663000008vacLine_TF.fits";
    std::string noisePath = "../test/data/LinemodelFitEstimateLeastSquareTestCase/10000663000008vacLine_ErrF.fits";
    std::string linecatalogPath = "../test/data/LinemodelFitEstimateLeastSquareTestCase/linecatalogamazedvacuum_B13D.txt";

    std::string opt_fittingmethod = "individual";
    Int32 lineTypeFilter = CRay::nType_Emission;
    Float64 initialVelocity = 100.0;
    Int32 forceFilter = CRay::nForce_Strong;

    Float64 z = 0.954114;

    //test the linemodel tplshape
    std::string opt_continuumcomponent = "fromspectrum";
    checkLeastSquareFast(spectrumPath, noisePath, linecatalogPath, opt_fittingmethod, opt_continuumcomponent, lineTypeFilter, forceFilter, initialVelocity, z);

    //test the fullmodel tplshape
    opt_continuumcomponent = "tplfit";
    checkLeastSquareFast(spectrumPath, noisePath, linecatalogPath, opt_fittingmethod, opt_continuumcomponent, lineTypeFilter, forceFilter, initialVelocity, z);

}

BOOST_AUTO_TEST_SUITE_END()
