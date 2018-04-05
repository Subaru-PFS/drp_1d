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

BOOST_AUTO_TEST_SUITE(LinemodelMultilineSupport)

/**
 * @brief checkSupport
 * @param linecatalogPath
*/
void checkSupport(std::string linecatalogPath,
                  std::string spcPath,
                  std::string nPath,
                  const TFloat64Range& lambdaRange,
                  Float64 redshift,
                  bool targetSupportSizeZero,
                  bool targetOustideLambdaRangeTrue)
{
    //some input params
    std::string spectrumPath = spcPath;//DATA_ROOT_DIR "LinemodelProfileTestCase/signalnoise_4lines_sig400_6000A_10000A.fits";       //unused, only needed for linemodel initialization
    std::string noisePath    = nPath;//DATA_ROOT_DIR "LinemodelProfileTestCase/noise_4lines_sig400_6000A_10000A.fits";     //unused, only needed for linemodel initialization
    Int32 lineTypeFilter = CRay::nType_Emission;
    Int32 forceFilter = -1;
    std::string opt_fittingmethod = "ones"; //all the elements amplitudes set to 1.0
    Float64 z = 0.0;

    // load spectrum
    CSpectrumIOFitsReader reader;
    CSpectrum spectrum;

    Bool retVal = reader.Read( spectrumPath.c_str(), std::shared_ptr<CSpectrum>(&spectrum));
    BOOST_CHECK_MESSAGE( retVal == true, "check load flux spectrum");
    CNoiseFromFile noise;
    retVal = noise.SetNoiseFilePath( noisePath.c_str() );
    BOOST_CHECK_MESSAGE( retVal == true, "check load noise spectrum");
    retVal = noise.AddNoise( spectrum ) ;
    BOOST_CHECK_MESSAGE( retVal == true, "check add noise spectrum");


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
    BOOST_CHECK_MESSAGE( rValue == true, "check load catalog");
    CRayCatalog::TRayVector lineList = lineCatalog.GetFilteredList(lineTypeFilter, forceFilter);
    BOOST_CHECK_MESSAGE( lineList.size()>0, "check loaded catalog not empty");


    std::string opt_continuumcomponent = "fromspectrum";
    std::string opt_lineWidthType = "velocitydriven";

    Float64 opt_resolution = 2350; //unused with velocity driven linewidth
    Float64 opt_velocityEmission = 300.0; //IMPORTANT value here in order to set the linewidth
    Float64 opt_velocityAbsorption = 300.0; //unused
    std::string opt_rules = "no";
    std::string opt_rigidity = "rules";
    std::string unused_calibrationPath="";


    //these tplcatalog related variables are unused here.
    CTemplateCatalog tplCatalog;
    Bool retValue = tplCatalog.Load( DATA_ROOT_DIR "templatecatalog/" );
    TStringList tplCategories;

    CLineModelElementList model(spectrum, spectrumContinuum, tplCatalog, tplCategories, unused_calibrationPath, lineList, opt_fittingmethod, opt_continuumcomponent, opt_lineWidthType, opt_resolution, opt_velocityEmission, opt_velocityAbsorption, opt_rules, opt_rigidity);


    //initialize the model spectrum
    const CSpectrumSpectralAxis& spectralAxis = model.m_SpectrumModel->GetSpectralAxis();
    //prepare model support
    model.initModelAtZ(redshift, lambdaRange, spectralAxis);

    BOOST_CHECK_MESSAGE( model.m_Elements.size() == 1, "Check -only 1 element- failed" );
    Int32 iElt = 0; //only 1 element
    BOOST_CHECK_MESSAGE( model.m_Elements[iElt]->GetRays().size() == 1, "Check -only 1 sub-element- failed" );
    Int32 lineIdx = 0;  //only 1 line in this model
    //TInt32Range support = model.m_Elements[iElt]->getTheoreticalSupportSubElt(lineIdx);
    TInt32Range support = model.m_Elements[iElt]->getSupportSubElt(lineIdx);
    bool IsOutsideLambdaRange = model.m_Elements[iElt]->IsOutsideLambdaRange(lineIdx);


    if( !targetSupportSizeZero )
    {
        BOOST_CHECK_MESSAGE( support.GetEnd()-support.GetBegin()+1 > 0, "Check support size>0 failed");

    }else{
        BOOST_CHECK_MESSAGE( support.GetEnd()-support.GetBegin()+1 <= 0, "Check support size==0 failed");
    }
    if( !targetOustideLambdaRangeTrue )
    {
        BOOST_CHECK_MESSAGE( IsOutsideLambdaRange==false, "Check .IsOutsideLambdaRange failed");

    }else{
        BOOST_CHECK_MESSAGE( IsOutsideLambdaRange==true, "Check .IsOutsideLambdaRange failed");
    }

}

/**
  *Test support 1EL fully included in the wavelength range: Redshift 0
  *
**/
BOOST_AUTO_TEST_CASE( Linemodel_multiline_support_1EL_center )
{
    std::string linecatalogPath = DATA_ROOT_DIR "LinemodelSupportTestCase/linecatalog_test_linemodel_support_1EL_15000.txt";

    std::string spectrumPath = DATA_ROOT_DIR "LinemodelSupportTestCase/spc_100_114181174_axespc_F.fits";
    std::string noisePath    = DATA_ROOT_DIR "LinemodelSupportTestCase/spc_100_114181174_axespc_ErrF.fits";

    TFloat64Range lambdaRange = TFloat64Range( 13000.0, 18000.0 );
    Float64 redshift = 0.0;

    bool targetSupportSizeZero = false;  //support should NOT be zero because the line is fully inside the lambda range
    bool targetOustideLambdaRangeTrue = false;
    checkSupport(linecatalogPath, spectrumPath, noisePath, lambdaRange, redshift, targetSupportSizeZero, targetOustideLambdaRangeTrue);
}

/**
  *Test support 1EL fully outside the wavelength range: Redshift 0
  *
**/
BOOST_AUTO_TEST_CASE( Linemodel_multiline_support_1EL_outside )
{
    std::string linecatalogPath = DATA_ROOT_DIR "LinemodelSupportTestCase/linecatalog_test_linemodel_support_1EL_15000.txt";

    std::string spectrumPath = DATA_ROOT_DIR "LinemodelSupportTestCase/spc_100_114181174_axespc_F.fits";
    std::string noisePath    = DATA_ROOT_DIR "LinemodelSupportTestCase/spc_100_114181174_axespc_ErrF.fits";

    TFloat64Range lambdaRange = TFloat64Range( 17000.0, 18000.0 );
    Float64 redshift = 0.0;

    bool targetSupportSizeZero = true; //support should be null because the line is outside lambda range
    bool targetOustideLambdaRangeTrue = true;
    checkSupport(linecatalogPath, spectrumPath, noisePath, lambdaRange, redshift, targetSupportSizeZero, targetOustideLambdaRangeTrue);
}

/**
  *Test support 1EL partially outside the wavelength range: 50% overlap, Redshift 0
  *
**/
BOOST_AUTO_TEST_CASE( Linemodel_multiline_support_1EL_borderOverlap50percent )
{
    std::string linecatalogPath = DATA_ROOT_DIR "LinemodelSupportTestCase/linecatalog_test_linemodel_support_1EL_15000.txt";

    std::string spectrumPath = DATA_ROOT_DIR "LinemodelSupportTestCase/spc_100_114181174_axespc_F.fits";
    std::string noisePath    = DATA_ROOT_DIR "LinemodelSupportTestCase/spc_100_114181174_axespc_ErrF.fits";

    TFloat64Range lambdaRange = TFloat64Range( 15000.0, 18000.0 );
    Float64 redshift = 0.0;

    bool targetSupportSizeZero = false; //support should NOT be zero because the line is 50% inside the lambda range
    bool targetOustideLambdaRangeTrue = false;
    checkSupport(linecatalogPath, spectrumPath, noisePath, lambdaRange, redshift, targetSupportSizeZero, targetOustideLambdaRangeTrue);
}

/**
  *Test support 1EL partially outside the wavelength range: 12% overlap, Redshift 0
  *
**/
BOOST_AUTO_TEST_CASE( Linemodel_multiline_support_1EL_borderOverlap12percent )
{
    std::string linecatalogPath = DATA_ROOT_DIR "LinemodelSupportTestCase/linecatalog_test_linemodel_support_1EL_15000.txt";

    std::string spectrumPath = DATA_ROOT_DIR "LinemodelSupportTestCase/spc_100_114181174_axespc_F.fits";
    std::string noisePath    = DATA_ROOT_DIR "LinemodelSupportTestCase/spc_100_114181174_axespc_ErrF.fits";

    TFloat64Range lambdaRange = TFloat64Range( 15045.6, 18000.0 );
    Float64 redshift = 0.0;

    bool targetSupportSizeZero = false; //support should NOT be zero because the line is 12% inside the lambda range
    bool targetOustideLambdaRangeTrue = true; //outsidelambdarange should be true because the line is 1% inside the lambda range only. The threshold is set to 33% in the code.
    checkSupport(linecatalogPath, spectrumPath, noisePath, lambdaRange, redshift, targetSupportSizeZero, targetOustideLambdaRangeTrue);
}

/**
  *Test support 1EL partially outside the wavelength range: 1% overlap, Redshift 0
  *
**/
BOOST_AUTO_TEST_CASE( Linemodel_multiline_support_1EL_borderOverlap1percent )
{
    std::string linecatalogPath = DATA_ROOT_DIR "LinemodelSupportTestCase/linecatalog_test_linemodel_support_1EL_15000.txt";

    std::string spectrumPath = DATA_ROOT_DIR "LinemodelSupportTestCase/spc_100_114181174_axespc_F.fits";
    std::string noisePath    = DATA_ROOT_DIR "LinemodelSupportTestCase/spc_100_114181174_axespc_ErrF.fits";

    TFloat64Range lambdaRange = TFloat64Range( 15058.8, 18000.0 );
    Float64 redshift = 0.0;

    bool targetSupportSizeZero = false; //support should NOT be zero
    bool targetOustideLambdaRangeTrue = true; //outsidelambdarange should be true because the line is 1% inside the lambda range only. The threshold is set to 33% in the code.
    checkSupport(linecatalogPath, spectrumPath, noisePath, lambdaRange, redshift, targetSupportSizeZero, targetOustideLambdaRangeTrue);
}


/**
  *Test support 1EL fully included in the wavelength range: Redshift 1.5102
  *
**/
BOOST_AUTO_TEST_CASE( Linemodel_multiline_support_1EL_borderOverlap_withNonZeroRedshift )
{
    std::string linecatalogPath = DATA_ROOT_DIR "LinemodelSupportTestCase/linecatalog_test_linemodel_support_1EL_OIIIb.txt";

    std::string spectrumPath = DATA_ROOT_DIR "LinemodelSupportTestCase/spc_100_114181174_axespc_F.fits";
    std::string noisePath    = DATA_ROOT_DIR "LinemodelSupportTestCase/spc_100_114181174_axespc_ErrF.fits";

    TFloat64Range lambdaRange = TFloat64Range( 0.0, 20000.0 );
    Float64 redshift = 1.5102;

    bool targetSupportSizeZero = false;  //support should NOT be zero
    bool targetOustideLambdaRangeTrue = false; //
    checkSupport(linecatalogPath, spectrumPath, noisePath, lambdaRange, redshift, targetSupportSizeZero, targetOustideLambdaRangeTrue);
}

BOOST_AUTO_TEST_SUITE_END()
