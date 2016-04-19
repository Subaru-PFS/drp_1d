

#include <epic/core/common/datatypes.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/continuum/irregularsamplingmedian.h>
#include <epic/redshift/spectrum/io/fitswriter.h>
#include <epic/core/log/log.h>
#include <epic/core/log/consolehandler.h>
#include <epic/redshift/spectrum/io/fitsreader.h>


#include <epic/redshift/continuum/indexes.h>

#include <math.h>
#include <boost/test/unit_test.hpp>

using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(Continuum)

BOOST_AUTO_TEST_CASE(Compute)
{
    // load continuum and associated simu ECN spectrum (Emission line + Continuum + Noise)
    CSpectrumIOFitsReader reader;
    CSpectrum s;
    CSpectrum s_continuumRef;

    Bool retVal = reader.Read( "../test/data/ContinuumTestCase/simu_ECN_continuum.fits", s_continuumRef );
    BOOST_CHECK( retVal == true );
    retVal = reader.Read( "../test/data/ContinuumTestCase/simu_ECN_all.fits", s );
    BOOST_CHECK( retVal == true );

    // Remove continuum
    CContinuumIrregularSamplingMedian continuum;
    CSpectrumFluxAxis fluxAxisWithoutContinuumCalc;

    retVal = continuum.RemoveContinuum( s, fluxAxisWithoutContinuumCalc );
    BOOST_CHECK( retVal == true );
    // test the extracted continuum to be lower than a threshold all over the lambda range
    Float64 threshold = 0.05;
    CSpectrumFluxAxis fluxAxis = s.GetFluxAxis();
    fluxAxis.Subtract(fluxAxisWithoutContinuumCalc);
    Float64 er = fluxAxis.ComputeRMSDiff(s_continuumRef.GetFluxAxis());
    BOOST_CHECK(er < threshold);
}

BOOST_AUTO_TEST_CASE(ContinuumIndexes)
{
    // load spectrum
    CSpectrumIOFitsReader reader;
    CSpectrum s;

    Bool retVal = reader.Read( "../test/data/ContinuumTestCase/spectrum_tpl_vvds-reddestdataExtensionData.txt_TF.fits", s );
    BOOST_CHECK( retVal == true );

    //compute continuum indexes
    CContinuumIndexes continuumIndexes;
    CContinuumIndexes::TContinuumIndexList indexesList = continuumIndexes.getIndexes( s, 0.0 );

    // test the extracted indexes
    //Lya index
    BOOST_CHECK_SMALL( abs(indexesList[0].Color - 0.34), 0.1 );
    BOOST_CHECK_SMALL( abs(indexesList[0].Break - 0.0), 0.1 );
    //OII index
    BOOST_CHECK_SMALL( abs(indexesList[1].Color - 0.9), 0.1 );
    BOOST_CHECK_SMALL( abs(indexesList[1].Break - -0.5), 0.1 );
    //OIII index
    BOOST_CHECK_SMALL( abs(indexesList[2].Color - 0.1), 0.1 );
    BOOST_CHECK_SMALL( abs(indexesList[2].Break - 0.4), 0.1 );
    //Halpha index
    BOOST_CHECK( std::isnan(indexesList[3].Color) );
    BOOST_CHECK( std::isnan(indexesList[3].Break) );


}


BOOST_AUTO_TEST_SUITE_END()




