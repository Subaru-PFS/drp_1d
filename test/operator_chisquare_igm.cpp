#include <epic/core/common/datatypes.h>
#include <epic/redshift/operator/chisquare2.h>
#include <epic/redshift/operator/chisquareresult.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/spectrum/io/fitsreader.h>
#include <epic/redshift/spectrum/io/genericreader.h>
#include <epic/redshift/extremum/extremum.h>
#include <epic/redshift/continuum/irregularsamplingmedian.h>
#include <epic/redshift/noise/fromfile.h>

#include <fstream>
#include <boost/math/special_functions.hpp>

#include <boost/test/unit_test.hpp>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(Operator_Chisquare_igm)


void UtilChisquareTestFit( const char* spectraPath, const char* noisePath, const char* tplPath, bool disableMask, const Float64 targetFittedAmplitude )
{
    Bool retVal;
    CSpectrum s;
    CTemplate t;

    Float64 z = 0.0;

    // Load spectrum and templates
    CSpectrumIOGenericReader reader;
    retVal = reader.Read( spectraPath, s );
    BOOST_CHECK( retVal );

    if( noisePath )
    {
        CNoiseFromFile noise;
        noise.SetNoiseFilePath( noisePath );
        noise.AddNoise( s );
    }

    retVal = reader.Read( tplPath, t );
    BOOST_CHECK( retVal );

    Float64 redshiftDelta = 0.0001;
    TFloat64List redshifts = TFloat64Range( 0.0, 0.0 ).SpreadOver( redshiftDelta );

    //building the mask
    Int32 sampleCount = s.GetFluxAxis().GetSamplesCount();
    std::vector<CMask> additional_spcMasks;
    CMask spcMask = Mask();
    spcMask.SetSize(sampleCount);

    //TBD: Remove all additional masks for this test...
    Float64 mean = -1.0;
    Float64 std = -1.0;
    s.GetMeanAndStdFluxInRange(TFloat64Range( 920, 9000 ), mean, std);
    Float64 thres = mean*0.5;

    for(Int32 km=0; km<spcMask.GetMasksCount(); km++)
    {
        if(!disableMask)
        {
            if(s.GetFluxAxis()[km]>thres)
            {
                spcMask[km]=0.0;
            }else
            {
                spcMask[km]=1.0;
            }
        }else{
            spcMask[km]=1.0;
        }
    }
    additional_spcMasks.push_back(spcMask);
    // end TBD


    std::string calibrationPath = "../test/data/Operator_Chisquare_igmTestCase/calibration";

    COperatorChiSquare2 chi(calibrationPath);
    auto r = std::dynamic_pointer_cast<CChisquareResult>( chi.Compute( s, t, TFloat64Range( 920, 9000 ), redshifts, 1.0, additional_spcMasks, "precomputedfinegrid", 0 ) );
    BOOST_CHECK( r != NULL );

    Float64 fit_amplitude = r->FitAmplitude[0];
    BOOST_CHECK_CLOSE_FRACTION( targetFittedAmplitude, fit_amplitude, 0.1 );
//    CChisquareResult referenceResult;

//    std::ifstream input( resultPath );
//    BOOST_CHECK( input.is_open() );

//    referenceResult.Load( input );

//    for( Int32 i=0; i<referenceResult.ChiSquare.size(); i++ )
//    {
//        if( boost::math::isnan( referenceResult.ChiSquare[i] ) )
//        {
//            BOOST_CHECK( boost::math::isnan( r->ChiSquare[i] ) );
//        }
//        else
//        {
//            BOOST_CHECK_CLOSE_FRACTION( referenceResult.ChiSquare[i], r->ChiSquare[i], 0.00001 );
//        }


//        BOOST_CHECK_CLOSE_FRACTION( referenceResult.Redshifts[i], r->Redshifts[i], 0.00001 );
//        BOOST_CHECK_CLOSE_FRACTION( referenceResult.Overlap[i], r->Overlap[i], 0.00001 );
//    }

}

BOOST_AUTO_TEST_CASE(ChisquareMaskTest)
{
    //reference test
    UtilChisquareTestFit( "../test/data/Operator_Chisquare_igmTestCase/fits_chisquare_igm_cstflux/spc_synth_constant1p0_meiksin-z3c5.fits",
                            "../test/data/Operator_Chisquare_igmTestCase/fits_chisquare_igm_cstflux/spc_synth_constant0p5_ErrF.fits",
                            "../test/data/Operator_Chisquare_igmTestCase/fits_chisquare_igm_cstflux/template_constantflux.dat", //constant template f=1.0 for all lambda
                            true,
                            2.0);
}


BOOST_AUTO_TEST_SUITE_END()
