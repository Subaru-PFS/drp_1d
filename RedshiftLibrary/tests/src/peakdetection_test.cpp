#include <RedshiftLibrary/operator/peakdetection.h>
#include <RedshiftLibrary/operator/peakdetectionresult.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/io/fitsreader.h>

#include <math.h>
#include <boost/test/unit_test.hpp>
#include "test-config.h"

using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(PeakDetection)

BOOST_AUTO_TEST_CASE(Compute)
{

    // load simu spectrum (emission lines + Noise)
    CSpectrumIOFitsReader reader;
    CSpectrum s;
	
    Bool retVal = reader.Read( DATA_ROOT_DIR "PeakDetectionTestCase/peakdetection_simu.fits", std::shared_ptr<CSpectrum>(&s) );

    BOOST_CHECK( retVal == true);

    TLambdaRange lambdaRange = s.GetLambdaRange();
    CPeakDetection detection(500.0, 15, 1, 0);
    auto peakDetectionResult = detection.Compute( s, lambdaRange); //using winsize=500 and cut=15 so that 3 only peaks are detected in the test signal for sure
    BOOST_CHECK( retVal == true );
    const TInt32RangeList& resPeaks = peakDetectionResult->PeakList;


    Float64 peakxpos[] = {1000, 5000, 8000};
    Float64 peaksigmas[] = {40.0, 10.0, 45.0};
    UInt32 n=4;

    // test number of peaks
    BOOST_CHECK(n == resPeaks.size());

    // test peak xpos
    Float64 toleranceXPos = 15.f;
    for(int i=0; i<n-1; i++){
        Float64 fwhm =  2*sqrt(2*log(2))*peaksigmas[i];
        Float64 infRef = peakxpos[i]-fwhm/2.0*1.5;
        infRef = max((Float64)infRef, 0.0);
        Float64 supRef = peakxpos[i]+fwhm/2.0*1.5;
        supRef = min((Float64)supRef, (Float64)s.GetFluxAxis().GetSamplesCount());

        Float64 infCalc = resPeaks[i].GetBegin();
        Float64 supCalc = resPeaks[i].GetEnd();
        BOOST_CHECK_CLOSE_FRACTION(infRef, infCalc, toleranceXPos);
        BOOST_CHECK_CLOSE_FRACTION(supRef, supCalc, toleranceXPos);
    }
}

BOOST_AUTO_TEST_SUITE_END()



