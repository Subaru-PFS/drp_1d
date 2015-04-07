#include "peakdetection.h"

#include <epic/redshift/peak/detection.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/spectrum/io/fitsreader.h>
#include <math.h>

using namespace NSEpicTest;
using namespace NSEpic;
using namespace std;

void CRedshiftPeakDetectionTestCase::setUp()
{
}

void CRedshiftPeakDetectionTestCase::tearDown()
{
}


void CRedshiftPeakDetectionTestCase::Compute()
{

    // load simu spectrum (emission lines + Noise)
    CSpectrumIOFitsReader reader;
    CSpectrum s;
	
    Bool retVal = reader.Read( "../test/data/PeakDetectionTestCase/peakdetection_simu.fits", s );


    TLambdaRange lambdaRange = s.GetLambdaRange();
    CPeakDetection detection;
    retVal = detection.Compute( s, lambdaRange );
    CPPUNIT_ASSERT_MESSAGE( "Failed to detect peak" , retVal == true );
    TInt32RangeList resPeaks = detection.GetResults();

    Float64 noiseSigma = 0.1f;
    Float64 peakxpos[] = {1000, 2500, 3500, 3700, 5000, 8000, 7000, 9950};
    Float64 peaksigmas[] = {40.0, 100.0, 40.0, 45.0, 10.0, 45.0, 30.0, 25.0};
    UInt32 n=8;

    // test number of peaks
    CPPUNIT_ASSERT(n = resPeaks.size());

    // test peak xpos
    Float64 toleranceXPos = 5.f; //todo: refine with noiseSigma... (noiseSigma)

    for(int i=0; i<n; i++){
        Float64 infRef = peakxpos[i]-peaksigmas[i];
        Float64 supRef = peakxpos[i]+peaksigmas[i];

        Float64 infCalc = resPeaks[i].GetBegin();
        Float64 supCalc = resPeaks[i].GetEnd();
        CPPUNIT_ASSERT_DOUBLES_EQUAL(infRef, infCalc, toleranceXPos);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(supRef, supCalc, toleranceXPos);
    }

}




