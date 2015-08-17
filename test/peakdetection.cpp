#include "peakdetection.h"

#include <epic/redshift/operator/peakdetection.h>
#include <epic/redshift/operator/peakdetectionresult.h>
#include <epic/core/common/datatypes.h>
#include <epic/core/common/constref.h>
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

    CPPUNIT_ASSERT_MESSAGE(  "load fits", retVal == true);

    TLambdaRange lambdaRange = s.GetLambdaRange();
    CPeakDetection detection(500.0, 15, 1, 0);
    CConstRef<CPeakDetectionResult> peakDetectionResult = detection.Compute( s, lambdaRange); //using winsize=500 and cut=15 so that 3 only peaks are detected in the test signal for sure
    CPPUNIT_ASSERT_MESSAGE( "compute detection" , retVal == true );
    const TInt32RangeList& resPeaks = peakDetectionResult->PeakList;


    Float64 peakxpos[] = {1000, 5000, 8000};
    Float64 peaksigmas[] = {40.0, 10.0, 45.0};
    UInt32 n=4;

    // test number of peaks
    CPPUNIT_ASSERT(n == resPeaks.size());

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
        CPPUNIT_ASSERT_DOUBLES_EQUAL(infRef, infCalc, toleranceXPos);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(supRef, supCalc, toleranceXPos);
    }

}




