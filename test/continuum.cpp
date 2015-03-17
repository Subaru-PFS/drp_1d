#include "continuum.h"

#include <epic/core/common/datatypes.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/continuum/median.h>
#include <epic/redshift/spectrum/io/fitswriter.h>

#include <epic/core/log/log.h>
#include <epic/core/log/consolehandler.h>
#include <math.h>

using namespace NSEpic;
using namespace std;

void CRedshiftContinuumTestCase::setUp()
{
}

void CRedshiftContinuumTestCase::tearDown()
{
}

#include <epic/redshift/spectrum/io/fitsreader.h>

void CRedshiftContinuumTestCase::Compute()
{

    // load continuum and associated simu ECN spectrum (Emission line + Continuum + Noise)
    CSpectrumIOFitsReader reader;
    CSpectrum s;
    CSpectrum s_continuumRef;
	
    Bool retVal = reader.Read( "../test/data/ContinuumTestCase/simu_ECN_continuum.fits", s_continuumRef );
    CPPUNIT_ASSERT( retVal == true );
    retVal = reader.Read( "../test/data/ContinuumTestCase/simu_ECN_all.fits", s );
    CPPUNIT_ASSERT( retVal == true );
    //Log.LogInfo("simu signals loaded");

    // Remove continuum 
    CContinuumMedian continuum;
    CSpectrumFluxAxis fluxAxisWithoutContinuumCalc;

    retVal = continuum.RemoveContinuum( s, fluxAxisWithoutContinuumCalc );
    CPPUNIT_ASSERT( retVal == true );
    //Log.LogInfo("continuum removed");

    // test the extracted continuum to be lower than a threshold all over the lambda range
    Float64 threshold = 0.05;
    CSpectrumFluxAxis fluxAxis = s.GetFluxAxis();
    fluxAxis.Subtract(fluxAxisWithoutContinuumCalc);
    //Log.LogInfo("subtraction done");
    Float64 er = fluxAxis.ComputeRMSDiff(s_continuumRef.GetFluxAxis());
    //Log.LogInfo("Continuum rms error: %f", er);
    CPPUNIT_ASSERT(er < threshold);
     

    CSpectrumFluxAxis& sfluxAxisPtr = s.GetFluxAxis();
    sfluxAxisPtr = fluxAxis;
    CSpectrumIOFitsWriter writer;
    retVal = writer.Write( "../test/data/ContinuumTestCase/simu_ECN_continuum_calc.fits", s );
}




