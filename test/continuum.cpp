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

    // load continuum and associated simu spectrum
    CSpectrumIOFitsReader reader;
    CSpectrum s;
    CSpectrum s_continuum;
	
    Bool retVal = reader.Read( "../test/data/simu_ECN_continuum_20150312B.fits", s_continuum );
    CPPUNIT_ASSERT( retVal == true );
    retVal = reader.Read( "../test/data/simu_ECN_all_20150312B.fits", s );
    CPPUNIT_ASSERT( retVal == true );

    // Remove continuum 
    CContinuumMedian continuum;
    CSpectrumFluxAxis fluxAxisWithoutContinuum;

    retVal = continuum.RemoveContinuum( s, fluxAxisWithoutContinuum );
    CSpectrumFluxAxis& fa = s.GetFluxAxis();
     

    //retVal = s.RemoveContinuum<CContinuumMedian>();
    CPPUNIT_ASSERT( retVal == true );

    // test the extracted continuum to be lower than a threshold all over the lambda range
    Float64 threshold = 0.05;
    CSpectrumFluxAxis refContinuumFluxAxis;
    refContinuumFluxAxis = s_continuum.GetFluxAxis();
    Int32 N = refContinuumFluxAxis.GetSamplesCount();
    Float64 weight = (Float64)N;
    Float64 er2=0.f;

    for( Int32 i=0; i<N; i++ )
    {
        er2 += (refContinuumFluxAxis[i] - (fa[i] - fluxAxisWithoutContinuum[i]))*(refContinuumFluxAxis[i] - (fa[i] - fluxAxisWithoutContinuum[i]))/weight;
    }
    Float64 er = sqrt(er2);
    Log.LogInfo("Continuum rms error: %.5f", er);
    CPPUNIT_ASSERT(er < threshold);
    
    fa = fluxAxisWithoutContinuum; 
    CSpectrumIOFitsWriter writer;
    retVal = writer.Write( "../test/data/simu_ECN_continuum_20150312B_calc.fits", s );
}




