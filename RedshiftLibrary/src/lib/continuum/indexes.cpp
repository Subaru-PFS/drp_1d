
#include "RedshiftLibrary/continuum/indexes.h"

#include <iostream>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/log/log.h"


using namespace NSEpic;
using namespace std;

/**
 * Calculates the continuum indexes for a given spectrum
 * See: Lilly et al. [arXiv:astro-ph/9507012]
 */
CContinuumIndexes::TContinuumIndexList CContinuumIndexes::getIndexes(const CSpectrum& spectrum , Float64 z)
{
    TFloat64RangeList restLambdaRanges_A;
    TFloat64RangeList restLambdaRanges_B;

    //add the ranges to be processed
    //usually corresponds to 0.858, 0.966 for A, and 1.073, 1.127 for B
    restLambdaRanges_A.push_back(TFloat64Range( 1043.0, 1174.0 )); //Lya, A
    restLambdaRanges_B.push_back(TFloat64Range( 1304.0, 1369.0 )); //Lya, B

    restLambdaRanges_A.push_back(TFloat64Range( 3200.0, 3600.0 )); //OII, A
    restLambdaRanges_B.push_back(TFloat64Range( 4000.0, 4200.0 )); //OII, B

    restLambdaRanges_A.push_back(TFloat64Range( 4290.0, 4830.0 )); //OIII, A
    restLambdaRanges_B.push_back(TFloat64Range( 5365.0, 5635.0 )); //OIII, B

    restLambdaRanges_A.push_back(TFloat64Range( 5632.0, 6341.0 )); //Halpha, A
    restLambdaRanges_B.push_back(TFloat64Range( 7043.0, 7397.6 )); //Halpha, B

    restLambdaRanges_A.push_back(TFloat64Range( 1329.0, 1497.0 )); //CIV, A
    restLambdaRanges_B.push_back(TFloat64Range( 1663.0, 1746.0 )); //CIV, B

    restLambdaRanges_A.push_back(TFloat64Range( 1637.0, 1843.0 )); //CIII, A
    restLambdaRanges_B.push_back(TFloat64Range( 2047.0, 2150.0 )); //CIII, B


    Int32 nIndexes = restLambdaRanges_A.size();
    CContinuumIndexes::TContinuumIndexList indexesList;
    for(Int32 i=0; i<nIndexes; i++)
    {
        TFloat64Range rangeA = TFloat64Range( restLambdaRanges_A[i].GetBegin()*(z+1), restLambdaRanges_A[i].GetEnd()*(z+1) );
        TFloat64Range rangeB = TFloat64Range( restLambdaRanges_B[i].GetBegin()*(z+1), restLambdaRanges_B[i].GetEnd()*(z+1) );
        Float64 Fa = NAN;
        Float64 std = NAN;
        bool retA = spectrum.GetMeanAndStdFluxInRange( rangeA, Fa, std );
        Float64 Fb = NAN;
        bool retB = spectrum.GetMeanAndStdFluxInRange( rangeB, Fb, std );

        Float64 a;
        Float64 b;
        bool retC = spectrum.GetLinearRegInRange( rangeA,  a, b);
        Float64 wlCenterB = (rangeB.GetEnd()+rangeB.GetBegin())/2.0;
        Float64 Fc = wlCenterB*a + b;

        SContinuumIndex sci;
        sci.Break = NAN;
        sci.Color = NAN;

        if(Fb > 0.0 && Fa > 0.0 && retA && retB)
        {
            sci.Color = -2.5*log10(Fa/Fb);
        }else if(Fb<=0.0 && Fa > 0.0 && retA && retB)
        {
            sci.Color = -6.0;
        }
        if(Fc > 0.0 && Fb > 0.0 && retB && retC)
        {
            sci.Break = -2.5*log10(Fb/Fc);
        }else if(Fc<=0.0 && Fb > 0.0 && retB && retC)
        {
            sci.Break = -6.0;
        }

        indexesList.push_back(sci);
    }

    return indexesList;
}


/**
 * Calculates the spectrum and continuum std for a given spectrum and continuum
 *
 */
CContinuumIndexes::SContinuumRelevance CContinuumIndexes::getRelevance(const CSpectrum& spectrum, const CSpectrum& continuum)
{
    Float64 Fa = NAN;
    Float64 stdS = NAN;
    Float64 stdC = NAN;

    spectrum.GetMeanAndStdFluxInRange( spectrum.GetLambdaRange(), Fa, stdS );
    continuum.GetMeanAndStdFluxInRange( continuum.GetLambdaRange(), Fa, stdC );

    SContinuumRelevance relevance;
    relevance.StdSpectrum = stdS;
    relevance.StdContinuum = stdC;
    return relevance;
}

