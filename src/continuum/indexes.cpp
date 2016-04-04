
#include <epic/redshift/continuum/indexes.h>

#include <iostream>

#include <epic/core/common/datatypes.h>
#include <epic/core/log/log.h>


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
    restLambdaRanges_A.push_back(TFloat64Range( 1043.0, 1174.0 )); //Lya, A
    restLambdaRanges_B.push_back(TFloat64Range( 1304.0, 1369.0 )); //Lya, B

    restLambdaRanges_A.push_back(TFloat64Range( 3200.0, 3600.0 )); //OII, A
    restLambdaRanges_B.push_back(TFloat64Range( 4000.0, 4200.0 )); //OII, B

    restLambdaRanges_A.push_back(TFloat64Range( 4290.0, 4830.0 )); //OIII, A
    restLambdaRanges_B.push_back(TFloat64Range( 5365.0, 5635.0 )); //OIII, B

    restLambdaRanges_A.push_back(TFloat64Range( 5632.0, 6341.0 )); //Halpha, A
    restLambdaRanges_B.push_back(TFloat64Range( 7043.0, 7397.6 )); //Halpha, B

    Int32 nIndexes = restLambdaRanges_A.size();
    CContinuumIndexes::TContinuumIndexList indexesList;
    for(Int32 i=0; i<nIndexes; i++)
    {
        TFloat64Range rangeA = TFloat64Range( restLambdaRanges_A[i].GetBegin()*(z+1), restLambdaRanges_A[i].GetEnd()*(z+1) );
        TFloat64Range rangeB = TFloat64Range( restLambdaRanges_B[i].GetBegin()*(z+1), restLambdaRanges_B[i].GetEnd()*(z+1) );
        Float64 Fa = NAN;
        bool retA = spectrum.GetMeanFluxInRange( rangeA, Fa );
        Float64 Fb = NAN;
        bool retB = spectrum.GetMeanFluxInRange( rangeB, Fb );
        SContinuumIndex sci;
        sci.Break = NAN;
        sci.Color = NAN;

        if(Fb != 0.0 && retA && retB)
        {
            sci.Color = -2.5*log10(Fa/Fb);
        }

        indexesList.push_back(sci);
    }

    return indexesList;
}
