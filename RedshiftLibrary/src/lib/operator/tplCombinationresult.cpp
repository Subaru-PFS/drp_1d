#include "RedshiftLibrary/operator/tplcombinationresult.h"
#include <cfloat>
using namespace NSEpic;

void CTplCombinationResult::Init(UInt32 n , Int32 EbmvListSize, Int32 MeiksinListSize, Int32 componentSize)
{
    ChiSquare.resize(n);
    FitEbmvCoeff.resize(n);
    FitMeiksinIdx.resize(n);
    FitCOV.resize(n);//covariance
    //LogPrior.resize(n);
    Redshifts.resize(n);
    Overlap.resize(n);
    Status.resize(n);
    SNR.resize(n);

    ChiSquareIntermediate.clear();
    IsmEbmvCoeffIntermediate.clear();
    IgmMeiksinIdxIntermediate.clear();
    
    std::vector<TFloat64List> _chi2ListList(EbmvListSize, TFloat64List(MeiksinListSize, DBL_MAX));
    std::vector<TFloat64List> _ismListList(EbmvListSize, TFloat64List(MeiksinListSize, NAN));
    std::vector<TInt32List>   _igmListList(EbmvListSize, TInt32List(MeiksinListSize, -1));

    ChiSquareIntermediate = std::vector<std::vector<TFloat64List>>(n, _chi2ListList);
    IsmEbmvCoeffIntermediate = std::vector<std::vector<TFloat64List>>(n, _ismListList);
    IgmMeiksinIdxIntermediate = std::vector<std::vector<TInt32List>>(n, _igmListList);

    FitAmplitude = std::vector<TFloat64List>(n,TFloat64List(componentSize, NAN));
    FitAmplitudeError = std::vector<TFloat64List>(n,TFloat64List(componentSize, NAN));
    FitAmplitudeSigma = std::vector<TFloat64List>(n,TFloat64List(componentSize, NAN));
}


