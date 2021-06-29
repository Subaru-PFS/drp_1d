#include "RedshiftLibrary/operator/tplcombinationresult.h"
#include <cfloat>
using namespace NSEpic;

void CTplCombinationResult::Init(UInt32 n , Int32 nISM, Int32 nIGM, Int32 componentSize)
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
    //Queue
    TFloat64List _ampList(componentSize, NAN);
    TFloat64List _chi2IGMList(nIGM, DBL_MAX);
    TFloat64List _EbmvCoeffIGMList(nIGM, NAN);
    TInt32List   _meiksinIdxIGMList(nIGM, -1);
    
    //fill for each ism
    std::vector<TFloat64List> _ChiSquareISMList;
    std::vector<TFloat64List> _IsmEbmvCoeffISMList;
    std::vector<TInt32List>   _IgmMeiksinIdxISMList;
    for(Int32 kism=0; kism<nISM; kism++)
    {
        _ChiSquareISMList.push_back(_chi2IGMList);
        _IsmEbmvCoeffISMList.push_back(_EbmvCoeffIGMList);
        _IgmMeiksinIdxISMList.push_back(_meiksinIdxIGMList);

    }

    //fill for each z
    ChiSquareIntermediate.clear();
    IsmEbmvCoeffIntermediate.clear();
    IgmMeiksinIdxIntermediate.clear();
    for(Int32 k=0; k<n; k++)
    {
        ChiSquareIntermediate.push_back(_ChiSquareISMList);
        IsmEbmvCoeffIntermediate.push_back(_IsmEbmvCoeffISMList);
        IgmMeiksinIdxIntermediate.push_back(_IgmMeiksinIdxISMList);

        //fit amplitude related data with NAN vetors
        FitAmplitude.push_back(_ampList);
        FitAmplitudeError.push_back(_ampList);
        FitAmplitudeSigma.push_back(_ampList);
    }
}


