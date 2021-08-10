#include <RedshiftLibrary/operator/templatefittingresult.h>
#include <RedshiftLibrary/extremum/extremum.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision
#include <boost/algorithm/string/predicate.hpp>

using namespace NSEpic;

void CTemplateFittingResult::Init(UInt32 n , Int32 nISM, Int32 nIGM)
{
    ChiSquare.resize( n );
    FitAmplitude.resize( n );
    FitAmplitudeError.resize( n );
    FitAmplitudeSigma.resize( n);
    FitEbmvCoeff.resize( n );
    FitMeiksinIdx.resize( n );
    FitDtM.resize( n );
    FitMtM.resize( n );
    LogPrior.resize( n );
    Redshifts.resize( n );
    Overlap.resize( n );
    Status.resize( n );
    SNR.resize( n );
    ChiSquareIntermediate.clear();
    IsmEbmvCoeffIntermediate.clear();
    IgmMeiksinIdxIntermediate.clear();
    
    std::vector<TFloat64List> _ChiSquareISMList;
    std::vector<TFloat64List> _IsmEbmvCoeffISMList;
    std::vector<TInt32List> _IgmMeiksinIdxISMList;

    TFloat64List _chi2IGMList(nIGM, DBL_MAX);
    TFloat64List _EbmvCoeffIGMList(nIGM, NAN);
    TInt32List _meiksinIdxIGMList(nIGM, -1);
    for(Int32 kism=0; kism<nISM; kism++)
    {
        _ChiSquareISMList.push_back(_chi2IGMList);
        _IsmEbmvCoeffISMList.push_back(_EbmvCoeffIGMList);
        _IgmMeiksinIdxISMList.push_back(_meiksinIdxIGMList);

    }
    for(Int32 k=0; k<n; k++)
    {
        ChiSquareIntermediate.push_back(_ChiSquareISMList);
        IsmEbmvCoeffIntermediate.push_back(_IsmEbmvCoeffISMList);
        IgmMeiksinIdxIntermediate.push_back(_IgmMeiksinIdxISMList);
    }
}


