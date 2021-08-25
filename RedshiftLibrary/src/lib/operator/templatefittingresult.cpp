#include <RedshiftLibrary/operator/templatefittingresult.h>
#include <RedshiftLibrary/extremum/extremum.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision
#include <boost/algorithm/string/predicate.hpp>

using namespace NSEpic;

void CTemplateFittingResult::Init(UInt32 n , Int32 EbmvListSize, Int32 MeiksinListSize)
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

    std::vector<TFloat64List> _chi2ListList(EbmvListSize, TFloat64List(MeiksinListSize, DBL_MAX));
    std::vector<TFloat64List> _ismListList(EbmvListSize, TFloat64List(MeiksinListSize, NAN));
    std::vector<TInt32List>   _igmListList(EbmvListSize, TInt32List(MeiksinListSize, -1));
        
    ChiSquareIntermediate = std::vector<std::vector<TFloat64List>>(n, _chi2ListList);
    IsmEbmvCoeffIntermediate = std::vector<std::vector<TFloat64List>>(n, _ismListList);
    IgmMeiksinIdxIntermediate = std::vector<std::vector<TInt32List>>(n, _igmListList);
}


