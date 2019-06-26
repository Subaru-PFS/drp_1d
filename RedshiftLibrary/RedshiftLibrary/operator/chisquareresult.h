#ifndef _REDSHIFT_OPERATOR_CHISQUARERESULT_
#define _REDSHIFT_OPERATOR_CHISQUARERESULT_

#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/operator/operator.h>

namespace NSEpic
{

class CChisquareResult : public COperatorResult
{

public:

    CChisquareResult();
    virtual ~CChisquareResult();

    void Init( UInt32 n, Int32 nISM, Int32 nIGM);

    void Save( const CDataStore& store, std::ostream& stream ) const;
    void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    inline Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const
    {
        return 1;
    }

    void Load( std::istream& stream );

    TFloat64List            Redshifts;

    //best fit results
    TFloat64List            ChiSquare;
    TFloat64List            FitAmplitude;
    TFloat64List            FitDustCoeff;
    TFloat64List            FitMeiksinIdx;
    TFloat64List            FitDtM;
    TFloat64List            FitMtM;
    TFloat64List            LogPrior;

    //intermediate chisquare results
    std::vector<std::vector<TFloat64List>> ChiSquareIntermediate; // chi2 for each intermediate results (for each config [z][Calzetti][Meiksin])
    std::vector<std::vector<TFloat64List>> IsmDustCoeffIntermediate; // calzetti dust coeff for each intermediate result (for each config [z][Calzetti][Meiksin])
    std::vector<std::vector<TInt32List>> IgmMeiksinIdxIntermediate; // meiksin idx for each intermediate result (for each config [z][Calzetti][Meiksin])
    //TODO: std::vector<std::vector<TFloat64List>> LogPriorIntermediate
    //TODO: std::vector<std::vector<TFloat64List>> AmpIntermediate //is needed for correct prior use in marg. mode tplmodel method

    Float64                 CstLog;
    TFloat64List            Overlap;
    TFloat64List            Extrema;
    COperator::TStatusList  Status;

};


}

#endif
