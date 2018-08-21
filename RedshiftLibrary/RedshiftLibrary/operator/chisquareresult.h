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
    TFloat64List            ChiSquare;
    TFloat64List            FitAmplitude;
    TFloat64List            FitDustCoeff;
    TFloat64List            FitMeiksinIdx;
    TFloat64List            FitDtM;
    TFloat64List            FitMtM;
    Float64                 CstLog;

    //intermediate chisquare results
    std::vector<std::vector<TFloat64List>> ChiSquareIntermediate; // full chi2 results (for each config [Calzetti, Meiksin])

    TFloat64List            Overlap;
    TFloat64List            Extrema;
    COperator::TStatusList  Status;

};


}

#endif
