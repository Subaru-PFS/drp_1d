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

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    Void Load( std::istream& stream );

    TFloat64List            Redshifts;
    TFloat64List            ChiSquare;
    TFloat64List            FitAmplitude;
    TFloat64List            FitDustCoeff;
    TFloat64List            FitMeiksinIdx;
    TFloat64List            FitDtM;
    TFloat64List            FitMtM;
    Float64                 CstLog;

    TFloat64List            Overlap;
    TFloat64List            Extrema;
    COperator::TStatusList  Status;

};


}

#endif
