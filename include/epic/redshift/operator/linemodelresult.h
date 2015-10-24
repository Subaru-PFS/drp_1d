#ifndef _REDSHIFT_OPERATOR_LINEMODELRESULT_
#define _REDSHIFT_OPERATOR_LINEMODELRESULT_

#include <epic/redshift/operator/result.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/operator/operator.h>

#include <epic/redshift/ray/catalog.h>

namespace NSEpic
{

class CLineModelResult : public COperatorResult
{

    DEFINE_MANAGED_OBJECT( CLineModelResult )

public:
    struct SLineModelSolution
    {
        std::vector<Float64> ElementId;     //id of the linemodel element it is part of
        std::vector<Float64> Amplitudes;
        std::vector<Float64> Errors;
        std::vector<Float64> Widths;
        std::vector<Bool> OutsideLambdaRange;
        std::vector<TInt32Range> fittingIndexRange;
        Int32 nDDL;
    };

    CLineModelResult();
    virtual ~CLineModelResult();

    Void Save( const COperatorResultStore& store, std::ostream& stream ) const;
    Void SaveLine( const COperatorResultStore& store, std::ostream& stream ) const;
    Void Load( std::istream& stream );

    Int32 GetNLinesOverCutThreshold(Int32 extremaIdx, Float64 cutThres);
    Float64 GetExtremaMerit(Int32 extremaIdx);

    TFloat64List            Redshifts;  // z axis
    TFloat64List            ChiSquare;  // chi2

    TFloat64List            Extrema;    // z extramas
    TFloat64List            LogArea;    // log area for each extrema
    TFloat64List            LogAreaCorrectedExtrema;    //corrected z for each extrema
    TFloat64List            SigmaZ; //sigmaz for each extrema
    std::vector<SLineModelSolution> LineModelSolutions; //linemodel for each extrema
    TFloat64List            bic;    // bayesian information criterion for each extrema

    COperator::TStatusList  Status;
    CRayCatalog::TRayVector restRayList;


};


}

#endif
