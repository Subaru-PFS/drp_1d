#ifndef _REDSHIFT_OPERATOR_LINEMODELRESULT_
#define _REDSHIFT_OPERATOR_LINEMODELRESULT_

#include <epic/redshift/processflow/result.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/operator/operator.h>

#include <epic/redshift/ray/catalog.h>
#include <epic/redshift/continuum/indexes.h>

namespace NSEpic
{

class CLineModelResult : public COperatorResult
{

public:

    struct SLineModelSolution
    {
        std::vector<Float64> ElementId;     //id of the linemodel element it is part of
        std::vector<Float64> Amplitudes;
        std::vector<CRay> Rays;
        std::vector<Float64> Errors;    //noise sigma
        std::vector<Float64> FittingError;    //ModelLeastSquare error under each line

        std::vector<Float64> Widths;
        std::vector<Bool> OutsideLambdaRange;
        std::vector<TInt32Range> fittingIndexRange;

        Float64 LyaWidthCoeff;
        Float64 LyaAlpha;
        Float64 LyaDelta;

        Int32 nDDL;
    };

    CLineModelResult();
    virtual ~CLineModelResult();

    Void ResizeExtremaResults(Int32 size);

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    Void Load( std::istream& stream );

    Int32 GetNLinesOverCutThreshold(Int32 extremaIdx, Float64 snrThres, Float64 fitThres) const;
    Float64 GetExtremaMerit(Int32 extremaIdx) const;

    TFloat64List            Redshifts;  // z axis
    TFloat64List            ChiSquare;  // chi2

    TFloat64List            Extrema;    // z extrema
    TFloat64List            ExtremaExtendedRedshifts;    // z range around extrema
    TFloat64List            Posterior;    // z extrema
    TFloat64List            StrongELSNR;
    TFloat64List            LogArea;    // log area for each extrema
    TFloat64List            LogAreaCorrectedExtrema;    //corrected z for each extrema
    TFloat64List            SigmaZ; //sigmaz for each extrema
    std::vector<SLineModelSolution> LineModelSolutions; //linemodel for each extrema
    TFloat64List            bic;    // bayesian information criterion for each extrema
    std::vector<CContinuumIndexes::TContinuumIndexList> ContinuumIndexes; //continuum indexes for each extrema
    std::vector<CMask>      OutsideLinesMask; //Mask with 0 under the lines and 1 anywhere else

    COperator::TStatusList  Status;
    CRayCatalog::TRayVector restRayList;
    Int32 nSpcSamples;
    Int32 dTransposeDNocontinuum;
    Int32 dTransposeD;


};


}

#endif
