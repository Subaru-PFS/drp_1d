#ifndef _REDSHIFT_OPERATOR_LINEMODELRESULT_
#define _REDSHIFT_OPERATOR_LINEMODELRESULT_

#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/operator/operator.h>

#include <RedshiftLibrary/ray/catalog.h>
#include <RedshiftLibrary/continuum/indexes.h>

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
        std::vector<Float64> CenterContinuumFlux;    //Continuum flux value at the center of each line
        std::vector<Float64> Sigmas;    //width for each line
        std::vector<Float64> Fluxs;    //Flux for each line
        std::vector<Float64> FluxErrors;    //Flux error for each line

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
    std::vector<bool> GetStrongLinesPresence( UInt32 filterType=0 ) const;
    Float64 GetExtremaMerit(Int32 extremaIdx) const;
    UInt32 GetExtremaIndex(UInt32 extremaIdx) const;

    Float64 GetMinChiSquare() const;
    Float64 GetMaxChiSquare() const;


    //Full Merit curve
    TFloat64List            Redshifts;  // z axis
    TFloat64List            ChiSquare;  // chi2

    //Extrema results
    TFloat64List            Extrema;    // z extrema
    TFloat64List            ExtremaMerit;    // extrema merit
    TFloat64List            DeltaZ;    // extrema delta z
    TFloat64List            mTransposeM;    // extrema model norm
    TFloat64List            ExtremaLastPass; //z extrema with more precision

    //Deprecated?
    TFloat64List            ExtremaExtendedRedshifts;    // z range around extrema
    TFloat64List            Posterior;    // z extrema
    TFloat64List            LogArea;    // log area for each extrema
    TFloat64List            LogAreaCorrectedExtrema;    //corrected z for each extrema
    TFloat64List            SigmaZ; //sigmaz for each extrema

    //
    TFloat64List            StrongELSNR;
    std::vector<SLineModelSolution> LineModelSolutions; //linemodel for each extrema
    TFloat64List            bic;    // bayesian information criterion for each extrema
    std::vector<CContinuumIndexes::TContinuumIndexList> ContinuumIndexes; //continuum indexes for each extrema
    std::vector<CMask>      OutsideLinesMask; //Mask with 0 under the lines and 1 anywhere else
    std::vector<std::string>      FittedTplName; //Name of the best template fitted for continuum
    TFloat64List      FittedTplAmplitude; //Amplitude for the best template fitted for continuum
    TFloat64List      FittedTplDustCoeff; //Calzetti dustcoeff for the best template fitted for continuum
    std::vector<Int32>      FittedTplMeiksinIdx; //Meiksin igm index for the best template fitted for continuum
    std::vector<std::string>      FittedTplcorrTplName; //Name of the best template fitted for tplcorr/tplshape

    //
    COperator::TStatusList  Status;
    CRayCatalog::TRayVector restRayList;
    Int32 nSpcSamples;
    Float64 dTransposeDNocontinuum;
    Float64 dTransposeD;


};


}

#endif
