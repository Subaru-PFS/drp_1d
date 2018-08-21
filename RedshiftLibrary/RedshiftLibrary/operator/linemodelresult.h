#ifndef _REDSHIFT_OPERATOR_LINEMODELRESULT_
#define _REDSHIFT_OPERATOR_LINEMODELRESULT_

#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/operator/operator.h>

#include <RedshiftLibrary/ray/catalog.h>
#include <RedshiftLibrary/continuum/indexes.h>
#include <RedshiftLibrary/linemodel/linemodelextremaresult.h>
#include <RedshiftLibrary/linemodel/linemodelsolution.h>

namespace NSEpic
{

class CLineModelResult : public COperatorResult
{
public:

    CLineModelResult();
    virtual ~CLineModelResult();

    Int32 Init(std::vector<Float64> redshifts, CRayCatalog::TRayVector restRays, Int32 nTplshapes, std::vector<Float64> tplshapesPriors);

    void Save( const CDataStore& store, std::ostream& stream ) const;
    void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    inline Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const
    {
        return 1;
    }

    void Load( std::istream& stream );

    Int32 GetNLinesOverCutThreshold(Int32 extremaIdx, Float64 snrThres, Float64 fitThres) const;
    std::vector<bool> GetStrongLinesPresence( UInt32 filterType, std::vector<CLineModelSolution> linemodelsols ) const;
    Float64 GetExtremaMerit(Int32 extremaIdx) const;
    UInt32 GetExtremaIndex(UInt32 extremaIdx) const;

    Float64 GetMinChiSquare() const;
    Float64 GetMaxChiSquare() const;

    Int32 ResizeChisquareTplShapes( Int32 nTplshapes, Int32 nRedshifts );
    Int32 SetChisquareTplshapeResult(Int32 index, TFloat64List chisquareTplshape, TFloat64List scaleMargCorrTplshape, std::vector<bool> strongEmissionLinePresentTplshape);
    TFloat64List GetChisquareTplshapeResult( Int32 index );
    TFloat64List GetScaleMargCorrTplshapeResult( Int32 index );
    std::vector<bool> GetStrongELPresentTplshapeResult( Int32 index );

    //Merit results
    TFloat64List            Redshifts;  // z axis
    TFloat64List            ChiSquare;  // min chi2
    TFloat64List            ScaleMargCorrection;  // margCorrection for min chi2

    std::vector<TFloat64List> ChiSquareTplshapes; // full chi2 results (for each tplshape)
    std::vector<Float64> PriorTplshapes; // model prior (for each tplshape)
    std::vector<TFloat64List> ScaleMargCorrectionTplshapes; // full scale marginalization correction results (for each tplshape)
    std::vector<std::vector<bool>> StrongELPresentTplshapes; // full strongELPresent results (for each tplshape)
    TFloat64List ChiSquareContinuum; // chi2 result for the continuum
    TFloat64List ScaleMargCorrectionContinuum; //  scale marginalization correction result for the continuum

    std::vector<CLineModelSolution> LineModelSolutions;

    //Extrema results
    CLineModelExtremaResult ExtremaResult;

    //
    COperator::TStatusList  Status;
    CRayCatalog::TRayVector restRayList;
    Int32 nSpcSamples;
    Float64 dTransposeDNocontinuum;
    Float64 dTransposeD;
    Float64 cstLog;
};


}

#endif
