#ifndef _REDSHIFT_OPERATOR_LINEMODELEXTREMARESULT_
#define _REDSHIFT_OPERATOR_LINEMODELEXTREMARESULT_

#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/operator/operator.h>

#include <RedshiftLibrary/ray/catalog.h>
#include <RedshiftLibrary/continuum/indexes.h>
#include <RedshiftLibrary/linemodel/linemodelsolution.h>

namespace NSEpic
{

class CLineModelExtremaResult : public COperatorResult
{

public:

    CLineModelExtremaResult();
    ~CLineModelExtremaResult();

    void Resize(Int32 size);

    void Save( const CDataStore& store, std::ostream& stream) const;
    void SaveJSON( const CDataStore& store, std::ostream& stream) const;

    void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    inline Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const
    {
        return 1;
    }
    //void Load( std::istream& stream );
    bool RemoveSecondPassCandidatebyIdx(Int32 idx);
    Int32 SetWOrderforFP_basedonSortedIDs(TInt32List& Rank_PDF, std::vector<std::string> ids) const;

    //Extrema results
    TFloat64List            ExtremaPDF;    // Ranks extrema
    std::vector<std::string> ExtremaIDs;    // Ranks extrema 

    TFloat64List            Extrema;    // z extrema
    TFloat64List            ExtremaMerit;    // extrema merit
    TFloat64List            ExtremaMeritContinuum; //extrema merit for continuum
    TFloat64List            DeltaZ;    // extrema delta z
    TFloat64List            mTransposeM;    // extrema model norm
    TFloat64List            CorrScaleMarg;    // extrema scale marg. correction
    std::vector<Int32>      NDof;   //non zero elements in the lambdarange
    TFloat64List            ExtremaLastPass; //z extrema with more precision
    TFloat64List            lmfitPass;// z found with lmfit
    TFloat64List            snrHa;
    TFloat64List            lfHa;
    TFloat64List            snrOII;
    TFloat64List            lfOII;

    //Deprecated?
    std::vector<TFloat64List> ExtremaExtendedRedshifts;    // z range around extrema
    TFloat64List            Posterior;    // z extrema
    TFloat64List            LogArea;    // log area for each extrema
    TFloat64List            LogAreaCorrectedExtrema;    //corrected z for each extrema
    TFloat64List            SigmaZ; //sigmaz for each extrema

    //
    TFloat64List            StrongELSNR;
    std::vector<std::vector<std::string>>            StrongELSNRAboveCut;
    TFloat64List            bic;    // bayesian information criterion for each extrema
    std::vector<CContinuumIndexes::TContinuumIndexList> ContinuumIndexes; //continuum indexes for each extrema
    std::vector<CMask>      OutsideLinesMask;   //Mask with 0 under the lines and 1 anywhere else
    TFloat64List            OutsideLinesSTDFlux;    //STD measured on the spectrum continuum substracted outside lines
    TFloat64List            OutsideLinesSTDError;   //STD measured on the error spectrum outside lines

    //line width
    TFloat64List      Elv;   //emission line width
    TFloat64List      Alv;   //absorption line width
    std::vector<TFloat64List>      GroupsELv;   //per fitting group line width , EL
    std::vector<TFloat64List>      GroupsALv;   //per fitting group line width , AL


    //template continuum
    std::vector<std::string>      FittedTplName;    //Name of the best template fitted for continuum
    TFloat64List            FittedTplAmplitude;     //Amplitude for the best template fitted for continuum
    TFloat64List            FittedTplAmplitudeError;     //Amplitude error for the best template fitted for continuum
    TFloat64List            FittedTplMerit;     //Chisquare for the best template fitted for continuum
    TFloat64List            FittedTplDustCoeff;     //Calzetti dustcoeff for the best template fitted for continuum
    std::vector<Int32>      FittedTplMeiksinIdx;    //Meiksin igm index for the best template fitted for continuum
    TFloat64List      FittedTplRedshift;    //Redshift for the best template fitted for continuum
    TFloat64List      FittedTplDtm;    //DTM for the best template fitted for continuum
    TFloat64List      FittedTplMtm;    //MTM for the best template fitted for continuum
    TFloat64List      FittedTplLogPrior;    //log prior for the best template fitted for continuum
    std::vector<TFloat64List>      FittedTplpCoeffs;    //poly coeffs for the best template fitted for continuum

    //template ratio
    std::vector<std::string>      FittedTplshapeName;   //Name of the best template fitted for tplcorr/tplshape
    TFloat64List      FittedTplshapeAmplitude;   //amp of the best template fitted for tplcorr/tplshape
    TFloat64List      FittedTplshapeDtm;   //dtm of the best template fitted for tplcorr/tplshape
    TFloat64List      FittedTplshapeMtm;   //mtm of the best template fitted for tplcorr/tplshape
    TFloat64List      FittedTplshapeIsmCoeff;   //IsmCoeff/EBMV of the best template fitted for tplcorr/tplshape

};


}

#endif
