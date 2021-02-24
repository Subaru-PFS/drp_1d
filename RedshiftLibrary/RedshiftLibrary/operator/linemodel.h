#ifndef _REDSHIFT_OPERATOR_LINEMODEL_
#define _REDSHIFT_OPERATOR_LINEMODEL_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/operator/operator.h>
#include <RedshiftLibrary/operator/linemodelresult.h>
#include <RedshiftLibrary/linemodel/elementlist.h>
#include <RedshiftLibrary/linemodel/multirollmodel.h>
#include <RedshiftLibrary/operator/modelspectrumresult.h>
#include <RedshiftLibrary/linemodel/modelfittingresult.h>
#include <RedshiftLibrary/operator/modelcontinuumfittingresult.h>
#include <RedshiftLibrary/linemodel/modelrulesresult.h>
#include <RedshiftLibrary/operator/spectraFluxResult.h>
#include <RedshiftLibrary/common/mask.h>
#include <RedshiftLibrary/processflow/resultstore.h>

#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/ray/catalog.h>

#include <RedshiftLibrary/linemodel/linemodelextremaresult.h>

namespace NSEpic
{

class COperatorLineModelExtremaResult 
{

public:

  COperatorLineModelExtremaResult() = default;
  COperatorLineModelExtremaResult(Int32 n)
  {
    Resize(n);
  }
  ~COperatorLineModelExtremaResult() = default;

  void Resize(Int32 size)
  {
    m_ranked_candidates.resize(size);
    FittedTplName.resize(size);
    FittedTplAmplitude.resize(size);
    FittedTplAmplitudeError.resize(size);
    FittedTplMerit.resize(size);
    FittedTplEbmvCoeff.resize(size);
    FittedTplMeiksinIdx.resize(size);
    FittedTplDtm.resize(size);
    FittedTplMtm.resize(size);
    FittedTplLogPrior.resize(size);
    FittedTplSNR.resize(size);
    
    m_savedModelSpectrumResults.resize(size);
    m_savedModelContinuumFittingResults.resize(size);
        MeritContinuum.resize(size);

    mTransposeM.resize(size);
    CorrScaleMarg.resize(size);
    NDof.resize(size);
    Redshift_lmfit.resize(size);
    snrHa.resize(size);
    lfHa.resize(size);
    snrOII.resize(size);
    lfOII.resize(size);

    ExtendedRedshifts.resize(size);
    NLinesOverThreshold.resize(size);
    LogArea.resize(size);
    LogAreaCorrectedExtrema.resize(size);
    SigmaZ.resize(size);

    StrongELSNR.resize(size);
    StrongELSNRAboveCut.resize(size);
    bic.resize(size);
    ContinuumIndexes.resize(size);
    OutsideLinesMask.resize(size);
    OutsideLinesSTDFlux.resize(size);
    OutsideLinesSTDError.resize(size);

    Elv.resize(size);
    Alv.resize(size);
    GroupsELv.resize(size);
    GroupsALv.resize(size);
    for(Int32 ke=0; ke<size; ke++)
    {
        GroupsELv[ke] = std::vector<Float64>(250, -1);   //WARNING: hardcoded ngroups max
        GroupsALv[ke] = std::vector<Float64>(250, -1);   //WARNING: hardcoded ngroups max
    }

    FittedTplRedshift.resize(size);
    FittedTplpCoeffs.resize(size);

    FittedTplratioName.resize(size);
    FittedTplratioAmplitude.resize(size);
    FittedTplratioDtm.resize(size);
    FittedTplratioMtm.resize(size);
    FittedTplratioIsmCoeff.resize(size);
    
    m_savedModelFittingResults.resize(size);
    m_savedModelRulesResults.resize(size);
    m_savedModelContinuumSpectrumResults.resize(size);
  }
  TFloat64List GetRedshifts() const
  {
    TFloat64List redshifts;
    redshifts.reserve(m_ranked_candidates.size());
    for (auto c: m_ranked_candidates) redshifts.push_back(c.second.Redshift);
    return redshifts;
  }

  TStringList GetIDs() const
  {
    TStringList ids;
    ids.reserve(m_ranked_candidates.size());
    for (auto c: m_ranked_candidates) ids.push_back(c.first);
    return ids;
}

   std::string ID(Int32 i) const {return m_ranked_candidates[i].first;}
 Float64 Redshift(Int32 i) const { return m_ranked_candidates[i].second.Redshift;}
 Float64 ValProba(Int32 i) const { return m_ranked_candidates[i].second.ValProba;}
 Float64 ValSumProba(Int32 i) const { return m_ranked_candidates[i].second.ValSumProba;}
 Float64 DeltaZ(Int32 i) const { return m_ranked_candidates[i].second.Deltaz;}

  Int32 size() const { return m_ranked_candidates.size();}
  
  Int32                       m_optMethod; //0: direct integration, 1:gaussian fit

  TCandidateZbyRank m_ranked_candidates;

  

    //template continuum
  TStringList       FittedTplName;    //Name of the best template fitted for continuum
  TFloat64List      FittedTplAmplitude;     //Amplitude for the best template fitted for continuum
  TFloat64List      FittedTplAmplitudeError;     //Amplitude error for the best template fitted for continuum
  TFloat64List      FittedTplMerit;     //Chisquare for the best template fitted for continuum
  TFloat64List      FittedTplEbmvCoeff;     //Calzetti ebmvcoeff for the best template fitted for continuum
  TInt32List        FittedTplMeiksinIdx;    //Meiksin igm index for the best template fitted for continuum
  TFloat64List      FittedTplDtm;    //DTM for the best template fitted for continuum
  TFloat64List      FittedTplMtm;    //MTM for the best template fitted for continuum
  TFloat64List      FittedTplLogPrior;    //log prior for the best template fitted for continuum
  TFloat64List      FittedTplSNR; 

  std::vector<std::shared_ptr<const CModelSpectrumResult>  > m_savedModelSpectrumResults;
  std::vector<std::shared_ptr<const CModelContinuumFittingResult>  > m_savedModelContinuumFittingResults;


    //Extrema results
    TFloat64List            MeritContinuum; //extrema merit for continuum

    TFloat64List            mTransposeM;    // extrema model norm
    TFloat64List            CorrScaleMarg;    // extrema scale marg. correction
    TInt32List              NDof;   //non zero elements in the lambdarange
    TFloat64List            Redshift_lmfit;// z found with lmfit
    TFloat64List            snrHa;
    TFloat64List            lfHa;
    TFloat64List            snrOII;
    TFloat64List            lfOII;

    std::vector<TFloat64List> ExtendedRedshifts;    // z range around extrema
    TFloat64List            NLinesOverThreshold;  
    TFloat64List            LogArea;   // log area for each extrema
    TFloat64List            LogAreaCorrectedExtrema;   // corrected z for each extrema
    TFloat64List            SigmaZ;    // sigmaz for each extrema

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

    //template continuum (+ base class)
    TFloat64List      FittedTplRedshift;    //Redshift for the best template fitted for continuum
    std::vector<TFloat64List>      FittedTplpCoeffs;    //poly coeffs for the best template fitted for continuum

    //template ratio
    std::vector<std::string>      FittedTplratioName;   //Name of the best template fitted for tplcorr/tplratio
    TFloat64List      FittedTplratioAmplitude;   //amp of the best template fitted for tplcorr/tplratio
    TFloat64List      FittedTplratioDtm;   //dtm of the best template fitted for tplcorr/tplratio
    TFloat64List      FittedTplratioMtm;   //mtm of the best template fitted for tplcorr/tplratio
    TFloat64List      FittedTplratioIsmCoeff;   //IsmCoeff/EBMV of the best template fitted for tplcorr/tplratio

    mutable std::map<int,TFloat64List> continuumIndexesColorCopy;
    mutable std::map<int,TFloat64List> continuumIndexesBreakCopy;
    
    std::vector<std::shared_ptr<const CModelFittingResult>  > m_savedModelFittingResults;
    std::vector<std::shared_ptr<const CModelRulesResult>  > m_savedModelRulesResults;
    std::vector<std::shared_ptr<const CSpectraFluxResult>  > m_savedModelContinuumSpectrumResults;

};


  
/**
 * \ingroup Redshift
 */
class COperatorLineModel
{

public:

    COperatorLineModel();
    virtual ~COperatorLineModel();


    Int32 Init( const CSpectrum& spectrum, 
                const TFloat64List& redshifts, 
                const std::string &opt_continuumcomponent,
                const Float64 nsigmasupport, 
                const Float64 halfwdwsize, 
                const Float64 radius);

    std::shared_ptr<COperatorResult> getResult();

    void PrecomputeContinuumFit(const CSpectrum &spectrum,
                                const CSpectrum &rebinnedSpectrum,
                                const CTemplateCatalog &tplCatalog,
                                const TStringList &tplCategoryList,
                                const std::string opt_calibrationPath,
                                const TFloat64Range &lambdaRange,
                                const TFloat64List& redshifts,
                                bool ignoreLinesSupport=false,
                                Int32 candidateIdx = -1);

    Int32 ComputeFirstPass(const CSpectrum& spectrum,
                           const CSpectrum& rebinnedSpc,
                           const CTemplateCatalog &tplCatalog,
                           const TStringList &tplCategoryList,
                           const std::string opt_calibrationPath,
                           const CRayCatalog& restraycatalog,
                           const std::string &opt_lineTypeFilter,
                           const std::string &opt_lineForceFilter,
                           const TFloat64Range& lambdaRange,
                           const std::string &opt_fittingmethod,
                           const std::string& opt_lineWidthType,
                           const Float64 opt_resolution,
                           const Float64 opt_velocityEmission,
                           const Float64 opt_velocityAbsorption,
                           const std::string &opt_continuumreest="no",
                           const std::string &opt_rules="all",
                           const std::string &opt_velocityFitting="no",
                           const UInt32 &opt_twosteplargegridstep_ratio=10,
                           const string &opt_twosteplargegridsampling="log",
                           const std::string &opt_rigidity="rules",
                           const string &opt_tplratioCatRelPath="",
                           const string &opt_offsetCatRelPath="");
    void CreateRedshiftLargeGrid(Int32 ratio, TFloat64List& largeGridRedshifts);
    Int32 SetFirstPassCandidates(const TCandidateZbyRank & candidatesz);

  Int32 Combine_firstpass_candidates(std::shared_ptr<const COperatorLineModelExtremaResult> firstpass_results_b);


    Int32 ComputeSecondPass(const CSpectrum& spectrum,
                            const CSpectrum &rebinnedSpectrum,
                            const CTemplateCatalog &tplCatalog,
                            const TStringList &tplCategoryList,
                            const std::string opt_calibrationPath,
                            const CRayCatalog& restraycatalog,
                            const std::string &opt_lineTypeFilter,
                            const std::string &opt_lineForceFilter,
                            const TFloat64Range& lambdaRange,
                            const std::string &opt_fittingmethod,
                            const std::string& opt_lineWidthType,
                            const Float64 opt_resolution,
                            const Float64 opt_velocityEmission,
                            const Float64 opt_velocityAbsorption,
                            const std::string &opt_continuumreest="no",
                            const std::string &opt_rules="all",
                            const std::string &opt_velocityFitting="no",
                            const std::string &opt_rigidity="rules",
                            const Float64 &opt_emvelocityfitmin=20.,
                            const Float64 &opt_emvelocityfitmax=500.,
                            const Float64 &opt_emvelocityfitstep=20.,
                            const Float64 &opt_absvelocityfitmin=150.,
                            const Float64 &opt_absvelocityfitmax=500.,
                            const Float64 &opt_absvelocityfitstep=20.,
                            const string &opt_continuumfit_method="fromfirstpass");

    Int32 EstimateSecondPassParameters(const CSpectrum &spectrum,
                                       const TFloat64Range &lambdaRange,
                                       const string &opt_continuumreest,
                                       const string &opt_fittingmethod,
                                       const std::string &opt_rigidity,
                                       const string &opt_velocityFitting,
                                       const Float64 &opt_emvelocityfitmin,
                                       const Float64 &opt_emvelocityfitmax,
                                       const Float64 &opt_emvelocityfitstep,
                                       const Float64 &opt_absvelocityfitmin,
                                       const Float64 &opt_absvelocityfitmax,
                                       const Float64 &opt_absvelocityfitstep);

    Int32 RecomputeAroundCandidates(const TFloat64Range &lambdaRange,
                                    const std::string &opt_continuumreest,
                                    const Int32 tplfit_option,
                                    const bool overrideRecomputeOnlyOnTheCandidate=false);

  std::shared_ptr<LineModelExtremaResult> SaveExtremaResults(const CSpectrum& spectrum,
                                                         const TFloat64Range& lambdaRange,
                                                         const TCandidateZbyRank & zCandidates,
                                                         const std::string &opt_continuumreest="no");

    void InitTplratioPriors();

    bool m_enableWidthFitByGroups = false;

    Float64 m_linesmodel_nsigmasupport;

    Int32 m_maxModelSaveCount;
    Float64 m_secondPass_halfwindowsize; // = 0.005;
    Float64 m_extremaRedshiftSeparation; 


    bool m_enableLoadContTemplate=false;
    Int32 m_iRollContaminated=-1;
    Float64 m_contLambdaOffset=0;
    std::shared_ptr<CTemplate> m_tplContaminant=NULL;
    Int32 initContaminant(std::shared_ptr<CModelSpectrumResult> contModelSpectrum, Int32 iRollContaminated, Float64 contLambdaOffset);
    std::shared_ptr<CModelSpectrumResult> GetContaminantSpectrumResult();
    std::shared_ptr<CModelSpectrumResult> m_savedContaminantSpectrumResult;

    bool m_opt_tplfit_fftprocessing = false; //we cant set it as the default since not taken into account when deiding on rebinning
    bool m_opt_tplfit_fftprocessing_secondpass = false;//true;
    Int32 m_opt_tplfit_dustFit = 1;
    Int32 m_opt_tplfit_extinction = 1;
    Int32 m_opt_fitcontinuum_maxN = 2;
    bool m_opt_tplfit_ignoreLinesSupport=false; //default: false, as ortho templates store makes this un-necessary
    Float64 m_opt_tplfit_continuumprior_betaA=1.0;
    Float64 m_opt_tplfit_continuumprior_betaTE=1.0;
    Float64 m_opt_tplfit_continuumprior_betaZ=1.0;
    std::string m_opt_tplfit_continuumprior_dirpath="";

    Int32 m_opt_tplratio_ismFit = 1;
    Int32 m_opt_firstpass_tplratio_ismFit=0;
    Int32 m_opt_firstpass_multiplecontinuumfit_disable=1;
    std::string m_opt_firstpass_fittingmethod;
    std::string m_opt_secondpasslcfittingmethod="-1";
    Float64 m_opt_tplratio_prior_betaA=1.0;
    Float64 m_opt_tplratio_prior_betaTE=1.0;
    Float64 m_opt_tplratio_prior_betaZ=1.0;
    std::string m_opt_tplratio_prior_dirpath="";
    std::string m_opt_continuumcomponent;
    Float64 m_opt_continuum_neg_amp_threshold = -INFINITY;
    std::string m_opt_enableLSF;

    std::string m_opt_lya_forcefit;
    std::string m_opt_lya_forcedisablefit;
    Float64 m_opt_lya_fit_asym_min;
    Float64 m_opt_lya_fit_asym_max;
    Float64 m_opt_lya_fit_asym_step;
    Float64 m_opt_lya_fit_width_min;
    Float64 m_opt_lya_fit_width_max;
    Float64 m_opt_lya_fit_width_step;
    Float64 m_opt_lya_fit_delta_min;
    Float64 m_opt_lya_fit_delta_max;
    Float64 m_opt_lya_fit_delta_step;

    std::string m_opt_enableImproveBalmerFit;
    Int32 m_continnuum_fit_option = 0;//default to "retryall" templates
    //candidates
  std::shared_ptr<COperatorLineModelExtremaResult> m_firstpass_extremaResult;
    COperatorLineModelExtremaResult m_secondpass_parameters_extremaResult;

private:

    std::shared_ptr<CLineModelResult> m_result;
    std::shared_ptr<CLineModelElementList> m_model;
    TFloat64List m_sortedRedshifts;
    std::shared_ptr<CTemplateCatalog> m_orthoTplCatalog;
    Int32 m_enableFastFitLargeGrid = 0;
    Int32 m_estimateLeastSquareFast = 0;
    Float64 m_extremaCount;
    Float64 m_Zlinemeasref;

    TFloat64List SpanRedshiftWindow(Float64 z) const;

    Float64 FitBayesWidth(const CSpectrumSpectralAxis& spectralAxis,const CSpectrumFluxAxis& fluxAxis, Float64 z, Int32 start, Int32 end);

    Bool AllAmplitudesAreZero(const TBoolList &amplitudesZero, Int32 nbZ);

    Int32 interpolateLargeGridOnFineGrid(TFloat64List redshiftsLargeGrid, TFloat64List redshiftsFineGrid, TFloat64List meritLargeGrid, TFloat64List &meritFineGrid);
    
    std::shared_ptr<COperatorTemplateFittingBase> templateFittingOperator;

    //lmfit

    bool mlmfit_modelInfoSave = false;
    std::vector<std::shared_ptr<CModelSpectrumResult>> mlmfit_savedModelSpectrumResults_lmfit;
    std::vector<std::shared_ptr<CModelFittingResult>> mlmfit_savedModelFittingResults_lmfit;
    std::vector<std::shared_ptr<CModelRulesResult>> mlmfit_savedModelRulesResults_lmfit;
    std::vector<std::shared_ptr<CSpectraFluxResult>> mlmfit_savedBaselineResult_lmfit;

    std::shared_ptr<CPriorHelper> m_phelperContinuum;
    std::shared_ptr<CTemplatesFitStore> m_tplfitStore_firstpass; 
    std::vector<std::shared_ptr<CTemplatesFitStore>> m_tplfitStore_secondpass;

};


}

#endif
