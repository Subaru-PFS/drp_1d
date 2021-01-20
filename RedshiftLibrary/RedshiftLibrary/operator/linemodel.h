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

namespace NSEpic
{

/**
 * \ingroup Redshift
 */
class COperatorLineModel
{

public:

    COperatorLineModel();
    virtual ~COperatorLineModel();

    std::shared_ptr<COperatorResult> Compute(const CSpectrum& spectrum,
                                             const CTemplateCatalog &tplCatalog,
                                             const TStringList &tplCategoryList,
                                             const std::string opt_calibrationPath,
                                             const CRayCatalog& restraycatalog,
                                             const std::string &opt_lineTypeFilter,
                                             const std::string &opt_lineForceFilter,
                                             const TFloat64Range& lambdaRange,
                                             const TFloat64List& redshifts ,
                                             const Int32 opt_extremacount,
                                             const std::string &opt_fittingmethod,
                                             const std::string &opt_continuumcomponent,
                                             const std::string& opt_lineWidthType,
                                             const Float64 opt_resolution,
                                             const Float64 opt_velocityEmission,
                                             const Float64 opt_velocityAbsorption,
                                             const std::string &opt_continuumreest="no",
                                             const std::string &opt_rules="all",
                                             const std::string &opt_velocityFitting="no",
                                             const Float64 &opt_twosteplargegridstep=0.001,
                                             const std::string &opt_twosteplargegridsampling="log",
                                             const std::string &opt_rigidity="rules",
                                             const string &opt_tplratioCatRelPath="",
                                             const string &opt_offsetCatRelPath="",
                                             const Float64 &opt_emvelocityfitmin=20.,
                                             const Float64 &opt_emvelocityfitmax=500.,
                                             const Float64 &opt_emvelocityfitstep=20.,
                                             const Float64 &opt_absvelocityfitmin=150.,
                                             const Float64 &opt_absvelocityfitmax=500.,
                                             const Float64 &opt_absvelocityfitstep=20.);

    Int32 Init( const CSpectrum& spectrum, 
                const TFloat64List& redshifts, 
                const std::string &opt_continuumcomponent,
                const Float64 nsigmasupport, 
                const Float64 halfwdwsize, 
                const Float64 radius);

    std::shared_ptr<COperatorResult> getResult();
    std::shared_ptr<CLineModelExtremaResult> GetFirstpassExtremaResult() const;

    void PrecomputeContinuumFit(const CSpectrum &spectrum,
                                const CTemplateCatalog &tplCatalog,
                                const TStringList &tplCategoryList,
                                const std::string opt_calibrationPath,
                                const TFloat64Range &lambdaRange,
                                const TFloat64List& redshifts,
                                bool ignoreLinesSupport=false,
                                Int32 candidateIdx = -1);

    Int32 ComputeFirstPass(const CSpectrum& spectrum,
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
                           const Float64 &opt_twosteplargegridstep=0.001,
                           const string &opt_twosteplargegridsampling="log",
                           const std::string &opt_rigidity="rules",
                           const string &opt_tplratioCatRelPath="",
                           const string &opt_offsetCatRelPath="");

    Int32 ComputeCandidates(const Int32 opt_extremacount,
                            const Int32 opt_sign,
                            const std::vector<Float64> floatValues,
                            const Float64 meritCut);
    Int32 Combine_firstpass_candidates(std::shared_ptr<CLineModelExtremaResult> firstpass_results_b);


    Int32 ComputeSecondPass(const CSpectrum& spectrum,
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

    Int32 RecomputeAroundCandidates(TPointList input_extremumList,
                                    const TFloat64Range &lambdaRange,
                                    const std::string &opt_continuumreest,
                                    const Int32 tplfit_option,
                                    const bool overrideRecomputeOnlyOnTheCandidate=false);

    std::shared_ptr<COperatorResult> computeWithUltimPass(const CSpectrum& spectrum,
                                                          const CTemplateCatalog& tplCatalog,
                                                          const TStringList& tplCategoryList,
                                                          const std::string opt_calibrationPath,
                                                          const CRayCatalog& restraycatalog,
                                                          const std::string& opt_lineTypeFilter,
                                                          const std::string& opt_lineForceFilter,
                                                          const TFloat64Range& lambdaRange,
                                                          const TFloat64List& redshifts,
                                                          const Int32 opt_extremacount,
                                                          const std::string& opt_fittingmethod,
                                                          const std::string& opt_continuumcomponent,
                                                          const std::string& opt_lineWidthType,
                                                          const Float64 opt_resolution,
                                                          const Float64 opt_velocityEmission,
                                                          const Float64 opt_velocityAbsorption,
                                                          const std::string& opt_continuumreest,
                                                          const std::string& opt_rules,
                                                          const std::string& opt_velocityFitting,
                                                          const Float64 &opt_twosteplargegridstep,
                                                          const string &opt_twosteplargegridsampling,
                                                          const std::string& opt_rigidity,
                                                          const string &opt_tplratioCatRelPath,
                                                          const string &opt_offsetCatRelPath,
                                                          const Float64 &opt_emvelocityfitmin,
                                                          const Float64 &opt_emvelocityfitmax,
                                                          const Float64 &opt_emvelocityfitstep,
                                                          const Float64 &opt_absvelocityfitmin,
                                                          const Float64 &opt_absvelocityfitmax,
                                                          const Float64 &opt_absvelocityfitstep);

    Int32 SaveResults(const CSpectrum& spectrum,
                      const TFloat64Range& lambdaRange,
                      const std::string &opt_continuumreest="no");

    void InitTplratioPriors();

    void storeGlobalModelResults( COperatorResultStore &resultStore );

    std::shared_ptr<CModelSpectrumResult> GetModelSpectrumResult(Int32 idx);
    std::shared_ptr<CSpectraFluxResult> GetModelSpectrumContinuumResult(Int32 idx);


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

    std::string m_opt_tplfit_method = "templatefittinglog"; //only for firstpass
    std::string m_opt_tplfit_method_secondpass = "templatefittinglog"; 
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
    TPointList m_firstpass_extremumList;
    std::vector<Int32> m_secondpass_indiceSortedCandidatesList;
    CLineModelExtremaResult m_firstpass_extremaResult;
    CLineModelExtremaResult m_secondpass_parameters_extremaResult;

private:

    std::shared_ptr<CLineModelResult> m_result;
    std::shared_ptr<CLineModelElementList> m_model;
    TFloat64List m_sortedRedshifts;
    std::shared_ptr<CTemplateCatalog> m_orthoTplCatalog;
    Int32 m_enableFastFitLargeGrid = 0;
    Int32 m_estimateLeastSquareFast = 0;
    Float64 m_extremaCount;
    Float64 m_Zlinemeasref;

    void ComputeArea1(CLineModelResult& results);
    void ComputeArea2(CLineModelResult& results);
    Float64 FitBayesWidth( CSpectrumSpectralAxis& spectralAxis, CSpectrumFluxAxis& fluxAxis, Float64 z, Int32 start, Int32 end);

    Bool AllAmplitudesAreZero(const TBoolList &amplitudesZero, Int32 nbZ);

    Int32 interpolateLargeGridOnFineGrid(TFloat64List redshiftsLargeGrid, TFloat64List redshiftsFineGrid, TFloat64List meritLargeGrid, TFloat64List &meritFineGrid);
    
    std::shared_ptr<COperatorTemplateFittingBase> templateFittingOperator;
    std::vector<std::shared_ptr<CModelSpectrumResult>  > m_savedModelSpectrumResults;
    std::vector<std::shared_ptr<CModelFittingResult>  > m_savedModelFittingResults;
    std::vector<std::shared_ptr<CModelContinuumFittingResult>  > m_savedModelContinuumFittingResults;
    std::vector<std::shared_ptr<CModelRulesResult>  > m_savedModelRulesResults;
    std::vector<std::shared_ptr<CSpectraFluxResult>  > m_savedModelContinuumSpectrumResults;

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
