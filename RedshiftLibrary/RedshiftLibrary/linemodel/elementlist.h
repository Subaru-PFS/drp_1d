#ifndef _REDSHIFT_LINEMODEL_ELEMENTLIST_
#define _REDSHIFT_LINEMODEL_ELEMENTLIST_

#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/common/datatypes.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <RedshiftLibrary/ray/catalog.h>
#include <RedshiftLibrary/ray/catalogsTplShape.h>
#include <RedshiftLibrary/ray/regulament.h>
#include <RedshiftLibrary/spectrum/spectrum.h>

#include <RedshiftLibrary/operator/chisquare2.h>

#include <RedshiftLibrary/operator/linemodelresult.h>
#include <RedshiftLibrary/operator/modelspectrumresult.h>
#include <RedshiftLibrary/linemodel/element.h>

#include <RedshiftLibrary/statistics/pdfz.h>

#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/linemodel/templatesfitstore.h>

#include <boost/shared_ptr.hpp>

#include <RedshiftLibrary/linemodel/lmfitcontroller.h>

#include <memory>
namespace NSEpic
{
  static Int32 defaultIdx = -1;

class CLineModelElementList
{

public:

    CLineModelElementList(const CSpectrum& spectrum,
                          const CTemplateCatalog& tplCatalog,
                          const CTemplateCatalog& orthoTplCatalog,
                          const TStringList& tplCategoryList,
                          const std::string calibrationPath,
                          const CRayCatalog::TRayVector& restRayList,
                          const std::string& opt_fittingmethod,
                          const std::string &opt_continuumcomponent,
                          const std::string& lineWidthType,
                          const Float64 nsigmasupport,
                          const Float64 resolution,
                          const Float64 velocityEmission,
                          const Float64 velocityAbsorption,
                          const std::string &opt_rules,
                          const std::string &opt_rigidity);

    ~CLineModelElementList();

    void LoadCatalog(const CRayCatalog::TRayVector& restRayList);
    void LoadCatalogOneMultiline(const CRayCatalog::TRayVector& restRayList);
    void LoadCatalogTwoMultilinesAE(const CRayCatalog::TRayVector& restRayList);

    void LogCatalogInfos();

    void PrepareContinuum();
    void EstimateSpectrumContinuum(Float64 opt_enhance_lines, const TFloat64Range &lambdaRange);

    void LoadFitContinuumOneTemplate(const TFloat64Range& lambdaRange, const CTemplate& tpl);
    void LoadFitContinuum(const TFloat64Range& lambdaRange, Int32 icontinuum, Int32 autoSelect);
    void setRedshift(Float64 redshift, bool reinterpolatedContinuum);
    Int32 ApplyContinuumOnGrid(const CTemplate& tpl, Float64 zcontinuum);
    Bool SolveContinuum(const CTemplate& tpl,
                        const TFloat64Range& lambdaRange,
                        const TFloat64List& redshifts,
                        Float64 overlapThreshold,
                        std::vector<CMask> maskList,
                        std::string opt_interp,
                        Int32 opt_extinction,
                        Int32 opt_dustFit,
                        Float64 &merit,
                        Float64& fitAmplitude,
                        Float64& fitAmplitudeError,
                        Bool& fitAmplitudeNegative,
                        Float64 &fitDustCoeff,
                        Int32 &fitMeiksinIdx,
                        Float64& fitDtM,
                        Float64& fitMtM,
                        Float64 &fitLogprior);
    const std::string & getFitContinuum_tplName();
    Float64 getFitContinuum_tplAmplitude();
    Float64 getFitContinuum_tplAmplitudeError();
    Float64 getFitContinuum_snr();
    Float64 getFitContinuum_tplMerit();
    void setFitContinuum_tplAmplitude(Float64 tplAmp, Float64 tplAmpErr, const std::vector<Float64> & polyCoeffs);
    Float64 getFitContinuum_tplIsmDustCoeff();
    Float64 getFitContinuum_tplIgmMeiksinIdx();
    void SetContinuumComponent(std::string component);
    Int32 SetFitContinuum_FitStore(const std::shared_ptr<const CTemplatesFitStore> & fitStore);
    Int32 SetFitContinuum_PriorHelper(const std::shared_ptr<const CPriorHelper> & priorhelper);
    void SetFitContinuum_SNRMax(Float64 snr_max);
    void SetFitContinuum_Option(Int32 opt);
    Int32 GetFitContinuum_Option();
    void SetFitContinuum_FitValues(std::string tplfit_name,
                                   Float64 tplfit_amp,
                                   Float64 tplfit_amperr,
                                   Float64 tplfit_chi2,
                                   Float64 tplfit_ebmv,
                                   Int32 tplfit_meiksinidx,
                                   Float64 tplfit_continuumredshift,
                                   Float64 tplfit_dtm,
                                   Float64 tplfit_mtm,
                                   Float64 tplfit_logprior,
                                   const std::vector<Float64> & polyCoeffs);

    Int32 LoadFitContaminantTemplate(const TFloat64Range& lambdaRange, const CTemplate& tpl);
    std::shared_ptr<CModelSpectrumResult> GetContaminantSpectrumResult();

    Bool initDtd(const TFloat64Range& lambdaRange);
    Float64 EstimateDTransposeD(const TFloat64Range& lambdaRange, const std::string & spcComponent);
    Float64 EstimateMTransposeM(const TFloat64Range& lambdaRange);
    Float64 EstimateLikelihoodCstLog(const TFloat64Range& lambdaRange);
    Float64 getDTransposeD(const TFloat64Range& lambdaRange);
    Float64 getLikelihood_cstLog(const TFloat64Range& lambdaRange);
    Int32 getMTransposeMCumulative(const TFloat64Range& lambdaRange, std::vector<Float64> & lbda, std::vector<Float64> & mtmCumul);

    const std::string & getTplshape_bestTplName();
    Float64 getTplshape_bestTplIsmCoeff();
    Float64 getTplshape_bestAmplitude();
    Float64 getTplshape_bestDtm();
    Float64 getTplshape_bestMtm();
    Int32 getTplshape_count();
    const std::vector<Float64> & getTplshape_priors();
    const std::vector<Float64> & GetChisquareTplshape();
    std::vector<Float64> GetPriorLinesTplshape();
    const std::vector<Float64> & GetScaleMargTplshape();
    const TBoolList & GetStrongELPresentTplshape();
    const std::vector<Int32> & GetNLinesAboveSNRTplshape();
    Int32 SetTplshape_PriorHelper(const std::shared_ptr<const CPriorHelper> & priorhelper);

    Int32 GetNElements();
    Int32 GetModelValidElementsNDdl();
    Int32 GetModelNonZeroElementsNDdl();
    std::vector<UInt32> GetModelValidElementsIndexes();
    bool IsElementIndexInDisabledList(Int32 index);
    void SetElementAmplitude(Int32 j, Float64 a, Float64 snr);
    Float64 GetElementAmplitude(Int32 j);
    void SetElementIndexesDisabledAuto();
    void ResetElementIndexesDisabled();

    void SetVelocityEmission(Float64 vel);
    void SetVelocityAbsorption(Float64 vel);
    void SetVelocityEmissionOneElement(Float64 vel, Int32 idxElt);
    void SetVelocityAbsorptionOneElement(Float64 vel, Int32 idxElt);
    Float64 GetVelocityEmission();
    Float64 GetVelocityAbsorption();
    Float64 GetVelocityInfFromInstrumentResolution();
    Int32 ApplyVelocityBound(Float64 inf, Float64 sup);
    void SetSourcesizeDispersion(Float64 sizeArcsec);
    std::vector<std::vector<Int32>> GetModelVelfitGroups(Int32 lineType );

    Bool initModelAtZ(Float64 redshift, const TFloat64Range& lambdaRange, const CSpectrumSpectralAxis &spectralAxis);

    Float64 fit(Float64 redshift,
                const TFloat64Range& lambdaRange,
                CLineModelSolution &modelSolution,
                CContinuumModelSolution &continuumModelSolution,
                Int32 contreest_iterations=0,
                bool enableLogging=0);

    Bool initTplratioCatalogs(std::string opt_tplratioCatRelPath, Int32 opt_tplratio_ismFit);
    void initLambdaOffsets(std::string offsetsCatalogsRelPath);

    Bool setTplshapeModel(Int32 itplshape, Bool enableSetVelocity=false);
    Bool setTplshapeAmplitude(std::vector<Float64> ampsElts, std::vector<Float64> errorsElts);

    std::vector<CLmfitController*> createLmfitControllers( const TFloat64Range& lambdaRange);
    void fitWithModelSelection(Float64 redshift, const TFloat64Range& lambdaRange, CLineModelSolution &modelSolution);
    void SetFittingMethod(std::string fitMethod);
    void SetSecondpassContinuumFitPrms(Int32 dustfit, Int32 meiksinfit, Int32 outsidelinemask, Int32 observedFrame);

    void SetAbsLinesLimit(Float64 limit);
    void SetLeastSquareFastEstimationEnabled(Int32 enabled);

    Float64 GetRedshift();

    void reinitModel();
    void refreshModel(Int32 lineTypeFilter=-1);
    void getModel(CSpectrumFluxAxis& modelfluxAxis, Int32 lineTypeFilter=-1);
    void reinitModelUnderElements(std::vector<UInt32> filterEltsIdx, Int32 lineIdx=-1);
    void refreshModelInitAllGrid();
    void refreshModelUnderElements(std::vector<UInt32> filterEltsIdx, Int32 lineIdx=-1);
    void refreshModelDerivVelUnderElements(std::vector<UInt32> filterEltsIdx);
    void refreshModelDerivVelAbsorptionUnderElements(std::vector<UInt32> filterEltsIdx);
    void refreshModelDerivVelEmissionUnderElements(std::vector<UInt32> filterEltsIdx);

    Bool addToSpectrumAmplitudeOffset(CSpectrumFluxAxis &modelfluxAxis);
    Int32 prepareAmplitudeOffset(const CSpectrumFluxAxis &spcFlux);


    void setModelSpcObservedOnSupportZeroOutside(const TFloat64Range &lambdaRange);
    CMask getOutsideLinesMask();
    Float64 getOutsideLinesSTD(Int32 which, TFloat64Range lambdarange);


    Int32 getSpcNSamples(const TFloat64Range& lambdaRange);
    Float64 getLeastSquareMeritFast(Int32 idxLine=-1);
    Float64 getLeastSquareContinuumMeritFast();
    Float64 getLeastSquareMerit(const TFloat64Range &lambdaRange);
    Float64 getLeastSquareContinuumMerit(const TFloat64Range &lambdaRange);
    Float64 getLeastSquareMeritUnderElements();
    Float64 getScaleMargCorrection(Int32 idxLine=-1);
    Float64 getContinuumScaleMargCorrection();
    Float64 getStrongerMultipleELAmpCoeff();
    std::vector<std::string> getLinesAboveSNR(Float64 snrcut=3.5);
    Float64 getCumulSNRStrongEL();
    Float64 getCumulSNROnRange( TInt32Range idxRange );
    bool GetModelStrongEmissionLinePresent();
    bool GetModelHaStrongest();
    Float64 getModelErrorUnderElement(UInt32 eltId);
    Float64 getContinuumMeanUnderElement(UInt32 eltId);
    Int32 LoadModelSolution(const CLineModelSolution&  modelSolution);
    CLineModelSolution GetModelSolution(Int32 opt_level=0);
    CContinuumModelSolution GetContinuumModelSolution();
    const CSpectrum&    GetModelSpectrum() const;
    CSpectrum GetSpectrumModelContinuum() const;
    const CSpectrum&    GetObservedSpectrumWithLinesRemoved(Int32 lineTypeFilter=-1);
    Float64 GetContinuumError(Int32 eIdx, Int32 subeIdx);
    Int32 GetFluxDirectIntegration(TInt32List eIdx_list,
                                   TInt32List subeIdx_list,
                                   Int32 opt_cont_substract_abslinesmodel,
                                   Float64& fluxdi,
                                   Float64& snrdi);
    const CSpectrumFluxAxis&    GetModelContinuum() const;
    Float64 getModelFluxVal(UInt32 idx) const;
    Float64 getModelFluxDerivEltVal(UInt32 DerivEltIdx, UInt32 idx) const;
    Float64 getModelFluxDerivContinuumAmpEltVal(UInt32 DerivEltIdx, UInt32 idx) const;
    Float64 getModelFluxDerivZContinuumVal(UInt32 idx)const;
    //void calculateUnscaleContinuumDerivZ();
    Float64 getModelFluxDerivZEltValNoContinuum(UInt32 DerivEltIdx, UInt32 idx) const;
    Float64 getModelFluxDerivZEltVal(UInt32 DerivEltIdx, UInt32 idx, Float64 continuumFluxDerivZ) const;
    Float64 getModelFluxDerivVelVal(UInt32 idx) const;
    Float64 getModelFluxDerivVelEmissionVal(UInt32 idx) const;
    Float64 getModelFluxDerivVelAbsorptionVal(UInt32 idx) const;
    Int32 estimateMeanSqFluxAndGradient(const Float64* varPack,
                                        const Float64 normFactor,
                                        std::vector<UInt32> filteredEltsIdx,
                                        std::vector<UInt32> xInds,
                                        Int32 lineType,
                                        Float64 *fluxdata,
                                        Float64* msqBuffer,
                                        Float64& f,
                                        Float64* g);


    std::vector<std::shared_ptr<CLineModelElement> > m_Elements;
    const CSpectrum & m_inputSpc;
    CSpectrum m_SpectrumModel;  //model
    Int32 FindElementIndex(Int32 LineCatalogIndex);
    Int32 FindElementIndex(std::string LineTagStr, Int32 linetype=-1, Int32& lineIdx=defaultIdx);
    std::vector<UInt32> getOverlappingElements(UInt32 ind , std::vector<UInt32> excludedInd, Float64 overlapThres=0.1);
    CRayCatalog::TRayVector m_RestRayList;
    CSpectrum m_SpcCorrectedUnderLines;  //observed spectrum corrected under the lines

    TStringList GetModelRulesLog();

    Int32 setPassMode(Int32 iPass);
    void SetForcedisableTplratioISMfit(bool opt);

    CRayCatalogsTplShape m_CatalogTplShape;
    std::vector<Float64> m_ChisquareTplshape;
    std::vector<std::vector<Float64>> m_FittedAmpTplshape;
    std::vector<std::vector<Float64>> m_FittedErrorTplshape;
    std::vector<std::vector<Float64>> m_MtmTplshape;
    std::vector<std::vector<Float64>> m_DtmTplshape;
    std::vector<std::vector<Float64>> m_LyaAsymCoeffTplshape;
    std::vector<std::vector<Float64>> m_LyaWidthCoeffTplshape;
    std::vector<std::vector<Float64>> m_LyaDeltaCoeffTplshape;
    std::vector<std::vector<Float64>> m_LinesLogPriorTplshape;

    bool m_enableAmplitudeOffsets;
    Float64 m_LambdaOffsetMin = -400.0;
    Float64 m_LambdaOffsetMax = 400.0;
    Float64 m_LambdaOffsetStep = 25.0;
    bool m_enableLambdaOffsetsFit;

    bool m_opt_lya_forcefit=false;
    bool m_opt_lya_forcedisablefit=false;
    Float64 m_opt_lya_fit_asym_min=0.0;
    Float64 m_opt_lya_fit_asym_max=4.0;
    Float64 m_opt_lya_fit_asym_step=1.0;
    Float64 m_opt_lya_fit_width_min=1.;
    Float64 m_opt_lya_fit_width_max=4.;
    Float64 m_opt_lya_fit_width_step=1.;
    Float64 m_opt_lya_fit_delta_min=0.;
    Float64 m_opt_lya_fit_delta_max=0.;
    Float64 m_opt_lya_fit_delta_step=1.;

    Int32 m_opt_fitcontinuum_maxCount = 2;
    bool m_opt_firstpass_forcedisableMultipleContinuumfit=true;
    bool m_opt_firstpass_forcedisableTplratioISMfit=true;
    std::string m_opt_firstpass_fittingmethod = "hybrid";
    std::string m_opt_secondpass_fittingmethod = "hybrid";

    bool m_opt_enable_improveBalmerFit=false;

private:

    Int32 fitAmplitudesHybrid(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& spcFluxAxisNoContinuum, const CSpectrumFluxAxis &continuumfluxAxis, Float64 redshift);
    void fitAmplitudesSimplex();
    Int32 fitAmplitudesLmfit( const CSpectrumFluxAxis& fluxAxis, CLmfitController * controller);
    Int32 fitAmplitudesLinSolve(std::vector<UInt32> EltsIdx, const CSpectrumSpectralAxis &spectralAxis, const CSpectrumFluxAxis &fluxAxis, const CSpectrumFluxAxis& continuumfluxAxis, std::vector<Float64> &ampsfitted, std::vector<Float64> &errorsfitted);
    Int32 fitAmplitudesLinSolveAndLambdaOffset(std::vector<UInt32> EltsIdx, const CSpectrumSpectralAxis &spectralAxis, const CSpectrumFluxAxis &fluxAxis, const CSpectrumFluxAxis& continuumfluxAxis, std::vector<Float64> &ampsfitted, std::vector<Float64> &errorsfitted, Bool enableOffsetFitting);

    Int32 fitAmplitudesLinesAndContinuumLinSolve(std::vector<UInt32> EltsIdx,
                                                 const TFloat64Range& lambdaRange,
                                                 const CSpectrumSpectralAxis& spectralAxis,
                                                 const CSpectrumFluxAxis& fluxAxis,
                                                 const CSpectrumFluxAxis& continuumfluxAxis,
                                                 std::vector<Float64>& ampsfitted,
                                                 std::vector<Float64>& errorsfitted,
                                                 Float64 &chisquare,
                                                 Int32 polyOrder=-1);

    bool m_forceDisableLyaFitting=false;
    bool m_forceLyaFitting=false;
    Int32 setLyaProfile( Float64 redshift, const CSpectrumSpectralAxis& spectralAxis );

    std::vector<UInt32> getSupportIndexes(std::vector<UInt32> EltsIdx);
    Float64 GetWeightingAnyLineCenterProximity(UInt32 sampleIndex, std::vector<UInt32> EltsIdx);
    std::vector<UInt32> getOverlappingElementsBySupport(UInt32 ind , Float64 overlapThres=0.1);
    std::vector<UInt32> ReestimateContinuumApprox(std::vector<UInt32> EltsIdx);
    std::vector<UInt32> ReestimateContinuumUnderLines(std::vector<UInt32> EltsIdx);
    void refreshModelAfterContReestimation(std::vector<UInt32> EltsIdx, CSpectrumFluxAxis& modelFluxAxis, CSpectrumFluxAxis& spcFluxAxisNoContinuum);

    std::vector<UInt32> findLineIdxInCatalog(const CRayCatalog::TRayVector& restRayList, std::string strTag, Int32 type);


    void addDoubleLine(const CRay &r1, const CRay &r2, Int32 index1, Int32 index2, Float64 nominalWidth, Float64 a1, Float64 a2);

    Int32 improveBalmerFit();
    void applyRules(bool enableLogs=false);
    CRegulament m_Regulament;


    std::vector<Float64> m_ScaleMargCorrTplshape;
    TBoolList m_StrongELPresentTplshape;
    std::vector<Int32> m_NLinesAboveSNRTplshape;

    Float64 m_Redshift;

    CSpectrumFluxAxis m_SpcFluxAxis;    //observed spectrum
    CSpectrumFluxAxis m_spcFluxAxisNoContinuum; //observed spectrum for line fitting
    CTemplate m_tplContaminantSpcRebin; //optionally used contaminant to be removed from observed spectrum
    TFloat64List& m_ErrorNoContinuum;
    CSpectrumFluxAxis m_SpcFluxAxisModelDerivVelEmi;
    CSpectrumFluxAxis m_SpcFluxAxisModelDerivVelAbs;
    Float64 m_dTransposeD; //the cached dtd (maximum chisquare value)
    TFloat64Range m_dTransposeDLambdaRange; //the lambdaRange used to computed cached dTransposeD values
    Float64 m_likelihood_cstLog; // constant term for the Likelihood calculation

    TAxisSampleList m_observeGridContinuumFlux;   //the continuum spectre without the amplitude coeff; m_ContinuumFLux = amp * m_observeGridContinuumFlux
    //Float64* m_unscaleContinuumFluxAxisDerivZ;
    CSpectrumFluxAxis m_ContinuumFluxAxis;  //rebined model continuum
    std::string m_ContinuumComponent;
    std::string m_LineWidthType;
    Float64 m_NSigmaSupport;
    Float64 m_resolution;
    Float64 m_velocityEmission;
    Float64 m_velocityAbsorption;
    Float64 m_velocityEmissionInit;
    Float64 m_velocityAbsorptionInit;

    Float64 m_nominalWidthDefaultEmission;
    Float64 m_nominalWidthDefaultAbsorption;

    std::string m_calibrationPath;
    std::string m_fittingmethod;
    std::vector<Int32> m_elementsDisabledIndexes;
    std::string m_rulesoption;
    std::string m_rigidity;
    bool m_forcedisableTplratioISMfit=false;

    CTemplateCatalog m_tplCatalog;
    CTemplateCatalog m_orthoTplCatalog;
    TStringList m_tplCategoryList;
    std::string m_tplshapeBestTplName;
    Float64 m_tplshapeBestTplIsmCoeff;
    Float64 m_tplshapeBestTplAmplitude;
    Float64 m_tplshapeBestTplDtm;
    Float64 m_tplshapeBestTplMtm;
    Int32 m_tplshapeLeastSquareFast = 0;    //for rigidity=tplshape: switch to use fast least square estimation
    std::shared_ptr<const CPriorHelper> m_tplshape_priorhelper;

    COperatorChiSquare2 m_chiSquareOperator;
    Int32 m_secondpass_fitContinuum_dustfit;
    Int32 m_secondpass_fitContinuum_igm;
    Int32 m_secondpass_fitContinuum_outsidelinesmask;
    Int32 m_secondpass_fitContinuum_observedFrame;

    std::shared_ptr<const CTemplatesFitStore> m_fitContinuum_tplfitStore;
    Int32 m_fitContinuum_option;
    std::string m_fitContinuum_tplName;
    Float64 m_fitContinuum_tplFitAmplitude=-1.0;
    Float64 m_fitContinuum_tplFitAmplitudeError;
    Float64 m_fitContinuum_tplFitMerit;
    Float64 m_fitContinuum_tplFitDustCoeff;
    Int32 m_fitContinuum_tplFitMeiksinIdx;
    Float64 m_fitContinuum_tplFitRedshift; // only used with m_fitContinuum_option==2 for now
    Float64 m_fitContinuum_tplFitDtM;
    Float64 m_fitContinuum_tplFitMtM;
    Float64 m_fitContinuum_tplFitLogprior;
    Float64 m_fitContinuum_tplFitSNRMax=0.0;
    std::vector<Float64> m_fitContinuum_tplFitPolyCoeffs;   // only used with m_fitContinuum_option==2 for now
    bool m_forcedisableMultipleContinuumfit=false;
    Float64 m_fitContinuum_tplFitAlpha=0.;
    std::shared_ptr<const CPriorHelper> m_fitContinuum_priorhelper;

    bool m_lmfit_noContinuumTemplate;
    bool m_lmfit_bestTemplate;
    bool m_lmfit_fitContinuum;
    bool m_lmfit_fitEmissionVelocity;
    bool m_lmfit_fitAbsorptionVelocity;

    std::vector<Float64> m_ampOffsetsX0;
    std::vector<Float64> m_ampOffsetsX1;
    std::vector<Float64> m_ampOffsetsX2;
    std::vector<Int32> m_ampOffsetsIdxStart;
    std::vector<Int32> m_ampOffsetsIdxStop;
};

}


#endif // _REDSHIFT_LINEMODEL_ELEMENTLIST_
