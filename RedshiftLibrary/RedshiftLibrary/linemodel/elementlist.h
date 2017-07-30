#ifndef ELEMENTLIST_H
#define ELEMENTLIST_H

#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/common/datatypes.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <RedshiftLibrary/ray/catalog.h>
#include <RedshiftLibrary/spectrum/spectrum.h>

#include <RedshiftLibrary/operator/chisquare2.h>

#include <RedshiftLibrary/operator/linemodelresult.h>
#include <RedshiftLibrary/linemodel/element.h>

#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/linemodel/templatesfitstore.h>

#include <boost/shared_ptr.hpp>

#include <RedshiftLibrary/linemodel/lmfitcontroller.h>

#include <memory>
namespace NSEpic
{
  static Int32 defaultIdx = -1;
  class CRegulament;
  class CRayCatalogsTplShape;

class CLineModelElementList
{

public:

    CLineModelElementList(const CSpectrum& spectrum,
                          const CSpectrum& spectrumNoContinuum,
                          const CTemplateCatalog& tplCatalog,
                          const TStringList& tplCategoryList,
                          const std::string calibrationPath,
                          const CRayCatalog::TRayVector& restRayList,
                          const std::string& opt_fittingmethod,
                          const std::string &opt_continuumcomponent,
                          const std::string& lineWidthType,
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

    void PrepareContinuum(Float64 z);
    void EstimateSpectrumContinuum(Float64 opt_enhance_lines=0);

    Int32 LoadFitContinuumOneTemplate(const TFloat64Range& lambdaRange, const CTemplate& tpl);
    Int32 LoadFitContinuum(const TFloat64Range& lambdaRange);
    void setRedshift(Float64 redshift, bool reinterpolatedContinuum);
    Int32 ApplyContinuumOnGrid(const CTemplate& tpl);
    Bool SolveContinuum(const CSpectrum& spectrum,
                        const CTemplate& tpl,
                        const TFloat64Range& lambdaRange,
                        const TFloat64List& redshifts,
                        Float64 overlapThreshold,
                        std::vector<CMask> maskList,
                        std::string opt_interp,
                        Int32 opt_extinction,
                        Int32 opt_dustFit,
                        Float64 &merit,
                        Float64& fitAmplitude,
                        Float64 &fitDustCoeff,
                        Int32 &fitMeiksinIdx,
                        Float64& fitDtM,
                        Float64& fitMtM);
    std::string getFitContinuum_tplName();
    Float64 getFitContinuum_tplAmplitude();
    Float64 getFitContinuum_tplMerit();
    void setFitContinuum_tplAmplitude(Float64 tplAmp);
    Float64 getFitContinuum_tplIsmDustCoeff();
    Float64 getFitContinuum_tplIgmMeiksinIdx();
    Float64* getPrecomputedGridContinuumFlux();
    void SetContinuumComponent(std::string component);
    Int32 SetFitContinuum_FitStore(CTemplatesFitStore* fitStore);

    Bool initDtd(const TFloat64Range& lambdaRange);
    Float64 EstimateDTransposeD(const TFloat64Range& lambdaRange, std::string spcComponent);
    Float64 EstimateMTransposeM(const TFloat64Range& lambdaRange);
    Float64 EstimateLikelihoodCstLog(const TFloat64Range& lambdaRange);
    Float64 getDTransposeD(const TFloat64Range& lambdaRange, std::string spcComponent);
    Float64 getLikelihood_cstLog(const TFloat64Range& lambdaRange);
    Int32 getMTransposeMCumulative(const TFloat64Range& lambdaRange, std::vector<Float64> lbda, std::vector<Float64> mtmCumul);

    std::string getTplshape_bestTplName();
    Int32 getTplshape_count();
    std::vector<Float64> GetChisquareTplshape();

    Int32 GetNElements();
    Int32 GetModelValidElementsNDdl();
    Int32 GetModelNonZeroElementsNDdl();
    std::vector<Int32> GetModelValidElementsIndexes();
    bool IsElementIndexInDisabledList(Int32 index);
    void SetElementAmplitude(Int32 j, Float64 a, Float64 snr);
    Float64 GetElementAmplitude(Int32 j);
    void SetElementIndexesDisabledAuto();
    void ResetElementIndexesDisabled();

    void SetVelocityEmission(Float64 vel);
    void SetVelocityAbsorption(Float64 vel);
    Float64 GetVelocityEmission();
    Float64 GetVelocityAbsorption();
    Float64 GetVelocityInfFromInstrumentResolution();
    Int32 ApplyVelocityBound(Float64 inf, Float64 sup);

    Bool initModelAtZ(Float64 redshift, const TFloat64Range& lambdaRange, const CSpectrumSpectralAxis &spectralAxis);

    Float64 fit(Float64 redshift, const TFloat64Range& lambdaRange, CLineModelSolution &modelSolution, Int32 contreest_iterations=0, bool enableLogging=0);

    std::vector<CLmfitController*> createLmfitControllers( const TFloat64Range& lambdaRange);
    void fitWithModelSelection(Float64 redshift, const TFloat64Range& lambdaRange, CLineModelSolution &modelSolution);
    void SetFittingMethod(std::string fitMethod);

    void SetAbsLinesLimit(Float64 limit);

    Float64 GetRedshift();

    void reinitModel();
    void refreshModel();
    void reinitModelUnderElements(std::vector<Int32> filterEltsIdx, Int32 lineIdx=-1 );
    void refreshModelInitAllGrid();
    void refreshModelUnderElements(std::vector<Int32> filterEltsIdx, Int32 lineIdx=-1 );
    void refreshModelDerivVelUnderElements(std::vector<Int32> filterEltsIdx);
    void refreshModelDerivVelAbsorptionUnderElements(std::vector<Int32> filterEltsIdx);
    void refreshModelDerivVelEmissionUnderElements(std::vector<Int32> filterEltsIdx);

    void setModelSpcObservedOnSupportZeroOutside(const TFloat64Range &lambdaRange);
    CMask getOutsideLinesMask();


    Int32 getSpcNSamples(const TFloat64Range& lambdaRange);
    Float64 getLeastSquareMeritFast(Int32 idxLine=-1);
    Float64 getLeastSquareContinuumMeritFast();
    Float64 getLeastSquareMerit(const TFloat64Range &lambdaRange);
    Float64 getLeastSquareMeritUnderElements();
    Float64 getScaleMargCorrection(Int32 idxLine=-1);
    Float64 getStrongerMultipleELAmpCoeff();
    Float64 getCumulSNRStrongEL();
    Float64 getCumulSNROnRange( TInt32Range idxRange );
    Float64 getModelErrorUnderElement(Int32 eltId);
    Float64 getContinuumMeanUnderElement(Int32 eltId);
    Int32 LoadModelSolution(const CLineModelSolution&  modelSolution);
    CLineModelSolution GetModelSolution();
    const CSpectrum&    GetModelSpectrum() const;
    const CSpectrum&    GetObservedSpectrumWithLinesRemoved() const;
    Float64 GetContinuumError(Int32 eIdx, Int32 subeIdx);
    Float64 GetFluxDirectIntegration(Int32 eIdx, Int32 subeIdx);
    const CSpectrumFluxAxis&    GetModelContinuum() const;
    Float64 getModelFluxVal(Int32 idx) const;
    Float64 getModelFluxDerivEltVal(Int32 DerivEltIdx, Int32 idx) const;
    Float64 getModelFluxDerivContinuumAmpEltVal(Int32 DerivEltIdx, Int32 idx) const;
    Float64 getModelFluxDerivZContinuumVal(Int32 idx)const;
    //void calculateUnscaleContinuumDerivZ();
    Float64 getModelFluxDerivZEltValNoContinuum(Int32 DerivEltIdx, Int32 idx) const;
    Float64 getModelFluxDerivZEltVal(Int32 DerivEltIdx, Int32 idx, Float64 continuumFluxDerivZ) const;
    Float64 getModelFluxDerivVelVal(Int32 idx) const;
    Float64 getModelFluxDerivVelEmissionVal(Int32 idx) const;
    Float64 getModelFluxDerivVelAbsorptionVal(Int32 idx) const;
    Int32 estimateMeanSqFluxAndGradient(const Float64* varPack,
                                        const Float64 normFactor,
                                        std::vector<Int32> filteredEltsIdx,
                                        std::vector<Int32> xInds,
                                        Int32 lineType,
                                        Float64 *fluxdata,
                                        Float64* msqBuffer,
                                        Float64& f,
                                        Float64* g);


    std::vector<boost::shared_ptr<CLineModelElement>  > m_Elements;
    std::shared_ptr<CSpectrum>  m_SpectrumModel;  //model
    CSpectrumFluxAxis m_SpcContinuumFluxAxis; //continuum spectrum used for the model
    Int32 FindElementIndex(Int32 LineCatalogIndex);
    Int32 FindElementIndex(std::string LineTagStr, Int32 linetype=-1, Int32& lineIdx=defaultIdx);
    std::vector<Int32> getOverlappingElements(Int32 ind , std::vector<Int32> excludedInd, Float64 overlapThres=0.1);
    CRayCatalog::TRayVector m_RestRayList;
    std::shared_ptr<CSpectrum>  m_SpcCorrectedUnderLines;  //observed spectrum corrected under the lines

    TStringList GetModelRulesLog();

    Int32 setPassMode(Int32 iPass);

private:

    Int32 fitAmplitudesHybrid(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& spcFluxAxisNoContinuum, const CSpectrumFluxAxis &continuumfluxAxis, Float64 redshift);
    void fitAmplitudesSimplex();
    Int32 fitAmplitudesLmfit( const CSpectrumFluxAxis& fluxAxis, CLmfitController * controller);
    Int32 fitAmplitudesLinSolve(std::vector<Int32> EltsIdx, const CSpectrumSpectralAxis &spectralAxis, const CSpectrumFluxAxis &fluxAxis, const CSpectrumFluxAxis& continuumfluxAxis, std::vector<Float64> &ampsfitted, std::vector<Float64> &errorsfitted);
    Int32 fitAmplitudesLBFGS(std::vector<Int32> filteredEltsIdx, const CSpectrumFluxAxis& fluxAxis, std::vector<Float64>& ampsfitted, Int32 lineType);

    bool m_forceDisableLyaFitting;
    Int32 setLyaProfile( Float64 redshift, const CSpectrumSpectralAxis& spectralAxis );

    std::vector<Int32> getSupportIndexes(std::vector<Int32> EltsIdx);
    std::vector<Int32> getOverlappingElementsBySupport(Int32 ind , Float64 overlapThres=0.1);
    std::vector<Int32> ReestimateContinuumApprox(std::vector<Int32> EltsIdx);
    std::vector<Int32> ReestimateContinuumUnderLines(std::vector<Int32> EltsIdx);
    void refreshModelAfterContReestimation(std::vector<Int32> EltsIdx, CSpectrumFluxAxis& modelFluxAxis, CSpectrumFluxAxis& spcFluxAxisNoContinuum);

    std::vector<Int32> findLineIdxInCatalog(const CRayCatalog::TRayVector& restRayList, std::string strTag, Int32 type);


    void addDoubleLine(const CRay &r1, const CRay &r2, Int32 index1, Int32 index2, Float64 nominalWidth, Float64 a1, Float64 a2);

    Int32 improveBalmerFit();
    void applyRules(bool enableLogs=false);
    CRegulament* m_Regulament;

    CRayCatalogsTplShape* m_CatalogTplShape;
    std::vector<Float64> m_ChisquareTplshape;

    Float64 m_Redshift;

    CSpectrumFluxAxis m_SpcFluxAxis;    //observed spectrum
    CSpectrumFluxAxis m_spcFluxAxisNoContinuum; //observed spectrum for line fitting
    Float64* m_ErrorNoContinuum;
    CSpectrumFluxAxis m_SpcFluxAxisModelDerivVelEmi;
    CSpectrumFluxAxis m_SpcFluxAxisModelDerivVelAbs;
    Float64 m_dTransposeDNocontinuum; //the cached dtd (maximum chisquare value)
    Float64 m_dTransposeDRaw; //the cached dtd (maximum chisquare value)
    TFloat64Range m_dTransposeDLambdaRange; //the lambdaRange used to computed cached dTransposeD values
    Float64 m_likelihood_cstLog; // constant term for the Likelihood calculation

    Float64*          m_observeGridContinuumFlux;   //the continuum spectre without the amplitude coeff; m_ContinuumFLux = amp * m_observeGridContinuumFlux
    Float64* m_unscaleContinuumFluxAxisDerivZ;
    CSpectrumFluxAxis m_ContinuumFluxAxis;  //rebined model continuum
    Float64 m_ContinuumWinsize;
    std::string m_ContinuumComponent;
    std::string m_LineWidthType;
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

    std::shared_ptr<CSpectrum> m_inputSpc;
    CTemplateCatalog m_tplCatalog;
    TStringList m_tplCategoryList;
    std::string m_tplshapeBestTplName;

    COperatorChiSquare2* m_chiSquareOperator;
    Int32 m_fitContinuum_dustfit;
    Int32 m_fitContinuum_igm;
    Int32 m_fitContinuum_outsidelinesmask;
    Int32 m_fitContinuum_observedFrame;

    CTemplatesFitStore* m_fitContinuum_tplfitStore;
    Int32 m_fitContinuum_option;
    std::string m_fitContinuum_tplName;
    Float64 m_fitContinuum_tplFitAmplitude;
    Float64 m_fitContinuum_tplFitMerit;
    Float64 m_fitContinuum_tplFitDustCoeff;
    Int32 m_fitContinuum_tplFitMeiksinIdx;
    Float64 m_fitContinuum_tplFitDtM;
    Float64 m_fitContinuum_tplFitMtM;

    bool m_lmfit_noContinuumTemplate;
    bool m_lmfit_bestTemplate;
    bool m_lmfit_fitContinuum;
    bool m_lmfit_fitEmissionVelocity;
    bool m_lmfit_fitAbsorptionVelocity;
};

}







#endif // ELEMENTLIST_H
