// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#ifndef _REDSHIFT_LINEMODEL_FITTING_
#define _REDSHIFT_LINEMODEL_FITTING_

#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/common/datatypes.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "RedshiftLibrary/ray/catalog.h"
#include "RedshiftLibrary/ray/catalogsTplShape.h"
#include "RedshiftLibrary/ray/regulament.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

#include "RedshiftLibrary/operator/templatefitting.h"

#include "RedshiftLibrary/operator/linemodelresult.h"
#include "RedshiftLibrary/operator/modelspectrumresult.h"
#include "RedshiftLibrary/linemodel/element.h"

#include "RedshiftLibrary/operator/pdfz.h"

#include "RedshiftLibrary/spectrum/template/catalog.h"
#include "RedshiftLibrary/linemodel/templatesfitstore.h"

#include <boost/shared_ptr.hpp>

#include "RedshiftLibrary/linemodel/elementlist.h"
#include "RedshiftLibrary/linemodel/lmfitcontroller.h"
#include "RedshiftLibrary/ray/linetags.h"

#include <memory>
namespace NSEpic
{

  //to avoid circular dependency
  class CLineProfileASYM;
  class CLineProfileASYMFIXED;

class CLineModelFitting
{

public:

    CLineModelFitting(  const CSpectrum& spectrum,
                        const TFloat64Range& lambdaRange,
                        const CTemplateCatalog& tplCatalog,
                        const TStringList& tplCategoryList,
                        const CRayCatalog::TRayVector& restRayList,
                        const std::string& opt_fittingmethod,
                        const std::string &opt_continuumcomponent,
                        const Float64 opt_continuum_neg_threshold,
                        const std::string& lineWidthType,
                        const Float64 nsigmasupport,
                        const Float64 velocityEmission,
                        const Float64 velocityAbsorption,
                        const std::string &opt_rules,
                        const std::string &opt_rigidity);

    void LoadCatalog(const CRayCatalog::TRayVector& restRayList);
    void LoadCatalogOneMultiline(const CRayCatalog::TRayVector& restRayList);
    void LoadCatalogTwoMultilinesAE(const CRayCatalog::TRayVector& restRayList);

    void LogCatalogInfos();

    void PrepareContinuum();
    void EstimateSpectrumContinuum(Float64 opt_enhance_lines, const TFloat64Range &lambdaRange);

    void LoadFitContinuumOneTemplate(const TFloat64Range& lambdaRange, const std::shared_ptr<const CTemplate> & tpl);
    void LoadFitContinuum(const TFloat64Range& lambdaRange, Int32 icontinuum, Int32 autoSelect);
    void setRedshift(Float64 redshift, bool reinterpolatedContinuum);
    Int32 ApplyContinuumOnGrid(const std::shared_ptr<const CTemplate>& tpl, Float64 zcontinuum);
    bool SolveContinuum(const std::shared_ptr<const CTemplate>& tpl,
                        const TFloat64List& redshifts,
                        Float64 overlapThreshold,
                        std::vector<CMask> maskList,
                        std::string opt_interp,
                        Int32 opt_extinction,
                        Int32 opt_dustFit,
                        Float64 &merit,
                        Float64& fitAmplitude,
                        Float64& fitAmplitudeError,
                        Float64& fitAmplitudeSigma,
                        Float64 &fitEbmvCoeff,
                        Int32 &fitMeiksinIdx,
                        Float64& fitDtM,
                        Float64& fitMtM,
                        Float64 &fitLogprior);
    const std::string & getFitContinuum_tplName() const;
    Float64 getFitContinuum_tplAmplitude() const;
    Float64 getFitContinuum_tplAmplitudeError() const;
    Float64 getFitContinuum_snr() const;
    Float64 getFitContinuum_tplMerit() const;
    Float64 getFitContinuum_tplMeritPhot() const;
    void setFitContinuum_tplAmplitude(Float64 tplAmp, Float64 tplAmpErr, const TFloat64List & polyCoeffs);
    Float64 getFitContinuum_tplIsmEbmvCoeff() const;
    Float64 getFitContinuum_tplIgmMeiksinIdx() const;
    void SetContinuumComponent(std::string component);
    Int32 SetFitContinuum_FitStore(const std::shared_ptr<const CTemplatesFitStore> & fitStore);
    const std::shared_ptr<const CTemplatesFitStore> & GetFitContinuum_FitStore() const;
    Int32 SetFitContinuum_PriorHelper(const std::shared_ptr<const CPriorHelper> & priorhelper);
    void SetFitContinuum_SNRMax(Float64 snr_max);
    void SetFitContinuum_Option(Int32 opt);
    Int32 GetFitContinuum_Option() const;
    void SetFitContinuum_FitValues(std::string tplfit_name,
                                   Float64 tplfit_amp,
                                   Float64 tplfit_amperr,
                                   Float64 tplfit_chi2,
                                   Float64 tplfit_chi2_phot,
                                   Float64 tplfit_ebmv,
                                   Int32 tplfit_meiksinidx,
                                   Float64 tplfit_continuumredshift,
                                   Float64 tplfit_dtm,
                                   Float64 tplfit_mtm,
                                   Float64 tplfit_logprior,
                                   const TFloat64List & polyCoeffs);

    Int32 LoadFitContaminantTemplate(const TFloat64Range& lambdaRange, const CTemplate& tpl);
    std::shared_ptr<CModelSpectrumResult> GetContaminantSpectrumResult() const;

    bool initDtd(const TFloat64Range& lambdaRange);
    Float64 EstimateDTransposeD(const TFloat64Range& lambdaRange, const std::string & spcComponent) const;
    Float64 EstimateMTransposeM(const TFloat64Range& lambdaRange) const;
    Float64 EstimateLikelihoodCstLog(const TFloat64Range& lambdaRange) const;
    Float64 getDTransposeD(const TFloat64Range& lambdaRange);
    Float64 getLikelihood_cstLog(const TFloat64Range& lambdaRange);
    Int32 getMTransposeMCumulative(const TFloat64Range& lambdaRange, TFloat64List & lbda, TFloat64List & mtmCumul) const;

    const std::string & getTplshape_bestTplName() const ;
    Float64 getTplshape_bestTplIsmCoeff() const ;
    Float64 getTplshape_bestAmplitude() const ;
    Float64 getTplshape_bestDtm() const ;
    Float64 getTplshape_bestMtm() const ;
    Int32 getTplshape_count() const;
    const TFloat64List & getTplshape_priors();
    const TFloat64List & GetChisquareTplshape() const ;
    TFloat64List GetPriorLinesTplshape() const ;
    const TFloat64List & GetScaleMargTplshape() const ;
    const TBoolList & GetStrongELPresentTplshape() const ;
    const TBoolList & getHaELPresentTplshape() const;
    const TInt32List & GetNLinesAboveSNRTplshape() const ;
    Int32 SetTplshape_PriorHelper(const std::shared_ptr<const CPriorHelper> & priorhelper);

    Int32 GetNElements() const;


    void SetVelocityEmission(Float64 vel);
    void SetVelocityAbsorption(Float64 vel);
    void SetVelocityEmissionOneElement(Float64 vel, Int32 idxElt);
    void SetVelocityAbsorptionOneElement(Float64 vel, Int32 idxElt);

  void setVelocity(Float64 vel,Int32 rayType);
  void setVelocity(Float64 vel, Int32 idxElt,Int32 rayType);

    Float64 GetVelocityEmission() const;
    Float64 GetVelocityAbsorption() const;
    Int32 ApplyVelocityBound(Float64 inf, Float64 sup);


    bool initModelAtZ(Float64 redshift, const TFloat64Range& lambdaRange, const CSpectrumSpectralAxis &spectralAxis);

    Float64 fit(Float64 redshift,
                const TFloat64Range& lambdaRange,
                CLineModelSolution &modelSolution,
                CContinuumModelSolution &continuumModelSolution,
                Int32 contreest_iterations=0,
                bool enableLogging=0);
    TFloat64Range& getLambdaRange(){return m_dTransposeDLambdaRange;};
    bool initTplratioCatalogs(Int32 opt_tplratio_ismFit);
   
    bool setTplshapeModel(Int32 itplshape, bool enableSetVelocity=false);
    bool setTplshapeAmplitude(const TFloat64List & ampsElts, const TFloat64List & errorsElts);

    std::vector<CLmfitController*> createLmfitControllers( const TFloat64Range& lambdaRange);
    void SetFittingMethod(const std::string &fitMethod);
    void SetSecondpassContinuumFitPrms(Int32 dustfit, Int32 meiksinfit, Int32 outsidelinemask, Int32 observedFrame);

    void SetAbsLinesLimit(Float64 limit);
    void SetLeastSquareFastEstimationEnabled(Int32 enabled);

    Float64 GetRedshift() const;

    void reinitModel();
    void refreshModel(Int32 lineTypeFilter=-1);
    CSpectrumFluxAxis getModel(Int32 lineTypeFilter=-1) const;
    void reinitModelUnderElements(const TInt32List & filterEltsIdx, Int32 lineIdx=-1);
    void refreshModelInitAllGrid();
    void refreshModelUnderElements(const TInt32List & filterEltsIdx, Int32 lineIdx=-1);
    void refreshModelDerivVelUnderElements(const TInt32List & filterEltsIdx);
    void refreshModelDerivVelAbsorptionUnderElements(const TInt32List & filterEltsIdx);
    void refreshModelDerivVelEmissionUnderElements(const TInt32List & filterEltsIdx);

   
    void setModelSpcObservedOnSupportZeroOutside(const TFloat64Range &lambdaRange);
    CMask getOutsideLinesMask() const;
    Float64 getOutsideLinesSTD(Int32 which, const TFloat64Range &lambdarange) const;


    Int32 getSpcNSamples(const TFloat64Range& lambdaRange) const;
    Float64 getLeastSquareMeritFast(Int32 idxLine=-1) const;
    Float64 getLeastSquareContinuumMeritFast() const;
    Float64 getLeastSquareMerit(const TFloat64Range &lambdaRange) const;
    Float64 getLeastSquareContinuumMerit(const TFloat64Range &lambdaRange) const;
    Float64 getLeastSquareMeritUnderElements() const;
    Float64 getScaleMargCorrection(Int32 idxLine=-1) const;
    Float64 getContinuumScaleMargCorrection() const;
    Float64 getStrongerMultipleELAmpCoeff() const;
    TStringList getLinesAboveSNR(Float64 snrcut=3.5) const;
    Float64 getCumulSNRStrongEL() const;
    Float64 getCumulSNROnRange( TInt32Range idxRange ) const;
    bool GetModelStrongEmissionLinePresent() const;
    bool GetModelHaStrongest() const;
   
    Float64 getContinuumMeanUnderElement(Int32 eltId) const;
    void LoadModelSolution(const CLineModelSolution&  modelSolution);
    CLineModelSolution GetModelSolution(Int32 opt_level=0);
    CContinuumModelSolution GetContinuumModelSolution() const;
    const CSpectrum&    GetModelSpectrum() const;
    CSpectrum GetSpectrumModelContinuum() const;
    const CSpectrum&    GetObservedSpectrumWithLinesRemoved(Int32 lineTypeFilter=-1);
    Float64 GetContinuumError(Int32 eIdx, Int32 subeIdx);
    Int32 GetFluxDirectIntegration(const TInt32List & eIdx_list,
                                   const TInt32List & subeIdx_list,
                                   Int32 opt_cont_substract_abslinesmodel,
                                   Float64& fluxdi,
                                   Float64& snrdi) const;
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
                                        const TInt32List & filteredEltsIdx,
                                        const TInt32List & xInds,
                                        Int32 lineType,
                                        Float64 *fluxdata,
                                        Float64* msqBuffer,
                                        Float64& f,
                                        Float64* g);

    CLineModelElementList m_Elements;
    const CSpectrum & m_inputSpc;
    CSpectrum m_SpectrumModel;  //model
    const CRayCatalog::TRayVector& m_RestRayList;
    CSpectrum m_SpcCorrectedUnderLines;  //observed spectrum corrected under the lines

    const TStringList & GetModelRulesLog() const;

    Int32 setPassMode(Int32 iPass);
    Int32 GetPassNumber() const;

    void SetForcedisableTplratioISMfit(bool opt);

    CRayCatalogsTplShape m_CatalogTplShape;
    TFloat64List m_ChisquareTplshape;
    std::vector<TFloat64List> m_FittedAmpTplshape;
    std::vector<TFloat64List> m_FittedErrorTplshape;
    std::vector<TFloat64List> m_MtmTplshape;
    std::vector<TFloat64List> m_DtmTplshape;
    std::vector<TFloat64List> m_LyaAsymCoeffTplshape;
    std::vector<TFloat64List> m_LyaWidthCoeffTplshape;
    std::vector<TFloat64List> m_LyaDeltaCoeffTplshape;
    std::vector<TFloat64List> m_LinesLogPriorTplshape;

    Int32 m_pass = 1;
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
    Float64 m_opt_fitcontinuum_neg_threshold = - INFINITY;
    bool m_opt_firstpass_forcedisableMultipleContinuumfit=true;
    bool m_opt_firstpass_forcedisableTplratioISMfit=true;
    std::string m_opt_firstpass_fittingmethod = "hybrid";
    std::string m_opt_secondpass_fittingmethod = "hybrid";

    bool m_opt_enable_improveBalmerFit=false;
    Float64 m_opt_haprior = -1.;
    void setHaPriorOption(Float64 opt){m_opt_haprior = opt;};
private:

    Int32 fitAmplitudesHybrid(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& spcFluxAxisNoContinuum, const CSpectrumFluxAxis &continuumfluxAxis, Float64 redshift);
    void fitAmplitudesSimplex();
    Int32 fitAmplitudesLmfit( const CSpectrumFluxAxis& fluxAxis, CLmfitController * controller);
    Int32 fitAmplitudesLinSolve(const TInt32List & EltsIdx, const CSpectrumSpectralAxis &spectralAxis, const CSpectrumFluxAxis &fluxAxis, const CSpectrumFluxAxis& continuumfluxAxis, TFloat64List &ampsfitted, TFloat64List &errorsfitted);
    Int32 fitAmplitudesLinSolveAndLambdaOffset(TInt32List EltsIdx, const CSpectrumSpectralAxis &spectralAxis, const CSpectrumFluxAxis &fluxAxis, const CSpectrumFluxAxis& continuumfluxAxis, TFloat64List &ampsfitted, TFloat64List &errorsfitted, bool enableOffsetFitting);

    Int32 fitAmplitudesLinesAndContinuumLinSolve(const TInt32List & EltsIdx,
                                                 const TFloat64Range& lambdaRange,
                                                 const CSpectrumSpectralAxis& spectralAxis,
                                                 const CSpectrumFluxAxis& fluxAxis,
                                                 const CSpectrumFluxAxis& continuumfluxAxis,
                                                 TFloat64List& ampsfitted,
                                                 TFloat64List& errorsfitted,
                                                 Float64 &chisquare,
                                                 Int32 polyOrder=-1);

    bool m_forceDisableLyaFitting=false;
    bool m_forceLyaFitting=false;

    bool SetMultilineNominalAmplitudes(Int32 iLine);
    bool SetMultilineNominalAmplitudesFast(Int32 iCatalog);
    Int32  setLyaProfile(Float64 redshift, 
                          const CRayCatalog::TRayVector& rayList,
                          bool tplratio=false);
    TAsymParams   FitAsymParameters(const Float64& redshift, 
                                    const Int32& idxLyaE,
                                    const TInt32List& filterEltsIdxLya, 
                                    const Int32& idxLineLyaE);
    void SetLSF();

    Float64 GetWeightingAnyLineCenterProximity(Int32 sampleIndex, const TInt32List & EltsIdx) const;
    TInt32List getOverlappingElementsBySupport(Int32 ind , Float64 overlapThres=0.1) const;
    TInt32List ReestimateContinuumApprox(const TInt32List & EltsIdx);
    TInt32List ReestimateContinuumUnderLines(const TInt32List & EltsIdx);
    void refreshModelAfterContReestimation(const TInt32List & EltsIdx, CSpectrumFluxAxis& modelFluxAxis, CSpectrumFluxAxis& spcFluxAxisNoContinuum) const;

    TInt32List findLineIdxInCatalog(const CRayCatalog::TRayVector& restRayList, const std::string &strTag, Int32 type) const;
    TPolynomCoeffs getPolynomCoeffs(Int32 eIdx) const;
    void           applyPolynomCoeffs(Int32 eIdx, const TPolynomCoeffs& polynom_coeffs);
    void addDoubleLine(const CRay &r1, const CRay &r2, Int32 index1, Int32 index2, Float64 nominalWidth, Float64 a1, Float64 a2);

    Int32 improveBalmerFit();
    void applyRules(bool enableLogs=false);
    CRegulament m_Regulament;


    TFloat64List m_ScaleMargCorrTplshape;
    TBoolList m_StrongELPresentTplshape;
    TBoolList m_StrongHalphaELPresentTplshape;
    TInt32List m_NLinesAboveSNRTplshape;

    Float64 m_Redshift;

    CSpectrumFluxAxis m_SpcFluxAxis;    //observed spectrum
    CSpectrumFluxAxis m_spcFluxAxisNoContinuum; //observed spectrum for line fitting
    std::shared_ptr<CTemplate> m_tplContaminantSpcRebin; //optionally used contaminant to be removed from observed spectrum
    CSpectrumNoiseAxis& m_ErrorNoContinuum;
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
    Float64 m_velocityEmission;
    Float64 m_velocityAbsorption;
    Float64 m_velocityEmissionInit;
    Float64 m_velocityAbsorptionInit;

    Float64 m_nominalWidthDefaultEmission;
    Float64 m_nominalWidthDefaultAbsorption;

    std::string m_fittingmethod;

    std::string m_rulesoption;
    std::string m_rigidity;
    bool m_forcedisableTplratioISMfit=false;

    CTemplateCatalog m_tplCatalog;
    TStringList m_tplCategoryList;
    std::string m_tplshapeBestTplName = "None";
    Float64 m_tplshapeBestTplIsmCoeff = NAN;
    Float64 m_tplshapeBestTplAmplitude = NAN;
    Float64 m_tplshapeBestTplDtm = NAN;
    Float64 m_tplshapeBestTplMtm = NAN;
    Int32 m_tplshapeLeastSquareFast = 0;    //for rigidity=tplshape: switch to use fast least square estimation
    std::shared_ptr<const CPriorHelper> m_tplshape_priorhelper;

    COperatorTemplateFitting m_templateFittingOperator;
    Int32 m_secondpass_fitContinuum_dustfit;
    Int32 m_secondpass_fitContinuum_igm;
    Int32 m_secondpass_fitContinuum_outsidelinesmask;
    Int32 m_secondpass_fitContinuum_observedFrame;

    std::shared_ptr<const CTemplatesFitStore> m_fitContinuum_tplfitStore;
    Int32 m_fitContinuum_option;
    std::string m_fitContinuum_tplName;
    Float64 m_fitContinuum_tplFitAmplitude = NAN;
    Float64 m_fitContinuum_tplFitAmplitudeError = NAN;
    Float64 m_fitContinuum_tplFitAmplitudeSigmaMAX = NAN;
    Float64 m_fitContinuum_tplFitMerit = NAN;
    Float64 m_fitContinuum_tplFitMerit_phot = NAN;
    Float64 m_fitContinuum_tplFitEbmvCoeff = NAN;
    Int32   m_fitContinuum_tplFitMeiksinIdx = -1;
    Float64 m_fitContinuum_tplFitRedshift = NAN; // only used with m_fitContinuum_option==2 for now
    Float64 m_fitContinuum_tplFitDtM = NAN;
    Float64 m_fitContinuum_tplFitMtM = NAN;
    Float64 m_fitContinuum_tplFitLogprior = NAN;
    Float64 m_fitContinuum_tplFitSNRMax=NAN;
    TFloat64List m_fitContinuum_tplFitPolyCoeffs;   // only used with m_fitContinuum_option==2 for now
    bool m_forcedisableMultipleContinuumfit=false;
    Float64 m_fitContinuum_tplFitAlpha=0.;
    std::shared_ptr<const CPriorHelper> m_fitContinuum_priorhelper;

    bool m_lmfit_noContinuumTemplate;
    bool m_lmfit_bestTemplate;
    bool m_lmfit_fitContinuum;
    bool m_lmfit_fitEmissionVelocity;
    bool m_lmfit_fitAbsorptionVelocity;

    const TFloat64Range m_lambdaRange;
    
    linetags ltags;
};

}


#endif // _REDSHIFT_LINEMODEL_ELEMENTLIST_
