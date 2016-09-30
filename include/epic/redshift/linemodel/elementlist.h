#ifndef ELEMENTLIST_H
#define ELEMENTLIST_H

#include <epic/core/common/range.h>
#include <epic/redshift/common/datatypes.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <epic/redshift/ray/catalog.h>
#include <epic/redshift/spectrum/spectrum.h>

#include <epic/redshift/operator/linemodelresult.h>
#include <epic/redshift/linemodel/element.h>
#include <epic/redshift/linemodel/singleline.h>

#include <epic/redshift/ray/catalogsTplShape.h>

#include <boost/shared_ptr.hpp>

#include <memory>


namespace NSEpic
{
  static Int32 defaultIdx = -1;
  class CRegulament;
  
class CLineModelElementList
{

public:

    CLineModelElementList(const CSpectrum& spectrum, const CSpectrum& spectrumNoContinuum, const CRayCatalog::TRayVector& restRayList, const std::string& opt_fittingmethod, const std::string &opt_continuumcomponent, const std::string& lineWidthType, const Float64 resolution, const Float64 velocityEmission, const Float64 velocityAbsorption, const std::string &opt_rules);
    ~CLineModelElementList();

    void LoadCatalog(const CRayCatalog::TRayVector& restRayList);
    void LoadCatalogSingleLines(const CRayCatalog::TRayVector& restRayList);
    void LogCatalogInfos();

    void LoadContinuum();
    void PrepareContinuum(Float64 z);

    void EstimateSpectrumContinuum();

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
    Float64 GetVelocitySup();
    Int32 ApplyVelocityBound();

    Float64 fit(Float64 redshift, const TFloat64Range& lambdaRange, CLineModelResult::SLineModelSolution &modelSolution, Int32 contreest_iterations=0, bool enableLogging=0);
    void fitWithModelSelection(Float64 redshift, const TFloat64Range& lambdaRange, CLineModelResult::SLineModelSolution &modelSolution);
    void SetFittingMethod(std::string fitMethod);

    Float64 GetRedshift();

    void reinitModel();
    void refreshModel();
    void reinitModelUnderElements(std::vector<Int32> filterEltsIdx, Int32 lineIdx=-1 );
    void refreshModelUnderElements(std::vector<Int32> filterEltsIdx, Int32 lineIdx=-1 );
    void refreshModelDerivSigmaUnderElements(std::vector<Int32> filterEltsIdx);

    void setModelSpcObservedOnSupportZeroOutside(const TFloat64Range &lambdaRange);

    Int32 getSpcNSamples(const TFloat64Range& lambdaRange);
    Float64 getLeastSquareMerit(const TFloat64Range &lambdaRange);
    Float64 getLeastSquareMeritUnderElements();
    Float64 getStrongerMultipleELAmpCoeff();
    Float64 getCumulSNRStrongEL();
    Float64 getCumulSNROnRange( TInt32Range idxRange );
    Float64 getModelErrorUnderElement(Int32 eltId);
    Float64 getContinuumMeanUnderElement(Int32 eltId);
    CLineModelResult::SLineModelSolution GetModelSolution();
    const CSpectrum&    GetModelSpectrum() const;
    const CSpectrum&    GetObservedSpectrumWithLinesRemoved() const;
    const CSpectrumFluxAxis&    GetModelContinuum() const;
    Float64 getModelFluxVal(Int32 idx) const;
    Float64 getModelFluxDerivEltVal(Int32 DerivEltIdx, Int32 idx) const;
    Float64 getModelFluxDerivSigmaVal(Int32 idx) const;

    std::vector<boost::shared_ptr<CLineModelElement>  > m_Elements;
    std::shared_ptr<CSpectrum>  m_SpectrumModel;  //model
    CSpectrumFluxAxis m_SpcContinuumFluxAxis; //continuum spectrum used for the model
    Int32 FindElementIndex(Int32 LineCatalogIndex);
    Int32 FindElementIndex(std::string LineTagStr, Int32 linetype=-1, Int32& lineIdx=defaultIdx);
    std::vector<Int32> getOverlappingElements(Int32 ind , std::vector<Int32> excludedInd, Float64 overlapThres=0.1);
    CRayCatalog::TRayVector m_RestRayList;
    std::shared_ptr<CSpectrum>  m_SpcCorrectedUnderLines;  //observed spectrum corrected under the lines

    TStringList GetModelRulesLog();

private:

    Int32 fitAmplitudesHybrid(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& spcFluxAxisNoContinuum, Float64 redshift);
    void fitAmplitudesSimplex();
    Int32 fitAmplitudesLmfit(std::vector<Int32> EltsIdx, const CSpectrumFluxAxis &fluxAxis, std::vector<Float64> &ampsfitted, Int32 lineType);
    Int32 fitAmplitudesLinSolve(std::vector<Int32> EltsIdx, const CSpectrumSpectralAxis &spectralAxis, const CSpectrumFluxAxis &fluxAxis, std::vector<Float64> &ampsfitted, std::vector<Float64> &errorsfitted);

    Int32 setLyaProfile( Float64 redshift, const CSpectrumSpectralAxis& spectralAxis );

    std::vector<Int32> getSupportIndexes(std::vector<Int32> EltsIdx);
    std::vector<Int32> getOverlappingElementsBySupport(Int32 ind , Float64 overlapThres=0.1);
    std::vector<Int32> ReestimateContinuumApprox(std::vector<Int32> EltsIdx);
    std::vector<Int32> ReestimateContinuumUnderLines(std::vector<Int32> EltsIdx);
    void refreshModelAfterContReestimation(std::vector<Int32> EltsIdx, CSpectrumFluxAxis& modelFluxAxis, CSpectrumFluxAxis& spcFluxAxisNoContinuum);

    std::vector<Int32> findLineIdxInCatalog(const CRayCatalog::TRayVector& restRayList, std::string strTag, Int32 type);

    void addSingleLine(const CRay &r, Int32 index, Float64 nominalWidth);
    void addDoubleLine(const CRay &r1, const CRay &r2, Int32 index1, Int32 index2, Float64 nominalWidth, Float64 a1, Float64 a2);

    void applyRules(bool enableLogs=false);
    CRegulament* m_Regulament;

    CRayCatalogsTplShape* m_CatalogTplShape;

    Float64 m_Redshift;

    CSpectrumFluxAxis m_SpcFluxAxis;    //observed spectrum
    CSpectrumFluxAxis m_spcFluxAxisNoContinuum; //observed spectrum for line fitting
    Float64* m_ErrorNoContinuum;
    CSpectrumFluxAxis m_SpcFluxAxisModelDerivSigma;

    Float64*          m_precomputedFineGridContinuumFlux;   //PFG buffer for model continuum
    CSpectrumFluxAxis m_ContinuumFluxAxis;  //rebined model continuum
    std::string m_ContinuumComponent;
    std::string m_LineWidthType;
    Float64 m_resolution;
    Float64 m_velocityEmission;
    Float64 m_velocityAbsorption;
    Float64 m_velocityEmissionInit;
    Float64 m_velocityAbsorptionInit;
    Float64 m_nominalWidthDefaultEmission;
    Float64 m_nominalWidthDefaultAbsorption;
    std::string m_fittingmethod;
    std::vector<Int32> m_elementsDisabledIndexes;
    std::string m_rulesoption;
};

}







#endif // ELEMENTLIST_H
