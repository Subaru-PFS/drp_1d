#ifndef MULTIMODEL_H
#define MULTIMODEL_H

#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/common/datatypes.h>

#include <math.h>
#include <float.h>

#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/ray/catalog.h>
#include <RedshiftLibrary/linemodel/templatesfitstore.h>
#include <RedshiftLibrary/operator/linemodelresult.h>
#include <RedshiftLibrary/linemodel/elementlist.h>
#include <RedshiftLibrary/linemodel/element.h>
#include <RedshiftLibrary/linemodel/multiline.h>



namespace NSEpic
{

class CMultiModel
{

public:

    CMultiModel(const CSpectrum& spectrum,
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

    ~CMultiModel();

    std::shared_ptr<CSpectrum> LoadRollSpectrum(std::string refSpcFullPath, Int32 iRoll);


    Int32 getTplshape_count();
    std::string getTplshape_bestTplName();
    std::vector<Float64> GetChisquareTplshape();
    std::vector<Float64> GetScaleMargTplshape();
    std::vector<bool> GetStrongELPresentTplshape();
    Float64 getLeastSquareContinuumMerit(const TFloat64Range& lambdaRange);
    Float64 getLeastSquareContinuumMeritFast();
    Float64 getContinuumScaleMargCorrection();

    Int32 getSpcNSamples(const TFloat64Range& lambdaRange);
    Float64 getDTransposeD(const TFloat64Range& lambdaRange, std::string spcComponent);
    Float64 getLikelihood_cstLog(const TFloat64Range& lambdaRange);
    Float64 getScaleMargCorrection(Int32 idxLine=-1);

    Float64 EstimateMTransposeM(const TFloat64Range& lambdaRange);
    Float64 getCumulSNRStrongEL();
    Int32 GetNElements();
    Int32 GetModelNonZeroElementsNDdl();
    CMask getOutsideLinesMask();

    std::string getFitContinuum_tplName();
    Float64 getFitContinuum_tplAmplitude();
    Float64 getFitContinuum_tplMerit();
    Float64 getFitContinuum_tplIsmDustCoeff();
    Float64 getFitContinuum_tplIgmMeiksinIdx();

    const CSpectrum& GetModelSpectrum() const;
    const CSpectrum GetSpectrumModelContinuum() const;
    const CSpectrumFluxAxis& GetModelContinuum() const;

    Float64 GetVelocityEmission();
    Float64 GetVelocityAbsorption();
    TStringList GetModelRulesLog();



    Int32 LoadModelSolution(const CLineModelSolution&  modelSolution);

    Int32 setPassMode(Int32 iPass);
    void SetAbsLinesLimit(Float64 limit);
    void SetLeastSquareFastEstimationEnabled(Int32 enabled);

    Int32 SetFitContinuum_FitStore(CTemplatesFitStore* fitStore);
    void SetFittingMethod(std::string fitMethod);
    void ResetElementIndexesDisabled();
    Int32 ApplyVelocityBound(Float64 inf, Float64 sup);
    void SetVelocityAbsorption(Float64 vel);
    void SetVelocityEmission(Float64 vel);



    Float64 fit(Float64 redshift, const TFloat64Range& lambdaRange, CLineModelSolution& modelSolution, Int32 contreest_iterations, bool enableLogging);

    std::vector<std::shared_ptr<CLineModelElementList>  > m_models;

private:

    std::string m_opt_rigidity;
    std::vector<Float64> m_chi2tplshape;
    Int32 mIndexExportModel = 0;

};

}







#endif // MULTIMODEL_H
