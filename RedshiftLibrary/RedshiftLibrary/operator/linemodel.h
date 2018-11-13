#ifndef _REDSHIFT_OPERATOR_LINEMODEL_
#define _REDSHIFT_OPERATOR_LINEMODEL_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/operator/operator.h>
#include <RedshiftLibrary/operator/linemodelresult.h>
#include <RedshiftLibrary/linemodel/elementlist.h>
#include <RedshiftLibrary/linemodel/multirollmodel.h>
#include <RedshiftLibrary/linemodel/modelspectrumresult.h>
#include <RedshiftLibrary/linemodel/modelfittingresult.h>
#include <RedshiftLibrary/linemodel/modelcontinuumfittingresult.h>
#include <RedshiftLibrary/linemodel/modelrulesresult.h>
#include <RedshiftLibrary/operator/spectraFluxResult.h>
#include <RedshiftLibrary/common/mask.h>

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

    std::shared_ptr<COperatorResult> Compute(CDataStore &dataStore,
                                              const CSpectrum& spectrum,
                                              const CSpectrum &spectrumContinuum,
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

    Int32 Init(const CSpectrum& spectrum, const TFloat64List& redshifts);
    std::shared_ptr<COperatorResult> getResult();

    void PrecomputeContinuumFit(const CSpectrum &spectrum,
                                const CSpectrum &spectrumContinuum,
                                const CTemplateCatalog &tplCatalog,
                                const TStringList &tplCategoryList,
                                const std::string opt_calibrationPath,
                                const TFloat64Range &lambdaRange,
                                const Float64 redshiftStep=0.00015,
                                const string zsampling="log",
                                bool ignoreLinesSupport=false);

    Int32 ComputeFirstPass(CDataStore &dataStore,
                                              const CSpectrum& spectrum,
                                              const CSpectrum &spectrumContinuum,
                                              const CTemplateCatalog &tplCatalog,
                                              const TStringList &tplCategoryList,
                                              const std::string opt_calibrationPath,
                                              const CRayCatalog& restraycatalog,
                                              const std::string &opt_lineTypeFilter,
                                              const std::string &opt_lineForceFilter,
                                              const TFloat64Range& lambdaRange,
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
                                              const string &opt_twosteplargegridsampling="log",
                                              const std::string &opt_rigidity="rules",
                                              const string &opt_tplratioCatRelPath="",
                                              const string &opt_offsetCatRelPath="");

    Int32 ComputeCandidates(const Int32 opt_extremacount,
                            const Int32 opt_sign,
                            const std::vector<Float64> floatValues);


    Int32 ComputeSecondPass(CDataStore &dataStore,
                            const CSpectrum& spectrum,
                            const CSpectrum &spectrumContinuum,
                            const CTemplateCatalog &tplCatalog,
                            const TStringList &tplCategoryList,
                            const std::string opt_calibrationPath,
                            const CRayCatalog& restraycatalog,
                            const std::string &opt_lineTypeFilter,
                            const std::string &opt_lineForceFilter,
                            const TFloat64Range& lambdaRange,
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
                            const std::string &opt_rigidity="rules",
                            const Float64 &opt_emvelocityfitmin=20.,
                            const Float64 &opt_emvelocityfitmax=500.,
                            const Float64 &opt_emvelocityfitstep=20.,
                            const Float64 &opt_absvelocityfitmin=150.,
                            const Float64 &opt_absvelocityfitmax=500.,
                            const Float64 &opt_absvelocityfitstep=20.);

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
                                    const std::string &opt_continuumreest);

    std::shared_ptr<COperatorResult> computeWithUltimPass(CDataStore &dataStore,
                                      const CSpectrum& spectrum,
                                      const CSpectrum& spectrumContinuum,
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

    void storeGlobalModelResults( CDataStore &dataStore );
    void storePerTemplateModelResults( CDataStore &dataStore, const CTemplate& tpl );
    std::shared_ptr<CModelSpectrumResult> GetModelSpectrumResult(Int32 idx);
    std::shared_ptr<CSpectraFluxResult> GetModelSpectrumContinuumResult(Int32 idx);


    bool m_enableWidthFitByGroups = false;

    Int32 m_maxModelSaveCount;
    Float64 m_secondPass_extensionradius = 0.005;
    Float64 m_secondPass_velfit_dzInfLim = -4e-4;
    Float64 m_secondPass_velfit_dzSupLim = 4e-4;
    Float64 m_secondPass_velfit_dzStep = 2e-4;

    bool m_enableLoadContTemplate=false;
    Int32 m_iRollContaminated=-1;
    Float64 m_contLambdaOffset=0;
    std::shared_ptr<CTemplate> m_tplContaminant=NULL;
    Int32 initContaminant(std::shared_ptr<CModelSpectrumResult> contModelSpectrum, Int32 iRollContaminated, Float64 contLambdaOffset);
    std::shared_ptr<CModelSpectrumResult> GetContaminantSpectrumResult();
    std::shared_ptr<CModelSpectrumResult> m_savedContaminantSpectrumResult;

    Int32 m_opt_tplfit_dustFit = 1;
    Int32 m_opt_tplfit_extinction = 1;
    Int32 m_opt_fitcontinuum_maxN = 2;
    bool m_opt_tplfit_ignoreLinesSupport=false; //default: false, as ortho templates store makes this un-necessary

    Int32 m_opt_tplratio_ismFit = 1;
    Int32 m_opt_firstpass_tplratio_ismFit=0;
    Int32 m_opt_firstpass_multiplecontinuumfit_disable=1;
    std::string m_opt_firstpass_fittingmethod;
    Int32 m_opt_secondpass_tplfit_fixfromfirstpass=1; //0: load fit continuum, 1 (default): use the best continuum from first pass
private:

    std::shared_ptr<CLineModelResult> m_result;
    std::shared_ptr<CLineModelElementList> m_model;
    TFloat64List m_sortedRedshifts;

    //candidates
    TPointList m_firstpass_extremumList;
    CLineModelExtremaResult m_firstpass_extremaResult;
    CLineModelExtremaResult m_secondpass_parameters_extremaResult;
    TPointList m_secondpass_recomputed_extremumList;
    std::vector<Int32> m_secondpass_indiceSortedCandidatesList;

    Int32 m_enableFastFitLargeGrid = 0;
    Int32 m_estimateLeastSquareFast = 0;


    void ComputeArea1(CLineModelResult& results);
    void ComputeArea2(CLineModelResult& results);
    Float64 FitBayesWidth( CSpectrumSpectralAxis& spectralAxis, CSpectrumFluxAxis& fluxAxis, Float64 z, Int32 start, Int32 end);

    Int32 interpolateLargeGridOnFineGrid(TFloat64List redshiftsLargeGrid, TFloat64List redshiftsFineGrid, TFloat64List meritLargeGrid, TFloat64List &meritFineGrid);

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

};


}

#endif
