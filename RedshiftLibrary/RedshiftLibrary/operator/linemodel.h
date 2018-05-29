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
                                              const std::string &opt_rigidity="rules",
                                              const Float64 &opt_emvelocityfitmin=20,
                                              const Float64 &opt_emvelocityfitmax=500,
                                             const Float64 &opt_absvelocityfitmin=150,
                                             const Float64 &opt_absvelocityfitmax=500);

    Int32 Init(const CSpectrum& spectrum, const TFloat64List& redshifts);
    std::shared_ptr<COperatorResult> getResult();

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
                                              const std::string &opt_rigidity="rules");

    Int32 ComputeCandidates(const Int32 opt_extremacount, const Int32 opt_sign, const std::vector<Float64> floatValues);


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
                                              const Float64 &opt_emvelocityfitmin=20,
                                              const Float64 &opt_emvelocityfitmax=500,
                                             const Float64 &opt_absvelocityfitmin=150,
                                             const Float64 &opt_absvelocityfitmax=500);

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
                                      const std::string& opt_rigidity,
                                      const Float64 &opt_velocityfitmin,
                                      const Float64 &opt_velocityfitmax);

    void storeGlobalModelResults( CDataStore &dataStore );
    void storePerTemplateModelResults( CDataStore &dataStore, const CTemplate& tpl );
    std::shared_ptr<CModelSpectrumResult> GetModelSpectrumResult(Int32 idx);
    std::shared_ptr<CSpectraFluxResult> GetModelSpectrumContinuumResult(Int32 idx);


    bool m_enableWidthFitByGroups;

    Int32 m_maxModelSaveCount;
    Float64 m_secondPass_extensionradius = 0.005;
    Float64 m_secondPass_velfit_dzInfLim = -4e-4;
    Float64 m_secondPass_velfit_dzSupLim = 4e-4;
    Float64 m_secondPass_velfit_dzStep = 2e-4;
    Float64 m_secondPass_velfit_vStep = 20.0;

    bool m_enableLoadContTemplate=false;
    Int32 m_iRollContaminated=-1;
    Float64 m_contLambdaOffset=0;
    std::shared_ptr<CTemplate> m_tplContaminant=NULL;
    Int32 initContaminant(std::shared_ptr<CModelSpectrumResult> contModelSpectrum, Int32 iRollContaminated, Float64 contLambdaOffset);
    std::shared_ptr<CModelSpectrumResult> GetContaminantSpectrumResult();
    std::shared_ptr<CModelSpectrumResult> m_savedContaminantSpectrumResult;

private:

    std::shared_ptr<CLineModelResult> m_result;
    std::shared_ptr<CLineModelElementList> m_model;
    TFloat64List m_sortedRedshifts;
    TPointList m_extremumList;

    Int32 m_enableFastFitLargeGrid;
    Int32 m_estimateLeastSquareFast;

    void ComputeArea1(CLineModelResult& results);
    void ComputeArea2(CLineModelResult& results);
    Float64 FitBayesWidth( CSpectrumSpectralAxis& spectralAxis, CSpectrumFluxAxis& fluxAxis, Float64 z, Int32 start, Int32 end);

    Int32 interpolateLargeGridOnFineGrid(TFloat64List redshiftsLargeGrid, TFloat64List redshiftsFineGrid, TFloat64List meritLargeGrid, TFloat64List &meritFineGrid);

    std::vector<std::shared_ptr<CModelSpectrumResult>  > m_savedModelSpectrumResults;
    std::vector<std::shared_ptr<CModelFittingResult>  > m_savedModelFittingResults;
    std::vector<std::shared_ptr<CModelRulesResult>  > m_savedModelRulesResults;
    std::vector<std::shared_ptr<CSpectraFluxResult>  > m_savedModelContinuumSpectrumResults;

};


}

#endif
