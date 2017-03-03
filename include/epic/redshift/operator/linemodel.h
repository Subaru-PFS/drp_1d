#ifndef _REDSHIFT_OPERATOR_LINEMODEL_
#define _REDSHIFT_OPERATOR_LINEMODEL_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>
#include <epic/redshift/operator/operator.h>
#include <epic/redshift/operator/linemodelresult.h>
#include <epic/redshift/linemodel/elementlist.h>
#include <epic/redshift/linemodel/modelspectrumresult.h>
#include <epic/redshift/linemodel/modelfittingresult.h>
#include <epic/redshift/linemodel/modelrulesresult.h>
#include <epic/redshift/common/mask.h>

#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/ray/catalog.h>

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
                                              const Float64 &opt_velocityfitmin=20,
                                              const Float64 &opt_velocityfitmax=800);

    void storeGlobalModelResults( CDataStore &dataStore );
    void storePerTemplateModelResults( CDataStore &dataStore, const CTemplate& tpl );

private:

    Void ModelFit(NSEpic::CLineModelElementList &model, const TFloat64Range& lambdaRange, Float64 redshift,
                  Float64& chiSquare, CLineModelResult::SLineModelSolution &modelSolution, Int32 contreest_iterations, bool enableLogging);

    void ComputeArea1(CLineModelResult& results);
    void ComputeArea2(CLineModelResult& results);
    Float64 FitBayesWidth( CSpectrumSpectralAxis& spectralAxis, CSpectrumFluxAxis& fluxAxis, Float64 z, Int32 start, Int32 end);

    Float64 PrecomputeLogErr(const CSpectrum& spectrum);

    Float64 mSumLogErr;

    std::vector<std::shared_ptr<CModelSpectrumResult>  > m_savedModelSpectrumResults;
    std::vector<std::shared_ptr<CModelFittingResult>  > m_savedModelFittingResults;
    std::vector<std::shared_ptr<CModelRulesResult>  > m_savedModelRulesResults;

};


}

#endif
