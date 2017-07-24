#ifndef _REDSHIFT_OPERATOR_LINEMODEL_
#define _REDSHIFT_OPERATOR_LINEMODEL_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/operator/operator.h>
#include <RedshiftLibrary/operator/linemodelresult.h>
#include <RedshiftLibrary/linemodel/elementlist.h>
#include <RedshiftLibrary/linemodel/modelspectrumresult.h>
#include <RedshiftLibrary/linemodel/modelfittingresult.h>
#include <RedshiftLibrary/linemodel/modelrulesresult.h>
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
                                              const Float64 &opt_velocityfitmin=20,
                                              const Float64 &opt_velocityfitmax=800);



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


    Int32 m_maxModelSaveCount;

private:

    Void ModelFit(NSEpic::CLineModelElementList &model, const TFloat64Range& lambdaRange, Float64 redshift,
                  Float64& chiSquare, CLineModelSolution &modelSolution, Int32 contreest_iterations, bool enableLogging);

    void ComputeArea1(CLineModelResult& results);
    void ComputeArea2(CLineModelResult& results);
    Float64 FitBayesWidth( CSpectrumSpectralAxis& spectralAxis, CSpectrumFluxAxis& fluxAxis, Float64 z, Int32 start, Int32 end);

    std::vector<std::shared_ptr<CModelSpectrumResult>  > m_savedModelSpectrumResults;
    std::vector<std::shared_ptr<CModelFittingResult>  > m_savedModelFittingResults;
    std::vector<std::shared_ptr<CModelRulesResult>  > m_savedModelRulesResults;

};


}

#endif
