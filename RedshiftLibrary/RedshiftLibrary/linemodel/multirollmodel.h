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
#ifndef _REDSHIFT_LINEMODEL_MULTIROLLMODEL_
#define _REDSHIFT_LINEMODEL_MULTIROLLMODEL_

#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/common/datatypes.h"

#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"
#include "RedshiftLibrary/ray/catalog.h"
#include "RedshiftLibrary/linemodel/templatesfitstore.h"
#include "RedshiftLibrary/operator/linemodelresult.h"
#include "RedshiftLibrary/linemodel/linemodelfitting.h"
#include "RedshiftLibrary/linemodel/element.h"

#include <cmath>
#include <cfloat>


namespace NSEpic
{

class CMultiRollModel
{

public:

    CMultiRollModel(const CSpectrum& spectrum,
                    const TFloat64Range& lambdaRange,
                    const CTemplateCatalog& tplCatalog,
                    const TStringList& tplCategoryList,
                    const std::string calibrationPath,
                    const CRayCatalog::TRayVector& restRayList,
                    const CRayCatalogsTplShape& tplRatioCatalog,
                    const std::string& opt_fittingmethod,
                    const std::string &opt_continuumcomponent,
                    const Float64 opt_continuum_neg_threshold,
                    const std::string& lineWidthType,
                    const Float64 velocityEmission,
                    const Float64 velocityAbsorption,
                    const std::string &opt_rules,
                    const std::string &opt_rigidity);

    ~CMultiRollModel();

    std::shared_ptr<CSpectrum> LoadRollSpectrum(std::string refSpcFullPath, Int32 iRoll, Int32 iRollOffset);
    Int32 LoadFitContaminantTemplate(Int32 iRoll, CTemplate& tpl, const TFloat64Range& lambdaRange);
    std::shared_ptr<CModelSpectrumResult> GetContaminantSpectrumResult(Int32 iRoll);

    Int32 getTplshape_count();
    std::vector<Float64> getTplshape_priors();
    std::string getTplshape_bestTplName();
    std::vector<Float64> GetChisquareTplshape();
    std::vector<Float64> GetScaleMargTplshape();
    TBoolList GetStrongELPresentTplshape();
    Float64 getLeastSquareContinuumMerit(const TFloat64Range& lambdaRange);
    Float64 getLeastSquareContinuumMeritFast();
    Float64 getContinuumScaleMargCorrection();
    bool initTplratioCatalogs(const CRayCatalogsTplShape& tplRatioCatalog, Int32 opt_tplratio_ismFit);

    Int32 getSpcNSamples(const TFloat64Range& lambdaRange);
    Float64 getDTransposeD(const TFloat64Range& lambdaRange);
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
    Float64 getFitContinuum_tplIsmEbmvCoeff();
    Float64 getFitContinuum_tplIgmMeiksinIdx();

    const CSpectrum& GetModelSpectrum() const;
    const CSpectrum GetSpectrumModelContinuum() const;
    const CSpectrumFluxAxis& GetModelContinuum() const;
    const CSpectrum& GetObservedSpectrumWithLinesRemoved(Int32 lineTypeFilter=-1);

    Float64 GetVelocityEmission();
    Float64 GetVelocityAbsorption();
    TStringList GetModelRulesLog();
    std::vector<std::vector<Int32>> GetModelVelfitGroups(Int32 lineType );
    void SetVelocityEmissionOneElement(Float64 vel, Int32 idxElt);
    void SetVelocityAbsorptionOneElement(Float64 vel, Int32 idxElt);


    Int32 LoadModelSolution(const CLineModelSolution& modelSolution);

    Int32 setPassMode(Int32 iPass);
    void SetAbsLinesLimit(Float64 limit);
    void SetLeastSquareFastEstimationEnabled(Int32 enabled);

    Int32 SetFitContinuum_FitStore(std::shared_ptr<const CTemplatesFitStore> & fitStore);
    void SetFittingMethod(std::string fitMethod);
    void ResetElementIndexesDisabled();
    Int32 ApplyVelocityBound(Float64 inf, Float64 sup);
    void SetVelocityAbsorption(Float64 vel);
    void SetVelocityEmission(Float64 vel);

    Float64 fit(Float64 redshift, const TFloat64Range& lambdaRange, CLineModelSolution& modelSolution, Int32 contreest_iterations, bool enableLogging);

    std::vector<std::shared_ptr<CLineModelFitting> > m_models;

private:

    std::string m_opt_rigidity;
    std::vector<Float64> m_chi2tplshape;
    Int32 mIndexExportModel = 0;

};

}


#endif // _REDSHIFT_LINEMODEL_MULTIROLLMODEL_
