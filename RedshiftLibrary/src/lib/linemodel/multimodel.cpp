#include <RedshiftLibrary/linemodel/multimodel.h>



using namespace NSEpic;

/**
 * \brief Constructor.
 **/
CMultiModel::CMultiModel(const CSpectrum& spectrum,
                          const CSpectrum &spectrumContinuum,
                          const CTemplateCatalog& tplCatalog,
                          const TStringList& tplCategoryList,
                          const std::string calibrationPath,
                          const CRayCatalog::TRayVector& restRayList,
                          const std::string& opt_fittingmethod,
                          const std::string& opt_continuumcomponent,
                          const std::string& widthType,
                          const Float64 resolution,
                          const Float64 velocityEmission,
                          const Float64 velocityAbsorption,
                          const std::string& opt_rules,
                          const std::string &opt_rigidity)
{

    Int32 nModels = 4;
    for(Int32 km=0; km<nModels; km++)
    {
        m_models.push_back(std::shared_ptr<CLineModelElementList> (new CLineModelElementList(spectrum,
                                                                                             spectrumContinuum,
                                                                                             tplCatalog,
                                                                                             tplCategoryList,
                                                                                             calibrationPath,
                                                                                             restRayList,
                                                                                             opt_fittingmethod,
                                                                                             opt_continuumcomponent,
                                                                                             widthType,
                                                                                             resolution,
                                                                                             velocityEmission,
                                                                                             velocityAbsorption,
                                                                                             opt_rules,
                                                                                             opt_rigidity
                                                                                             )));
    }

}

/**
 * \brief Empty destructor.
 **/
CMultiModel::~CMultiModel()
{
}

Int32 CMultiModel::setPassMode(Int32 iPass)
{
    return 0;
}



Int32 CMultiModel::getTplshape_count()
{
    if(m_models.size()>0)
    {
        return m_models[0]->getTplshape_count();
    }
    else
    {
        return -1;
    }
}

Int32 CMultiModel::getSpcNSamples(const TFloat64Range& lambdaRange)
{
    if(m_models.size()>0)
    {
        return m_models[0]->getSpcNSamples(lambdaRange);
    }
    else
    {
        return -1;
    }
}

Int32 CMultiModel::SetFitContinuum_FitStore(CTemplatesFitStore* fitStore)
{
    Int32 ret=-1;
    for(Int32 km=0; km<m_models.size(); km++)
    {
        ret = m_models[km]->SetFitContinuum_FitStore(fitStore);
    }
    return ret;
}

Float64 CMultiModel::getDTransposeD(const TFloat64Range& lambdaRange, std::string spcComponent)
{
    Float64 valf=0.0;
    for(Int32 km=0; km<m_models.size(); km++)
    {
        valf += m_models[km]->getDTransposeD(lambdaRange, spcComponent);
    }
    return valf;
}

//TODO
Float64 CMultiModel::getLikelihood_cstLog(const TFloat64Range& lambdaRange)
{
    Float64 valf=0.0;
    for(Int32 km=0; km<m_models.size(); km++)
    {
        valf += 0.0;//m_models[km]->getLikelihood_cstLog(lambdaRange);
        //TODO: this is a log value to be combined !
    }
    return valf;
}

void CMultiModel::SetAbsLinesLimit(Float64 limit)
{
    for(Int32 km=0; km<m_models.size(); km++)
    {
        m_models[km]->SetAbsLinesLimit(limit);
    }
}

void CMultiModel::SetLeastSquareFastEstimationEnabled(Int32 enabled)
{
    for(Int32 km=0; km<m_models.size(); km++)
    {
        m_models[km]->SetLeastSquareFastEstimationEnabled(enabled);
    }
}

//TODO: this fitting function needs to be fine tuned for this multimodel use case.
//TBD First-order 0: only use individual fitting method: get the mtm, and dtm values for each model, then estimate the amplitude accordingly for all models
Float64 CMultiModel::fit(Float64 redshift, const TFloat64Range& lambdaRange, CLineModelSolution& modelSolution, Int32 contreest_iterations, bool enableLogging)
{
    Float64 valf=0.0;
    for(Int32 km=0; km<m_models.size(); km++)
    {
        valf += m_models[km]->fit(redshift, lambdaRange, modelSolution, contreest_iterations, enableLogging);
    }
    return valf;
}

//TODO
Float64 CMultiModel::getScaleMargCorrection(Int32 idxLine)
{
    Float64 valf=0.0;
    for(Int32 km=0; km<m_models.size(); km++)
    {
        valf += 0.0;//m_models[km]->getScaleMargCorrection(idxLine);
        //TODO: this is a log value to be combined !
    }
    return valf;
}

std::vector<Float64> CMultiModel::GetChisquareTplshape()
{
    std::vector<Float64> chi2tplshape;
    if(m_models.size()>0)
    {
        chi2tplshape = m_models[0]->GetChisquareTplshape();
    }
    for(Int32 km=1; km<m_models.size(); km++)
    {
        std::vector<Float64> _chi2tplshape = m_models[km]->GetChisquareTplshape();
        for(Int32 ktpl=0; ktpl<_chi2tplshape.size(); ktpl++)
        {
            chi2tplshape[ktpl] += _chi2tplshape[ktpl];
        }
    }

    return chi2tplshape;
}

std::vector<Float64> CMultiModel::GetScaleMargTplshape()
{
    std::vector<Float64> scaleMargtplshape;
    if(m_models.size()>0)
    {
        scaleMargtplshape = m_models[0]->GetScaleMargTplshape();
    }
    for(Int32 km=1; km<m_models.size(); km++)
    {
        std::vector<Float64> _scaleMargtplshape = m_models[km]->GetChisquareTplshape();
        for(Int32 ktpl=0; ktpl<_scaleMargtplshape.size(); ktpl++)
        {
            scaleMargtplshape[ktpl] += _scaleMargtplshape[ktpl];
        }
    }

    return scaleMargtplshape;
}

Float64 CMultiModel::getLeastSquareContinuumMerit(const TFloat64Range& lambdaRange)
{
    Float64 valf=0.0;
    for(Int32 km=0; km<m_models.size(); km++)
    {
        valf += m_models[km]->getLeastSquareContinuumMerit(lambdaRange);
    }
    return valf;
}

Float64 CMultiModel::getLeastSquareContinuumMeritFast()
{
    Float64 valf=0.0;
    for(Int32 km=0; km<m_models.size(); km++)
    {
        valf += m_models[km]->getLeastSquareContinuumMeritFast();
    }
    return valf;
}

Float64 CMultiModel::getContinuumScaleMargCorrection()
{
    Float64 valf=0.0;
    for(Int32 km=0; km<m_models.size(); km++)
    {
        valf += m_models[km]->getContinuumScaleMargCorrection();
    }
    return valf;
}

//todo: check if all values are shared between all the models
Int32 CMultiModel::LoadModelSolution(const CLineModelSolution&  modelSolution)
{
    for(Int32 km=0; km<m_models.size(); km++)
    {
        m_models[km]->LoadModelSolution(modelSolution);
    }
    return 0;
}

void CMultiModel::SetFittingMethod(std::string fitMethod)
{
    for(Int32 km=0; km<m_models.size(); km++)
    {
        m_models[km]->SetFittingMethod(fitMethod);
    }
}

//todo: tbd: which model should be returned in this multimodel case ?
const CSpectrum& CMultiModel::GetModelSpectrum() const
{
    if(m_models.size()>0)
    {
        return m_models[0]->GetModelSpectrum();
    }
    else
    {
        std::shared_ptr<CSpectrum> spcDumb=0;
        return *spcDumb;
    }
}

const CSpectrumFluxAxis& CMultiModel::GetModelContinuum() const
{
    if(m_models.size()>0)
    {
        return m_models[0]->GetModelContinuum();
    }
    else
    {
        std::shared_ptr<CSpectrumFluxAxis> spcDumb=0;
        return *spcDumb;
    }
}

//todo: tbd: which one should be returned in this multimodel case ?
Float64 CMultiModel::GetVelocityEmission()
{
    if(m_models.size()>0)
    {
        return m_models[0]->GetVelocityEmission();
    }
    else
    {
        return -1;
    }
}

//todo: tbd: which one should be returned in this multimodel case ?
Float64 CMultiModel::GetVelocityAbsorption()
{
    if(m_models.size()>0)
    {
        return m_models[0]->GetVelocityAbsorption();
    }
    else
    {
        return -1;
    }
}

//todo: add model number tag in front of the log-string items ?
TStringList CMultiModel::GetModelRulesLog()
{
    TStringList tstr;
    for(Int32 km=0; km<m_models.size(); km++)
    {
        TStringList _tstr = m_models[km]->GetModelRulesLog();
        for(Int32 klist=0; klist<_tstr.size(); klist++)
        {
            tstr.push_back(_tstr[klist]);
        }
    }
    return tstr;

}

void CMultiModel::ResetElementIndexesDisabled()
{
    for(Int32 km=0; km<m_models.size(); km++)
    {
        m_models[km]->ResetElementIndexesDisabled();
    }
}

Int32 CMultiModel::ApplyVelocityBound(Float64 inf, Float64 sup)
{
    for(Int32 km=0; km<m_models.size(); km++)
    {
        m_models[km]->ApplyVelocityBound(inf, sup);
    }
    return 0;
}

void CMultiModel::SetVelocityAbsorption(Float64 vel)
{
    for(Int32 km=0; km<m_models.size(); km++)
    {
        m_models[km]->SetVelocityAbsorption(vel);
    }
}

void CMultiModel::SetVelocityEmission(Float64 vel)
{
    for(Int32 km=0; km<m_models.size(); km++)
    {
        m_models[km]->SetVelocityEmission(vel);
    }
}

Float64 CMultiModel::EstimateMTransposeM(const TFloat64Range& lambdaRange)
{
    Float64 valf=0.0;
    for(Int32 km=0; km<m_models.size(); km++)
    {
        valf += m_models[km]->EstimateMTransposeM(lambdaRange);
    }
    return valf;
}

Float64 CMultiModel::getCumulSNRStrongEL()
{
    Float64 valf=0.0;
    for(Int32 km=0; km<m_models.size(); km++)
    {
        valf += m_models[km]->getCumulSNRStrongEL();
    }
    return valf;
}

Int32 CMultiModel::GetNElements()
{
    if(m_models.size()>0)
    {
        return m_models[0]->GetNElements();
    }
    else
    {
        return -1;
    }
}

Int32 CMultiModel::GetModelNonZeroElementsNDdl()
{
    if(m_models.size()>0)
    {
        return m_models[0]->GetModelNonZeroElementsNDdl();
    }
    else
    {
        return -1;
    }
}

//todo: make the returned mask the concatenation of the mask for all models
CMask CMultiModel::getOutsideLinesMask()
{
    if(m_models.size()>0)
    {
        return m_models[0]->getOutsideLinesMask();
    }
    else
    {
        return -1;
    }
}

//todo: tbd: which one should be returned in this multimodel case ?
std::string CMultiModel::getFitContinuum_tplName()
{
    if(m_models.size()>0)
    {
        return m_models[0]->getFitContinuum_tplName();
    }
    else
    {
        return "";
    }
}

//todo: tbd: which one should be returned in this multimodel case ?
Float64 CMultiModel::getFitContinuum_tplAmplitude()
{
    if(m_models.size()>0)
    {
        return m_models[0]->getFitContinuum_tplAmplitude();
    }
    else
    {
        return -1.;
    }
}

//todo: tbd: which one should be returned in this multimodel case ?
Float64 CMultiModel::getFitContinuum_tplMerit()
{
    if(m_models.size()>0)
    {
        return m_models[0]->getFitContinuum_tplMerit();
    }
    else
    {
        return -1.;
    }
}

//todo: tbd: which one should be returned in this multimodel case ?
Float64 CMultiModel::getFitContinuum_tplIsmDustCoeff()
{
    if(m_models.size()>0)
    {
        return m_models[0]->getFitContinuum_tplIsmDustCoeff();
    }
    else
    {
        return -1.;
    }
}

//todo: tbd: which one should be returned in this multimodel case ?
Float64 CMultiModel::getFitContinuum_tplIgmMeiksinIdx()
{
    if(m_models.size()>0)
    {
        return m_models[0]->getFitContinuum_tplIgmMeiksinIdx();
    }
    else
    {
        return -1.;
    }
}

//todo: tbd: which one should be returned in this multimodel case ?
std::string CMultiModel::getTplshape_bestTplName()
{
    if(m_models.size()>0)
    {
        return m_models[0]->getTplshape_bestTplName();
    }
    else
    {
        return "";
    }
}
