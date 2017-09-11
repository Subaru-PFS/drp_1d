#include <RedshiftLibrary/linemodel/multimodel.h>

#include <boost/filesystem.hpp>
#include <RedshiftLibrary/spectrum/io/genericreader.h>
#include <RedshiftLibrary/noise/fromfile.h>

namespace bfs = boost::filesystem;
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

    m_opt_rigidity = opt_rigidity;

    Int32 nModels = 2;
    for(Int32 km=0; km<nModels; km++)
    {
        std::shared_ptr<CSpectrum> spcRoll = LoadRollSpectrum(spectrum.GetFullPath(), km+1);
        m_models.push_back(std::shared_ptr<CLineModelElementList> (new CLineModelElementList(*spcRoll,
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
        if(km==0){
            m_models[km]->SetSourcesizeDispersion(0.1);
        }else if(km==1){
            m_models[km]->SetSourcesizeDispersion(0.5);
        }
    }

}

/**
 * \brief Empty destructor.
 **/
CMultiModel::~CMultiModel()
{
}

//temporary hack: probably, all the rolls should be read by the client instead
std::shared_ptr<CSpectrum> CMultiModel::LoadRollSpectrum(std::string refSpcFullPath, Int32 iRoll)
{
    bool verbose=false;
    std::shared_ptr<CSpectrum> spc = std::shared_ptr<CSpectrum>( new CSpectrum() );

    if(verbose)
    {
        Log.LogInfo( "Multimodel: load roll #%d: ref spc full path = %s", iRoll, refSpcFullPath.c_str() );
    }
    std::string spcName = bfs::path( refSpcFullPath ).stem().string().c_str();
    std::string spcPath = bfs::path( refSpcFullPath ).parent_path().string().c_str();
    if(verbose)
    {
        Log.LogInfo( "Multimodel: load roll #%d: ref spc name = %s", iRoll, spcName.c_str() );
        Log.LogInfo( "Multimodel: load roll #%d: ref spc path = %s", iRoll, spcPath.c_str() );
    }
    Int32 substring_start = 0;
    Int32 substring_n;
    std::string strTag = "_roll1_F";
    std::size_t foundstra = spcName.find(strTag.c_str());
    if (foundstra!=std::string::npos){
        substring_n = (Int32)foundstra;
    }else{
        Log.LogWarning( "Multimodel - load roll hack - unable to find strTag=%s", strTag.c_str());
        return spc;
    }

    std::string newSpcRollName = boost::str(boost::format("%s_roll%d_F.fits") % spcName.substr(substring_start, substring_n).append("") % iRoll);
    std::string newNoiseRollName = boost::str(boost::format("%s_roll%d_ErrF.fits") % spcName.substr(substring_start, substring_n).append("") % iRoll);

    bfs::path newSpcRollPath = bfs::path(spcPath)/bfs::path(newSpcRollName.c_str());
    bfs::path newNoiseRollPath = bfs::path(spcPath)/bfs::path(newNoiseRollName.c_str());


    Log.LogInfo( "Multimodel: load roll #%d: new spc roll path = %s", iRoll, newSpcRollPath.c_str() );
    Log.LogInfo( "Multimodel: load roll #%d: new noise roll path = %s", iRoll, newNoiseRollPath.c_str() );
    //

    //Read the fits data
    CSpectrumIOGenericReader reader;
    Bool rValue = reader.Read( newSpcRollPath.c_str(), *spc );
    if( !rValue )
    {
        Log.LogError("Failed to read input spectrum file: (%s)", newSpcRollPath.c_str() );
        spc = NULL;
        return spc;
    }

    //add noise
    {
        CNoiseFromFile noise;
        if( ! noise.SetNoiseFilePath( newNoiseRollPath.c_str() ) )
        {
            Log.LogError("Failed to load noise spectrum");
            return spc;
        }

        if( ! noise.AddNoise( *spc ) )
        {
            Log.LogError( "Failed to apply noise from spectrum: %s", newNoiseRollPath.c_str() );
            return spc;
        }
    }

    return spc;
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
Float64 CMultiModel::
fit(Float64 redshift, const TFloat64Range& lambdaRange, CLineModelSolution& modelSolution, Int32 contreest_iterations, bool enableLogging)
{
    //first individual fitting: get the amps, dtm, mtm calculated
    Float64 merit=0.0;
    for(Int32 km=0; km<m_models.size(); km++)
    {
        CLineModelSolution _modelSolution;
        m_models[km]->fit(redshift, lambdaRange, _modelSolution, contreest_iterations, enableLogging);
    }

    if(m_opt_rigidity != "tplshape")
    {
        /*
    //set amps from ref model
    Int32 irefModel = 0;
    std::vector<Float64> amps;
    for(Int32 k=0; k<m_models[irefModel]->m_Elements.size(); k++)
    {
        Float64 _amp = m_models[irefModel]->m_Elements[k]->GetElementAmplitude();
        amps.push_back(_amp);
    }
    //*/

        //*
        //estimate error weighted average amps over models
        std::vector<Float64> amps;
        for(Int32 k=0; k<m_models[0]->m_Elements.size(); k++)
        {
            Float64 weightSum = 0;
            amps.push_back(0.0);
            for(Int32 km=0; km<m_models.size(); km++)
            {
                Float64 err = m_models[km]->m_Elements[k]->GetElementError();
                if(err>0.0)
                {
                    Float64 weight = 1/(err*err);
                    amps[k] += m_models[km]->m_Elements[k]->GetElementAmplitude()*weight;
                    weightSum += weight;
                }
            }
            if(weightSum>0.0)
            {
                amps[k]/=weightSum;
            }
        }
        //*/

        /*
        //set amps from combined chi2 calculation: work in progress...
        std::vector<Float64> dtm_combined;
        std::vector<Float64> mtm_combined;
        for(Int32 k=0; k<m_models[0]->m_Elements.size(); k++)
        {
            dtm_combined.push_back(0.0);
            mtm_combined.push_back(0.0);
            for(Int32 km=0; km<m_models.size(); km++)
            {
                dtm_combined[k] += m_models[km]->m_Elements[k]->GetDtmFree();
                mtm_combined[k] += m_models[km]->m_Elements[k]->GetSumGauss();
            }
        }
        std::vector<Float64> amps;
        for(Int32 k=0; k<m_models[0]->m_Elements.size(); k++)
        {
            Float64 dtm = std::max(0.0, dtm_combined[k]);
            Float64 _amp = 0.0;
            if(mtm_combined[k]>0.0)
            {
                dtm/mtm_combined[k];
            }
            amps.push_back(_amp);
        }
        //*/

        //*
        for(Int32 km=0; km<m_models.size(); km++)
        {
            for(Int32 k=0; k<m_models[km]->m_Elements.size(); k++)
            {
                m_models[km]->m_Elements[k]->SetFittedAmplitude(amps[k], 0.0);
            }
            m_models[km]->refreshModel();
        }
        //Get updated merit
        Float64 valf=0.0;
        for(Int32 km=0; km<m_models.size(); km++)
        {
            valf += m_models[km]->getLeastSquareMerit(lambdaRange);
        }
        merit=valf;
        //*/
    }else{
        Int32 nTplshape = m_models[0]->getTplshape_count();
        std::vector<Float64> chi2tplshape(nTplshape, DBL_MAX);
        Float64 minChi2Tplshape = DBL_MAX;

        for(Int32 kts=0; kts<nTplshape; kts++)
        {

            //estimate the error weighted average amps over models
            std::vector<Float64> amps;
            for(Int32 k=0; k<m_models[0]->m_Elements.size(); k++)
            {                
//                //set the tplshape
//                m_models[km].initTplshapeModel(kts, false);
//                //set the tplshape amplitude
//                m_models[km].setTplshapeAmplitude( ampsElts, errorsElts);

                Float64 weightSum = 0;
                amps.push_back(0.0);
                for(Int32 km=0; km<m_models.size(); km++)
                {
                    Float64 err = m_models[km]->m_FittedErrorTplshape[kts][k];
                    if(err>0.0)
                    {
                        Float64 weight = 1/(err*err);
                        amps[k] += m_models[km]->m_FittedAmpTplshape[kts][k]*weight;
                        weightSum += weight;
                    }
                }
                if(weightSum>0.0)
                {
                    amps[k]/=weightSum;
                }
            }


            //re-compute the lst-square and store it for current tplshape
            for(Int32 km=0; km<m_models.size(); km++)
            {

                //set the tplshape
                m_models[km]->initTplshapeModel(kts, false);
                //set the tplshape amplitude
                //m_models[km]->setTplshapeAmplitude( ampsElts, errorsElts);

                for(Int32 k=0; k<m_models[km]->m_Elements.size(); k++)
                {
                    m_models[km]->m_Elements[k]->SetFittedAmplitude(amps[k], 0.0);
                }
                m_models[km]->refreshModel();
            }
            //Get updated merit
            Float64 valf=0.0;
            for(Int32 km=0; km<m_models.size(); km++)
            {
                valf += m_models[km]->getLeastSquareMerit(lambdaRange);
            }
            chi2tplshape[kts] = valf;
            if(minChi2Tplshape>valf)
            {
                minChi2Tplshape=valf;
                modelSolution=m_models[mIndexExportModel]->GetModelSolution();
            }
        }
        merit=minChi2Tplshape;
    }


    return merit;
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

//todo: tbd, are these booleans to be combined by OR or AND ?
std::vector<bool> CMultiModel::GetStrongELPresentTplshape()
{
    std::vector<bool> strongElPresentplshape;
    if(m_models.size()>0)
    {
        strongElPresentplshape = m_models[0]->GetStrongELPresentTplshape();
    }
    for(Int32 km=1; km<m_models.size(); km++)
    {
        std::vector<bool> _strongElPresentplshape = m_models[km]->GetStrongELPresentTplshape();
        for(Int32 ktpl=0; ktpl<_strongElPresentplshape.size(); ktpl++)
        {
            strongElPresentplshape[ktpl] = strongElPresentplshape[ktpl] || _strongElPresentplshape[ktpl];
        }
    }

    return strongElPresentplshape;
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
    if(m_models.size()>mIndexExportModel)
    {
        return m_models[mIndexExportModel]->GetModelSpectrum();
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
