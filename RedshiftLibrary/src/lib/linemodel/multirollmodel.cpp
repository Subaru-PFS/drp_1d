#include <RedshiftLibrary/linemodel/multirollmodel.h>

#include <boost/filesystem.hpp>
#include <RedshiftLibrary/spectrum/io/genericreader.h>
#include <RedshiftLibrary/noise/fromfile.h>
#include <RedshiftLibrary/spectrum/combination.h>

#include <RedshiftLibrary/continuum/median.h>
#include <RedshiftLibrary/continuum/waveletsdf.h>
#include <RedshiftLibrary/continuum/irregularsamplingmedian.h>

namespace bfs = boost::filesystem;
using namespace NSEpic;
using namespace std;

/**
 * \brief Constructor.
 **/
CMultiRollModel::CMultiRollModel(const CSpectrum& spectrum,
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

    Int32 nModels = 3;
    Int32 irollOffset = 0; //0 if using roll0_F, roll1_F etc..., 1 if using roll1_F, roll2_F, etc...
    //also, in this hack, iRollOffset determines which roll is in the spectrumlist.
    std::vector<std::shared_ptr<CSpectrum>> spcRolls;
    for(Int32 km=0; km<nModels; km++)
    {
        std::shared_ptr<CSpectrum> spcRoll = LoadRollSpectrum(spectrum.GetFullPath(), km+irollOffset, irollOffset);
        spcRolls.push_back(spcRoll);
    }
    CSpectrum spcContinuumForMultimodel = CSpectrum(spectrumContinuum);

    Bool enableOverrideContinuumFromCombined=true;
    //auto-deactivate override of continuum if not estimated from spectrum
    if(opt_continuumcomponent!="fromspectrum")
    {
        enableOverrideContinuumFromCombined = false;
    }
    if(enableOverrideContinuumFromCombined)
    {
        CSpectrumCombination spcCombination;
        CSpectrum spcCombined = CSpectrum(*spcRolls[0]);
        Int32 retComb = spcCombination.Combine(spcRolls, spcCombined);
        if( retComb !=0 )
        {
            Log.LogError( "    multirollmodel: Unable to combine rolls for continuum estimate override");
        }

        // Estimate continuum spectrum
        TFloat64List dumb_mask;
        for(Int32 km=0; km<spcCombined.GetSampleCount(); km++)
        {
            dumb_mask.push_back(1);
        }
        std::shared_ptr<CSpectrum> spectrumWithoutContinuum = std::shared_ptr<CSpectrum>( new CSpectrum(spcCombined, dumb_mask) );
        //*spectrumWithoutContinuum = *spcCombined;
        std::shared_ptr<CSpectrum> spectrumForContinuumEstimation = std::shared_ptr<CSpectrum>( new CSpectrum(spcCombined, dumb_mask) );
        //*spectrumForContinuumEstimation = *spcCombined; //NB: could be set to an individual roll instead.

        std::string medianRemovalMethod = spectrumContinuum.GetContinuumEstimationMethod();

        //const char*     nameBaseline;         // baseline filename
        Log.LogInfo( "    multirollmodel: Continuum estimation: using %s", medianRemovalMethod.c_str() );
        if( medianRemovalMethod== "IrregularSamplingMedian")
        {
            // nameBaseline = "preprocess/baselineISMedian";
            CContinuumIrregularSamplingMedian continuum;
            Float64 opt_medianKernelWidth=spectrumContinuum.GetMedianWinsize();
            continuum.SetMedianKernelWidth(opt_medianKernelWidth);
            continuum.SetMeanKernelWidth(opt_medianKernelWidth);
            spectrumWithoutContinuum->RemoveContinuum( continuum );
            spectrumWithoutContinuum->SetMedianWinsize(opt_medianKernelWidth);
        }else if( medianRemovalMethod== "Median")
        {
            // nameBaseline = "preprocess/baselineMedian";
            CContinuumMedian continuum;
            Float64 opt_medianKernelWidth=spectrumContinuum.GetMedianWinsize();
            continuum.SetMedianKernelWidth(opt_medianKernelWidth);
            spectrumWithoutContinuum->RemoveContinuum( continuum );
            spectrumWithoutContinuum->SetMedianWinsize(opt_medianKernelWidth);

        }else if( medianRemovalMethod== "waveletsDF")
        {
            // nameBaseline = "preprocess/baselineDF";
            Int64 nscales=spectrumContinuum.GetDecompScales();
            std::string dfBinPath=spectrumContinuum.GetWaveletsDFBinPath();
            CContinuumDF continuum(dfBinPath);
            spectrumWithoutContinuum->SetDecompScales(nscales);
            bool ret = spectrumWithoutContinuum->RemoveContinuum( continuum );
            if( !ret ) //doesn't seem to work. TODO: check that the df errors lead to a ret=false value
            {
                Log.LogError( "    multirollmodel: Failed to apply continuum substraction for multimodel-combined spectrum" );
                return;
            }
        }else if( medianRemovalMethod== "raw")
        {
            // nameBaseline = "preprocess/baselineRAW";
            CSpectrumFluxAxis& spcFluxAxis = spectrumWithoutContinuum->GetFluxAxis();
            spcFluxAxis.SetSize( spectrumForContinuumEstimation->GetSampleCount() );
            CSpectrumSpectralAxis& spcSpectralAxis = spectrumWithoutContinuum->GetSpectralAxis();
            spcSpectralAxis.SetSize( spectrumForContinuumEstimation->GetSampleCount()  );


            for(Int32 k=0; k<spectrumWithoutContinuum->GetSampleCount(); k++)
            {
                spcFluxAxis[k] = 0.0;
            }
        }else if( medianRemovalMethod== "zero")
        {
            // nameBaseline = "preprocess/baselineZERO";
            CSpectrumFluxAxis& spcFluxAxis = spectrumWithoutContinuum->GetFluxAxis();
            spcFluxAxis.SetSize( spectrumForContinuumEstimation->GetSampleCount() );
            CSpectrumSpectralAxis& spcSpectralAxis = spectrumWithoutContinuum->GetSpectralAxis();
            spcSpectralAxis.SetSize( spectrumForContinuumEstimation->GetSampleCount()  );


            for(Int32 k=0; k<spectrumWithoutContinuum->GetSampleCount(); k++)
            {
                spcFluxAxis[k] = spectrumForContinuumEstimation->GetFluxAxis()[k];
            }
        }
        spectrumWithoutContinuum->SetContinuumEstimationMethod(medianRemovalMethod);

        spcContinuumForMultimodel = *spectrumForContinuumEstimation;
        CSpectrumFluxAxis spcfluxAxis = spcContinuumForMultimodel.GetFluxAxis();
        spcfluxAxis.Subtract( spectrumWithoutContinuum->GetFluxAxis() );
        CSpectrumFluxAxis& sfluxAxisPtr = spcContinuumForMultimodel.GetFluxAxis();
        sfluxAxisPtr = spcfluxAxis;
        spcContinuumForMultimodel.SetMedianWinsize(spectrumWithoutContinuum->GetMedianWinsize());
        spcContinuumForMultimodel.SetDecompScales(spectrumWithoutContinuum->GetDecompScales());
        spcContinuumForMultimodel.SetContinuumEstimationMethod(spectrumWithoutContinuum->GetContinuumEstimationMethod());
        spcContinuumForMultimodel.SetWaveletsDFBinPath(spectrumWithoutContinuum->GetWaveletsDFBinPath());
        Log.LogInfo("    multirollmodel: ===============================================");


    }

    for(Int32 km=0; km<nModels; km++)
    {
        m_models.push_back(std::shared_ptr<CLineModelElementList> (new CLineModelElementList(*spcRolls[km],
                                                                                             spcContinuumForMultimodel,
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



        //hardcoded source size definition for each roll
        if(km==0){
            m_models[km]->SetSourcesizeDispersion(0.35);
        }else if(km==1){
            m_models[km]->SetSourcesizeDispersion(0.35);
        }else if(km==2){
            m_models[km]->SetSourcesizeDispersion(0.35);
        }else if(km==3){
            m_models[km]->SetSourcesizeDispersion(0.35);
        }
    }


}

/**
 * \brief Empty destructor.
 **/
CMultiRollModel::~CMultiRollModel()
{
}

//temporary hack: probably, all the rolls should be read by the client instead
std::shared_ptr<CSpectrum> CMultiRollModel::LoadRollSpectrum(std::string refSpcFullPath, Int32 iRoll, Int32 iRollOffset)
{
    bool verbose=false;
    std::shared_ptr<CSpectrum> spc = std::shared_ptr<CSpectrum>( new CSpectrum() );

    if(verbose)
    {
        Log.LogInfo( "    multirollmodel: load roll #%d: ref spc full path = %s", iRoll, refSpcFullPath.c_str() );
    }
    std::string spcName = bfs::path( refSpcFullPath ).stem().string().c_str();
    std::string spcPath = bfs::path( refSpcFullPath ).parent_path().string().c_str();
    if(verbose)
    {
        Log.LogInfo( "    multirollmodel: load roll #%d: ref spc name = %s", iRoll, spcName.c_str() );
        Log.LogInfo( "    multirollmodel: load roll #%d: ref spc path = %s", iRoll, spcPath.c_str() );
    }
    Int32 substring_start = 0;
    Int32 substring_n;
    std::string strTag = boost::str(boost::format("_roll%d_F") % iRollOffset);
    std::size_t foundstra = spcName.find(strTag.c_str());
    if (foundstra!=std::string::npos){
        substring_n = (Int32)foundstra;
    }else{
        Log.LogWarning( "    multirollmodel: load roll hack - unable to find strTag=%s", strTag.c_str());
        return spc;
    }

    std::string newSpcRollName = boost::str(boost::format("%s_roll%d_F.fits") % spcName.substr(substring_start, substring_n).append("") % iRoll);
    std::string newNoiseRollName = boost::str(boost::format("%s_roll%d_ErrF.fits") % spcName.substr(substring_start, substring_n).append("") % iRoll);

    bfs::path newSpcRollPath = bfs::path(spcPath)/bfs::path(newSpcRollName.c_str());
    bfs::path newNoiseRollPath = bfs::path(spcPath)/bfs::path(newNoiseRollName.c_str());


    Log.LogInfo( "    multirollmodel: load roll #%d: new spc roll path = %s", iRoll, newSpcRollPath.c_str() );
    Log.LogInfo( "    multirollmodel: load roll #%d: new noise roll path = %s", iRoll, newNoiseRollPath.c_str() );
    //

    //Read the fits data
    CSpectrumIOGenericReader reader;

    reader.Read( newSpcRollPath.c_str(), *spc );
    CNoiseFromFile noise;
    noise.SetNoiseFilePath( newNoiseRollPath.c_str(), reader );
    noise.AddNoise( *spc );

    return spc;
}

Int32 CMultiRollModel::LoadFitContaminantTemplate(Int32 iRoll, CTemplate &tpl, const TFloat64Range& lambdaRange)
{
    if(m_models.size()>iRoll)
    {
        return m_models[iRoll]->LoadFitContaminantTemplate(lambdaRange, tpl);
    }
    return 0;
}

std::shared_ptr<CModelSpectrumResult> CMultiRollModel::GetContaminantSpectrumResult(Int32 iRoll)
{
    if(m_models.size()>iRoll)
    {
        return m_models[iRoll]->GetContaminantSpectrumResult();
    }
    return 0;
}

Int32 CMultiRollModel::setPassMode(Int32 iPass)
{
    return 0;
}


Bool CMultiRollModel::initTplratioCatalogs(std::string opt_tplratioCatRelPath, Int32 opt_tplratio_ismFit)
{
    Bool ret=-1;
    for(Int32 km=0; km<m_models.size(); km++)
    {
        ret = m_models[km]->initTplratioCatalogs(opt_tplratioCatRelPath, opt_tplratio_ismFit);
    }

    //
    if(m_models.size()>0 && m_opt_rigidity=="tplshape")
    {
        Int32 nTplshape = m_models[0]->getTplshape_count();
        m_chi2tplshape.resize(nTplshape);
    }

    return ret;
}

Bool CMultiRollModel::initLambdaOffsets(std::string offsetsCatalogsRelPath)
{
    Bool ret=-1;
    for(Int32 km=0; km<m_models.size(); km++)
    {
        m_models[km]->initLambdaOffsets(offsetsCatalogsRelPath);
        ret = true;
    }
    return ret;
}

Int32 CMultiRollModel::getTplshape_count()
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

std::vector<Float64> CMultiRollModel::getTplshape_priors()
{
    if(m_models.size()>0)
    {
        return m_models[0]->getTplshape_priors();
    }
    else
    {
        std::vector<Float64> dumb;
        return dumb;
    }
}

Int32 CMultiRollModel::getSpcNSamples(const TFloat64Range& lambdaRange)
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

Int32 CMultiRollModel::SetFitContinuum_FitStore(CTemplatesFitStore* fitStore)
{
    Int32 ret=-1;
    for(Int32 km=0; km<m_models.size(); km++)
    {
        ret = m_models[km]->SetFitContinuum_FitStore(fitStore);
    }
    return ret;
}

Float64 CMultiRollModel::getDTransposeD(const TFloat64Range& lambdaRange, std::string spcComponent)
{
    Float64 valf=0.0;
    for(Int32 km=0; km<m_models.size(); km++)
    {
        valf += m_models[km]->getDTransposeD(lambdaRange, spcComponent);
    }
    return valf;
}

//TODO
Float64 CMultiRollModel::getLikelihood_cstLog(const TFloat64Range& lambdaRange)
{
    Float64 valf=0.0;
    for(Int32 km=0; km<m_models.size(); km++)
    {
        valf += 0.0;//m_models[km]->getLikelihood_cstLog(lambdaRange);
        //TODO: this is a log value to be combined !
    }
    return valf;
}

void CMultiRollModel::SetAbsLinesLimit(Float64 limit)
{
    for(Int32 km=0; km<m_models.size(); km++)
    {
        m_models[km]->SetAbsLinesLimit(limit);
    }
}

void CMultiRollModel::SetLeastSquareFastEstimationEnabled(Int32 enabled)
{
    for(Int32 km=0; km<m_models.size(); km++)
    {
        m_models[km]->SetLeastSquareFastEstimationEnabled(enabled);
    }
}

//TODO: this fitting function needs to be fine tuned for this multimodel use case.
//TBD First-order 0: only use individual fitting method: get the mtm, and dtm values for each model, then estimate the amplitude accordingly for all models
Float64 CMultiRollModel::fit(Float64 redshift,
                             const TFloat64Range& lambdaRange,
                             CLineModelSolution& modelSolution,
                             Int32 contreest_iterations,
                             bool enableLogging)
{
    //first individual fitting: get the amps, dtm, mtm calculated
    Float64 merit=0.0;
    for(Int32 km=0; km<m_models.size(); km++)
    {
        CLineModelSolution _modelSolution;
        CContinuumModelSolution continuumModelSolution;
        merit += m_models[km]->fit(redshift,
                                   lambdaRange,
                                   _modelSolution,
                                   continuumModelSolution,
                                   contreest_iterations,
                                   enableLogging);
        modelSolution = _modelSolution;
    }

    bool enableOverrideMonolinemodelFit = true;

    if(enableOverrideMonolinemodelFit){
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
                        Float64 weight = 1.;//1/(err*err);
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
            //std::vector<Float64> chi2tplshape(nTplshape, DBL_MAX);
            Float64 minChi2Tplshape = DBL_MAX;
            Int32 iBestTplshape = -1;


            std::vector<std::vector<Float64>> multifit_amps;
            for(Int32 kts=0; kts<nTplshape; kts++)
            {

                /*
                //estimate the error weighted average amps over models
                for(Int32 k=0; k<m_models[0]->m_Elements.size(); k++)
                {
                    amps.push_back(0.0);

                    Float64 weightSum = 0;
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

                    if(enableLogging)
                    {
                        Log.LogInfo( "Multifit: for tplshape=%d, for kElt=%d found A=%f", kts, k, amps[k] );
                    }
                }
                //*/

                //*
                //set amps from cumulated dtm, mtm
                std::vector<Float64> amps(m_models.size(), 0.0);
                std::vector<Float64> dtm_combined;
                std::vector<Float64> mtm_combined;
                for(Int32 k=0; k<m_models[0]->m_Elements.size(); k++)
                {
                    dtm_combined.push_back(0.0);
                    mtm_combined.push_back(0.0);
                    for(Int32 km=0; km<m_models.size(); km++)
                    {
                        dtm_combined[k] += m_models[km]->m_DtmTplshape[kts][k];
                        mtm_combined[k] += m_models[km]->m_MtmTplshape[kts][k];
                    }
                }
                for(Int32 k=0; k<m_models[0]->m_Elements.size(); k++)
                {
                    Float64 dtm = std::max(0.0, dtm_combined[k]);
                    Float64 _amp = 0.0;
                    if(mtm_combined[k]>0.0)
                    {
                        _amp = dtm/mtm_combined[k];
                    }
                    amps[k]=_amp;

                    if(enableLogging)
                    {
                        Log.LogDetail( "    multirollmodel: Multifit: for tplshape=%d, for kElt=%d found A=%f", kts, k, amps[k] );
                    }
                }
                //*/



                multifit_amps.push_back(amps);

                //re-compute the lst-square and store it for current tplshape
                for(Int32 km=0; km<m_models.size(); km++)
                {

                    //set the tplshape
                    m_models[km]->setTplshapeModel(kts, false);
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

                if(enableLogging)
                {
                    Log.LogDetail( "    multirollmodel: Multifit: for tplshape=%d, found lst-sq=%f", kts, valf );
                }
                m_chi2tplshape[kts] = valf;
                if(minChi2Tplshape>valf)
                {
                    minChi2Tplshape=valf;
                    iBestTplshape = kts;
                    modelSolution=m_models[mIndexExportModel]->GetModelSolution();
                }
            }
            merit=minChi2Tplshape;

            //set the model to the min chi2 model, for export
            if(enableLogging)
            {
                for(Int32 km=0; km<m_models.size(); km++)
                {
                    m_models[km]->setTplshapeModel(iBestTplshape, false);
                    for(Int32 k=0; k<m_models[km]->m_Elements.size(); k++)
                    {
                        m_models[km]->m_Elements[k]->SetFittedAmplitude(multifit_amps[iBestTplshape][k], 0.0);
                    }
                    m_models[km]->refreshModel();
                }
            }

        }
    }

    return merit;
}

//TODO
Float64 CMultiRollModel::getScaleMargCorrection(Int32 idxLine)
{
    Float64 valf=0.0;
    for(Int32 km=0; km<m_models.size(); km++)
    {
        valf += 0.0;//m_models[km]->getScaleMargCorrection(idxLine);
        //TODO: this is a log value to be combined !
    }
    return valf;
}

std::vector<Float64> CMultiRollModel::GetChisquareTplshape()
{
    /*
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
    //*/

    return m_chi2tplshape;
}

std::vector<Float64> CMultiRollModel::GetScaleMargTplshape()
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
std::vector<bool> CMultiRollModel::GetStrongELPresentTplshape()
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

Float64 CMultiRollModel::getLeastSquareContinuumMerit(const TFloat64Range& lambdaRange)
{
    Float64 valf=0.0;
    for(Int32 km=0; km<m_models.size(); km++)
    {
        valf += m_models[km]->getLeastSquareContinuumMerit(lambdaRange);
    }
    return valf;
}

Float64 CMultiRollModel::getLeastSquareContinuumMeritFast()
{
    Float64 valf=0.0;
    for(Int32 km=0; km<m_models.size(); km++)
    {
        valf += m_models[km]->getLeastSquareContinuumMeritFast();
    }
    return valf;
}

Float64 CMultiRollModel::getContinuumScaleMargCorrection()
{
    Float64 valf=0.0;
    for(Int32 km=0; km<m_models.size(); km++)
    {
        valf += m_models[km]->getContinuumScaleMargCorrection();
    }
    return valf;
}

//todo: check if all values are shared between all the models
Int32 CMultiRollModel::LoadModelSolution(const CLineModelSolution&  modelSolution)
{
    for(Int32 km=0; km<m_models.size(); km++)
    {
        m_models[km]->LoadModelSolution(modelSolution);
    }
    return 0;
}

void CMultiRollModel::SetFittingMethod(std::string fitMethod)
{
    for(Int32 km=0; km<m_models.size(); km++)
    {
        m_models[km]->SetFittingMethod(fitMethod);
    }
}

const CSpectrum& CMultiRollModel::GetModelSpectrum() const
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

const CSpectrum& CMultiRollModel::GetObservedSpectrumWithLinesRemoved(Int32 lineTypeFilter)
{
    if(m_models.size()>mIndexExportModel)
    {
        return m_models[mIndexExportModel]->GetObservedSpectrumWithLinesRemoved(lineTypeFilter);
    }
    else
    {
        std::shared_ptr<CSpectrum> spcDumb=0;
        return *spcDumb;
    }
}

const CSpectrum CMultiRollModel::GetSpectrumModelContinuum() const
{
    if(m_models.size()>mIndexExportModel)
    {
        return m_models[mIndexExportModel]->GetSpectrumModelContinuum();
    }
    else
    {
        std::shared_ptr<CSpectrum> spcDumb=0;
        return *spcDumb;
    }
}

const CSpectrumFluxAxis& CMultiRollModel::GetModelContinuum() const
{
    if(m_models.size()>mIndexExportModel)
    {
        return m_models[mIndexExportModel]->GetModelContinuum();
    }
    else
    {
        std::shared_ptr<CSpectrumFluxAxis> spcDumb=0;
        return *spcDumb;
    }
}

//todo: tbd: which one should be returned in this multimodel case ?
Float64 CMultiRollModel::GetVelocityEmission()
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
Float64 CMultiRollModel::GetVelocityAbsorption()
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

std::vector<std::vector<Int32>> CMultiRollModel::GetModelVelfitGroups( Int32 lineType )
{
    if(m_models.size()>0)
    {
        return m_models[0]->GetModelVelfitGroups(lineType);
    }
    else
    {
        std::vector<std::vector<Int32>> dumb;
        return dumb;
    }
}

void CMultiRollModel::SetVelocityEmissionOneElement(Float64 vel, Int32 idxElt)
{
    for(Int32 km=0; km<m_models.size(); km++)
    {
        m_models[km]->SetVelocityEmissionOneElement(vel, idxElt);
    }
}

void CMultiRollModel::SetVelocityAbsorptionOneElement(Float64 vel, Int32 idxElt)
{
    for(Int32 km=0; km<m_models.size(); km++)
    {
        m_models[km]->SetVelocityAbsorptionOneElement(vel, idxElt);
    }
}

//todo: add model number tag in front of the log-string items ?
TStringList CMultiRollModel::GetModelRulesLog()
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

void CMultiRollModel::ResetElementIndexesDisabled()
{
    for(Int32 km=0; km<m_models.size(); km++)
    {
        m_models[km]->ResetElementIndexesDisabled();
    }
}

Int32 CMultiRollModel::ApplyVelocityBound(Float64 inf, Float64 sup)
{
    for(Int32 km=0; km<m_models.size(); km++)
    {
        m_models[km]->ApplyVelocityBound(inf, sup);
    }
    return 0;
}

void CMultiRollModel::SetVelocityAbsorption(Float64 vel)
{
    for(Int32 km=0; km<m_models.size(); km++)
    {
        m_models[km]->SetVelocityAbsorption(vel);
    }
}

void CMultiRollModel::SetVelocityEmission(Float64 vel)
{
    for(Int32 km=0; km<m_models.size(); km++)
    {
        m_models[km]->SetVelocityEmission(vel);
    }
}

Float64 CMultiRollModel::EstimateMTransposeM(const TFloat64Range& lambdaRange)
{
    Float64 valf=0.0;
    for(Int32 km=0; km<m_models.size(); km++)
    {
        valf += m_models[km]->EstimateMTransposeM(lambdaRange);
    }
    return valf;
}

Float64 CMultiRollModel::getCumulSNRStrongEL()
{
    Float64 valf=0.0;
    for(Int32 km=0; km<m_models.size(); km++)
    {
        valf += m_models[km]->getCumulSNRStrongEL();
    }
    return valf;
}

Int32 CMultiRollModel::GetNElements()
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

Int32 CMultiRollModel::GetModelNonZeroElementsNDdl()
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
CMask CMultiRollModel::getOutsideLinesMask()
{
    if(m_models.size()>0)
    {
        return m_models[0]->getOutsideLinesMask();
    }
    else
    {
      throw runtime_error("getOutsideLinesMask: Invalid size");
    }
}

//todo: tbd: which one should be returned in this multimodel case ?
std::string CMultiRollModel::getFitContinuum_tplName()
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
Float64 CMultiRollModel::getFitContinuum_tplAmplitude()
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
Float64 CMultiRollModel::getFitContinuum_tplMerit()
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
Float64 CMultiRollModel::getFitContinuum_tplIsmDustCoeff()
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
Float64 CMultiRollModel::getFitContinuum_tplIgmMeiksinIdx()
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
std::string CMultiRollModel::getTplshape_bestTplName()
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
