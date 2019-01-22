#include <RedshiftLibrary/linemodel/elementlist.h>
#include <RedshiftLibrary/linemodel/multiline.h>
#include <RedshiftLibrary/linemodel/modelfittingresult.h>
#include <RedshiftLibrary/ray/regulament.h>
#include <RedshiftLibrary/ray/catalogsTplShape.h>
#include <RedshiftLibrary/ray/catalogsOffsets.h>
#include <RedshiftLibrary/ray/linetags.h>
#include <RedshiftLibrary/ray/ray.h>

#include <gsl/gsl_multifit.h>
#include <RedshiftLibrary/spectrum/io/genericreader.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <RedshiftLibrary/continuum/waveletsdf.h>
#include <RedshiftLibrary/continuum/irregularsamplingmedian.h>
#include <RedshiftLibrary/spectrum/io/fitswriter.h>
#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/log/log.h>

#include <math.h>
#include <boost/numeric/conversion/bounds.hpp>
#include <boost/format.hpp>
#include <boost/chrono/thread_clock.hpp>
#include <boost/algorithm/string.hpp>
#include <algorithm>

#include <stdlib.h>
#include <stdio.h>
//#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include "lmfitfunctions.c"

#include <numeric>




using namespace NSEpic;
using namespace std;

/**
 * \brief Prepares the state for Linemodel operation.
 * Loads the catalog.
 * Sets many state variables.
 * Sets the continuum either as a nocontinuum or a fromspectrum.
 **/
CLineModelElementList::CLineModelElementList(const CSpectrum& spectrum,
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
                          const std::string &opt_rigidity):
  m_ErrorNoContinuum(m_spcFluxAxisNoContinuum.GetError())
{
    //members for chi2 continuum fitting
    m_tplCatalog = tplCatalog;
    m_tplCategoryList = tplCategoryList;
    m_inputSpc = std::shared_ptr<CSpectrum>( new CSpectrum(spectrum) );

    m_ContinuumComponent = opt_continuumcomponent;
    m_LineWidthType = widthType;
    m_resolution = resolution;
    m_velocityEmission = velocityEmission;
    m_velocityAbsorption = velocityAbsorption;
    m_velocityEmissionInit = m_velocityEmission;
    m_velocityAbsorptionInit = m_velocityAbsorption;
    m_fittingmethod = opt_fittingmethod;
    m_rulesoption = opt_rules;
    m_rigidity = opt_rigidity;
    m_calibrationPath = calibrationPath;

    //tplfit continuum option: warning, these options not used when using the precomputed continuum fit store (which is recommended)
    m_fitContinuum_dustfit = -10;
    m_fitContinuum_igm = 1;
    m_fitContinuum_outsidelinesmask = 0;
    m_fitContinuum_observedFrame = 0;

    //m_nominalWidthDefaultEmission = 1.15;// suited to new pfs simulations
    m_nominalWidthDefaultEmission = 13.4;// euclid 1 px
    m_nominalWidthDefaultAbsorption = m_nominalWidthDefaultEmission;

    m_enableAmplitudeOffsets = false; //this is highly experimental for now.
    m_enableLambdaOffsetsFit = true; //enable lambdaOffsetFit. Once enabled, the offset fixed value or the fitting on/off switch is done through the offset calibration file.

    m_SpectrumModel = std::shared_ptr<CSpectrum>( new CSpectrum(spectrum) );
    m_SpcCorrectedUnderLines = std::shared_ptr<CSpectrum>( new CSpectrum(spectrum) );
    UInt32 spectrumSampleCount = spectrum.GetSampleCount();
    m_SpcFluxAxis.SetSize( spectrumSampleCount );
    m_SpcContinuumFluxAxis = spectrumContinuum.GetFluxAxis();
    m_ContinuumWinsize = spectrumContinuum.GetMedianWinsize();
    Log.LogDetail("    model: Continuum winsize found is %.2f A", m_ContinuumWinsize);

    m_ContinuumFluxAxis.SetSize( spectrumSampleCount );
    m_SpcFluxAxisModelDerivVelEmi.SetSize( spectrumSampleCount );
    m_SpcFluxAxisModelDerivVelAbs.SetSize( spectrumSampleCount );
    CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();
    const CSpectrumFluxAxis& spectrumFluxAxis = spectrum.GetFluxAxis();
    m_spcFluxAxisNoContinuum.SetSize( spectrumSampleCount );

    const TFloat64List& error = spectrumFluxAxis.GetError();
    TFloat64List& errorSpc = m_SpcFluxAxis.GetError();
    TFloat64List& errorSpcContinuum = m_SpcContinuumFluxAxis.GetError();
    // sets the error vectors
    for( UInt32 i=0; i<spectrumSampleCount; i++ )
    {
        m_ErrorNoContinuum[i] = error[i];
        errorSpc[i] = error[i];
        errorSpcContinuum[i] = error[i];
    }

    if( m_ContinuumComponent=="nocontinuum" )
    {
        //the continuum is set to zero and the observed spectrum is the spectrum without continuum
        for( UInt32 i=0; i<modelFluxAxis.GetSamplesCount(); i++ )
        {
            modelFluxAxis[i] = 0.0;
            m_ContinuumFluxAxis[i] = 0.0;
            m_spcFluxAxisNoContinuum[i] = spectrumFluxAxis[i]-m_SpcContinuumFluxAxis[i];
            m_SpcFluxAxis[i] = m_spcFluxAxisNoContinuum[i];
        }
    }
    if( m_ContinuumComponent == "fromspectrum" || m_ContinuumComponent == "tplfit")
    {
        //the continuum is set to the SpcContinuum and the observed spectrum is the raw spectrum
        CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();
        for(UInt32 i=0; i<modelFluxAxis.GetSamplesCount(); i++)
        {
            m_ContinuumFluxAxis[i] = m_SpcContinuumFluxAxis[i];
            m_spcFluxAxisNoContinuum[i] = spectrumFluxAxis[i]-m_SpcContinuumFluxAxis[i];

            modelFluxAxis[i] = m_ContinuumFluxAxis[i];
            m_SpcFluxAxis[i] = spectrumFluxAxis[i];
        }
    }
    m_observeGridContinuumFlux = NULL;
    m_unscaleContinuumFluxAxisDerivZ =NULL;
    m_chiSquareOperator = NULL;

    //NB: fitContinuum_option: this is the initialization (default value), eventually overriden in SetFitContinuum_FitStore() when a fitStore gets available
    m_fitContinuum_option = 0; //0=interactive fitting, 1=use precomputed fit store, 2=use fixed values (typical use for second pass recompute)
    m_fitContinuum_tplfitStore = NULL;

    if(m_ContinuumComponent == "tplfit" ||opt_fittingmethod == "lmfit" )
    {
        m_chiSquareOperator = new COperatorChiSquare2(calibrationPath);
        m_observeGridContinuumFlux = new Float64[modelFluxAxis.GetSamplesCount()]();
        if(m_observeGridContinuumFlux == NULL){
          throw runtime_error("unable to allocate m_observeGridContinuumFlux");
        }
        if(0)
        {
            Log.LogInfo( "Elementlist: fitContinuum_dustfit = %d", m_fitContinuum_dustfit );
            Log.LogInfo( "Elementlist: fitContinuum_igm = %d", m_fitContinuum_igm );
            Log.LogInfo( "Elementlist: fitContinuum_outsidelinesmask = %d", m_fitContinuum_outsidelinesmask );
            Log.LogInfo( "Elementlist: fitContinuum_observedFrame = %d", m_fitContinuum_observedFrame );
        }
    }

    if(opt_fittingmethod=="lmfit"){
      m_lmfit_noContinuumTemplate = (m_ContinuumComponent == "fromspectrum");
      m_lmfit_bestTemplate = true;
      m_lmfit_fitContinuum = true;
      m_lmfit_fitEmissionVelocity = true;
      m_lmfit_fitAbsorptionVelocity = true;
      Log.LogInfo("Lmfit parameters : tplAlreadyFit %d ; dobestTemplate %d, fitContinuum %d, fitEmissionVelocity %d, fitAbsorptionVelocity %d",m_lmfit_noContinuumTemplate, m_lmfit_bestTemplate, m_lmfit_fitContinuum, m_lmfit_fitEmissionVelocity, m_lmfit_fitAbsorptionVelocity );
    }
    //*/


    // "New style" rules initialization:
    m_Regulament = new CRegulament ( );
    m_Regulament->CreateRulesFromJSONFiles( );
    m_Regulament->EnableRulesAccordingToParameters ( m_rulesoption );

    // Load the line catalog
    m_RestRayList = restRayList;
    Log.LogDebug( "About to load catalog." );
    if(m_rigidity != "tplshape")
    {
        //load the regular catalog
        LoadCatalog(restRayList);
    }else{
        //load the tplshape catalog with only 1 element for all lines
        //LoadCatalogOneMultiline(restRayList);
        //load the tplshape catalog with 2 elements: 1 for the Em lines + 1 for the Abs lines
        LoadCatalogTwoMultilinesAE(restRayList);
    }
    LogCatalogInfos();

    /*
    //check the continuum flux axis
    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    for(Int32 j=0; j<m_ContinuumFluxAxis.GetSamplesCount(); j++)
    {
        if(isnan(m_ContinuumFluxAxis[j]))
        {
            Log.LogError( "CLineModelElementList(): NaN value found for the ContinuumFluxAxis at lambda=%f", spectralAxis[j] );
            throw runtime_error("CLineModelElementList() NaN value found");
            break;
        }
    }
    */
}

/**
 * \brief Empty destructor.
 **/
CLineModelElementList::~CLineModelElementList()
{
    //Log.LogInfo("Linemodel: Elementlist destructor call");
    if(m_observeGridContinuumFlux)
    {
        delete[] m_observeGridContinuumFlux;
    }
    if(m_chiSquareOperator)
    {
        delete m_chiSquareOperator;
    }
}

void CLineModelElementList::initLambdaOffsets(std::string offsetsCatalogsRelPath)
{
  CLineCatalogsOffsets* ctlgOffsets = new CLineCatalogsOffsets();
  ctlgOffsets->Init(m_calibrationPath, offsetsCatalogsRelPath);
  // load static offset catalog, idx=0
  ctlgOffsets->SetLinesOffsets( *this, 0);

  // load auto stack, hack from reference catalog
  //std::string spcName = m_inputSpc->GetName();
  //ctlgOffsets->SetLinesOffsetsAutoSelectStack(*this, spcName);
}


Bool CLineModelElementList::initTplratioCatalogs(std::string opt_tplratioCatRelPath, Int32 opt_tplratio_ismFit)
{
    //tplshape catalog initialization : used for rigidities tplcorr and tplshape
    m_CatalogTplShape = new CRayCatalogsTplShape();

    bool ret = m_CatalogTplShape->Init(m_calibrationPath, opt_tplratioCatRelPath, opt_tplratio_ismFit);
    if(!ret)
    {
        Log.LogError("Unable to initialize the the tpl-shape catalogs. aborting...");
        return false;
    }
    m_CatalogTplShape->InitLineCorrespondingAmplitudes(*this);
    m_CatalogTplShape->SetMultilineNominalAmplitudesFast( *this, 0 );
    //m_CatalogTplShape->SetMultilineNominalAmplitudes( *this, 0 );
    //m_RestRayList = m_CatalogTplShape->GetRestLinesList(0);
    //LoadCatalog(m_RestRayList);
    //LogCatalogInfos();

    //Resize tplshape buffers
    m_ChisquareTplshape.resize(m_CatalogTplShape->GetCatalogsCount());
    m_ScaleMargCorrTplshape.resize(m_CatalogTplShape->GetCatalogsCount());
    m_StrongELPresentTplshape.resize(m_CatalogTplShape->GetCatalogsCount());
    m_FittedAmpTplshape.resize(m_CatalogTplShape->GetCatalogsCount());
    m_FittedErrorTplshape.resize(m_CatalogTplShape->GetCatalogsCount());
    m_MtmTplshape.resize(m_CatalogTplShape->GetCatalogsCount());
    m_DtmTplshape.resize(m_CatalogTplShape->GetCatalogsCount());
    for(Int32 ktplshape=0; ktplshape<m_CatalogTplShape->GetCatalogsCount(); ktplshape++)
    {
        m_FittedAmpTplshape[ktplshape].resize(m_Elements.size());
        m_FittedErrorTplshape[ktplshape].resize(m_Elements.size());
        m_MtmTplshape[ktplshape].resize(m_Elements.size());
        m_DtmTplshape[ktplshape].resize(m_Elements.size());
    }

    m_tplshapeLeastSquareFast = false;
    return true;
}

/**
 * @brief setPassMode
 * @param iPass
 * set the fitting parameters according the the iPass argument.
 * @return
 */
Int32 CLineModelElementList::setPassMode(Int32 iPass)
{
    if(iPass==1)
    {
        m_forceDisableLyaFitting = true;
        m_forcedisableTplratioISMfit = m_opt_firstpass_forcedisableTplratioISMfit;
        m_forcedisableMultipleContinuumfit = m_opt_firstpass_forcedisableMultipleContinuumfit;

        m_fittingmethod = m_opt_firstpass_fittingmethod;
    }
    if(iPass==2)
    {
        m_forceDisableLyaFitting = false;
        m_forcedisableTplratioISMfit = false;
        m_forcedisableMultipleContinuumfit=false;

        m_fittingmethod = m_opt_secondpass_fittingmethod;
    }


    return true;
}

/**
 * \brief Returns a pointer to m_SpectrumModel.
 **/
const CSpectrum& CLineModelElementList::GetModelSpectrum() const
{
    return *m_SpectrumModel;
}

CSpectrum CLineModelElementList::GetSpectrumModelContinuum() const
{
    //std::shared_ptr<CSpectrum> spc = std::shared_ptr<CSpectrum>( new CSpectrum(*m_SpectrumModel) );
    CSpectrum spc(*m_SpectrumModel);
    for(UInt32 k=0; k<spc.GetSampleCount(); k++)
    {
        spc.GetFluxAxis()[k]=m_ContinuumFluxAxis[k];
    }
    return spc;
}

/**
 * \brief Returns a pointer to a spectrum containing the observed spectrum with the fitted lines subtracted
 **/
const CSpectrum& CLineModelElementList::GetObservedSpectrumWithLinesRemoved(Int32 lineTypeFilter)
{
    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();

    //
    std::vector<Float64> fluxAndContinuum(spectralAxis.GetSamplesCount(), 0.0);
    if(lineTypeFilter==CRay::nType_Emission)
    {
        refreshModel(CRay::nType_Absorption);
        for( UInt32 t=0;t<spectralAxis.GetSamplesCount();t++)
        {
            fluxAndContinuum[t] = m_SpectrumModel->GetFluxAxis()[t];
        }
        refreshModel(CRay::nType_Emission);
    }else if(lineTypeFilter==CRay::nType_Absorption)
    {
        refreshModel(CRay::nType_Emission);
        for( UInt32 t=0;t<spectralAxis.GetSamplesCount();t++)
        {
            fluxAndContinuum[t] = m_SpectrumModel->GetFluxAxis()[t];
        }
        refreshModel(CRay::nType_Absorption);
    }

    const CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();

    //create new spectrum, which is corrected under the lines
    //std::shared_ptr<CSpectrum> spcCorrectedUnderLines= std::shared_ptr<CSpectrum>( new CSpectrum(*m_SpectrumModel));
    //CSpectrum spcCorrectedUnderLines(*m_SpectrumModel);
    CSpectrumFluxAxis& fluxAxisNothingUnderLines = m_SpcCorrectedUnderLines->GetFluxAxis();
    Float64* Y = fluxAxisNothingUnderLines.GetSamples();

    for( UInt32 t=0;t<spectralAxis.GetSamplesCount();t++)
    {
        Y[t] = m_SpcFluxAxis[t]-modelFluxAxis[t]+m_ContinuumFluxAxis[t];
    }

    bool enableSmoothBlendContinuUnderLines=true;
    if(enableSmoothBlendContinuUnderLines)
    {
        Float64 alphaMax = 0.9; //alpha blend = 0: only lineSubtractedFlux, alpha=1: only continuum

        std::vector<UInt32> validEltsIdx = GetModelValidElementsIndexes();
        std::vector<UInt32> nonZeroValidEltsIdx;
        for(UInt32 j=0; j<validEltsIdx.size(); j++)
        {
            Int32 lineType = m_Elements[validEltsIdx[j]]->m_Rays[0].GetType();
            if(!(lineTypeFilter==-1 || lineTypeFilter==lineType))
            {
                continue;
            }

            if(m_Elements[validEltsIdx[j]]->GetElementAmplitude()>0.0)
            {
                nonZeroValidEltsIdx.push_back(validEltsIdx[j]);
            }
        }
        if(nonZeroValidEltsIdx.size()>0)
        {
            std::vector<UInt32> supportIdxes = getSupportIndexes( nonZeroValidEltsIdx );
            if(supportIdxes.size()>0)
            {
                for( UInt32 i=1; i<supportIdxes.size(); i++ )
                {
                    Float64 weighting = GetWeightingAnyLineCenterProximity(supportIdxes[i], nonZeroValidEltsIdx);
                    Float64 alpha = alphaMax*weighting;
                    Y[supportIdxes[i]] = (1.-alpha)*Y[supportIdxes[i]] + alpha*(fluxAndContinuum[supportIdxes[i]]);
                }
            }
        }
    }

    return *m_SpcCorrectedUnderLines;
}



/**
 * @brief GetContinuumError
 * Estimate the error on the continuum in a given window around the line center.
 * 1. calculate the observed spectrum flux with lines subtracted (fitted line model)
 * 2. estimate the std in a window corresponding to the sliding median window size (default/hardcoded=150A)
 * @param subeIdx
 * @param spectralAxis
 * @param spectrumfluxAxis
 * @param redshift
 * @param continuumfluxAxis
 * @return -1: zero samples found for error estimation, -9: not enough samples found for error estimation, -99: bad continuum winsize
 */
Float64 CLineModelElementList::GetContinuumError(Int32 eIdx, Int32 subeIdx)
{
//    if(m_ContinuumWinsize<=0)
//    {
//        return -99;
//    }
    UInt32 nMinValueForErrorEstimation=10;

    const CSpectrum& noLinesSpectrum = GetObservedSpectrumWithLinesRemoved();
    const CSpectrumFluxAxis& noLinesFluxAxis = noLinesSpectrum.GetFluxAxis();
    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    TFloat64Range lambdaRange = spectralAxis.GetLambdaRange(); //using the full wavelength range for this error estimation
    Float64 winsizeAngstrom = 150;

    TInt32Range indexRange = m_Elements[eIdx]->EstimateIndexRange(subeIdx,
                                                                    spectralAxis,
                                                                    m_Redshift,
                                                                    lambdaRange,
                                                                    winsizeAngstrom);

    //estimate sum square error between continuum and nolines-spectrum
    Float64 sum=0.0;
    UInt32 nsum = 0;
    for( Int32 t=indexRange.GetBegin();t<indexRange.GetEnd();t++)
    {
        Float64 diff = noLinesFluxAxis[t]-m_ContinuumFluxAxis[t];
        sum += diff*diff;
        nsum++;
    }

    Float64 error = -1.;
    if(nsum<nMinValueForErrorEstimation || nsum<1)
    {
        return -9.;
    }else{
        error = sqrt(sum/nsum);
    }

    return error;
}

/**
 * \brief Returns a pointer to the (re-)estimated continuum flux.
 **/
const CSpectrumFluxAxis &CLineModelElementList::GetModelContinuum() const
{
    return m_ContinuumFluxAxis;
}

/**
 * @brief CLineModelElementList::GetFluxDirectIntegration
 * Integrates the flux (F-continuum) in a lbda range around the center observed wavelength of the line.
 * The wavelength range is defined by the instrument resolution and a hardcoded nsigma factor
 * @param eIdx
 * @param subeIdx
 * @return
 */
Int32 CLineModelElementList::GetFluxDirectIntegration(TInt32List eIdx_list, TInt32List subeIdx_list, Float64& fluxdi, Float64& snrdi)
{
    Int32 ret=0;
    Float64 nsigma = 6.; //total range: ie. range will be mu-nsigma/2; mu+nsigma/2

    Int32 nlines = eIdx_list.size();
    if(nlines!=subeIdx_list.size())
    {
        return -1;
    }

    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    TFloat64Range lambdaRange = spectralAxis.GetLambdaRange(); //using the full wavelength range for this error estimation

    TInt32List indexes;
    for(Int32 kl=0; kl<nlines; kl++)
    {
        Int32 eIdx=eIdx_list[kl];
        Int32 subeIdx=subeIdx_list[kl];

        Float64 mu = m_Elements[eIdx]->GetObservedPosition(subeIdx, m_Redshift);
        Float64 instrumentSigma = mu/m_resolution;
        Float64 winsizeAngstrom = instrumentSigma*nsigma;
        TInt32Range indexRange = m_Elements[eIdx]->EstimateIndexRange(subeIdx,
                                                                      spectralAxis,
                                                                      m_Redshift,
                                                                      lambdaRange,
                                                                      winsizeAngstrom);


        for( Int32 t=indexRange.GetBegin();t<indexRange.GetEnd()-1;t++)
        {
            indexes.push_back(t);
        }
    }
    //remove duplicate indexes
    std::sort(indexes.begin(), indexes.end());
    indexes.erase( std::unique( indexes.begin(), indexes.end() ), indexes.end() );

    Int32 n_indexes = indexes.size();

    //estimate the integrated flux between obs. spectrum and continuum: trapezoidal intg
    Float64 sumFlux=0.0;
    Float64 sumErr=0.0;
    UInt32 nsum = 0;
    Int32 t;
    for( Int32 kt=0;kt<n_indexes;kt++)
    {
        t=indexes[kt];
        //trapez
        Float64 fa = m_SpcFluxAxis[t]-m_ContinuumFluxAxis[t];
        Float64 fb = m_SpcFluxAxis[t+1]-m_ContinuumFluxAxis[t+1];
        Float64 diffFlux = (spectralAxis[t+1]-spectralAxis[t])*(fb+fa)*0.5;
        sumFlux += diffFlux;


        Float64 ea = m_ErrorNoContinuum[t]*m_ErrorNoContinuum[t];
        Float64 eb = m_ErrorNoContinuum[t+1]*m_ErrorNoContinuum[t+1];
        Float64 diffError = (spectralAxis[t+1]-spectralAxis[t])*(eb+ea)*0.5;
        sumErr += diffError;

        //direct - missing dlambda
        //sumFlux += m_SpcFluxAxis[t]-m_ContinuumFluxAxis[t];
        //sumErr += m_ErrorNoContinuum[t]*m_ErrorNoContinuum[t];

        nsum++;
    }

    fluxdi = 0.;
    if(nsum>0)
    {
        fluxdi = sumFlux;
    }else{
        ret=-1;
    }

    snrdi = -1;
    if(sumErr>0.0)
    {
        snrdi = fluxdi/sqrt(sumErr);
    }else{
        ret=-1;
    }

    return ret;

}

Float64 CLineModelElementList::getModelFluxVal(UInt32 idx) const
{
    CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();
    if(idx<modelFluxAxis.GetSamplesCount()){
        return modelFluxAxis[idx];
    }

    return -1.0;
}


Float64 CLineModelElementList::getModelFluxDerivEltVal(UInt32 DerivEltIdx, UInt32 idx) const
{

    CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    Int32 iElts=DerivEltIdx;
    if(idx<spectralAxis.GetSamplesCount())
    {
        Float64 derivateVal = m_Elements[iElts]->GetModelDerivAmplitudeAtLambda(spectralAxis[idx], m_Redshift, m_ContinuumFluxAxis[idx]);
        return derivateVal;
    }

    return -1.0;
}

Float64 CLineModelElementList::getModelFluxDerivContinuumAmpEltVal(UInt32 DerivEltIdx, UInt32 idx) const
{

    CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    Int32 iElts=DerivEltIdx;
    if(idx<spectralAxis.GetSamplesCount())
    {
        Float64 derivateVal = m_Elements[iElts]->GetModelDerivContinuumAmpAtLambda(spectralAxis[idx], m_Redshift, m_observeGridContinuumFlux[idx]);
        return derivateVal;
    }

    return -1.0;
}

Float64 CLineModelElementList::getModelFluxDerivZContinuumVal(UInt32 idx)const{
  CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
  if(idx<spectralAxis.GetSamplesCount())
  {
      Float64 derivateVal = m_ContinuumFluxAxis[idx] * (-1) * spectralAxis[idx] / (1+m_Redshift) / (1+m_Redshift);
      return derivateVal;
  }
  // if(idx< spectralAxis.GetSamplesCount()){
  //   return m_unscaleContinuumFluxAxisDerivZ[idx] * m_fitContinuum_tplFitAmplitude;
  // }
  return -1.0;
}

// void CLineModelElementList::calculateUnscaleContinuumDerivZ(){
//
//   CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
//
//   if(m_unscaleContinuumFluxAxisDerivZ == NULL){
//     m_unscaleContinuumFluxAxisDerivZ = new Float64[spectralAxis.GetSamplesCount()]();
//     if(m_unscaleContinuumFluxAxisDerivZ == NULL){
//       Log.LogError("unable to allocate m_unscaleContinuumFluxAxisDerivZ");
//       return;
//     }
//   }
//   for (int idx = 0;idx<spectralAxis.GetSamplesCount(); idx ++ ){
//     m_unscaleContinuumFluxAxisDerivZ[idx] = m_observeGridContinuumFlux[idx] * spectralAxis[idx] / (1+m_Redshift) / (1+m_Redshift);
//   }
//
// }

/*
 * Calculate the partial deriv of the flux by z, when the continuum is not a variable of z
 */
Float64 CLineModelElementList::getModelFluxDerivZEltValNoContinuum(UInt32 DerivEltIdx, UInt32 idx) const{
  CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
  Int32 iElts=DerivEltIdx;
  if(idx<spectralAxis.GetSamplesCount())
  {
      Float64 derivateVal = m_Elements[iElts]->GetModelDerivZAtLambdaNoContinuum(spectralAxis[idx], m_Redshift, m_ContinuumFluxAxis[idx]);
      return derivateVal;
  }

  return -1.0;
}

/*
 * Calculate the partial deriv of the flux by z, when the continuum is a variable of z, the partial deriv of z of the continuum at this lambda is given by continuumFluxDerivZ
 */
Float64 CLineModelElementList::getModelFluxDerivZEltVal(UInt32 DerivEltIdx, UInt32 idx, Float64 continuumFluxDerivZ) const{
  CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
  Int32 iElts=DerivEltIdx;
  if(idx<spectralAxis.GetSamplesCount())
  {
      Float64 derivateVal = m_Elements[iElts]->GetModelDerivZAtLambda(spectralAxis[idx], m_Redshift, m_ContinuumFluxAxis[idx], continuumFluxDerivZ);
      return derivateVal;
  }

  return -1.0;
}

/*
Change on this function, m_SpcFluxAxisModelDerivVel has been remove,
replace by m_SpcFluxAxisModelDerivVelEmi and m_SpcFluxAxisModelDerivVelAbs
*/
Float64 CLineModelElementList::getModelFluxDerivVelVal(UInt32 idx) const
{
    if(idx<m_SpcFluxAxisModelDerivVelEmi.GetSamplesCount()){
        return m_SpcFluxAxisModelDerivVelEmi[idx]+ m_SpcFluxAxisModelDerivVelAbs[idx];
    }

    return -1.0;
}

Float64 CLineModelElementList::getModelFluxDerivVelEmissionVal(UInt32 idx) const
{
    if(idx<m_SpcFluxAxisModelDerivVelEmi.GetSamplesCount()){
        return m_SpcFluxAxisModelDerivVelEmi[idx];
    }

    return -1.0;
}

Float64 CLineModelElementList::getModelFluxDerivVelAbsorptionVal(UInt32 idx) const
{
    if(idx<m_SpcFluxAxisModelDerivVelAbs.GetSamplesCount()){
        return m_SpcFluxAxisModelDerivVelAbs[idx];
    }

    return -1.0;
}



/**
 * \brief For each ray in each group of the argument, finds the associated line in the catalog and saves this information to m_Elements.
 * Converts the argument restRayList to a group list.
 * For each entry in this list:
 *   For each ray in this entry:
 *     Finds the index in the catalog from the ray name and type.
 *     Saves the line, the catalog index and the nominal amplitude for the line thusly associated to this ray.
 * If at least one line was found, save this result in m_Elements.
 **/
void CLineModelElementList::LoadCatalog(const CRayCatalog::TRayVector& restRayList)
{
    CRayCatalog crctlg;
    std::vector<CRayCatalog::TRayVector> groupList = crctlg.ConvertToGroupList(restRayList);
    for(UInt32 ig=0; ig<groupList.size(); ig++)
    {
        std::vector<CRay> lines;
        std::vector<Float64> amps;
        std::vector<UInt32> inds;
        for(UInt32 i=0; i<groupList[ig].size(); i++)
        {
            std::vector<UInt32> idx = findLineIdxInCatalog( restRayList, groupList[ig][i].GetName(), groupList[ig][i].GetType());
            inds.push_back(idx[0]);
            amps.push_back(groupList[ig][i].GetNominalAmplitude());
            lines.push_back(groupList[ig][i]);
        }
        if(lines.size()>0)
        {
            m_Elements.push_back(boost::shared_ptr<CLineModelElement> (new CMultiLine(lines, m_LineWidthType, m_resolution, m_velocityEmission, m_velocityAbsorption, amps, m_nominalWidthDefaultAbsorption, inds)));
        }
    }
}

void CLineModelElementList::LoadCatalogOneMultiline(const CRayCatalog::TRayVector& restRayList)
{

    std::vector<CRay> lines;
    std::vector<Float64> amps;
    std::vector<UInt32> inds;
    for(UInt32 ir=0; ir<restRayList.size(); ir++)
    {
        inds.push_back(ir);
        amps.push_back(restRayList[ir].GetNominalAmplitude());
        lines.push_back(restRayList[ir]);
    }

    if(lines.size()>0)
    {
        m_Elements.push_back(boost::shared_ptr<CLineModelElement> (new CMultiLine(lines, m_LineWidthType, m_resolution, m_velocityEmission, m_velocityAbsorption, amps, m_nominalWidthDefaultAbsorption, inds)));
    }
}

void CLineModelElementList::LoadCatalogTwoMultilinesAE(const CRayCatalog::TRayVector& restRayList)
{
    std::vector<CRay::EType> types = { CRay::nType_Absorption, CRay::nType_Emission };

    for(Int32 iType=0; iType<2; iType++)
    {
        std::vector<CRay> lines;
        std::vector<Float64> amps;
        std::vector<UInt32> inds;
        for(UInt32 ir=0; ir<restRayList.size(); ir++)
        {
            if( restRayList[ir].GetType() == types[iType]){
                inds.push_back(ir);
                amps.push_back(restRayList[ir].GetNominalAmplitude());
                lines.push_back(restRayList[ir]);
            }
        }

        if(lines.size()>0)
        {
            m_Elements.push_back(boost::shared_ptr<CLineModelElement> (new CMultiLine(lines, m_LineWidthType, m_resolution, m_velocityEmission, m_velocityAbsorption, amps, m_nominalWidthDefaultAbsorption, inds)));
        }
    }
}



/**
 * \brief LogDetail the number of lines for each element, and their nominal amplitudes.
 **/
void CLineModelElementList::LogCatalogInfos()
{
    Log.LogDetail( "\n");
    Log.LogDetail( "LineModel Infos: %d elements", m_Elements.size());
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        Int32 nRays = m_Elements[iElts]->GetSize();
        if(nRays<1)
        {
            Log.LogDetail( "LineModel ctlg: elt %d (%s): no lines", iElts, m_Elements[iElts]->GetElementTypeTag().c_str());

        }
        for(UInt32 j=0; j<nRays; j++){
            std::string nominalAmpStr = "";
            if(nRays>0){
                nominalAmpStr = boost::str(boost::format("(nominal amp = %.4e)") % m_Elements[iElts]->GetNominalAmplitude(j));
            }
            Log.LogDetail( "LineModel ctlg: elt %d (%s): line %d = %s %s", iElts, m_Elements[iElts]->GetElementTypeTag().c_str(), j, m_Elements[iElts]->GetRayName(j).c_str(), nominalAmpStr.c_str());
        }
    }
    Log.LogDetail( "\n");
}


void CLineModelElementList::LoadFitContinuumOneTemplate(const TFloat64Range& lambdaRange, const CTemplate& tpl){
  Float64 merit = DBL_MAX;

  Float64 fitContinuumAmplitude = -1.0;
  Float64 fitDustCoeff = -1.0;
  Int32 fitMeiksinIdx = -1;
  Float64 fitDtM = -1.0;
  Float64 fitMtM = -1.0;
  Float64 overlapThreshold = 1.0;

  bool ignoreLinesSupport=m_fitContinuum_outsidelinesmask;
  std::vector<CMask> maskList;
  if(ignoreLinesSupport){
      maskList.resize(1);
      maskList[0]=getOutsideLinesMask();
  }
  std::vector<Float64> redshifts(1, m_Redshift);
  std::string opt_interp = "precomputedfinegrid"; //"lin";
  Int32 opt_extinction = m_fitContinuum_igm;
  Int32 opt_dustFit = m_fitContinuum_dustfit;

  if(m_observeGridContinuumFlux == NULL)
  {
    throw runtime_error("Elementlist, cannot SolveContinuum without m_observeGridContinuumFlux");
  }
  Bool ret = SolveContinuum( *m_inputSpc,
                             tpl,
                             lambdaRange,
                             redshifts,
                             overlapThreshold,
                             maskList,
                             opt_interp,
                             opt_extinction,
                             opt_dustFit,
                             merit,
                             fitContinuumAmplitude,
                             fitDustCoeff,
                             fitMeiksinIdx,
                             fitDtM,
                             fitMtM);
  /*//debug
  // export for debug
  FILE* fspc = fopen( "Continuum.txt", "w+" );
  for( Int32 t=0;t<m_ContinuumFluxAxis.GetSamplesCount();t++)
  {
      fprintf( fspc, "%f %f %f\n", spectralAxis[t],m_ContinuumFluxAxis[t]);
  }
  fclose( fspc );
  //*/
  Log.LogInfo("LMfit : Continuum Amp set Init at %0.00f", fitContinuumAmplitude);
  Log.LogInfo("Lmfit: solve succes %s", (ret? "t": "f"));
  ApplyContinuumOnGrid(tpl, m_Redshift);
}

/**
 * \brief Generates a continuum from the fitting with a set of templates : uses the chisquare2 operator
 **/
void CLineModelElementList::LoadFitContinuum(const TFloat64Range& lambdaRange, Int32 icontinuum)
{
    Log.LogDebug("Elementlist, m_fitContinuum_option=%d", m_fitContinuum_option);
    if(m_observeGridContinuumFlux == NULL)
    {
      throw runtime_error("Elementlist, cannot loadfitcontinuum without precomputedGridTplFlux");
    }

    Float64 bestMerit = DBL_MAX;
    Float64 bestFitAmplitude = -1.0;
    Float64 bestFitDustCoeff = -1.0;
    Int32 bestFitMeiksinIdx = -1;
    Float64 bestFitRedshift = m_Redshift;
    Float64 bestFitDtM = -1.0;
    Float64 bestFitMtM = -1.0;
    std::string bestTplName="";
    std::vector<Float64> bestFitPolyCoeffs;

    if(m_fitContinuum_option==0){
        //hardcoded parameters
        std::string opt_interp = "lin"; //"precomputedfinegrid"; //
        Int32 opt_extinction = m_fitContinuum_igm;
        Int32 opt_dustFit = m_fitContinuum_dustfit;
        Float64 overlapThreshold = 1.0;

        bool ignoreLinesSupport=m_fitContinuum_outsidelinesmask;
        std::vector<CMask> maskList;
        if(ignoreLinesSupport){
            maskList.resize(1);
            maskList[0]=getOutsideLinesMask();
        }

        std::vector<Float64> redshifts(1, m_Redshift);
        if(m_fitContinuum_observedFrame){
            redshifts[0] = 0.0;
        }

        for( UInt32 i=0; i<m_tplCategoryList.size(); i++ )
        {
            std::string category = m_tplCategoryList[i];

            for( UInt32 j=0; j<m_tplCatalog.GetTemplateCount( category ); j++ )
            {
                const CTemplate& tpl = m_tplCatalog.GetTemplate( category, j );

                Float64 merit = DBL_MAX;
                Float64 fitAmplitude = -1.0;
                Float64 fitDustCoeff = -1.0;
                Int32 fitMeiksinIdx = -1;
                Float64 fitDtM = -1.0;
                Float64 fitMtM = -1.0;
                Bool ret = SolveContinuum( *m_inputSpc, tpl, lambdaRange,
					   redshifts, overlapThreshold,
					   maskList, opt_interp, opt_extinction,
					   opt_dustFit, merit, fitAmplitude,
					   fitDustCoeff, fitMeiksinIdx, fitDtM, fitMtM);

                if(ret)
                {
                    if(merit<bestMerit)
                    {
                        bestMerit = merit;
                        bestFitAmplitude = fitAmplitude;
                        bestFitDustCoeff = fitDustCoeff;
                        bestFitMeiksinIdx = fitMeiksinIdx;
                        bestFitDtM = fitDtM;
                        bestFitMtM = fitMtM;
                        bestTplName = tpl.GetName();
                    }
                }else{
                    Log.LogError("Failed to load-fit continuum");
                    throw std::runtime_error("Failed to load-fit continuum");
                }
            }
        }
    }else if(m_fitContinuum_option==1){
        CTemplatesFitStore::TemplateFitValues fitValues = m_fitContinuum_tplfitStore->GetFitValues(m_Redshift, icontinuum);
        if(fitValues.tplName=="")
        {
            throw runtime_error("Empty template name");
        }
        bestMerit = fitValues.merit;
        bestFitAmplitude = fitValues.fitAmplitude;
        bestFitDustCoeff = fitValues.ismDustCoeff;
        bestFitMeiksinIdx = fitValues.igmMeiksinIdx;
        bestFitDtM = fitValues.fitDtM;
        bestFitMtM = fitValues.fitMtM;
        bestTplName = fitValues.tplName;
    }else if(m_fitContinuum_option==2){
        //values unmodified
        bestTplName = m_fitContinuum_tplName;
        bestMerit = m_fitContinuum_tplFitMerit;
        bestFitAmplitude = m_fitContinuum_tplFitAmplitude;
        bestFitDustCoeff = m_fitContinuum_tplFitDustCoeff;
        bestFitMeiksinIdx = m_fitContinuum_tplFitMeiksinIdx;
        bestFitDtM = m_fitContinuum_tplFitDtM;
        bestFitMtM = m_fitContinuum_tplFitMtM;
        bestFitRedshift = m_fitContinuum_tplFitRedshift;
        bestFitPolyCoeffs = m_fitContinuum_tplFitPolyCoeffs;
    }else if(m_fitContinuum_option==3){
        //redo the fit only for the current template/ISM/IGM continuum

        //hardcoded parameters
        std::string opt_interp = "lin"; //"precomputedfinegrid"; //
        Int32 opt_extinction = m_fitContinuum_igm;
        Int32 opt_dustFit = m_fitContinuum_dustfit;
        Float64 overlapThreshold = 1.0;

        bool ignoreLinesSupport=m_fitContinuum_outsidelinesmask;
        std::vector<CMask> maskList;
        if(ignoreLinesSupport){
            maskList.resize(1);
            maskList[0]=getOutsideLinesMask();
        }

        std::vector<Float64> redshifts(1, m_Redshift);
        if(m_fitContinuum_observedFrame){
            redshifts[0] = 0.0;
        }

        for( UInt32 i=0; i<m_tplCategoryList.size(); i++ )
        {
            std::string category = m_tplCategoryList[i];

            for( UInt32 j=0; j<m_tplCatalog.GetTemplateCount( category ); j++ )
            {
                const CTemplate& tpl = m_tplCatalog.GetTemplate( category, j );

                if(tpl.GetName()==m_fitContinuum_tplName)
                {
                    Float64 merit = DBL_MAX;
                    Float64 fitAmplitude = -1.0;
                    Float64 fitDustCoeff = -1.0;
                    Int32 fitMeiksinIdx = -1;
                    Float64 fitDtM = -1.0;
                    Float64 fitMtM = -1.0;
                    Bool ret = SolveContinuum( *m_inputSpc,
                                               tpl,
                                               lambdaRange,
                                               redshifts,
                                               overlapThreshold,
                                               maskList,
                                               opt_interp,
                                               opt_extinction,
                                               opt_dustFit,
                                               merit,
                                               fitAmplitude,
                                               fitDustCoeff,
                                               fitMeiksinIdx,
                                               fitDtM,
                                               fitMtM);

                    if(ret)
                    {
                        bestMerit = merit;
                        bestFitAmplitude = fitAmplitude;
                        bestFitDustCoeff = fitDustCoeff;
                        bestFitMeiksinIdx = fitMeiksinIdx;
                        bestFitDtM = fitDtM;
                        bestFitMtM = fitMtM;
                        bestTplName = tpl.GetName();
                    }
                    break;
                }
            }
        }
    }else{
        throw runtime_error("Elementlist, cannot parse fitContinuum_option");
    }

    if(bestTplName!="")
    {
        //Log.LogInfo( "For z=%.5f : Best continuum tpl found: %s", m_Redshift, bestTplName.c_str());
        //
        //Retrieve the best template
        for( UInt32 i=0; i<m_tplCategoryList.size(); i++ )
        {
            std::string category = m_tplCategoryList[i];

            for( UInt32 j=0; j<m_tplCatalog.GetTemplateCount( category ); j++ )
            {
                const CTemplate& tpl = m_tplCatalog.GetTemplate( category, j );

                if(tpl.GetName()==bestTplName)
                {
                    m_fitContinuum_tplFitAmplitude = bestFitAmplitude;
                    m_fitContinuum_tplFitMerit = bestMerit;
                    m_fitContinuum_tplFitDustCoeff = bestFitDustCoeff;
                    m_fitContinuum_tplFitMeiksinIdx = bestFitMeiksinIdx;
                    m_fitContinuum_tplFitRedshift = bestFitRedshift;

                    ApplyContinuumOnGrid(tpl, m_fitContinuum_tplFitRedshift);

                    m_fitContinuum_tplFitDtM = bestFitDtM;
                    m_fitContinuum_tplFitMtM = bestFitMtM;
                    m_fitContinuum_tplFitPolyCoeffs = bestFitPolyCoeffs;
                    setFitContinuum_tplAmplitude(bestFitAmplitude, bestFitPolyCoeffs);

                    Log.LogDebug( "    model : LoadFitContinuum, loaded: %s", bestTplName.c_str());
                    break;
                }
            }
        }
    }else{
      throw runtime_error( "Failed to load and fit continuum");
    }
}

void CLineModelElementList::setFitContinuum_tplAmplitude(Float64 tplAmp, std::vector<Float64> polyCoeffs){
    const CSpectrumSpectralAxis& spcSpectralAxis = m_SpectrumModel->GetSpectralAxis();

    //Float64 alpha = 0.5; //alpha blend = 1: only m_SpcContinuumFluxAxis, alpha=0: only tplfit
    m_fitContinuum_tplFitAmplitude = tplAmp;
    m_fitContinuum_tplFitPolyCoeffs = polyCoeffs;
    for (UInt32 k=0; k<m_ContinuumFluxAxis.GetSamplesCount(); k++){
        //m_ContinuumFluxAxis[k] = (1.-alpha)*m_observeGridContinuumFlux[k]*tplAmp + (alpha)*m_SpcContinuumFluxAxis[k];
        m_ContinuumFluxAxis[k] = m_observeGridContinuumFlux[k]*tplAmp;
        for(Int32 kCoeff=0; kCoeff<polyCoeffs.size(); kCoeff++)
        {
            Float64 lbdaTerm=1.0;
            for(Int32 kt=0; kt<kCoeff; kt++)
            {
                lbdaTerm *=  spcSpectralAxis[k];
            }
            m_ContinuumFluxAxis[k] += polyCoeffs[kCoeff]*lbdaTerm;
        }
        m_spcFluxAxisNoContinuum[k] = m_SpcFluxAxis[k]-m_ContinuumFluxAxis[k];
    }
}



/*
Change the actaul value of redshift.
the continuum can be reinterpolate.
*/
void CLineModelElementList::setRedshift(Float64 redshift, bool reinterpolatedContinuum){
  m_Redshift = redshift;

  if(reinterpolatedContinuum){
    for( UInt32 i=0; i<m_tplCategoryList.size(); i++ )
    {
        std::string category = m_tplCategoryList[i];

        for( UInt32 j=0; j<m_tplCatalog.GetTemplateCount( category ); j++ )
        {
            const CTemplate& tpl = m_tplCatalog.GetTemplate( category, j );

            if(tpl.GetName()==m_fitContinuum_tplName)
            {
              ApplyContinuumOnGrid(tpl, redshift);
            }
          }
        }
  }

}

/**
Apply the template continuum by interpolating the grid as define in Init COntinuum
**/
Int32 CLineModelElementList::ApplyContinuumOnGrid(const CTemplate& tpl, Float64 zcontinuum){
    m_fitContinuum_tplName = tpl.GetName();

    Int32 n = tpl.GetSampleCount();
    CSpectrumFluxAxis tplFluxAxis = tpl.GetFluxAxis();
    const CSpectrumSpectralAxis& tplSpectralAxis = tpl.GetSpectralAxis();

    //inialize and allocate the gsl objects
    Float64* Ysrc = tplFluxAxis.GetSamples();
    const Float64* Xsrc = tplSpectralAxis.GetSamples();
    //apply dust attenuation
    const Float64* dustCoeffArray = m_chiSquareOperator->getDustCoeff(m_fitContinuum_tplFitDustCoeff, Xsrc[n-1]);
    //apply igm meiksin extinction
    const Float64* meiksinCoeffArray = m_chiSquareOperator->getMeiksinCoeff(m_fitContinuum_tplFitMeiksinIdx, m_Redshift, Xsrc[n-1]);

    if(dustCoeffArray!=0)
    {
        Float64 lambda = 0.0;
        for(Int32 ktpl=0; ktpl<n; ktpl++)
        {
            lambda = Xsrc[ktpl];
            Ysrc[ktpl]*=dustCoeffArray[Int32(lambda)]; //dust coeff is rounded at the nearest 1 angstrom value
        }
        delete dustCoeffArray;
    }
    if(meiksinCoeffArray!=0)
    {
        Float64 lambda = 0.0;
        for(Int32 ktpl=0; ktpl<n; ktpl++)
        {
            lambda = Xsrc[ktpl];
            Ysrc[ktpl]*=meiksinCoeffArray[Int32(lambda)]; //igm meiksin coeff is rounded at the nearest 1 angstrom value
        }
        delete meiksinCoeffArray;
    }

    //spline
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n);
    gsl_spline_init (spline, Xsrc, Ysrc, n);
    gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();
    Int32 k = 0;
    Float64 x = 0.0;
    const CSpectrumSpectralAxis& spcSpectralAxis = m_SpectrumModel->GetSpectralAxis();

    for(k=0; k<spcSpectralAxis.GetSamplesCount(); k++){
        x = spcSpectralAxis[k]/(1+zcontinuum);
        if(x < tplSpectralAxis[0] || x > tplSpectralAxis[n-1]){
            m_observeGridContinuumFlux[k] = 0.0;
        }else{
            m_observeGridContinuumFlux[k] = gsl_spline_eval (spline, x, accelerator);//m_fitContinuum_tplFitAmplitude*

        }
        /*//debug:
      // save reflex data
      FILE* f = fopen( "observegrid.txt", "w+" );
      for( Int32 t=0;t<spcSpectralAxis.GetSamplesCount();t++)
      {
          fprintf( f, "%f %f\n", spcSpectralAxis[t], m_observeGridContinuumFlux[t]);
      }
      fclose( f );
      //*/

    }



  gsl_spline_free (spline);
  gsl_interp_accel_free (accelerator);
  return 0;
}

Bool CLineModelElementList::SolveContinuum(const CSpectrum& spectrum,
                                           const CTemplate& tpl,
                                           const TFloat64Range& lambdaRange,
                                           const TFloat64List& redshifts,
                                           Float64 overlapThreshold,
                                           std::vector<CMask> maskList,
                                           std::string opt_interp,
                                           Int32 opt_extinction,
                                           Int32 opt_dustFit,
                                           Float64& merit,
                                           Float64& fitAmplitude,
                                           Float64& fitDustCoeff,
                                           Int32& fitMeiksinIdx,
                                           Float64& fitDtM,
                                           Float64& fitMtM)
{
    // Compute merit function
    //Log.LogInfo("Solving continuum for %s at z=%.4e", tpl.GetName().c_str(), redshifts[0]);
    //CRef<CChisquareResult>  chisquareResult = (CChisquareResult*)chiSquare.ExportChi2versusAZ( _spc, _tpl, lambdaRange, redshifts, overlapThreshold );
    auto  chisquareResult = std::dynamic_pointer_cast<CChisquareResult>( m_chiSquareOperator->Compute( spectrum,
                                                                                                       tpl,
                                                                                                       lambdaRange,
                                                                                                       redshifts,
                                                                                                       overlapThreshold,
                                                                                                       maskList,
                                                                                                       opt_interp,
                                                                                                       opt_extinction,
                                                                                                       opt_dustFit ) );
    if( !chisquareResult )
    {

        Log.LogInfo( "Failed to compute chi square value");
        return false;
    }else{
        // Store results
        merit = chisquareResult->ChiSquare[0];
        fitAmplitude = chisquareResult->FitAmplitude[0];
        fitDustCoeff = chisquareResult->FitDustCoeff[0];
        fitMeiksinIdx = chisquareResult->FitMeiksinIdx[0];
        fitDtM = chisquareResult->FitDtM[0];
        fitMtM = chisquareResult->FitMtM[0];
        return true;
    }

}

Int32 CLineModelElementList::LoadFitContaminantTemplate(const TFloat64Range& lambdaRange, const CTemplate& tpl){
    Float64 fitAmplitude = -1.0;
    Float64 overlapThreshold = 1.0;

    bool ignoreLinesSupport=false;
    std::vector<CMask> maskList;
    if(ignoreLinesSupport){
        maskList.resize(1);
        maskList[0]=getOutsideLinesMask();
    }
    std::vector<Float64> redshifts(1, 0.0); //fitting an already redshifted model
    std::string opt_interp = "lin";
    Int32 opt_extinction = 0;
    Int32 opt_dustFit = -1;

    //prepare observed spectrum without continuum
    const CSpectrumFluxAxis& inputSpectrumFluxAxis = m_inputSpc->GetFluxAxis();
    std::shared_ptr<CSpectrum> _spc = std::shared_ptr<CSpectrum>( new CSpectrum(*m_inputSpc) );
    CSpectrumFluxAxis _spcFluxAxis = _spc->GetFluxAxis();
    Float64* Y_spc = _spcFluxAxis.GetSamples();
    for(UInt32 i=0; i<inputSpectrumFluxAxis.GetSamplesCount(); i++)
    {
        Y_spc[i] = inputSpectrumFluxAxis[i]-m_ContinuumFluxAxis[i];
    }
    //prepare tpl contaminant
    const CSpectrumSpectralAxis& spcSpectralAxis = m_SpectrumModel->GetSpectralAxis();
    const std::string& category = "emission";
    m_tplContaminantSpcRebin = std::shared_ptr<CTemplate>( new CTemplate( "contaminantrebin", category ) );
    CSpectrumFluxAxis& tplContaminantRebinFluxAxis = m_tplContaminantSpcRebin->GetFluxAxis();
    CSpectrumSpectralAxis& tplContaminantRebinSpcAxis = m_tplContaminantSpcRebin->GetSpectralAxis();
    tplContaminantRebinFluxAxis.SetSize(spcSpectralAxis.GetSamplesCount());
    tplContaminantRebinSpcAxis.SetSize(spcSpectralAxis.GetSamplesCount());

    //prepare INTERPOLATED contaminant
    Int32 n = tpl.GetSampleCount();
    CSpectrumFluxAxis tplFluxAxis = tpl.GetFluxAxis();
    const CSpectrumSpectralAxis& tplSpectralAxis = tpl.GetSpectralAxis();
    Float64* Ysrc = tplFluxAxis.GetSamples();
    const Float64* Xsrc = tplSpectralAxis.GetSamples();
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n);
    gsl_spline_init (spline, Xsrc, Ysrc, n);
    gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();
    Int32 k = 0;
    Float64 x = 0.0;
    for(k=0; k<spcSpectralAxis.GetSamplesCount(); k++){
        x = spcSpectralAxis[k];
        tplContaminantRebinSpcAxis[k]=spcSpectralAxis[k];
        if(x < tplSpectralAxis[0] || x > tplSpectralAxis[n-1]){
            tplContaminantRebinFluxAxis[k] = 0.0;
        }else{
            tplContaminantRebinFluxAxis[k] = gsl_spline_eval (spline, x, accelerator);

        }
    }
    gsl_spline_free (spline);
    gsl_interp_accel_free (accelerator);


    //*
    //Fit contaminant template AMPLITUDE
    if(m_chiSquareOperator == NULL)
    {
        m_chiSquareOperator = new COperatorChiSquare2(m_calibrationPath);
    }
    auto  chisquareResult = std::dynamic_pointer_cast<CChisquareResult>( m_chiSquareOperator->Compute( *_spc,
                                                                                                       *m_tplContaminantSpcRebin,
                                                                                                       lambdaRange,
                                                                                                       redshifts,
                                                                                                       overlapThreshold,
                                                                                                       maskList,
                                                                                                       opt_interp,
                                                                                                       opt_extinction,
                                                                                                       opt_dustFit ) );
    if( !chisquareResult )
    {

        Log.LogError( "    model: contaminant - Failed to compute chi square value");
        return false;
    }else{
        //
        fitAmplitude = chisquareResult->FitAmplitude[0];
    }
    Log.LogInfo( "    model: contaminant raw fit amplitude is %f.", fitAmplitude);
    //*/
    if(fitAmplitude > 1.0)
    {
        Log.LogInfo( "    model: contaminant raw fit amplitude was %f. NOW SET TO 1.0", fitAmplitude);
        fitAmplitude = 1.0;
    }
    //*/
    //fitAmplitude = 1.0;
    Log.LogInfo( "    model: contaminant fit amplitude = %f", fitAmplitude);
    for(k=0; k<spcSpectralAxis.GetSamplesCount(); k++)
    {
        tplContaminantRebinFluxAxis[k] = fitAmplitude*tplContaminantRebinFluxAxis[k];
    }

    /*//debug:
    FILE* f = fopen( "contaminantFittedInterpolated.txt", "w+" );
    for(k=0; k<m_SpcFluxAxis.GetSamplesCount(); k++)
    {
        fprintf( f, "%f\t%e\n", spcSpectralAxis[k], tplContaminantRebinFluxAxis[k]);
    }
    fclose( f );
    //*/

    if( m_ContinuumComponent=="nocontinuum" )
    {
        for(k=0; k<m_SpcFluxAxis.GetSamplesCount(); k++){
            m_SpcFluxAxis[k] = inputSpectrumFluxAxis[k]-m_ContinuumFluxAxis[k]-tplContaminantRebinFluxAxis[k];
        }
    }else if( m_ContinuumComponent == "fromspectrum" || m_ContinuumComponent == "tplfit")
    {
        for(k=0; k<m_SpcFluxAxis.GetSamplesCount(); k++){
            m_SpcFluxAxis[k] = inputSpectrumFluxAxis[k]-tplContaminantRebinFluxAxis[k];
        }
    }

    return 1;
}

std::shared_ptr<CModelSpectrumResult> CLineModelElementList::GetContaminantSpectrumResult()
{
    std::shared_ptr<CModelSpectrumResult>  resultcont = std::shared_ptr<CModelSpectrumResult>( new CModelSpectrumResult(*m_tplContaminantSpcRebin) );
    return resultcont;
}

std::string CLineModelElementList::getFitContinuum_tplName()
{
    return m_fitContinuum_tplName;
}

Float64 CLineModelElementList::getFitContinuum_tplAmplitude()
{
    return m_fitContinuum_tplFitAmplitude;
}

Float64 CLineModelElementList::getFitContinuum_tplMerit()
{
    return m_fitContinuum_tplFitMerit;
}

Float64 CLineModelElementList::getFitContinuum_tplIsmDustCoeff()
{
    return m_fitContinuum_tplFitDustCoeff;
}

Float64 CLineModelElementList::getFitContinuum_tplIgmMeiksinIdx()
{
    return m_fitContinuum_tplFitMeiksinIdx;
}

void CLineModelElementList::SetContinuumComponent(std::string component)
{
    m_ContinuumComponent = component;
}

Int32 CLineModelElementList::SetFitContinuum_FitStore(CTemplatesFitStore* fitStore)
{
    m_fitContinuum_option = 1; //enable use of the fit store
    Log.LogInfo( "Elementlist: enabling fitContinuum store.");
    m_fitContinuum_tplfitStore = fitStore;
    return 1;
}

void CLineModelElementList::SetFitContinuum_Option(Int32 opt)
{
    m_fitContinuum_option = opt;
}

Int32 CLineModelElementList::GetFitContinuum_Option()
{
    return m_fitContinuum_option;
}

void CLineModelElementList::SetFitContinuum_FitValues(std::string tplfit_name,
                                                      Float64 tplfit_amp,
                                                      Float64 tplfit_chi2,
                                                      Float64 tplfit_ebmv,
                                                      Int32 tplfit_meiksinidx,
                                                      Float64 tplfit_continuumredshift,
                                                      Float64 tplfit_dtm,
                                                      Float64 tplfit_mtm,
                                                      std::vector<Float64> polyCoeffs)
{
    m_fitContinuum_tplName = tplfit_name;
    m_fitContinuum_tplFitAmplitude = tplfit_amp;
    m_fitContinuum_tplFitMerit = tplfit_chi2;
    m_fitContinuum_tplFitDustCoeff = tplfit_ebmv;
    m_fitContinuum_tplFitMeiksinIdx = tplfit_meiksinidx;
    m_fitContinuum_tplFitRedshift = tplfit_continuumredshift;


    m_fitContinuum_tplFitDtM = tplfit_dtm;
    m_fitContinuum_tplFitMtM = tplfit_mtm;
    m_fitContinuum_tplFitPolyCoeffs = polyCoeffs;
}

/**
 * \brief This function prepares the continuum for use in the fit with the line elements.
 * Rebin with PFG buffer
 * Find and apply amplitude factor from previously fitted tpl
 **/
void CLineModelElementList::PrepareContinuum(Float64 z)
{
    const CSpectrumSpectralAxis& targetSpectralAxis = m_SpectrumModel->GetSpectralAxis();
    Float64* Yrebin = m_ContinuumFluxAxis.GetSamples();

    if(m_observeGridContinuumFlux == NULL){
        UInt32 nSamples = targetSpectralAxis.GetSamplesCount();
        for ( UInt32 i = 0; i<nSamples; i++)
        {
            Yrebin[i] = m_SpcContinuumFluxAxis[i];
        }
        return;
    }


    std::vector<Float64> polyCoeffs=m_fitContinuum_tplFitPolyCoeffs;
    setFitContinuum_tplAmplitude(m_fitContinuum_tplFitAmplitude, polyCoeffs);
    return;
}

std::string CLineModelElementList::getTplshape_bestTplName()
{
    return m_tplshapeBestTplName;
}
Float64 CLineModelElementList::getTplshape_bestTplIsmCoeff()
{
    return m_tplshapeBestTplIsmCoeff;
}

Float64 CLineModelElementList::getTplshape_bestAmplitude()
{
    return m_tplshapeBestTplAmplitude;
}


Int32 CLineModelElementList::getTplshape_count()
{
    if(m_rigidity!="tplshape")
    {
        return 0;
    }
    return m_CatalogTplShape->GetCatalogsCount();
}

std::vector<Float64> CLineModelElementList::getTplshape_priors()
{
    if(m_rigidity!="tplshape")
    {
        std::vector<Float64> dumb;
        return dumb;
    }
    return m_CatalogTplShape->getCatalogsPriors();
}

std::vector<Float64> CLineModelElementList::GetChisquareTplshape()
{
    return m_ChisquareTplshape;
}

std::vector<Float64> CLineModelElementList::GetScaleMargTplshape()
{
    return m_ScaleMargCorrTplshape;
}

std::vector<bool> CLineModelElementList::GetStrongELPresentTplshape()
{
    return m_StrongELPresentTplshape;
}

Bool CLineModelElementList::initModelAtZ(Float64 redshift, const TFloat64Range& lambdaRange, const CSpectrumSpectralAxis &spectralAxis)
{
    m_Redshift = redshift;

    //prepare the elements support
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        m_Elements[iElts]->prepareSupport(spectralAxis, redshift, lambdaRange);
    }

    return true;
}

Bool CLineModelElementList::setTplshapeModel(Int32 itplshape, Bool enableSetVelocity)
{
    m_CatalogTplShape->SetLyaProfile(*this, itplshape);

    //m_CatalogTplShape->SetMultilineNominalAmplitudes( *this, ifitting );
    m_CatalogTplShape->SetMultilineNominalAmplitudesFast( *this, itplshape );

    if(enableSetVelocity)
    {
        //Set the velocities from templates: todo auto switch when velfit is ON
        m_CatalogTplShape->GetCatalogVelocities(itplshape, m_velocityEmission, m_velocityAbsorption);
    }

    Log.LogDebug( "    model : setTplshapeModel, loaded: %d = %s", itplshape, m_CatalogTplShape->GetCatalogName(itplshape).c_str());
    return true;
}

Bool CLineModelElementList::setTplshapeAmplitude(std::vector<Float64> ampsElts, std::vector<Float64> errorsElts)
{
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        m_Elements[iElts]->SetFittedAmplitude(ampsElts[iElts], errorsElts[iElts]);
    }
    return true;
}



Bool CLineModelElementList::initDtd(const TFloat64Range& lambdaRange)
{
    m_dTransposeDLambdaRange = lambdaRange;
    m_dTransposeDNocontinuum = EstimateDTransposeD(lambdaRange, "nocontinuum");
    m_dTransposeDRaw = EstimateDTransposeD(lambdaRange, "raw");
    m_likelihood_cstLog = EstimateLikelihoodCstLog(lambdaRange);
    return true;
}

/**
 * \brief Prepares the context and fits the Linemodel to the spectrum, returning the merit of the fit.
 * Prepare the continuum.
 * Initialize the model spectrum.
 * Prepare the elements.
 * Fit the amplitudes of each element independently.
 * Fit the amplitude of all elements together with iterative solver: Nelder Mead Simplex.
 * Fit the amplitude of all elements together with linear solver: gsl_multifit_wlinear.
 * Fit the amplitudes of each element independently, unless there is overlap.
 * Apply a continuum iterative re-estimation with lines removed from the initial spectrum.
 * Apply rules.
 * Create spectrum model.
 * Return merit.
 **/
Float64 CLineModelElementList::fit(Float64 redshift,
                                   const TFloat64Range& lambdaRange,
                                   CLineModelSolution& modelSolution,
                                   CContinuumModelSolution &continuumModelSolution,
                                   Int32 contreest_iterations,
                                   bool enableLogging)
{
    //initialize the model spectrum
    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();

    if(m_fittingmethod == "lmfit"){
      SetVelocityEmission(m_velocityEmissionInit);
      SetVelocityAbsorption(m_velocityAbsorptionInit);
    }

    initModelAtZ(redshift, lambdaRange, spectralAxis);

    if(! (m_dTransposeDLambdaRange.GetBegin()==lambdaRange.GetBegin() && m_dTransposeDLambdaRange.GetEnd()==lambdaRange.GetEnd() ) )
    {
        initDtd(lambdaRange);
    }


    Int32 nfitting=1; //multiple fitting steps for rigidity=tplshape/tplratio
    if(m_rigidity=="tplshape")
    {
        nfitting=m_CatalogTplShape->GetCatalogsCount();
    }
    Int32 savedIdxFitted=-1; //for rigidity=tplshape

    Int32 ncontinuumfitting=1;
    Int32 savedIdxContinuumFitted=-1; //for continuum tplfit
    if(m_ContinuumComponent == "tplfit" && !m_forcedisableMultipleContinuumfit)
    {
        ncontinuumfitting=m_opt_fitcontinuum_maxCount;
    }


    Float64 merit = DBL_MAX; //initializing 'on the fly' best-merit
    std::vector<Float64> meritTplratio(nfitting, DBL_MAX); //initializing 'on the fly' best-merit per tplratio
    for(Int32 icontfitting=0; icontfitting<ncontinuumfitting; icontfitting++)
    {
        if(m_ContinuumComponent == "tplfit") //the support has to be already computed when LoadFitContinuum() is called
        {
            LoadFitContinuum(lambdaRange, icontfitting);
        }
        if(m_ContinuumComponent != "nocontinuum"){
            //prepare the continuum
            if(m_ContinuumComponent == "tplfit" && m_fitContinuum_observedFrame){
                PrepareContinuum(0.0);
            }else{
                PrepareContinuum(redshift);
            }
        }
        //EstimateSpectrumContinuum();

        if(m_ContinuumComponent == "nocontinuum")
        {
        }
        else
        {
            for(UInt32 i=0; i<modelFluxAxis.GetSamplesCount(); i++)
            {
                modelFluxAxis[i] = m_ContinuumFluxAxis[i];
                m_spcFluxAxisNoContinuum[i] = m_SpcFluxAxis[i]-m_ContinuumFluxAxis[i];
            }
        }

        if(m_enableAmplitudeOffsets)
        {
            prepareAmplitudeOffset(m_spcFluxAxisNoContinuum);
        }

        for(Int32 ifitting=0; ifitting<nfitting; ifitting++)
        {
            if(m_rigidity!="tplshape")
            {
                if(!m_forceDisableLyaFitting)
                {
                    //prepare the Lya width and asym coefficients if the asymfit profile option is met
                    setLyaProfile(redshift, spectralAxis);
                }
            }else{

                if(m_forcedisableTplratioISMfit)
                {
                    if(m_CatalogTplShape->GetIsmCoeff(ifitting)>0)
                    {
                        //copy the values for ebmv=ebmv_fixed (=0) here
                        m_ChisquareTplshape[ifitting] = m_ChisquareTplshape[ifitting-1];
                        m_ScaleMargCorrTplshape[ifitting] = m_ScaleMargCorrTplshape[ifitting-1];
                        m_StrongELPresentTplshape[ifitting] = m_StrongELPresentTplshape[ifitting-1];
                        meritTplratio[ifitting] = meritTplratio[ifitting-1];
                        for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
                        {
                            m_FittedAmpTplshape[ifitting][iElts] = m_FittedAmpTplshape[ifitting-1][iElts];
                            m_FittedErrorTplshape[ifitting][iElts] = m_FittedErrorTplshape[ifitting-1][iElts];
                            m_DtmTplshape[ifitting][iElts] = m_DtmTplshape[ifitting-1][iElts];
                            m_MtmTplshape[ifitting][iElts] = m_MtmTplshape[ifitting-1][iElts];
                        }
                        continue;
                    }
                }
                setTplshapeModel(ifitting, false);
                //prepare the Lya width and asym coefficients if the asymfit profile option is met
                //INFO: tpl-shape are always ASYMFIXED for the lyaE profile, as of 2016-01-11
                setLyaProfile(redshift, spectralAxis);
            }

            //generate random amplitudes
            if(m_fittingmethod=="random")
            {
                srand(time(0));
                Float64 randNumFloat = (Float64) rand() / (Float64) (RAND_MAX);

                Float64 coeffAmpEmission = pow(10.0, randNumFloat*3.0-1.0);
                randNumFloat = (Float64) rand() / (Float64) (RAND_MAX);
                Float64 coeffAmpAbsorption = pow(10.0, randNumFloat*1.0-1.0);
                Log.LogInfo( "\nLineModel simulation: coeffAmpEmission = %.2f", coeffAmpEmission);
                Log.LogInfo( "LineModel simulation: coeffAmpAbsorption = %.2f", coeffAmpAbsorption);
                //fit the model amplitudes individually
                for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
                {
                    Float64 meanContinuum = getContinuumMeanUnderElement(iElts);
                    Float64 err = 1e-22;
                    Float64 amax = meanContinuum;
                    if(m_RestRayList[m_Elements[iElts]->m_LineCatalogIndexes[0]].GetType() == CRay::nType_Absorption)
                    {
                        amax = meanContinuum*0.5*coeffAmpAbsorption;
                    }else{
                        amax = meanContinuum*coeffAmpEmission;
                    }
                    randNumFloat = (Float64) rand() / (Float64) (RAND_MAX);
                    Float64 a = randNumFloat*amax;
                    if(a<0.0){
                        a=0.0;
                    }
                    //get the max nominal amplitude
                    Int32 nRays = m_Elements[iElts]->GetSize();
                    Float64 maxNominalAmp = -1.0;
                    for(UInt32 j=0; j<nRays; j++){
                        if(maxNominalAmp<m_Elements[iElts]->GetNominalAmplitude(j))
                        {
                            maxNominalAmp = m_Elements[iElts]->GetNominalAmplitude(j);
                        }
                    }

                    SetElementAmplitude(iElts, a/maxNominalAmp, err);
                }
            }

            //generate amplitudes from existing linemodel fit solution file .csv
            if(m_fittingmethod=="fromfile")
            {
                CModelFittingResult result;
                std::string linemodelFitResultsPath = "/home/aschmitt/gitlab/cpf-redshift/tools/simulation"; //this is the simulation folder absolute path !
                linemodelFitResultsPath.append("/amazed/linecatalogs/linemodelsolve.linemodel_fit_extrema_0.csv");

                Log.LogInfo("\nLineModel fitting from File: loading: %s", linemodelFitResultsPath.c_str());
                result.Load(linemodelFitResultsPath.c_str());
                if(result.GetLineModelSolution().Amplitudes.size()<2)
                {
                    Log.LogError( "\nLineModel fitting from File: unable to load from file: %s", linemodelFitResultsPath.c_str());
                }

                for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
                {
                    Float64 meanContinuum = getContinuumMeanUnderElement(iElts);

                    //get the max nominal amplitude
                    Int32 nRays = m_Elements[iElts]->GetSize();
                    Float64 nominalAmp = -1.0;
                    Float64 err=1.5*1e-20; //not used
                    for(UInt32 j=0; j<nRays; j++){
                        Int32 lineIndex = m_Elements[iElts]->m_LineCatalogIndexes[j];
                        nominalAmp = m_Elements[iElts]->GetNominalAmplitude(j);
                        Float64 a = result.GetLineModelSolution().Amplitudes[lineIndex];

                        if(m_RestRayList[m_Elements[iElts]->m_LineCatalogIndexes[0]].GetType() == CRay::nType_Absorption)
                        {
                            Float64 ampFitted = meanContinuum*a/nominalAmp;
                            SetElementAmplitude(iElts, ampFitted, err);
                        }else{
                            Float64 ampFitted = a/nominalAmp;
                            SetElementAmplitude(iElts, ampFitted, err);
                        }

                    }


                }

            }


            //fit the amplitudes of each element independently
            if(m_fittingmethod=="individual")
            {
                //fit the model amplitudes individually
                for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
                {
                    m_Elements[iElts]->fitAmplitude(spectralAxis, m_spcFluxAxisNoContinuum, m_ContinuumFluxAxis, redshift);

                    //optionnally log some fitting details
                    if(0)
                    {
                        Log.LogDebug("    model: elt #%d individual fit", iElts);
                    Float64 dbg_amp = m_Elements[iElts]->GetElementAmplitude();
                    Log.LogDebug("    model:     fitted elt amp = %f", dbg_amp);
                    Int32 nSubE = m_Elements[iElts]->GetSize();
                    for(Int32 iSubElts=0; iSubElts<nSubE; iSubElts++)
                    {
                        Log.LogDebug("    model:     sub #%d - fitted amp = %f", iSubElts, m_Elements[iElts]->GetFittedAmplitude(iSubElts));
                        Log.LogDebug("    model:     sub #%d - outside range = %d", iSubElts, m_Elements[iElts]->IsOutsideLambdaRange(iSubElts));
                    }
                    Log.LogDebug("    model:     dtm = %e", m_Elements[iElts]->GetSumCross());
                    Log.LogDebug("    model:     dtd = %e", m_Elements[iElts]->GetSumGauss());
                    }
                }
            }


            //fit the amplitude of all elements together (but Emission or Absorption separately) with iterative   solver: lmfit
            if(m_fittingmethod=="lmfit")
            {
                //Log.LogInfo( "Linemodel: Fitting method lmfit, z= %f", m_Redshift);
                SetVelocityEmission(m_velocityEmissionInit);
                SetVelocityAbsorption(m_velocityAbsorptionInit);
                m_Redshift = redshift;
                Float64 bestMerit = DBL_MAX;
                std::shared_ptr<CLmfitController> bestController;// at each run of lmfit we create a controller, we keep the best evaluation
                std::vector<CLmfitController*> controllers = createLmfitControllers(lambdaRange);
                if(controllers.size() == 0){
                    Log.LogError("LMfit : No Controller created");
                }
                for( UInt32 i=0; i<controllers.size(); i++ )
                {
                  CLmfitController *controller = controllers[i];
                  //Log.LogInfo("Continuum Template use : %s", controller->getTemplate()->GetName().c_str());
                  if(!controller->isNoContinuum() && !controller->isContinuumLoaded()){
                    LoadFitContinuumOneTemplate(lambdaRange, *controller->getTemplate());
                  }
                  // adding element base on configuration
                  std::vector<UInt32> validEltsIdx = GetModelValidElementsIndexes();
                  for (UInt32 iElt = 0; iElt < validEltsIdx.size(); iElt++)
                  {
                    // if(controller->isLineTypeVelocityFitted(m_RestRayList[m_Elements[validEltsIdx[iElt]]->m_LineCatalogIndexes[0]].GetType()))
                    // {
                        controller->addElement(validEltsIdx[iElt]);
                        m_Elements[validEltsIdx[iElt]]->fitAmplitude(spectralAxis, m_spcFluxAxisNoContinuum,m_ContinuumFluxAxis , redshift);
                    // }
                  }

                  controller->calculatedIndices();
                  bool NegVal = true;
                  while( NegVal)
                  {
                    Int32 retVal = fitAmplitudesLmfit(m_SpcFluxAxis, controller);
                    //Log.LogInfo( "LineModel Infos: Lmfit retVal = %d", retVal);
                    NegVal = controller->removeNegAmpLine();
                    //Log.LogInfo( "LineModel Infos: Lmfit NegVal = %s", (NegVal?"t":"f") );
                    //Log.LogInfo( "LineModel Infos: Lmfit merit = %.0f ; bestMerit : %.0f ", controller->getMerit(), bestMerit );
                    if(retVal==0 && !NegVal && controller->getMerit()< bestMerit ){

                      bestMerit = controller->getMerit();
                      bestController = std::shared_ptr<CLmfitController>( controller);
                      //Log.LogInfo("Pointeur best Controller setted");
                    }
                    //NegVal = false;//no while
                  }
                }


                //We set the best model found
                if(bestController){
                  if(bestController->isContinuumFitted()){
                    if(bestController->isRedshiftFitted()){
                        m_Redshift = bestController->getRedshift();
                    }
                    const CTemplate* tpl = bestController->getTemplate();
                    ApplyContinuumOnGrid(*tpl, m_Redshift);
                    std::vector<Float64> polyCoeffs;
                    setFitContinuum_tplAmplitude(bestController->getContinuumAmp(), polyCoeffs);
                      //setFitContinuum_tplAmplitude(  bestController->getContinuumAmp());
                  }else if (bestController-> isRedshiftFitted()){
                    setRedshift(bestController-> getRedshift(), bestController->isContinuumLoaded());

                  }
                  std::vector<UInt32> bestFilteredEltsIdx = bestController->getFilteredIdx();

                  for (UInt32 iElt = 0; iElt < bestFilteredEltsIdx.size(); iElt++)
                  {
                      Float64 amp = bestController->getLineAmp(iElt);
                      Float64 ampErr = bestController->getLineAmpErr(iElt);
                      //Log.LogInfo("Lmfit : amp [%d] = %.1f, z = %f", iElt, amp, m_Redshift);
                      SetElementAmplitude(bestFilteredEltsIdx[iElt], amp, ampErr);
                  }
                  if(bestController->isEmissionVelocityFitted()){
                    SetVelocityEmission(bestController->getEmissionVelocity());
                    //Log.LogInfo("Lmfit : VelEmission = %.1f, z = %f", bestController->getEmissionVelocity(), m_Redshift);
                  }
                  if(bestController->isAbsorptionVelocityFitted()){
                    SetVelocityAbsorption(bestController->getAbsorptionVelocity());
                    //Log.LogInfo("Lmfit : VelAbsorption= %.1f, z = %f", bestController->getAbsorptionVelocity(), m_Redshift);
                   }
                   //prepare the elements support
                   for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
                   {
                       m_Elements[iElts]->prepareSupport(spectralAxis, redshift, lambdaRange);
                   }
                   refreshModelInitAllGrid();
                   //modelSolution = GetModelSolution();
                }else{
                  Log.LogError( "LineModel Lmfit: not able to fit values at z %f", m_Redshift);
                  //continue;
                  return DBL_MAX;
                }
            }

            //fit the amplitude of all elements together with linear solver: gsl_multifit_wlinear
            if(m_fittingmethod=="svd")
            {
                std::vector<UInt32> validEltsIdx = GetModelValidElementsIndexes();
                std::vector<Float64> ampsfitted;
                std::vector<Float64> errorsfitted;
                fitAmplitudesLinSolveAndLambdaOffset(validEltsIdx,
                                                     spectralAxis,
                                                     m_spcFluxAxisNoContinuum,
                                                     m_ContinuumFluxAxis,
                                                     ampsfitted,
                                                     errorsfitted,
                                                     m_enableLambdaOffsetsFit);
            }

            //fit the amplitude of all elements AND continuum amplitude together with linear solver: gsl_multifit_wlinear
            if(m_fittingmethod=="svdlc" || m_fittingmethod=="svdlcp2")
            {
                //1. fit only the current continuum
                //prepare continuum on the observed grid
                Log.LogDebug( "    model: fitting svdlc, with continuum-tpl=%s", m_fitContinuum_tplName.c_str());

                //re-interpolate the continuum on the grid
                for( UInt32 i=0; i<m_tplCategoryList.size(); i++ )
                {
                    std::string category = m_tplCategoryList[i];

                    for( UInt32 j=0; j<m_tplCatalog.GetTemplateCount( category ); j++ )
                    {
                        const CTemplate& tpl = m_tplCatalog.GetTemplate( category, j );

                        if(tpl.GetName()==m_fitContinuum_tplName)
                        {
                            ApplyContinuumOnGrid(tpl, m_fitContinuum_tplFitRedshift);
                        }
                    }
                }

                m_fitContinuum_tplFitAmplitude = 1.0;
                std::vector<Float64> polyCoeffs_unused;
                setFitContinuum_tplAmplitude(m_fitContinuum_tplFitAmplitude, polyCoeffs_unused);

                std::vector<UInt32> validEltsIdx = GetModelValidElementsIndexes();
                std::vector<Float64> ampsfitted;
                std::vector<Float64> errorsfitted;
                Float64 chi2_cl = DBL_MAX;
                Int32 fitc_polyOrder = -1;
                if(m_fittingmethod=="svdlcp2")
                {
                    fitc_polyOrder = 2;
                }
                fitAmplitudesLinesAndContinuumLinSolve(validEltsIdx,
                                                       lambdaRange,
                                                       spectralAxis,
                                                       m_SpcFluxAxis,
                                                       m_ContinuumFluxAxis,
                                                       ampsfitted,
                                                       errorsfitted,
                                                       chi2_cl,
                                                       fitc_polyOrder);
                Log.LogDebug( "    model: fitting svdlc done");


                m_fitContinuum_tplFitAmplitude = ampsfitted[validEltsIdx.size()];
                std::vector<Float64> polyCoeffs;
                for(Int32 kpoly=validEltsIdx.size()+1; kpoly<ampsfitted.size();kpoly++)
                {
                    polyCoeffs.push_back(ampsfitted[kpoly]);
                }
                setFitContinuum_tplAmplitude(m_fitContinuum_tplFitAmplitude, polyCoeffs);
                //PrepareContinuum(m_Redshift);
                for(UInt32 i=0; i<modelFluxAxis.GetSamplesCount(); i++)
                {
                    modelFluxAxis[i] = m_ContinuumFluxAxis[i];
                    m_spcFluxAxisNoContinuum[i] = m_SpcFluxAxis[i]-m_ContinuumFluxAxis[i];
                }

                //2. loop on the continuum templates from catalog
                /*
                bool stop_loop_tpl = false;
                Float64 bestChi2 = DBL_MAX;
                Float64 bestAmp = 1.0;
                std::string bestTplName = "-1";
                for( UInt32 i=0; i<m_tplCategoryList.size(); i++ )
                {
                    if(stop_loop_tpl)
                    {
                        break;
                    }
                    std::string category = m_tplCategoryList[i];

                    for( UInt32 j=0; j<m_tplCatalog.GetTemplateCount( category ); j++ )
                    {
                        if(stop_loop_tpl)
                        {
                            break;
                        }
                        const CTemplate& tpl = m_tplCatalog.GetTemplate( category, j );

                        //prepare continuum on the observed grid
                        m_fitContinuum_tplFitDustCoeff = 0;
                        m_fitContinuum_tplFitMeiksinIdx = 0;
                        m_fitContinuum_tplName = tpl.GetName();
                        Log.LogDebug( "    model: fitting svdlc, with continuum-tpl=%s", m_fitContinuum_tplName.c_str());
                        setRedshift(m_Redshift, true);
                        m_fitContinuum_tplFitAmplitude = 1.0;
                        setFitContinuum_tplAmplitude(m_fitContinuum_tplFitAmplitude);

                        std::vector<UInt32> validEltsIdx = GetModelValidElementsIndexes();
                        std::vector<Float64> ampsfitted;
                        std::vector<Float64> errorsfitted;
                        Float64 chi2_cl = DBL_MAX;
                        fitAmplitudesLinesAndContinuumLinSolve(validEltsIdx,
                                                               lambdaRange,
                                                               spectralAxis,
                                                               m_SpcFluxAxis,
                                                               m_ContinuumFluxAxis,
                                                               ampsfitted,
                                                               errorsfitted,
                                                               chi2_cl);
                        Log.LogDebug( "    model: fitting svdlc done");

                        if(chi2_cl<bestChi2)
                        {
                            m_fitContinuum_tplFitAmplitude = ampsfitted[ampsfitted.size()-1];
                            if(1) //not needed here ?
                            {
                                Log.LogDebug( "    model: fitting svdlc, amp found=%e", m_fitContinuum_tplFitAmplitude);
                                Log.LogDebug( "    model: fitting svdlc, chi2 found=%e", chi2_cl);
                                m_fitContinuum_tplFitMerit = -1;
                                //m_fitContinuum_tplFitDustCoeff = -1;
                                //m_fitContinuum_tplFitMeiksinIdx = -1;
                                m_fitContinuum_tplFitDtM = -1;
                                m_fitContinuum_tplFitMtM = -1;
                                setFitContinuum_tplAmplitude(m_fitContinuum_tplFitAmplitude);
                                //PrepareContinuum(m_Redshift);

                                for(UInt32 i=0; i<modelFluxAxis.GetSamplesCount(); i++)
                                {
                                    modelFluxAxis[i] = m_ContinuumFluxAxis[i];
                                    m_spcFluxAxisNoContinuum[i] = m_SpcFluxAxis[i]-m_ContinuumFluxAxis[i];
                                }
                                //refreshModel();
                            }
                            bestChi2 = chi2_cl;
                            bestAmp = m_fitContinuum_tplFitAmplitude;
                            bestTplName = m_fitContinuum_tplName;
                        }
                        //stop_loop_tpl = true; //only process one template
                    }
                }

                m_fitContinuum_tplName = bestTplName;
                setRedshift(m_Redshift, true);
                m_fitContinuum_tplFitAmplitude = bestAmp;
                setFitContinuum_tplAmplitude(m_fitContinuum_tplFitAmplitude);
                //PrepareContinuum(m_Redshift);
                for(UInt32 i=0; i<modelFluxAxis.GetSamplesCount(); i++)
                {
                    modelFluxAxis[i] = m_ContinuumFluxAxis[i];
                    m_spcFluxAxisNoContinuum[i] = m_SpcFluxAxis[i]-m_ContinuumFluxAxis[i];
                }
                //*/

            }

            //fit the amplitudes of each element independently, unless there is overlap
            if(m_fittingmethod=="hybrid")
            {
                fitAmplitudesHybrid(spectralAxis, m_spcFluxAxisNoContinuum, m_ContinuumFluxAxis, redshift);

                //apply a continuum iterative re-estimation with lines removed from the initial spectrum
                Int32 nIt = contreest_iterations;
                Int32 it=0;
                while(it<nIt){
                    applyRules();

                    //*
                    //iterative continuum estimation :: RAW SLOW METHOD
                    refreshModel();
                    Float64 enhanceLines = 0;
                    //*
                    if(nIt>2*it && nIt>3.0 && it<=3){
                        enhanceLines = 2.0-((Float64)it*0.33);
                    }

                    //*/
                    /*
                    if(it==0 && nIt>1){
                        enhanceLines = 1.5;
                    }
                    */
                    EstimateSpectrumContinuum(enhanceLines, lambdaRange);

                    for(UInt32 i=0; i<modelFluxAxis.GetSamplesCount(); i++){
                        modelFluxAxis[i] = m_ContinuumFluxAxis[i];
                        m_spcFluxAxisNoContinuum[i] = m_SpcFluxAxis[i]-m_ContinuumFluxAxis[i];
                    }
                    //*/

                    /*
                //iterative continuum estimation approx: APPROX. METHOD
                std::vector<Int32> validEltsIdx = GetModelValidElementsIndexes();
                std::vector<Int32> refreshIdxs = ReestimateContinuumApprox(validEltsIdx);
                refreshModelAfterContReestimation(refreshIdxs, modelFluxAxis, spcFluxAxisNoContinuum);
                //*/

                    /*
                //iterative continuum estimation approx: FAST METHOD
                std::vector<Int32> validEltsIdx = GetModelValidElementsIndexes();
                std::vector<Int32> highSNRvalidEltsIdx;
                Float64 snrthres = 5.0;
                for( UInt32 i=0; i<validEltsIdx.size(); i++ )
                {
                    Int32 eltIdx = validEltsIdx[i];
                    Bool isSnrHigh = false;
                    Int32 nrays = m_Elements[eltIdx]->GetSize();
                    for(Int32 iray=0; iray<nrays; iray++)
                    {
                        Float64 A = m_Elements[eltIdx]->GetFittedAmplitude(iray);
                        Float64 Sigma = m_Elements[eltIdx]->GetFittedAmplitudeErrorSigma(iray);

                        Float64 snr = A/Sigma;
                        if(snr > snrthres){
                            isSnrHigh = true;
                            break;
                        }
                    }
                    if(isSnrHigh){
                        highSNRvalidEltsIdx.push_back(eltIdx);
                    }
                }
                std::vector<Int32> refreshIdxs = ReestimateContinuumUnderLines(highSNRvalidEltsIdx);
                refreshModelAfterContReestimation(refreshIdxs, modelFluxAxis, spcFluxAxisNoContinuum);
                //*/


                    fitAmplitudesHybrid(spectralAxis, m_spcFluxAxisNoContinuum, m_ContinuumFluxAxis, redshift);
                    it++;
                }
            }

            //set all the amplitudes to 1.0
            if(m_fittingmethod=="ones")
            {
                for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
                {
                    m_Elements[iElts]->SetFittedAmplitude(1.0, 1.0);
                }
            }

            if(m_rigidity=="rules")
            {
                //Apply rules,
                applyRules( enableLogging );

                refreshModel();
                //create spectrum model
                Int32 modelSolutionLevel = Int32(enableLogging);
                modelSolution = GetModelSolution(modelSolutionLevel);
                continuumModelSolution = GetContinuumModelSolution();
                m_tplshapeBestTplName = "None";

                merit = getLeastSquareMerit(lambdaRange);
            }


            //correct lines amplitude with tplshapePrior (tpl-corr): Warning: Rules must all be deactivated
            if(m_rigidity=="tplcorr")
            {
                refreshModel();
                //create spectrum model
                modelSolution = GetModelSolution();
                continuumModelSolution = GetContinuumModelSolution();
                //Log.LogInfo( "LineModel Infos: TPLCORR");
                std::vector<Float64> correctedAmplitudes;
                correctedAmplitudes.resize(modelSolution.Amplitudes.size());
                std::string bestTplName = "";

                m_CatalogTplShape->GetBestFit( modelSolution.Rays, modelSolution.Amplitudes, modelSolution.Errors, correctedAmplitudes, bestTplName);
                for( UInt32 iRestRay=0; iRestRay<m_RestRayList.size(); iRestRay++ )
                {
                    Int32 eIdx = FindElementIndex(iRestRay);
                    Int32 subeIdx = m_Elements[eIdx]->FindElementIndex(iRestRay);
                    Float64 er = m_Elements[eIdx]->GetFittedAmplitudeErrorSigma(subeIdx); //not modifying the fitting error for now
                    Float64 nominalAmp = m_Elements[eIdx]->GetNominalAmplitude(subeIdx);
                    m_Elements[eIdx]->SetFittedAmplitude(correctedAmplitudes[iRestRay]/nominalAmp, er);
                }
                refreshModel();
                modelSolution = GetModelSolution();
                continuumModelSolution = GetContinuumModelSolution();
                m_tplshapeBestTplName = bestTplName;

                merit = getLeastSquareMerit(lambdaRange);
            }

            if(m_rigidity=="tplshape")
            {
                Float64 _merit;
                //if(enableLogging || m_ContinuumComponent == "tplfit")
                if(enableLogging)
                {
                    refreshModel();
                    _merit = getLeastSquareMerit(lambdaRange);
                    //_merit = getLeastSquareMeritFast();
                }else{
                    if(m_tplshapeLeastSquareFast)
                    {
                        _merit = getLeastSquareMeritFast();
                    }else
                    {
                        refreshModel();
                        _merit = getLeastSquareMerit(lambdaRange);
                    }
                }

                if(_merit<meritTplratio[ifitting])
                {
                    meritTplratio[ifitting] = _merit;
                    m_ChisquareTplshape[ifitting] = _merit;
                    m_ScaleMargCorrTplshape[ifitting] = getScaleMargCorrection();
                    m_StrongELPresentTplshape[ifitting] = GetModelStrongEmissionLinePresent();

                    //Saving the model A, errorA, and dtm, mtm, ... (for all tplratios, needed ?)
                    //NB: this is only needed for the index=savedIdxFitted ultimately
                    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
                    {
                        bool savedAmp=false;

                        UInt32 nRays = m_Elements[iElts]->GetSize();
                        for(UInt32 j=0; j<nRays; j++){
                            if(savedAmp)
                            {
                                break;
                            }
                            Float64 amp = m_Elements[iElts]->GetFittedAmplitude(j);
                            if(amp>0.0 && !m_Elements[iElts]->IsOutsideLambdaRange(j))
                            {

                                Float64 amp_error = m_Elements[iElts]->GetFittedAmplitudeErrorSigma(j);
                                Float64 nominal_amp = m_Elements[iElts]->GetNominalAmplitude(j);
                                m_FittedAmpTplshape[ifitting][iElts] = amp/nominal_amp;
                                m_FittedErrorTplshape[ifitting][iElts] = amp_error/nominal_amp;
                                m_DtmTplshape[ifitting][iElts] = m_Elements[iElts]->GetSumCross();
                                m_MtmTplshape[ifitting][iElts] = m_Elements[iElts]->GetSumGauss();

                                savedAmp=true;
                                break;
                            }
                        }
                    }
                }

                if(merit>_merit)
                {
                    merit=_merit;

                    savedIdxContinuumFitted = icontfitting;
                    savedIdxFitted = ifitting;


                    modelSolution = GetModelSolution();
                    continuumModelSolution = GetContinuumModelSolution();
                    m_tplshapeBestTplName = m_CatalogTplShape->GetCatalogName(savedIdxFitted);
                    m_tplshapeBestTplIsmCoeff = m_CatalogTplShape->GetIsmCoeff(savedIdxFitted);
                    m_tplshapeBestTplAmplitude = m_FittedAmpTplshape[savedIdxFitted][0]; //Should be only 1 elt in tpl ratio mode...
                }
            }

            //if(m_rigidity=="tplcorr")
            //{
            //  Float64 tplshapePriorCoeff = m_CatalogTplShape->GetBestFit( modelSolution.Rays, modelSolution.Amplitudes,  );
            //  Log.LogDebug( "Linemodel: tplshapePriorCoeff = %f", tplshapePriorCoeff);
            //  merit = sqrt(merit*merit/2.0-log(tplshapePriorCoeff));
            //}

            if(m_ContinuumComponent == "nocontinuum"){
                reinitModel();
            }
        }
    }

    if(enableLogging)
    {
        if(m_ContinuumComponent == "tplfit")
        {
            if(ncontinuumfitting>1 && m_fittingmethod!="svdlc")
            {
                LoadFitContinuum(lambdaRange, savedIdxContinuumFitted);
            }
            if(m_fittingmethod=="svdlc")
            {
                //
            }
            Log.LogInfo( "    model - Linemodel: fitcontinuum = %d (%s, with ebmv=%.3f), and A=%e",
                         savedIdxContinuumFitted,
                         m_fitContinuum_tplName.c_str(),
                         m_fitContinuum_tplFitDustCoeff,
                         m_fitContinuum_tplFitAmplitude);
        }

        if(m_rigidity=="tplshape")
        {
            //m_CatalogTplShape->SetMultilineNominalAmplitudes( *this, savedIdxFitted );
            bool retSetMultiAmplFast = m_CatalogTplShape->SetMultilineNominalAmplitudesFast( *this, savedIdxFitted );
            if( !retSetMultiAmplFast ){
                Log.LogError( "Linemodel: tplshape, Unable to set Multiline NominalAmplitudes from Tplshape !");
            }


            //Set the velocities from templates: todo auto switch when velfit is ON
            //m_CatalogTplShape->GetCatalogVelocities(savedIdxFitted, m_velocityEmission, m_velocityAbsorption);
            for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
            {
                Log.LogInfo( "    model - Linemodel: tplratio = %d (%s, with ebmv=%.3f), and A=%e", savedIdxFitted, m_tplshapeBestTplName.c_str(), m_tplshapeBestTplIsmCoeff, m_FittedAmpTplshape[savedIdxFitted][iElts]);
                m_Elements[iElts]->SetFittedAmplitude(m_FittedAmpTplshape[savedIdxFitted][iElts], m_FittedErrorTplshape[savedIdxFitted][iElts]);
                m_Elements[iElts]->SetSumCross(m_DtmTplshape[savedIdxFitted][iElts]);
                m_Elements[iElts]->SetSumGauss(m_MtmTplshape[savedIdxFitted][iElts]);
            }
            //Lya
            bool retLyaProfile = m_CatalogTplShape->SetLyaProfile(*this, savedIdxFitted);
            if( !retLyaProfile ){
                Log.LogError( "Linemodel: tplshape, Unable to retrieve Lya Profile from Tplshape !");
            }

            //prepare the Lya width and asym coefficients if the asymfit profile option is met
            setLyaProfile(redshift, spectralAxis);

            refreshModel();

            Int32 modelSolutionLevel = Int32(enableLogging);
            modelSolution = GetModelSolution(modelSolutionLevel);
            continuumModelSolution = GetContinuumModelSolution();
        }
    }
    return merit;
}


std::vector<CLmfitController*> CLineModelElementList::createLmfitControllers( const TFloat64Range& lambdaRange){
  std::vector<CLmfitController*> useLmfitControllers;
  if(m_lmfit_noContinuumTemplate){
    useLmfitControllers.push_back(new CLmfitController(m_lmfit_fitEmissionVelocity,  m_lmfit_fitAbsorptionVelocity));
  }else{
    if(m_lmfit_bestTemplate){
      LoadFitContinuum(lambdaRange, -1);
      for( UInt32 i=0; i<m_tplCategoryList.size(); i++ )
      {
        std::string category = m_tplCategoryList[i];

        for( UInt32 j=0; j<m_tplCatalog.GetTemplateCount( category ) ; j++ )
        {
            const CTemplate& tpl = m_tplCatalog.GetTemplate( category, j );
            if(m_fitContinuum_tplName == tpl.GetName()){
              bool continumLoaded = true;
              useLmfitControllers.push_back(new CLmfitController(tpl, continumLoaded, m_lmfit_fitContinuum,m_lmfit_fitEmissionVelocity,  m_lmfit_fitAbsorptionVelocity));
            }
        }
      }
    }else{
      for( UInt32 i=0; i<m_tplCategoryList.size(); i++ )
      {
        std::string category = m_tplCategoryList[i];

        for( UInt32 j=0; j<m_tplCatalog.GetTemplateCount( category ) ; j++ )
        {
            const CTemplate& tpl = m_tplCatalog.GetTemplate( category, j );
              bool continumLoaded = false;
              useLmfitControllers.push_back(new CLmfitController(tpl, continumLoaded, m_lmfit_fitContinuum,m_lmfit_fitEmissionVelocity,  m_lmfit_fitAbsorptionVelocity));
            }
        }
      }
    }
    return useLmfitControllers;
  }





void CLineModelElementList::SetFittingMethod(std::string fitMethod)
{
    m_fittingmethod = fitMethod;

    if(fitMethod == "lmfit" ) // todo: add the tplfit continuum component option check
    {
      if(m_observeGridContinuumFlux == NULL){
        CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();
        m_chiSquareOperator = new COperatorChiSquare2(m_calibrationPath);
        m_observeGridContinuumFlux = new Float64[modelFluxAxis.GetSamplesCount()]();
        if(m_observeGridContinuumFlux == NULL){
          Log.LogError("unable to allocate m_observeGridContinuumFlux");
          return;
        }
        if(0)
        {
            Log.LogInfo( "Elementlist: fitContinuum_dustfit = %d", m_fitContinuum_dustfit );
            Log.LogInfo( "Elementlist: fitContinuum_igm = %d", m_fitContinuum_igm );
            Log.LogInfo( "Elementlist: fitContinuum_outsidelinesmask = %d", m_fitContinuum_outsidelinesmask );
            Log.LogInfo( "Elementlist: fitContinuum_observedFrame = %d", m_fitContinuum_observedFrame );
        }
      }
   }
}

void CLineModelElementList::SetLeastSquareFastEstimationEnabled(Int32 enabled)
{
    m_tplshapeLeastSquareFast = enabled;
}

void CLineModelElementList::SetAbsLinesLimit(Float64 limit)
{
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        m_Elements[iElts]->SetAbsLinesLimit(limit);
    }
}

/**
 * \brief Init the whole spectrum model with continuum.
 **/
void CLineModelElementList::reinitModel()
{
    CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();

    /*
    //check the continuum flux axis
    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    for(Int32 j=0; j<m_ContinuumFluxAxis.GetSamplesCount(); j++)
    {
        if(isnan(m_ContinuumFluxAxis[j]))
        {
            Log.LogError( "CLineModelElementList::reinitModel: NaN value found for the ContinuumFluxAxis at lambda=%f", spectralAxis[j] );
            throw runtime_error("CLineModelElementList::reinitModel: NaN value found");
            break;
        }
    }
    */

    //init spectrum model with continuum
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        m_Elements[iElts]->initSpectrumModel(modelFluxAxis, m_ContinuumFluxAxis);
    }
}

/**
 * \brief Init the argument elements from the spectrum model with continuum.
 **/
void CLineModelElementList::reinitModelUnderElements(std::vector<UInt32>  filterEltsIdx, Int32 lineIdx )
{
    Int32 iElts;
    CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();
    //init spectrum model with continuum
    for( UInt32 i=0; i<filterEltsIdx.size(); i++ )
    {
        iElts = filterEltsIdx[i];
        m_Elements[iElts]->initSpectrumModel(modelFluxAxis, m_ContinuumFluxAxis, lineIdx);
    }
}

/**
 * \brief Adds a new model to each m_Elements entry.
 * Calls reinitModel.
 * For each entry in m_Elements, addToSpectrumModel using the reinitModel output as arguments.
 **/
void CLineModelElementList::refreshModel(Int32 lineTypeFilter)
{
    reinitModel();    

    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();

    Bool enableMinMaxLog = false;
    if(enableMinMaxLog)
    {
        //check the model min/max values
        //Int32 imin = spectralAxis.GetIndexAtWaveLength(12500.);
        //Int32 imax = spectralAxis.GetIndexAtWaveLength(18500.);
        Int32 imin = 0;
        Int32 imax = modelFluxAxis.GetSamplesCount();
        Float64 fmin = DBL_MAX;
        Float64 fmax = -DBL_MAX;
        for(Int32 j=imin; j<imax; j++)
        {
            if(modelFluxAxis[j]<fmin)
            {
                fmin = modelFluxAxis[j];
            }
            if(modelFluxAxis[j]>fmax)
            {
                fmax = modelFluxAxis[j];
            }
        }
        Log.LogDebug( "CLineModelElementList::refreshModel AFTER REINIT: model min=%e and model max=%e", fmin, fmax );
    }

    /*
    //check the model reinited for nan values
    Int32 imin = spectralAxis.GetIndexAtWaveLength(12500.);
    Int32 imax = spectralAxis.GetIndexAtWaveLength(18500.);
    //Int32 imin = 0;
    //Int32 imax = modelFluxAxis.GetSamplesCount();

    for(Int32 j=imin; j<imax; j++)
    {
        if(isnan(modelFluxAxis[j]))
        {
            Log.LogError( "CLineModelElementList::refreshModel: NaN value found for the reinited model spectrum at lambda=%f", spectralAxis[j] );
            throw runtime_error("CLineModelElementList::refreshModel: NaN value found");
            break;
        }
    }
    //*/

    /*
    CSpectrumFluxAxis contFluxAxisWithAmpOffset = CSpectrumFluxAxis(m_ContinuumFluxAxis.GetSamplesCount());
    for( Int32 i=0; i<m_ContinuumFluxAxis.GetSamplesCount(); i++ )
    {
        contFluxAxisWithAmpOffset[i] = m_ContinuumFluxAxis[i];
    }
    //*/
    if(m_enableAmplitudeOffsets)
    {
        //add amplitude offsets
        addToSpectrumAmplitudeOffset(modelFluxAxis);
        //addToSpectrumAmplitudeOffset(contFluxAxisWithAmpOffset);
    }

    //create spectrum model
    UInt32 nElements = m_Elements.size();
    for( UInt32 iElts=0; iElts<nElements; iElts++ )
    {
        Int32 lineType = m_Elements[iElts]->m_Rays[0].GetType();
        if(lineTypeFilter==-1 || lineTypeFilter==lineType)
        {
            m_Elements[iElts]->addToSpectrumModel(spectralAxis, modelFluxAxis, m_ContinuumFluxAxis, m_Redshift);
        }

        /*
        //check the model for nan values
        Int32 imin = spectralAxis.GetIndexAtWaveLength(12500.);
        Int32 imax = spectralAxis.GetIndexAtWaveLength(18500.);
        //Int32 imin = 0;
        //Int32 imax = modelFluxAxis.GetSamplesCount();

        for(Int32 j=imin; j<imax; j++)
        {
            if(isnan(modelFluxAxis[j]))
            {
                Log.LogError( "CLineModelElementList::refreshModel: NaN value found for the model spectrum at lambda=%f", spectralAxis[j] );
                Log.LogError( "CLineModelElementList::refreshModel: NaN value found for the model spectrum when added line %d (nElts=%d)", iElts, m_Elements.size() );
                throw runtime_error("CLineModelElementList::refreshModel: NaN value found");
                break;
            }
        }
        //*/
    }

    if(enableMinMaxLog)
    {
        //check the model min/max values
        Int32 imin = spectralAxis.GetIndexAtWaveLength(12500.);
        Int32 imax = spectralAxis.GetIndexAtWaveLength(18500.);
        //Int32 imin = 0;
        //Int32 imax = modelFluxAxis.GetSamplesCount();
        Float64 fmin = DBL_MAX;
        Float64 fmax = -DBL_MAX;
        for(Int32 j=imin; j<imax; j++)
        {
            if(modelFluxAxis[j]<fmin)
            {
                fmin = modelFluxAxis[j];
            }
            if(modelFluxAxis[j]>fmax)
            {
                fmax = modelFluxAxis[j];
            }
        }
        Log.LogDebug( "CLineModelElementList::refreshModel AFTER addToSpectrumModel: model min=%e and model max=%e", fmin, fmax );
    }

}


Bool CLineModelElementList::addToSpectrumAmplitudeOffset( CSpectrumFluxAxis& modelfluxAxis )
{
    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();

    Log.LogInfo( "Elementlist: Adding n=%d ampOffsets", m_ampOffsetsIdxStart.size());
    for( UInt32 i=0; i<m_ampOffsetsIdxStart.size(); i++ )
    {
        for( UInt32 k=m_ampOffsetsIdxStart[i]; k<=m_ampOffsetsIdxStop[i]; k++ )
        {
            modelfluxAxis[k] += m_ampOffsetsX0[i] + m_ampOffsetsX1[i]*spectralAxis[k]+ m_ampOffsetsX2[i]*spectralAxis[k]*spectralAxis[k];
        }

    }
    return true;
}

Int32 CLineModelElementList::prepareAmplitudeOffset(const CSpectrumFluxAxis& spcFlux)
{
    m_ampOffsetsX0.clear();
    m_ampOffsetsX1.clear();
    m_ampOffsetsX2.clear();
    m_ampOffsetsIdxStart.clear();
    m_ampOffsetsIdxStop.clear();

    std::vector<UInt32> validEltsIdx = GetModelValidElementsIndexes();
    if(validEltsIdx.size()<1)
    {
        return -1;
    }
    std::vector<UInt32> supportIdxes = getSupportIndexes( validEltsIdx );
    if(supportIdxes.size()<1)
    {
        return -1;
    }

    Int32 idxPrevious = supportIdxes[0];
    m_ampOffsetsIdxStart.push_back(supportIdxes[0]);
    for( UInt32 i=1; i<supportIdxes.size(); i++ )
    {
        UInt32 idxCurrent = supportIdxes[i];
        if(idxCurrent>idxPrevious+1)
        {
            m_ampOffsetsIdxStop.push_back(idxPrevious);
            m_ampOffsetsIdxStart.push_back(idxCurrent);
        }
        idxPrevious = idxCurrent;
    }
    m_ampOffsetsIdxStop.push_back(supportIdxes[supportIdxes.size()-1]);

    /*
    //estimate mean fluxNoCOnt in each support
    for( Int32 i=0; i<m_ampOffsetsIdxStart.size(); i++ )
    {
        Float64 sum = 0.0;
        Int32 count = 0;

        for( Int32 k=m_ampOffsetsIdxStart[i]; k<=m_ampOffsetsIdxStop[i]; k++ )
        {
            sum +=spcFlux[k];
            count +=1;
        }
        Float64 mean = 0.0;
        if(count>0)
        {
            mean = sum/(Float64)count;
        }

        m_ampOffsetsA.push_back(mean);
    }
    //*/
    for( UInt32 i=0; i<m_ampOffsetsIdxStart.size(); i++ )
    {
        m_ampOffsetsX0.push_back(0.0);
        m_ampOffsetsX1.push_back(0.0);
        m_ampOffsetsX2.push_back(0.0);
    }

    if(1)
    {
        for( UInt32 i=0; i<m_ampOffsetsIdxStart.size(); i++ )
        {
            Log.LogInfo( "Elementlist: i=%d, m_ampOffsetsIdxStart: %d, m_ampOffsetsIdxStop: %d", i, m_ampOffsetsIdxStart[i], m_ampOffsetsIdxStop[i] );
        }
    }
    return 0;
}

/**
 * \brief refreshing all the grid and Adds a new model to each m_Elements entry .
 * Calls iterate on model flux to set it equal to continuum.
 * For each entry in m_Elements, addToSpectrumModel using the reinitModel output as arguments.
 **/
void CLineModelElementList::refreshModelInitAllGrid()
{

    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();

    for(UInt32 k= 0;k<modelFluxAxis.GetSamplesCount();k++){
      modelFluxAxis[k] = m_ContinuumFluxAxis[k];
    }

    //create spectrum model
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        m_Elements[iElts]->addToSpectrumModel(spectralAxis, modelFluxAxis, m_ContinuumFluxAxis, m_Redshift);
    }
}

/**
 * \brief Adds a new model to each m_Elements entry specified on the argument.
 * Works as refreshModel.
 **/
void CLineModelElementList::refreshModelUnderElements(std::vector<UInt32> filterEltsIdx, Int32 lineIdx )
{
    reinitModelUnderElements(filterEltsIdx, lineIdx);
    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();
    //create spectrum model
    Int32 iElts;
    for( UInt32 i=0; i<filterEltsIdx.size(); i++ )
    {
        iElts = filterEltsIdx[i];
        m_Elements[iElts]->addToSpectrumModel(spectralAxis, modelFluxAxis, m_ContinuumFluxAxis, m_Redshift, lineIdx);
    }
}

void CLineModelElementList::refreshModelDerivVelEmissionUnderElements(std::vector<UInt32> filterEltsIdx)
{
    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    std::vector<UInt32> supportIdxes = getSupportIndexes( filterEltsIdx );

    for( UInt32 i=0; i<supportIdxes.size(); i++ )
    {
        m_SpcFluxAxisModelDerivVelEmi[supportIdxes[i]]=0.0;
    }

    //create spectrum model partial derivate vs sigma
    Int32 iElts;
    for( UInt32 i=0; i<filterEltsIdx.size(); i++ )
    {
        iElts = filterEltsIdx[i];
        m_Elements[iElts]->addToSpectrumModelDerivVel(spectralAxis, m_SpcFluxAxisModelDerivVelEmi, m_ContinuumFluxAxis,  m_Redshift, true);
    }
}

void CLineModelElementList::refreshModelDerivVelAbsorptionUnderElements(std::vector<UInt32> filterEltsIdx)
{
    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    std::vector<UInt32> supportIdxes = getSupportIndexes( filterEltsIdx );

    for( UInt32 i=0; i<supportIdxes.size(); i++ )
    {
        m_SpcFluxAxisModelDerivVelAbs[supportIdxes[i]]=0.0;
    }

    //create spectrum model partial derivate vs sigma
    Int32 iElts;
    for( UInt32 i=0; i<filterEltsIdx.size(); i++ )
    {
        iElts = filterEltsIdx[i];
        m_Elements[iElts]->addToSpectrumModelDerivVel(spectralAxis, m_SpcFluxAxisModelDerivVelAbs, m_ContinuumFluxAxis,  m_Redshift, false);
    }
}


void CLineModelElementList::refreshModelDerivVelUnderElements(std::vector<UInt32> filterEltsIdx)
{
    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    std::vector<UInt32> supportIdxes = getSupportIndexes( filterEltsIdx );

    for( UInt32 i=0; i<supportIdxes.size(); i++ )
    {
        m_SpcFluxAxisModelDerivVelEmi[supportIdxes[i]]=0.0;
        m_SpcFluxAxisModelDerivVelAbs[supportIdxes[i]]=0.0;
    }

    //create spectrum model partial derivate vs sigma
    Int32 iElts;
    for( UInt32 i=0; i<filterEltsIdx.size(); i++ )
    {
        iElts = filterEltsIdx[i];
        m_Elements[iElts]->addToSpectrumModelDerivVel(spectralAxis, m_SpcFluxAxisModelDerivVelEmi, m_ContinuumFluxAxis,  m_Redshift, true);
        m_Elements[iElts]->addToSpectrumModelDerivVel(spectralAxis, m_SpcFluxAxisModelDerivVelAbs, m_ContinuumFluxAxis,  m_Redshift, false);
    }
}

void CLineModelElementList::setModelSpcObservedOnSupportZeroOutside(  const TFloat64Range& lambdaRange )
{
    m_Redshift = 0.0;

    //initialize the model spectrum
    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();

    //prepare the elements
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        m_Elements[iElts]->prepareSupport(spectralAxis, m_Redshift, lambdaRange);
    }

    std::vector<UInt32> validEltsIdx = GetModelValidElementsIndexes();
    std::vector<UInt32> supportIdxes = getSupportIndexes( validEltsIdx );
    for( UInt32 i=0; i<spectralAxis.GetSamplesCount(); i++ )
    {
        modelFluxAxis[i]=0.0;
    }
    for( UInt32 i=0; i<supportIdxes.size(); i++ )
    {
        modelFluxAxis[supportIdxes[i]] = m_SpcFluxAxis[supportIdxes[i]];
    }
}


/**
 * \brief Creates and returns a Mask with 0 in the lines support, 1 under the lines
 **/
CMask CLineModelElementList::getOutsideLinesMask()
{
    CMask _mask;
    //initialize the model spectrum
    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    _mask.SetSize(spectralAxis.GetSamplesCount());

    std::vector<UInt32> validEltsIdx = GetModelValidElementsIndexes();
    std::vector<UInt32> supportIdxes = getSupportIndexes( validEltsIdx );
    for( UInt32 i=0; i<spectralAxis.GetSamplesCount(); i++ )
    {
        _mask[i]=1;
    }
    for( UInt32 i=0; i<supportIdxes.size(); i++ )
    {
        _mask[supportIdxes[i]] = 0;
    }
    return _mask;
}


/**
 * \brief Estimates the STD outside the lines for the observed-model spectrum
 * NB: supposes the spectrum whithout continuum has a null mean value
 * input: which = 1: uses the spectrum flux continuum subtracted to compute STD
 * input: which = 2: uses the spectrum error to compute STD
 **/
Float64 CLineModelElementList::getOutsideLinesSTD( Int32 which, TFloat64Range lambdarange)
{
    if(which!=1 && which!=2)
    {
        Log.LogError( "    model: getOutsideLinesSTD - Failed to parse input argument, which" );
        return -1;
    }
    Float64 estimated_std = -1.;
    CMask _mask = getOutsideLinesMask();

    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    Float64 sum2 = 0.0;
    Int32 nsum=0;
    Float64 imin = spectralAxis.GetIndexAtWaveLength(lambdarange.GetBegin());
    Float64 imax = spectralAxis.GetIndexAtWaveLength(lambdarange.GetEnd());
    for( UInt32 i=imin; i<imax; i++ )
    {
        if(_mask[i])
        {
            if(which==1)
            {
                sum2 += m_spcFluxAxisNoContinuum[i]*m_spcFluxAxisNoContinuum[i];
                nsum++;
            }else if(which==2)
            {
                sum2 += m_ErrorNoContinuum[i]*m_ErrorNoContinuum[i];
                nsum++;
            }
        }
    }
    if(nsum>0)
    {
        estimated_std = sqrt(sum2/nsum);
    }else{
        return -1.;
    }

    return estimated_std;
}

/**
 * \brief Tries to fit subelements considering their overlap.
 * For each entry in GetModelValidElementsIndexes:
 *   If subelement in the entry already fitted, go for the next entry.
 *   getOverlappingElements for the fitted subelements.
 *   If the overlap is smaller than 2, call fitAmplitude on the entry.
 *   If the overlap is greater than or equal to 2:
 *     Call fitAmplitudeLinSolve with the subelements as argument.
 *     Store all non-negative fits.
 *     Set to 0.0 all negative fits.
 *     If the size of non-negative fits is 1, call the entry's fitAmplitude.
 *     If the size of non-negative fits is not 1:
 *       If the size of non-negative fits is greater than 1:
 *         Call fitAmplitudesLinSolve with the indexes of the non-negative subelements.
 *         If the above call return is different than 1:
 *           For each non-negative subelement, if the amplitude fitted is greater than 0, call fitAmplitude on its entry. Else, SetElementAmplitude to 0.
 *   Update the index of already-fitted subelements.
 **/
Int32 CLineModelElementList::fitAmplitudesHybrid(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& spcFluxAxisNoContinuum, const CSpectrumFluxAxis &continuumfluxAxis, Float64 redshift)
{
  bool verbose=false;
  std::vector<UInt32> validEltsIdx = GetModelValidElementsIndexes();
  std::vector<UInt32> indexesFitted;
  for( UInt32 iValidElts=0; iValidElts<validEltsIdx.size(); iValidElts++ )
  {
      UInt32 iElts = validEltsIdx[iValidElts];
      //skip if already fitted
      bool alreadyfitted=false;
      for(UInt32 i=0; i<indexesFitted.size(); i++)
      {
          if(iElts == indexesFitted[i])
          {
              alreadyfitted=true;
              break;
          }
      }
      if(alreadyfitted)
      {
          continue;
      }
      //do the fit on the ovelapping elements
      Float64 overlapThres = 0.15; //15% seemed necessary for Ha/SII complex when lines are very wide (either because of PSF or source size)
      std::vector<UInt32> overlappingInds = getOverlappingElements(iElts, indexesFitted, overlapThres);

      //setting the fitting group info
      for(UInt32 ifit=0; ifit<overlappingInds.size(); ifit++)
      {
          std::string fitGroupTag = boost::str(boost::format("hy%d") % iValidElts);
          m_Elements[overlappingInds[ifit]]->m_fittingGroupInfo=fitGroupTag;
      }


      if(verbose)
      {
          //Log.LogDebug( "Redshift: %f", m_Redshift);
          Log.LogDebug( "    model: hybrid fit: #%d - N overlapping=%d", iValidElts, overlappingInds.size());
          for(UInt32 ifit=0; ifit<overlappingInds.size(); ifit++)
          {
              Log.LogDebug( "    model: hybrid fit:     overlapping #%d - eltIdx=%d", ifit, overlappingInds[ifit]);
          }
      }
      if(!m_enableAmplitudeOffsets && overlappingInds.size()<2)
      {
          if(verbose)
          {
            Log.LogDebug( "    model: hybrid fit:     Individual fit");
          }
          m_Elements[iElts]->fitAmplitudeAndLambdaOffset(spectralAxis, spcFluxAxisNoContinuum, continuumfluxAxis, redshift, -1, m_enableLambdaOffsetsFit, m_LambdaOffsetStep, m_LambdaOffsetMin, m_LambdaOffsetMax);
      }
      else
      {
          if(verbose)
          {
            Log.LogDebug( "    model: hybrid fit:     SVD fit");
          }
          //fit individually: mainly for mtm and dtm estimation
          for(UInt32 ifit=0; ifit<overlappingInds.size(); ifit++)
          {
              m_Elements[overlappingInds[ifit]]->fitAmplitudeAndLambdaOffset(spectralAxis, spcFluxAxisNoContinuum, continuumfluxAxis, redshift, -1, m_enableLambdaOffsetsFit, m_LambdaOffsetStep, m_LambdaOffsetMin, m_LambdaOffsetMax);
          }
          std::vector<Float64> ampsfitted;
          std::vector<Float64> errorsfitted;
          Int32 retVal = -1;
          retVal = fitAmplitudesLinSolveAndLambdaOffset(overlappingInds, spectralAxis, spcFluxAxisNoContinuum, continuumfluxAxis, ampsfitted, errorsfitted, m_enableLambdaOffsetsFit);
          // if all the amplitudes fitted don't have the same sign, do it separately
          std::vector<UInt32> overlappingIndsSameSign;
          if(retVal!=1 && ampsfitted.size()>0)
          {
              for(UInt32 ifit=0; ifit<overlappingInds.size(); ifit++)
              {
                  if(ampsfitted[ifit]>0)
                  {
                      overlappingIndsSameSign.push_back(overlappingInds[ifit]);
                      //m_Elements[overlappingInds[ifit]]->fitAmplitude(spectralAxis, spcFluxAxisNoContinuum, redshift);
                  }
                  else
                  {
                      SetElementAmplitude(overlappingInds[ifit], 0.0, errorsfitted[ifit]);
                  }
              }
              //fit the rest of the overlapping elements (same sign) together
              if(!m_enableAmplitudeOffsets && overlappingIndsSameSign.size()==1)
              {
                  m_Elements[overlappingIndsSameSign[0]]->fitAmplitudeAndLambdaOffset(spectralAxis, spcFluxAxisNoContinuum, continuumfluxAxis, redshift, -1, m_enableLambdaOffsetsFit, m_LambdaOffsetStep, m_LambdaOffsetMin, m_LambdaOffsetMax);
              }else if(overlappingIndsSameSign.size()>0){
                  Int32 retVal2=-1;
                  retVal2 = fitAmplitudesLinSolveAndLambdaOffset(overlappingIndsSameSign, spectralAxis, spcFluxAxisNoContinuum, continuumfluxAxis, ampsfitted, errorsfitted, m_enableLambdaOffsetsFit);

                  if(retVal2!=1){
                      for(UInt32 ifit=0; ifit<overlappingIndsSameSign.size(); ifit++)
                      {
                          if(ampsfitted[ifit]>0){
                              m_Elements[overlappingIndsSameSign[ifit]]->fitAmplitudeAndLambdaOffset(spectralAxis, spcFluxAxisNoContinuum, continuumfluxAxis, redshift, -1, m_enableLambdaOffsetsFit, m_LambdaOffsetStep, m_LambdaOffsetMin, m_LambdaOffsetMax);
                          }else{
                              SetElementAmplitude(overlappingIndsSameSign[ifit], 0.0, errorsfitted[ifit]);
                          }
                      }
                  }
              }

          }
      }

      //update the already fitted list
      for(UInt32 i=0; i<overlappingInds.size(); i++)
      {
          indexesFitted.push_back(overlappingInds[i]);
      }

  }

  improveBalmerFit();

  return 0;
}

/**
 * @brief CLineModelElementList::estimateMeanSqFluxAndGradient
 * @param varPack: the variables of the model being fitted
 * @param normFactor
 * @param filteredEltsIdx: index of the elements included in the fit, also sets the size of varPack (n=nElts+1)
 * @param xInds: indexes of the samples where the data of the model is fitted
 * @param lineType: E or A
 * @param fluxdata: data to be fitted (already reshaped, no need to use xInds for this vector)
 * @param msqBuffer: buffer for fast computing of the meansquare
 * @param f: output meansquare residual
 * @param g: output meansquare gradient residual
 * @return
 */
Int32 CLineModelElementList::estimateMeanSqFluxAndGradient(const Float64* varPack,
                                                           const Float64 normFactor,
                                                           std::vector<UInt32> filteredEltsIdx,
                                                           std::vector<UInt32> xInds,
                                                           Int32 lineType,
                                                           Float64* fluxdata,
                                                           Float64* msqBuffer,
                                                           Float64& f,
                                                           Float64* g)
{
    // update the linemodel amplitudes
    for (UInt32 iElt = 0; iElt < filteredEltsIdx.size(); iElt++)
    {
        Float64 amp = varPack[iElt]/normFactor;
        SetElementAmplitude(filteredEltsIdx[iElt], amp, 0.0);
    }
    // update the linemodel velocity/linewidth
    Int32 idxVelocity = filteredEltsIdx.size();
    Float64 velocity = varPack[idxVelocity];
    if(lineType==CRay::nType_Emission)
    {
        SetVelocityEmission(velocity);
    }else
    {
        SetVelocityAbsorption(velocity);
    }

    Int32 nsamples = xInds.size();
    // retrieve the model
    refreshModelUnderElements(filteredEltsIdx);
    f = 0.0;
    for (Int32 i = 0; i < nsamples; i++)
    {
        Float64 Yi = getModelFluxVal(xInds[i])*normFactor;
        msqBuffer[i] = Yi - fluxdata[i];
        f += msqBuffer[i]*msqBuffer[i];
    }

    refreshModelDerivVelUnderElements(filteredEltsIdx);
    for (UInt32 iElt = 0; iElt < filteredEltsIdx.size(); iElt++)
    {
        g[iElt]=0.0;
    }
    for (Int32 i = 0; i < nsamples; i++)
    {
        for (UInt32 iElt = 0; iElt < filteredEltsIdx.size(); iElt++)
        {
            Float64 dm = getModelFluxDerivEltVal(filteredEltsIdx[iElt], xInds[i]);
            Float64 grad = 2*dm*msqBuffer[i];
            g[iElt] += grad;
        }
        //*
        Int32 iElt = filteredEltsIdx.size();
        //getModelFluxDerivVelVal change for LMFit, emssion and aborption have there own deriv sigma value
        Float64 dm = getModelFluxDerivVelVal(xInds[i])*normFactor;
        Float64 grad = 2*dm*msqBuffer[i];
        g[iElt] += grad;
        //*/
    }

    return 0;
}


//temporary stuff, dev of the lmfit method: int

void print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
  printf ("iter: %3lu x = % 15.8f % 15.8f % 15.8f "
          "|f(x)| = %g\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2),
          gsl_blas_dnrm2 (s->f));
}


Int32 CLineModelElementList::fitAmplitudesLmfit( const CSpectrumFluxAxis& fluxAxis, CLmfitController* controller)
//std::vector<Float64>& ampsfitted, std::vector<Float64>& ampsErrfitted, Float64& velocityFitted, Float64& continuumAmpFitted, Float64& merit, Int32 lineType
{
    //http://www.gnu.org/software/gsl/manual/html_node/Example-programs-for-Nonlinear-Least_002dSquares-Fitting.html
    Bool verbose = true;

    if(verbose)
    {
        //Log.LogInfo( "fitAmplitudesLmfit");
    }
    std::vector<UInt32> filteredEltsIdx= controller->getFilteredIdx();
    UInt32 nddl = filteredEltsIdx.size();
    if(nddl<1){
        Log.LogError("Lmfit :No ray amplitude to fit ");
        return -1;
    }

    std::vector<UInt32> xInds;
    if(! controller->isContinuumFitted()){
        xInds = getSupportIndexes( filteredEltsIdx );
    }else{
        xInds=  std::vector<UInt32>(fluxAxis.GetSamplesCount());
        //boost::push_back(xInds, fluxAxis.GetSamplesCount());
        std::iota(std::begin(xInds), std::end(xInds), 0);
    }

    const Float64* flux = fluxAxis.GetSamples();
    Float64 amplitudeContinuumInit = m_fitContinuum_tplFitAmplitude;

    //create gsl solver object
    const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
    gsl_multifit_fdfsolver *s;
    int status, info;
    size_t i;
    size_t n; //n samples on the support, /* number of data points to fit */
    size_t p = controller->getNumberParameters(); //DOF = n amplitudes to fit (1 for each element) + 1 (EL velocity) +1 Absorption VEl +1 Amplitude Template

    n = xInds.size();
    if(n<nddl){
        controller->resizeAmpsLine();
        Log.LogError( "LineModel Infos: LMfit not enough smaples on support");
        for (Int32 iddl = 0; iddl < nddl; iddl++)
        {
            controller->setAmpLine(iddl, 0.0,0.0);
        }
        return -1;
    }

    //set the norm factor
    Int32 idx = 0;
    Float64 maxabsval = DBL_MIN;
    Float64 normFactor;
    // if(controller->needCalculateNormFactor()){
      for (i = 0; i < n; i++)
      {
          idx = xInds[i];

          if(maxabsval<std::abs(flux[idx]))
          {
              maxabsval=std::abs(flux[idx]);
          }
      }
      normFactor = 1.0/maxabsval;
      controller->setNormFactor(normFactor);
    // }else{
    //   // the norm factor is the same as the last call
    //   normFactor=controller->getNormFactor();
    // }

    if(verbose)
    {
        //Log.LogInfo("normFactor = '%.3e'", normFactor);
    }

    gsl_matrix *J = gsl_matrix_alloc(n, p);
    gsl_matrix *covar = gsl_matrix_alloc (p, p);
    double y[n], weights[n];
    struct lmfitdata d = {n,y,this, xInds, m_observeGridContinuumFlux, controller};
    gsl_multifit_function_fdf f;

    Float64* x_init = (Float64*) calloc( p, sizeof( Float64 ) );
    //initialize lmfit with previously estimated individual/hybrid fit method
    Float64 bestAmpLine = -1;
    for(Int32 kp=0; kp<nddl; kp++)
    {
        Float64 amp =  m_Elements[filteredEltsIdx[kp]]->GetElementAmplitude();
        if( amp>bestAmpLine){
            bestAmpLine = amp;
        }
    }
    controller->setNormAmpLine(1.0);
//    if(bestAmpLine>0.){
//        controller->setNormAmpLine(1/bestAmpLine);//setNormAmpLine(controller->lineAmp_LmToModel(bestAmpLine));
//        if(verbose){
//            Log.LogInfo("Lmfit : normAmpLine : %f", controller->getNormAmpLine());
//        }
//    }
    for(Int32 kp=0; kp<nddl; kp++)
    {
        Float64 ampInitGuess = m_Elements[filteredEltsIdx[kp]]->GetElementAmplitude();//std::max(m_Elements[filteredEltsIdx[kp]]->GetElementAmplitude() *1.5, bestAmpLine*0.001) ;
        if(ampInitGuess <0.){
            Log.LogInfo("Amp negative %d , set to 0",ampInitGuess );
            ampInitGuess = 0.;
        }
        //Float64 ampInitGuess = 0;

        x_init[kp] = controller->lineAmp_ModelToLm(ampInitGuess);//*normFactor;
        if(verbose)
        {
            Log.LogInfo( "LineModel LMfit: set init guess amp [%d] Model = %f lmfit: %f", filteredEltsIdx[kp], ampInitGuess, x_init[kp]);
        }
        x_init[kp] = ampInitGuess*normFactor;
    }

    if(controller->isEmissionVelocityFitted())
    {
        Float64 normEmiFactor = 1.0;//bestAmpLine/GetVelocityEmission();//1.0;//
        controller->setNormEmiFactor(normEmiFactor);
        x_init[controller->getIndEmissionVel()] = controller->emiVel_ModelToLm(GetVelocityEmission());//bestAmpLine;
        if(verbose){
            Log.LogInfo( "LineModel LMfit: normEmiFactor = %f", normEmiFactor);
             Log.LogInfo( "LineModel LMfit: set init guess Emission velocity = %f ", GetVelocityEmission());
        }
    }
    if(controller->isAbsorptionVelocityFitted())
    {
        Float64 normAbsFactor =1.0;// bestAmpLine/GetVelocityAbsorption();//1.0f;//

        controller->setNormAbsFactor(normAbsFactor);
        x_init[controller->getIndAbsorptionVel()] = controller->absVel_ModelToLm(GetVelocityAbsorption());//bestAmpLine;
        if(verbose){
             Log.LogInfo( "LineModel LMfit: normAbsFactor = %f", normAbsFactor);
             Log.LogInfo( "LineModel LMfit: set init guess Absorption velocity = %f ", GetVelocityAbsorption());
        }

    }
    if(controller->isContinuumFitted()){
      x_init[controller->getIndContinuumAmp()]= controller->continuumAmp_ModelToLm(getFitContinuum_tplAmplitude());//*normFactor;
      if(verbose){
           Log.LogInfo( "LineModel LMfit: set init guess continuum amp = %f ", getFitContinuum_tplAmplitude());
      }
    }

    if(controller->isRedshiftFitted()){
      x_init[controller->getIndRedshift()]= m_Redshift;
      if(verbose){
           Log.LogInfo( "LineModel LMfit: set init guess redshift= %f ", m_Redshift);
      }
    }


    gsl_vector_view x = gsl_vector_view_array (x_init, p);
    gsl_vector_view w = gsl_vector_view_array(weights, n);
    //const gsl_rng_type * type;
    //gsl_rng * r;
    gsl_vector *res_f;
    double chi, chi0;

    const double xtol = 1e-8;
    const double gtol = 1e-8;
    const double ftol =1-8;
    Int32 maxIterations = 250;

    //gsl_rng_env_setup();

    //type = gsl_rng_default;
    //r = gsl_rng_alloc (type);

    f.f = &lmfit_f;
    f.df = &lmfit_df;
    f.n = n;
    f.p = p;
    f.params = &d;

    // This is the data to be fitted
    Float64 ei;
    for (i = 0; i < n; i++)
    {
        idx = xInds[i];
        ei = m_ErrorNoContinuum[idx]*normFactor;
        weights[i] = 1.0 / (ei * ei);
        y[i] = flux[idx]*normFactor;
    }

    s = gsl_multifit_fdfsolver_alloc (T, n, p);

    /* initialize solver with starting point and weights */
    gsl_multifit_fdfsolver_wset (s, &f, &x.vector, &w.vector);

    /* compute initial residual norm */
    res_f = gsl_multifit_fdfsolver_residual(s);
    chi0 = gsl_blas_dnrm2(res_f);

    /* solve the system with a maximum of maxIterations iterations */
    status = gsl_multifit_fdfsolver_driver(s, maxIterations, xtol, gtol, ftol, &info);

    gsl_multifit_fdfsolver_jac(s, J);
    gsl_multifit_covar (J, 0.0, covar);

    /* compute final residual norm */
    chi = gsl_blas_dnrm2(res_f);

    if(verbose)
    {
#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

        Log.LogInfo( "summary from method '%s'",
                gsl_multifit_fdfsolver_name(s));
        Log.LogInfo( "number of iterations: %zu",
                gsl_multifit_fdfsolver_niter(s));
        Log.LogInfo("function evaluations: %zu", f.nevalf);
        Log.LogInfo("Jacobian evaluations: %zu", f.nevaldf);
        Log.LogInfo( "reason for stopping: %s",
                (info == 1) ? "small step size" :(info==2)? "small gradient": "small change in f");
        Log.LogInfo( "initial |f(x)| = %g", chi0);
        Log.LogInfo("final   |f(x)| = %g", chi);

        {
            double dof = n - p;
            //double c = GSL_MAX_DBL(1, chi / sqrt(dof));

            Log.LogInfo( "chisq/dof = %g",  pow(chi, 2.0) / dof);

            for(Int32 k=0; k<p; k++)
            {
                if(FIT(k)<1e-3)
                {
                    //Log.LogInfo ("A %d     = %.3e +/- %.8f", k, FIT(k), c*ERR(k));
                }else{
                    //Log.LogInfo ( "A %d     = %.5f +/- %.8f", k, FIT(k), c*ERR(k));
                }

            }
        }
        Log.LogInfo ( "status = %s (%d)", gsl_strerror (status), status);
    }

    //=======================Result treament ====================
    Int32 ampPos = 1;
    controller->resizeAmpsLine();
    Float64 dof = n - p;
    Float64 c = GSL_MAX_DBL(1, chi / sqrt(dof));
    for (Int32 iddl = 0; iddl < nddl; iddl++)
    {
        Float64 a = controller->lineAmp_LmToModel(gsl_vector_get(s->x,iddl));///normFactor;
        Float64 sigma = controller->lineAmp_LmToModel(c*sqrt(gsl_matrix_get(covar,iddl,iddl)));///normFactor;
        if(a<0){
            ampPos = 0;
        }
        controller->setAmpLine(iddl, a,sigma);
        if(verbose){
          Log.LogInfo( "LineModel Infos: Lmfit set final amp [%] = %f with err = %f", iddl, a, sigma);
        }
    }
    //
    // // Analyse if result is good or not
    // if(controller->isContinuumFitted()){
    //   Float64 continuumTplAmp = sqrt(gsl_vector_get(s->x,controller->getIndContinuumAmp()));///normFactor;
    //   if(continuumTplAmp<0){
    //     ampPos =0;
    //   }
    // }
    //
    // if(controller->isEmissionVelocityFitted()){
    //   Float64 vel = gsl_vector_get(s->x,controller->getIndEmissionVel())/controller->getNormEmiFactor();
    //   if(vel<0){
    //     ampPos =0;
    //   }
    // }
    // if(controller->isAbsorptionVelocityFitted()){
    //   Float64 vel = gsl_vector_get(s->x,controller->getIndAbsorptionVel())/controller->getNormAbsFactor();
    //   if(vel<0){
    //     ampPos =0;
    //   }
    // }

    Int32 outputStatus=0;


    // finally populate the fitting results to the output
    if(ampPos && status==0)
    {
        controller->setMerit(chi/normFactor/normFactor);
        if(verbose){
            Log.LogInfo("LineModel Lmfit: Resultat accepted");
        }
        if(controller->isContinuumFitted()){
          Int32 continuumId = controller->getIndContinuumAmp();
          Float64 continuumTplAmp = controller->continuumAmp_LmToModel(gsl_vector_get(s->x,continuumId));///normFactor;
          Float64 continuumTplAmpErr = controller->continuumAmp_LmToModel(c*sqrt(gsl_matrix_get(covar,continuumId,continuumId)));///normFactor;
          controller->setContinummAmp(continuumTplAmp, continuumTplAmpErr);
          if(verbose){
            Log.LogInfo("LineModel Lmfit:  continuumTplAmp calculated =%.3f err =%.3f", continuumTplAmp, continuumTplAmpErr);
          }
        }

        if(controller->isEmissionVelocityFitted()){
          Int32 id = controller->getIndEmissionVel();
          Float64 vel = controller->emiVel_LmToModel(gsl_vector_get(s->x,id));
          Float64 errVel = controller->emiVel_LmToModel(c*sqrt(gsl_matrix_get(covar,id,id)));
          controller->setVelocityEmission(vel, errVel);
          if(verbose) {
            Log.LogInfo( "LineModel Infos: Lmfit Emission velocity found = %.3f with err = %.3f", vel, errVel);
          }
        }
        //TODO ajouter des condition sur les err?
        if(controller->isAbsorptionVelocityFitted()){
          Int32 id = controller->getIndAbsorptionVel();
          Float64 vel = controller->absVel_LmToModel(gsl_vector_get(s->x,id));
          Float64 errVel = controller->absVel_LmToModel(c*sqrt(gsl_matrix_get(covar,id,id)));
          controller->setVelocityAbsorption(vel, errVel);
          if(verbose) {
            Log.LogInfo( "LineModel Infos: Lmfit Absorption velocity found = %.3f with err = %.3f", vel, errVel);
          }
        }

        if(controller->isRedshiftFitted()){
          Int32 id = controller->getIndRedshift();
          Float64 vel =gsl_vector_get(s->x,id);
          Float64 errVel = c*sqrt(gsl_matrix_get(covar,id,id));
          controller->setRedshift(vel, errVel);
          if(verbose) {
            Log.LogInfo( "LineModel Infos: Lmfit Redshift found = %f with err = %f", vel, errVel);
          }
        }

        outputStatus = 0;
    }else{
        Log.LogInfo("LineModel Lmfit: Resultat dump");
        controller->setMerit(DBL_MAX);

        if(controller->isEmissionVelocityFitted()){
          controller->setVelocityEmission(m_velocityEmissionInit, 0);
        }

        if(controller->isAbsorptionVelocityFitted()){
          controller->setVelocityEmission(m_velocityAbsorptionInit, 0);
        }
        if(controller->isContinuumFitted()){
          controller->setContinummAmp(amplitudeContinuumInit, 0);
        }
        outputStatus = 1;
    }

    gsl_multifit_fdfsolver_free (s);
    gsl_matrix_free (covar);
    gsl_matrix_free (J);
    //gsl_rng_free (r);

    return outputStatus;
}

/**
 * \brief Returns a sorted set of samples indices present in the supports of the argument.
 * For each EltsIdx entry, if the entry is not outside lambda range, get the support of each subelement.
 * For each selected support, get the sample index. Sort this list and remove multiple entries. Return this clean list.
 **/
std::vector<UInt32> CLineModelElementList::getSupportIndexes( std::vector<UInt32> EltsIdx )
{
    std::vector<UInt32> indexes;

    TInt32RangeList support;
    for( UInt32 i=0; i<EltsIdx.size(); i++ )
    {
        Int32 iElts = EltsIdx[i];

        if(m_Elements[iElts]->IsOutsideLambdaRange()){
            continue;
        }
        TInt32RangeList s = m_Elements[iElts]->getSupport();
        for( UInt32 iS=0; iS<s.size(); iS++ )
        {
            support.push_back(s[iS]);
        }
    }

    for( UInt32 iS=0; iS<support.size(); iS++ )
    {
        for( UInt32 j=support[iS].GetBegin(); j<support[iS].GetEnd(); j++ )
        {
            indexes.push_back(j);
        }
    }

    std::sort(indexes.begin(), indexes.end());
    indexes.erase( std::unique( indexes.begin(), indexes.end() ), indexes.end() );

    return indexes;
}

Float64 CLineModelElementList::GetWeightingAnyLineCenterProximity(UInt32 sampleIndex, std::vector<UInt32> EltsIdx)
{
    Float64 maxWeight=0.0;
    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    Float64 currentLbda = spectralAxis[sampleIndex];

    for( UInt32 i=0; i<EltsIdx.size(); i++ )
    {
        Int32 iElts = EltsIdx[i];



        TInt32RangeList s = m_Elements[iElts]->getTheoreticalSupport();
        for( UInt32 iS=0; iS<s.size(); iS++ )
        {
            Float64 weight = 0;
            TInt32Range range = s[iS];
            if(sampleIndex>range.GetBegin() && sampleIndex<range.GetEnd())
            {
                Float64 minRangeLbda = spectralAxis[range.GetBegin()];
                Float64 maxRangeLbda = spectralAxis[range.GetEnd()];

                Float64 fullInterval = maxRangeLbda-minRangeLbda;
                if(fullInterval>0.0)
                {
                    Float64 leftInterval = currentLbda-minRangeLbda;
                    Float64 rightInterval = maxRangeLbda-currentLbda;
                    Float64 coeff = std::abs(leftInterval - rightInterval);

                    weight = 1.-(Float64)coeff/(Float64)fullInterval;
                }
            }

            if(maxWeight<weight)
            {
                maxWeight=weight;
            }
        }
    }

    return maxWeight;
}


/**
 * \brief Returns a sorted set of line indices present in the supports of the argument. 
 * Create a vector named indexes.
 * If the argument ind is an index to m_Elements that IsOutSideLambdaRange, return indexes.
 * For each entry in m_Elements:
 *   If the entry has a different linetype than the line corresponding to ind, go to the next entry.
 *   If the entry IsOutsideLambdaRange, go to the enxt entry.
 *   For each subentry in the support of entry:
 *     For each subsubentry in the support of ind:
 *       If the overlap in the spectralAxis is smaller than -1 * overlapThres * winsize, add entry to indexes.
 * Sort indexes, remove duplicates from indexes, and return indexes.
 **/
std::vector<UInt32> CLineModelElementList::getOverlappingElementsBySupport( UInt32 ind, Float64 overlapThres )
{
    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    std::vector<UInt32> indexes;

    if(m_Elements[ind]->IsOutsideLambdaRange())
      {
        indexes.push_back(ind);
        return indexes;
      }
    TInt32RangeList refsupport = m_Elements[ind]->getSupport();
    CRay ray = m_RestRayList[m_Elements[ind]->m_LineCatalogIndexes[0]];
    Int32 linetype = ray.GetType();
    Float64 mu = ray.GetPosition()*(1+m_Redshift);
    CRay::TProfile profile = ray.GetProfile();
    Float64 c = m_Elements[ind]->GetLineWidth(mu, m_Redshift, ray.GetIsEmission(), profile);
    Float64 winsize = m_Elements[ind]->GetNSigmaSupport(profile)*c;
    Float64 overlapThresholdMin = winsize*overlapThres;
    //overlapThresholdMin = 0.0;

    Int32 x1=0;
    Int32 y1=0;
    Int32 x2=0;
    Int32 y2=0;
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        if(m_RestRayList[m_Elements[iElts]->m_LineCatalogIndexes[0]].GetType() != linetype){
            continue;
        }

        if(m_Elements[iElts]->IsOutsideLambdaRange()){
            continue;
        }
        TInt32RangeList s = m_Elements[iElts]->getSupport();
        for( UInt32 iS=0; iS<s.size(); iS++ )
        {
            for( UInt32 iRefS=0; iRefS<refsupport.size(); iRefS++ )
            {
                x1 = refsupport[iRefS].GetBegin();
                x2 = refsupport[iRefS].GetEnd();
                y1 = s[iS].GetBegin();
                y2 = s[iS].GetEnd();

                //Log.LogInfo( "hybrid fit: iRefS=%d - support=%d,%d", iRefS, x1, x2);
                //Log.LogInfo( "hybrid fit: iS=%d - support=%d,%d", iS, y1, y2);

//                if( std::max(x1,y1) < std::min(x2,y2) ){
//                    indexes.push_back(iElts);
//                    break;
//                }

                Float64 max = spectralAxis[std::max(x1,y1)];
                Float64 min = spectralAxis[std::min(x2,y2)];
                if( max-min < -overlapThresholdMin )
		  {
                    indexes.push_back(iElts);
                    break;
                }
            }
        }
    }

    std::sort(indexes.begin(), indexes.end());
    indexes.erase( std::unique( indexes.begin(), indexes.end() ), indexes.end() );

    return indexes;
}

/**
 * \brief Returns a sorted, de-duplicated list of indices of lines whose support overlap ind's support and are not listed in the argument excludedInd.
 **/
std::vector<UInt32> CLineModelElementList::getOverlappingElements(UInt32 ind, std::vector<UInt32> excludedInd, Float64 overlapThres)
{
    std::vector<UInt32> indexes;

    if(m_Elements[ind]->IsOutsideLambdaRange()){
        indexes.push_back(ind);
        return indexes;
    }

    std::vector<CRay> raysRef = m_Elements[ind]->GetRays();
    Int32 linetypeRef = raysRef[0].GetType();

    Int32 xinf=0;
    Int32 yinf=0;
    Int32 xsup=0;
    Int32 ysup=0;
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        //check linetype
        if(m_RestRayList[m_Elements[iElts]->m_LineCatalogIndexes[0]].GetType() != linetypeRef){
            continue;
        }

        //check if outside lambdarange
        if(m_Elements[iElts]->IsOutsideLambdaRange()){
            continue;
        }

        //check if in exclusion list
        bool excluded=false;
        for( UInt32 iexcl=0; iexcl<excludedInd.size(); iexcl++ )
        {
            if(iElts == excludedInd[iexcl]){
               excluded = true;
               break;
            }
        }
        if(excluded){
            continue;
        }

        std::vector<CRay> raysElt = m_Elements[iElts]->GetRays();

        for( UInt32 iRayElt=0; iRayElt<raysElt.size(); iRayElt++ )
        {
            for( UInt32 iRayRef=0; iRayRef<raysRef.size(); iRayRef++ )
            {
                Float64 muRef = raysRef[iRayRef].GetPosition()*(1+m_Redshift);
                CRay::TProfile profileRef = raysRef[iRayRef].GetProfile();
                Float64 cRef = m_Elements[ind]->GetLineWidth(muRef, m_Redshift, raysRef[iRayRef].GetIsEmission(), profileRef);
                Float64 winsizeRef = m_Elements[ind]->GetNSigmaSupport(profileRef)*cRef;
                Float64 overlapSizeMin = winsizeRef*overlapThres;
                xinf = muRef-winsizeRef/2.0;
                xsup = muRef+winsizeRef/2.0;

                Float64 muElt = raysElt[iRayElt].GetPosition()*(1+m_Redshift);
                CRay::TProfile profileElt = raysElt[iRayElt].GetProfile();
                Float64 cElt = m_Elements[iElts]->GetLineWidth(muElt, m_Redshift, raysElt[iRayElt].GetIsEmission(), profileElt);
                Float64 winsizeElt = m_Elements[iElts]->GetNSigmaSupport(profileElt)*cElt;
                yinf = muElt-winsizeElt/2.0;
                ysup = muElt+winsizeElt/2.0;

                Float64 max = std::max(xinf,yinf);
                Float64 min = std::min(xsup,ysup);
                if( max-min < -overlapSizeMin ){
                    indexes.push_back(iElts);
                    break;
                }
            }
        }
    }

    std::sort(indexes.begin(), indexes.end());
    indexes.erase( std::unique( indexes.begin(), indexes.end() ), indexes.end() );

    return indexes;
}

Int32 CLineModelElementList::fitAmplitudesLinSolveAndLambdaOffset(std::vector<UInt32> EltsIdx,
                                                                   const CSpectrumSpectralAxis& spectralAxis,
                                                                   const CSpectrumFluxAxis& fluxAxis,
                                                                   const CSpectrumFluxAxis& continuumfluxAxis,
                                                                   std::vector<Float64>& ampsfitted,
                                                                   std::vector<Float64>& errorsfitted,
                                                                   Bool enableOffsetFitting)
{
    Int32 ret=-1;
    Int32 nSteps = int((m_LambdaOffsetMax-m_LambdaOffsetMin)/m_LambdaOffsetStep+0.5);

    bool atLeastOneOffsetToFit = false;
    if(enableOffsetFitting)
    {
        for(Int32 iE=0; iE<EltsIdx.size(); iE++)
        {
            Float64 nRays = m_Elements[iE]->m_Rays.size();
            for(Int32 iR=0; iR<nRays; iR++)
            {
                //check if the line is to be fitted
                if(m_Elements[iE]->m_Rays[iR].GetOffsetFitEnabled())
                {
                    atLeastOneOffsetToFit = true;
                    break;
                }
            }
        }
    }

    if(!atLeastOneOffsetToFit)
    {
        nSteps = 1;
    }

    Float64 bestMerit = DBL_MAX;
    Int32 idxBestMerit = -1;
    for(Int32 iO=0; iO<nSteps; iO++)
    {
        //set offset value
        if(atLeastOneOffsetToFit)
        {
            Float64 offset = m_LambdaOffsetMin+m_LambdaOffsetStep*iO;
            for(Int32 iE=0; iE<EltsIdx.size(); iE++)
            {
                Float64 nRays = m_Elements[iE]->m_Rays.size();


                for(Int32 iR=0; iR<nRays; iR++)
                {
                    if(m_Elements[iE]->m_Rays[iR].GetOffsetFitEnabled())
                    {
                        m_Elements[iE]->m_Rays[iR].SetOffset(offset);
                    }
                }
            }
        }

        //fit for this offset
        ret = fitAmplitudesLinSolve( EltsIdx, spectralAxis, fluxAxis, continuumfluxAxis, ampsfitted, errorsfitted);

        //check fitting
        if(atLeastOneOffsetToFit)
        {
            Float64 sumFit = 0.0;
            refreshModelUnderElements(EltsIdx);

            //todo: replace lambdarange using elements limits for speed
            for(Int32 iE=0; iE<EltsIdx.size(); iE++)
            {
                Float64 _fit = getModelErrorUnderElement(iE);
                sumFit += _fit;
            }
            if(sumFit<bestMerit)
            {
                bestMerit = sumFit;
                idxBestMerit = iO;
            }
        }

    }

    if(idxBestMerit>=0 && atLeastOneOffsetToFit)
    {
        //set offset value
        if(atLeastOneOffsetToFit)
        {
            Float64 offset = m_LambdaOffsetMin+m_LambdaOffsetStep*idxBestMerit;
            for(Int32 iE=0; iE<EltsIdx.size(); iE++)
            {
                Float64 nRays = m_Elements[iE]->m_Rays.size();
                for(Int32 iR=0; iR<nRays; iR++)
                {
                    if(m_Elements[iE]->m_Rays[iR].GetOffsetFitEnabled())
                    {
                        m_Elements[iE]->m_Rays[iR].SetOffset(offset);
                    }
                }
            }
        }
        //fit again for this offset
        ret = fitAmplitudesLinSolve( EltsIdx, spectralAxis, fluxAxis, continuumfluxAxis, ampsfitted, errorsfitted);
    }


    return ret;
}

/**
 * \brief Use GSL to fit linearly the elements listed in argument EltsIdx.
 * If size of argument EltsIdx is less than 1 return -1.
 **/
Int32 CLineModelElementList::fitAmplitudesLinSolve( std::vector<UInt32> EltsIdx,
                                                    const CSpectrumSpectralAxis& spectralAxis,
                                                    const CSpectrumFluxAxis& fluxAxis,
                                                    const CSpectrumFluxAxis& continuumfluxAxis,
                                                    std::vector<Float64>& ampsfitted,
                                                    std::vector<Float64>& errorsfitted)
{
  //boost::chrono::thread_clock::time_point start_prep = boost::chrono::thread_clock::now();

    bool verbose = false;
    bool useAmpOffset = m_enableAmplitudeOffsets;
    Int32 idxAmpOffset = -1;

    Int32 idx = 0;


    Int32 nddl = EltsIdx.size();
    if(nddl<1){
        return -1;
    }
    std::vector<UInt32> xInds = getSupportIndexes( EltsIdx );
    if(xInds.size()<1)
    {
        return -1;
    }

    if(useAmpOffset)
    {
        nddl +=3;
        //find the amplitudeOffset Support that corresponds to these elts
        for( Int32 i=0; i<m_ampOffsetsIdxStart.size(); i++ )
        {
            if(xInds[0]>= m_ampOffsetsIdxStart[i] && xInds[0]<=m_ampOffsetsIdxStop[i])
            {
                idxAmpOffset=i;
                break;
            }
        }
        if(verbose)
        {
            Log.LogInfo("AmplitudeOffset enabled: idx offset = %d", idxAmpOffset);
        }
    }

    for (Int32 iddl = 0; iddl < EltsIdx.size(); iddl++)
    {
        SetElementAmplitude(EltsIdx[iddl], 1.0, 0.0);
    }

    const Float64* spectral = spectralAxis.GetSamples();
    const Float64* flux = fluxAxis.GetSamples();

    //Linear fit
    int i, n;
    Float64 fval;
    double chisq;
    gsl_matrix *X, *cov;
    gsl_vector *y, *w, *c;

    n = xInds.size();
    if(n<nddl){
        ampsfitted.resize(nddl);
        errorsfitted.resize(nddl);
        for (Int32 iddl = 0; iddl < nddl; iddl++)
        {
            ampsfitted[iddl] = 0.0;
            errorsfitted[iddl] = 1e12;//some high number
        }
        return -1;
    }

    X = gsl_matrix_alloc (n, nddl);
    y = gsl_vector_alloc (n);
    w = gsl_vector_alloc (n);
    c = gsl_vector_alloc (nddl);
    cov = gsl_matrix_alloc (nddl, nddl);

    // Normalize
    Float64 maxabsval = DBL_MIN;
    for (i = 0; i < n; i++)
    {
        idx = xInds[i];
        if(maxabsval<std::abs(flux[idx]))
        {
            maxabsval=std::abs(flux[idx]);
        }
    }
    Float64 normFactor = 1.0/maxabsval;
    if(verbose)
    {
        Log.LogInfo("normFactor = '%.3e'\n", normFactor);
    }

    // Prepare the fit data
    for (i = 0; i < n; i++)
    {
        double xi, yi, ei;
        idx = xInds[i];
        xi = spectral[idx];
        yi = flux[idx]*normFactor;
        ei = m_ErrorNoContinuum[idx]*normFactor;

        for (Int32 iddl = 0; iddl < EltsIdx.size(); iddl++)
        {
            fval =  m_Elements[EltsIdx[iddl]]->getModelAtLambda(xi, m_Redshift, continuumfluxAxis[idx]);
            gsl_matrix_set (X, i, iddl, fval);

            if(verbose)
            {
                fprintf(stderr, "fval = '%.3e'\n", fval);
            }
        }

        if(useAmpOffset)
        {
            gsl_matrix_set (X, i, EltsIdx.size(), 1.0);
            gsl_matrix_set (X, i, EltsIdx.size()+1, xi);
            gsl_matrix_set (X, i, EltsIdx.size()+2, xi*xi);
        }

        gsl_vector_set (y, i, yi);
        gsl_vector_set (w, i, 1.0/(ei*ei));
    }

    //
    // boost::chrono::thread_clock::time_point stop_prep = boost::chrono::thread_clock::now();
    //Float64 duration_prep = boost::chrono::duration_cast<boost::chrono::microseconds>(stop_prep - start_prep).count();
    // boost::chrono::thread_clock::time_point start_fit = boost::chrono::thread_clock::now();

    {
      gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, nddl);
      gsl_multifit_wlinear (X, w, y, c, cov, &chisq, work);
      gsl_multifit_linear_free (work);
    }

    //
    //boost::chrono::thread_clock::time_point stop_fit = boost::chrono::thread_clock::now();
    //Float64 duration_fit = boost::chrono::duration_cast<boost::chrono::microseconds>(stop_fit - start_fit).count();
    //Log.LogInfo( "LineModel linear fit: prep = %.3f - fit = %.3f", duration_prep, duration_fit);

    if(verbose)
    {
#define C(i) (gsl_vector_get(c,(i)))
#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))
        if(1){
            Log.LogInfo("# best fit: Y = %g X1 + %g X2 ...", C(0), C(1));
            Log.LogInfo("# covariance matrix:");
            Log.LogInfo("[");
            Log.LogInfo("  %+.5e, %+.5e", COV(0,0), COV(0,1));
            Log.LogInfo("  %+.5e, %+.5e", COV(1,0), COV(1,1));

            //        Log.LogInfo("[ %+.5e, %+.5e, %+.5e  \n", COV(0,0), COV(0,1), COV(0,2));
            //        Log.LogInfo("  %+.5e, %+.5e, %+.5e  \n", COV(1,0), COV(1,1), COV(1,2));
            //        Log.LogInfo("  %+.5e, %+.5e, %+.5e ]\n", COV(2,0), COV(2,1), COV(2,2));

            Log.LogInfo("]");
            Log.LogInfo("# chisq/n = %g", chisq/n);
        }

        for (Int32 iddl = 0; iddl < nddl; iddl++)
        {
            Float64 a = gsl_vector_get(c,iddl)/normFactor;
            if(verbose)
            {
                Log.LogInfo("# Found amplitude %d: %+.5e", iddl, a);
            }
        }
    }

    Int32 sameSign = 1;
    Float64 a0 = gsl_vector_get(c,0)/normFactor;
    for (Int32 iddl = 1; iddl < EltsIdx.size(); iddl++)
    {
        Float64 a = gsl_vector_get(c,iddl)/normFactor;
        Float64 product = a0*a;
        if(product<0){
            sameSign = 0;
        }
    }

    if(verbose)
    {
        Log.LogInfo("# Found amplitudes with sameSign=%d", sameSign);
    }
    if(sameSign){
        for (Int32 iddl = 0; iddl < EltsIdx.size(); iddl++)
        {
            Float64 a = gsl_vector_get(c,iddl)/normFactor;
            Float64 cova = gsl_matrix_get(cov,iddl,iddl);
            Float64 sigma = sqrt(cova)/normFactor;
            SetElementAmplitude(EltsIdx[iddl], a, sigma);
        }
        //refreshModel();
    }else{
        ampsfitted.resize(EltsIdx.size());
        errorsfitted.resize(EltsIdx.size());
        for (Int32 iddl = 0; iddl < EltsIdx.size(); iddl++)
        {
            Float64 a = gsl_vector_get(c,iddl)/normFactor;
            Float64 cova = gsl_matrix_get(cov,iddl,iddl);
            Float64 sigma = sqrt(cova)/normFactor;
            SetElementAmplitude(EltsIdx[iddl], a, sigma);
            ampsfitted[iddl] = (a);
            errorsfitted[iddl] = (sigma);
        }
    }
    if(useAmpOffset)
    {
        Float64 x0 = gsl_vector_get(c,EltsIdx.size())/normFactor;
        m_ampOffsetsX0[idxAmpOffset] = x0;
        Float64 x1 = gsl_vector_get(c,EltsIdx.size()+1)/normFactor;
        m_ampOffsetsX1[idxAmpOffset] = x1;
        Float64 x2 = gsl_vector_get(c,EltsIdx.size()+2)/normFactor;
        m_ampOffsetsX2[idxAmpOffset] = x2;
    }

    gsl_matrix_free (X);
    gsl_vector_free (y);
    gsl_vector_free (w);
    gsl_vector_free (c);
    gsl_matrix_free (cov);

    return sameSign;
}


/**
 * @brief CLineModelElementList::fitAmplitudesLinesAndContinuumLinSolve
 * @param EltsIdx:  elements to be fitted
 * @param lambdaRange
 * @param spectralAxis
 * @param fluxAxis: must contain the observed spectrum
 * @param continuumfluxAxis: must contain the continuum to be amp-fitted
 * @param ampsfitted: output
 * @param errorsfitted: output
 * @param polyOrder: order of the polynom to be fitted along with the continuum (-1 = disabled)
 * @return
 *
 * WARNING: not sure about fitting abs. lines with this method...
 */
Int32 CLineModelElementList::fitAmplitudesLinesAndContinuumLinSolve( std::vector<UInt32> EltsIdx,
                                                                     const TFloat64Range& lambdaRange,
                                                                     const CSpectrumSpectralAxis& spectralAxis,
                                                                     const CSpectrumFluxAxis& fluxAxis,
                                                                     const CSpectrumFluxAxis& continuumfluxAxis,
                                                                     std::vector<Float64>& ampsfitted,
                                                                     std::vector<Float64>& errorsfitted,
                                                                     Float64& chisquare,
                                                                     Int32 polyOrder)
{
    //boost::chrono::thread_clock::time_point start_prep = boost::chrono::thread_clock::now();

    bool verbose = false;
    Int32 idx = 0;

    Int32 nddl = EltsIdx.size()+1; //number of param to be fitted=nlines+continuum
    if(nddl<2){
        return -1;
    }

    if(polyOrder>=0)
    {
        nddl += polyOrder+1;
    }

    std::vector<UInt32> xInds;
    for(Int32 i=0; i<spectralAxis.GetSamplesCount(); i++)
    {
        if(spectralAxis[i]>=lambdaRange.GetBegin() && spectralAxis[i]<=lambdaRange.GetEnd())
        {
            xInds.push_back(i);
        }
    }
    if(xInds.size()<1)
    {
        return -1;
    }

    for (Int32 iddl = 0; iddl < EltsIdx.size(); iddl++)
    {
        SetElementAmplitude(EltsIdx[iddl], 1.0, 0.0);
    }

    const Float64* spectral = spectralAxis.GetSamples();
    const Float64* flux = fluxAxis.GetSamples();

    //Linear fit
    int i, n;
    Float64 fval;
    double chisq;
    gsl_matrix *X, *cov;
    gsl_vector *y, *w, *c;

    n = xInds.size();
    if(n<nddl){
        ampsfitted.resize(nddl);
        errorsfitted.resize(nddl);
        for (Int32 iddl = 0; iddl < nddl; iddl++)
        {
            ampsfitted[iddl] = 0.0;
            errorsfitted[iddl] = 1e12;//some high number
        }
        return -1;
    }

    X = gsl_matrix_alloc (n, nddl);
    y = gsl_vector_alloc (n);
    w = gsl_vector_alloc (n);
    c = gsl_vector_alloc (nddl);
    cov = gsl_matrix_alloc (nddl, nddl);

    // Normalize
    Float64 maxabsval = DBL_MIN;
    for (i = 0; i < n; i++)
    {
        idx = xInds[i];
        if(maxabsval<std::abs(flux[idx]))
        {
            maxabsval=std::abs(flux[idx]);
        }
    }
    Float64 normFactor = 1.0/maxabsval;
    if(verbose)
    {
        Log.LogInfo("normFactor = '%.3e'\n", normFactor);
    }

    // Prepare the fit data
    for (i = 0; i < n; i++)
    {
        double xi, yi, ei, ci;
        idx = xInds[i];
        xi = spectral[idx];
        yi = flux[idx]*normFactor;
        ei = m_ErrorNoContinuum[idx]*normFactor;
        ci = continuumfluxAxis[idx];

        for (Int32 iddl = 0; iddl < EltsIdx.size(); iddl++)
        {
            fval =  m_Elements[EltsIdx[iddl]]->getModelAtLambda(xi, m_Redshift, ci);
            gsl_matrix_set (X, i, iddl, fval);

            if(false && verbose)
            {
                fprintf(stderr, "fval = '%.3e'\n", fval);
            }
        }

        gsl_matrix_set (X, i, EltsIdx.size(), ci);


        if(polyOrder>=0)
        {
            for(Int32 kCoeff=0; kCoeff<polyOrder+1; kCoeff++)
            {
                Float64 vect=1.0;
                for(Int32 kt=0; kt<kCoeff; kt++)
                {
                    vect *= xi;
                }
                gsl_matrix_set (X, i, EltsIdx.size()+1+kCoeff, vect);
            }
        }

        gsl_vector_set (y, i, yi);
        gsl_vector_set (w, i, 1.0/(ei*ei));
    }

    //
    // boost::chrono::thread_clock::time_point stop_prep = boost::chrono::thread_clock::now();
    //Float64 duration_prep = boost::chrono::duration_cast<boost::chrono::microseconds>(stop_prep - start_prep).count();
    // boost::chrono::thread_clock::time_point start_fit = boost::chrono::thread_clock::now();


    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, nddl);
    gsl_multifit_wlinear (X, w, y, c, cov, &chisq, work);
    gsl_multifit_linear_free (work);


    //
    //boost::chrono::thread_clock::time_point stop_fit = boost::chrono::thread_clock::now();
    //Float64 duration_fit = boost::chrono::duration_cast<boost::chrono::microseconds>(stop_fit - start_fit).count();
    //Log.LogInfo( "LineModel linear fit: prep = %.3f - fit = %.3f", duration_prep, duration_fit);

    if(verbose)
    {
#define C(i) (gsl_vector_get(c,(i)))
#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))
        if(1){
            Log.LogInfo("# best fit: Y = %g X1 + %g X2 ...", C(0), C(1));
            Log.LogInfo("# covariance matrix:");
            Log.LogInfo("[");
            Log.LogInfo("  %+.5e, %+.5e", COV(0,0), COV(0,1));
            Log.LogInfo("  %+.5e, %+.5e", COV(1,0), COV(1,1));

            //        Log.LogInfo("[ %+.5e, %+.5e, %+.5e  \n", COV(0,0), COV(0,1), COV(0,2));
            //        Log.LogInfo("  %+.5e, %+.5e, %+.5e  \n", COV(1,0), COV(1,1), COV(1,2));
            //        Log.LogInfo("  %+.5e, %+.5e, %+.5e ]\n", COV(2,0), COV(2,1), COV(2,2));

            Log.LogInfo("]");
            Log.LogInfo("# chisq = %g", chisq);
            Log.LogInfo("# chisq/n = %g", chisq/n);
        }

        for (Int32 iddl = 0; iddl < nddl; iddl++)
        {
            Float64 a = gsl_vector_get(c,iddl)/normFactor;
            if(verbose)
            {
                Log.LogInfo("# Found amplitude %d: %+.5e", iddl, a);
            }
        }
    }

    Int32 sameSign = 1;
    Float64 a0 = gsl_vector_get(c,0)/normFactor;
    for (Int32 iddl = 1; iddl < EltsIdx.size(); iddl++)
    {
        Float64 a = gsl_vector_get(c,iddl)/normFactor;
        Float64 product = a0*a;
        if(product<0){
            sameSign = 0;
        }
    }

    if(verbose)
    {
        Log.LogInfo("# Found n=%d amplitudes with sameSign=%d", EltsIdx.size(), sameSign);
    }

    ampsfitted.resize(nddl);
    errorsfitted.resize(nddl);
    for (Int32 iddl = 0; iddl < nddl; iddl++)
    {
        Float64 a = gsl_vector_get(c,iddl)/normFactor;
        Float64 cova = gsl_matrix_get(cov,iddl,iddl);
        Float64 sigma = sqrt(cova)/normFactor;
        if(iddl<EltsIdx.size())
        {
            SetElementAmplitude(EltsIdx[iddl], a, sigma);
        }
        ampsfitted[iddl] = (a);
        errorsfitted[iddl] = (sigma);
    }


    if(verbose && polyOrder>=0)
    {
        for(Int32 kCoeff=0; kCoeff<polyOrder+1; kCoeff++)
        {
            Float64 p = gsl_vector_get(c, EltsIdx.size()+1+kCoeff)/normFactor;
            Log.LogInfo("# Found p%d poly amplitude = %+.5e", kCoeff, p);
        }
    }

    if(verbose)
    {
        Log.LogInfo("# Returning (L+C) n=%d amplitudes", ampsfitted.size());
    }

    gsl_matrix_free (X);
    gsl_vector_free (y);
    gsl_vector_free (w);
    gsl_vector_free (c);
    gsl_matrix_free (cov);

    chisquare = chisq;
    return sameSign;
}


/**
* @brief CLineModelElementList::setLyaProfile
* If a Lya line is present with ASYMFIT profile, fit the width and asymmetry parameters
* If a Lya line is present with ASYMFIXED profile, set the width and asymmetry parameters according to profile parameters given in the string
* @param redshift, spectralAxis
* @return 1 if successfully fitted, 0 if error, 2 if Lya not present, 3 if Lya not configured to be fitted in the catalog
*/
Int32 CLineModelElementList::setLyaProfile(Float64 redshift, const CSpectrumSpectralAxis &spectralAxis)
{
    Int32 verbose=0;
    if(verbose)
    {
        Log.LogInfo("Setting Lya Profile");
    }

    //1. retrieve the Lya index
    linetags ltags;
    std::string lyaTag = ltags.lya_em;
    Int32 idxLineLyaE = -1;
    Int32 idxLyaE = FindElementIndex(lyaTag, -1, idxLineLyaE);
    if( idxLyaE<0 || idxLineLyaE<0 )
    {
        return 2; //Lya alpha not found
    }
    std::vector<UInt32> filterEltsIdxLya;
    filterEltsIdxLya.push_back(idxLyaE);

    //2. check if the profile has to be fitted
    bool asimfitProfileFound = false;
    bool asimfixedProfileFound = false;

    Int32 nrays = m_Elements[idxLyaE]->GetSize();
    for(Int32 iray=0; iray<nrays; iray++)
    {
        Int32 lineIndex = m_Elements[idxLyaE]->m_LineCatalogIndexes[iray];
        CRay::TProfile profile = m_RestRayList[lineIndex].GetProfile();
        bool lineOutsideLambdaRange = m_Elements[idxLyaE]->IsOutsideLambdaRange();
        if( profile==CRay::ASYMFIT && !lineOutsideLambdaRange)
        {
            asimfitProfileFound = true;
            break;
        }else if(profile==CRay::ASYMFIXED)
        {
            //use the manual fixed profile parameters from catalog profile string
            asimfixedProfileFound = true;
            TAsymParams asym_params = m_RestRayList[lineIndex].GetAsymParams();
            m_Elements[idxLyaE]->SetAsymfitWidthCoeff(asym_params.width);
            m_Elements[idxLyaE]->SetAsymfitAlphaCoeff(asym_params.alpha);
            m_Elements[idxLyaE]->SetAsymfitDelta(asym_params.delta);
            break;
        }
    }
    if( !asimfitProfileFound && !asimfixedProfileFound)
    {
        return 3; //no profile found with parameters to be fitted
    }

    //FIT the profile parameters
    if(asimfitProfileFound)
    {
        if(verbose)
        {
            Log.LogInfo("Fitting Lya Profile: width, asym, delta");
        }
        //3. find the best width and asym coeff. parameters
        Float64 widthCoeffStep = 1.0;
        Float64 widthCoeffMin = 1.0;
        Float64 widthCoeffMax = 4.0;
        Int32 nWidthSteps = int((widthCoeffMax-widthCoeffMin)/widthCoeffStep+0.5);
        Float64 asymCoeffStep = 0.5;
        Float64 asymCoeffMin = 0.0;
        Float64 asymCoeffMax = 2.5;
        Int32 nAsymSteps = int((asymCoeffMax-asymCoeffMin)/asymCoeffStep+0.5);
        Float64 deltaStep = 0.5;
        Float64 deltaMin = 0.0;
        Float64 deltaMax = 0.0;//4.0;
        Int32 nDeltaSteps = int((deltaMax-deltaMin)/deltaStep+0.5);
        nDeltaSteps = 1;//disable delta fit here: will be parameterized through the offset object

        Float64 bestWidth = widthCoeffMin;
        Float64 bestAlpha = asymCoeffMin;
        Float64 bestDelta = deltaMin;
        Float64 meritMin = boost::numeric::bounds<float>::highest();

        for(Int32 iDelta=0; iDelta<nDeltaSteps; iDelta++)
        {
            Float64 delta = deltaMin + deltaStep*iDelta;
            for(Int32 iWidth=0; iWidth<nWidthSteps; iWidth++)
            {
                Float64 asymWidthCoeff = widthCoeffMin + widthCoeffStep*iWidth;
                for(Int32 iAsym=0; iAsym<nAsymSteps; iAsym++)
                {
                    Float64 asymAlphaCoeff = asymCoeffMin + asymCoeffStep*iAsym;

                    m_Elements[idxLyaE]->SetAsymfitDelta(delta);
                    m_Elements[idxLyaE]->SetAsymfitWidthCoeff(asymWidthCoeff);
                    m_Elements[idxLyaE]->SetAsymfitAlphaCoeff(asymAlphaCoeff);

                    idxLineLyaE = -1;
                    m_Elements[idxLyaE]->fitAmplitude(spectralAxis, m_spcFluxAxisNoContinuum, m_ContinuumFluxAxis, redshift, idxLineLyaE);
                    Float64 m=m_dTransposeDNocontinuum;
                    if(0)
                    {
                        refreshModelUnderElements(filterEltsIdxLya, idxLineLyaE);
                        m = getModelErrorUnderElement(idxLyaE);
                    }else{
                        m = getLeastSquareMeritFast(idxLineLyaE);
                    }
                    if( m<meritMin )
                    {
                        meritMin = m;
                        bestWidth = m_Elements[idxLyaE]->GetAsymfitWidthCoeff();
                        bestAlpha = m_Elements[idxLyaE]->GetAsymfitAlphaCoeff();
                        bestDelta = m_Elements[idxLyaE]->GetAsymfitDelta();
                    }
                }
            }
        }

        //4. set the associated Lya members in the element definition
        if(verbose)
        {
            Log.LogInfo("Lya Profile found: width=%f, asym=%f, delta=%f", bestWidth, bestAlpha, bestDelta);
        }
        m_Elements[idxLyaE]->SetAsymfitWidthCoeff(bestWidth);
        m_Elements[idxLyaE]->SetAsymfitAlphaCoeff(bestAlpha);
        m_Elements[idxLyaE]->SetAsymfitDelta(bestDelta);
    }

    return 1;
}


/**
* @brief CLineModelElementList::ReestimateContinuumUnderLines
* For each line, reestimate the continuum using the original median routines on a sub segment of the original spectrum
* - the subsegment is a contiguous segment around the support of all the elements
* @param EltsIdx
* @return
*/
std::vector<UInt32> CLineModelElementList::ReestimateContinuumUnderLines(std::vector<UInt32> EltsIdx)
{
    if(EltsIdx.size()<1){
        std::vector<UInt32> empty;
        return empty;
    }

    //smoothing factor in continuum median filter
    Float64 smoof = 150;

    //modify m_ContinuumFluxAxis
    CSpectrumFluxAxis& fluxAxisModified = m_ContinuumFluxAxis;
    Float64* Ycont = fluxAxisModified.GetSamples();
    CSpectrumSpectralAxis spectralAxis = m_SpectrumModel->GetSpectralAxis();


    std::vector<UInt32> xInds = getSupportIndexes( EltsIdx );
    Int32 minInd = xInds[0];
    Int32 maxInd = xInds[xInds.size()-1];
    //
    UInt32 iminMerge = spectralAxis.GetIndexAtWaveLength(spectralAxis[minInd]-2*smoof);
    if(iminMerge==0){
        iminMerge = 1;
    }
    UInt32 imaxMerge = spectralAxis.GetIndexAtWaveLength(spectralAxis[maxInd]+2*smoof);
    //
    //
    UInt32 iminMerge2 = spectralAxis.GetIndexAtWaveLength(spectralAxis[minInd]-smoof);
    if(iminMerge2==0){
        iminMerge2 = 1;
    }
    UInt32 imaxMerge2 = spectralAxis.GetIndexAtWaveLength(spectralAxis[maxInd]+smoof);
    //
    UInt32 imin = spectralAxis.GetIndexAtWaveLength(spectralAxis[minInd]-3*smoof);
    if(imin==0){
        imin = 1;
    }
    Int32 imax = spectralAxis.GetIndexAtWaveLength(spectralAxis[maxInd]+3*smoof);
    Int32 spcSize = imax-imin+1;


    //prepare the line model with no continuum;
    CSpectrumFluxAxis modelFluxAxisTmp = m_SpectrumModel->GetFluxAxis();
    for(UInt32 i=0; i<spcSize; i++){
        modelFluxAxisTmp[i+imin]=0.0;
    }
    for(UInt32 idx=0; idx<EltsIdx.size(); idx++){
        UInt32 eltIdx = EltsIdx[idx];
        m_Elements[eltIdx]->addToSpectrumModel(spectralAxis, modelFluxAxisTmp, m_ContinuumFluxAxis, m_Redshift);

    }

    //gather the modified indexes
    std::vector<UInt32> modifiedIdxs;

    //create the spcBuffer with spectrum minus the lines
    CSpectrum spcBuffer;
    CSpectrumSpectralAxis *_SpectralAxis = new CSpectrumSpectralAxis(spcSize, false);
    CSpectrumFluxAxis *_FluxAxis = new CSpectrumFluxAxis(spcSize);

    for(Int32 i=0; i<spcSize; i++){
        (*_SpectralAxis)[i] = spectralAxis[i+imin];
        (*_FluxAxis)[i] = m_SpcFluxAxis[i+imin] - modelFluxAxisTmp[i+imin];
        //                if( error!= NULL ){
        //                    (*_FluxAxis).GetError()[i] = tmpError[i];
        //                }
    }
    spcBuffer.GetSpectralAxis() = *_SpectralAxis;
    spcBuffer.GetFluxAxis() = *_FluxAxis;

    /*
    // export for debug
    FILE* fspc = fopen( "ReestimateContinuumUnderLines_correctedSpc_dbg.txt", "w+" );
    Float64 coeffSaveSpc = 1e16;
    for( Int32 t=0;t<spcBuffer.GetSampleCount();t++)
    {
        fprintf( fspc, "%f %f %f\n", t, spcBuffer.GetSpectralAxis()[t], (m_SpcFluxAxis[t+imin])*coeffSaveSpc, (spcBuffer.GetFluxAxis()[t])*coeffSaveSpc);
    }
    fclose( fspc );
    //*/

    //apply continuum routine on this spcbuffer
    CContinuumIrregularSamplingMedian continuum;
    CSpectrumFluxAxis fluxAxisWithoutContinuumCalc;
    continuum.RemoveContinuum( spcBuffer, fluxAxisWithoutContinuumCalc );


    /*
    // export for debug
    FILE* f = fopen( "continuum_reestimated_underlines_dbg.txt", "w+" );
    Float64 coeffSave = 1e16;
    for( Int32 t=iminMerge;t<imaxMerge;t++)
    {
        fprintf( f, "%f %f %f\n", t, spectralAxis[t], (m_SpcContinuumFluxAxis[t])*coeffSave, (spcBuffer.GetFluxAxis()[t-imin] - fluxAxisWithoutContinuumCalc[t-imin])*coeffSave);
    }
    fclose( f );
    //*/

    Float64 modified = 0.0;
    Float64 coeff=0.0;
    //merge raw continuum free with the newly calculated cont. under the line, (todo: with cross-fade on the borders)
    for(UInt32 i=iminMerge; i<imaxMerge; i++){
        modified = spcBuffer.GetFluxAxis()[i-imin] - fluxAxisWithoutContinuumCalc[i-imin];
        coeff = 1.0;

        if(i<=iminMerge2){
            coeff = (Float64(i-iminMerge)/Float64(iminMerge2-iminMerge));
        }else if(i>=imaxMerge2){
            coeff = 1.0-(Float64(i-imaxMerge2)/Float64(imaxMerge-imaxMerge2));
        }

        Ycont[i] = coeff*modified + (1-coeff)*Ycont[i];
        modifiedIdxs.push_back(i);
    }

    std::sort(modifiedIdxs.begin(), modifiedIdxs.end());
    modifiedIdxs.erase( std::unique( modifiedIdxs.begin(), modifiedIdxs.end() ), modifiedIdxs.end() );
    return modifiedIdxs;
}

/**
 * \brief Modifies the fluxAxis where SNR is good enough, and returns the list of modified indexes of fluxAxis.
 **/
std::vector<UInt32> CLineModelElementList::ReestimateContinuumApprox(std::vector<UInt32> EltsIdx)
{
    //smoothing factor in continuum median filter
    Float64 smoof = 150;

    //
    std::vector<UInt32> modifiedIdxs;


    //modify m_ContinuumFluxAxis
    CSpectrumFluxAxis& fluxAxisModified = m_ContinuumFluxAxis;
    Float64* Ycont = fluxAxisModified.GetSamples();
    CSpectrumSpectralAxis spectralAxis = m_SpectrumModel->GetSpectralAxis();

    for(Int32 idx=0; idx<EltsIdx.size(); idx++){
        Int32 eltIdx = EltsIdx[idx];

        Int32 nrays = m_Elements[eltIdx]->GetSize();
        for(Int32 iray=0; iray<nrays; iray++)
        {
            TInt32Range s = m_Elements[eltIdx]->getSupportSubElt(iray);
            Int32 smin = s.GetBegin();
            Int32 smax = s.GetEnd();
            if(smax - smin < 2){
                continue;
            }
            Int32 imin = spectralAxis.GetIndexAtWaveLength(spectralAxis[smin]-smoof);
            if(imin==0){
                imin = 1;
            }
            Int32 imax = spectralAxis.GetIndexAtWaveLength(spectralAxis[smax]+smoof);
            Int32 sSize = imax-imin+1;

            Float64 A = m_Elements[eltIdx]->GetFittedAmplitude(iray);
            Float64 Sigma = m_Elements[eltIdx]->GetFittedAmplitudeErrorSigma(iray);

            if(A<=0 || std::abs(Sigma)>std::abs(A)){ //todo: check this error sigma rule, should we add a sigma thres. ?
                continue;
            }
            A*= m_Elements[eltIdx]->GetSignFactor(iray);

            Float64 integratedA = A*m_Elements[eltIdx]->GetWidth(iray, m_Redshift)*sqrt(2*M_PI);
            Float64 coeffA = integratedA/(Float64)sSize;

            Float64 term=0.0;
            for(Int32 i=imin; i<imax; i++){
	      Float64 dx=spectralAxis[imin]-spectralAxis[imin-1];
                if(i>smin && i<smax){
                    term = coeffA/dx;
                }
                else if(i>imin && i<=smin){
                    term = ((Float64(i-imin)/Float64(smin-imin))) * coeffA/dx;
                }else if(i>=smax && i<imax){
                    term = (1.0-(Float64(i-smax)/Float64(imax-smax))) * coeffA/dx;
                }else{
                    term=0.0;
                }
                Ycont[i] -= term;

                modifiedIdxs.push_back(i);
            }
        }
    }

    std::sort(modifiedIdxs.begin(), modifiedIdxs.end());
    modifiedIdxs.erase( std::unique( modifiedIdxs.begin(), modifiedIdxs.end() ), modifiedIdxs.end() );
    return modifiedIdxs;
}

/**
 * \brief Copies the continuum flux to the model flux, and the no continuum flux receives the value of spectrum flux minus the continuum flux.
 **/
void CLineModelElementList::refreshModelAfterContReestimation(std::vector<UInt32> xInds, CSpectrumFluxAxis& modelFluxAxis, CSpectrumFluxAxis& spcFluxAxisNoContinuum)
{
    Int32 n = xInds.size();

    Int32 idx = 0;
    for (Int32 i = 0; i < n; i++)
    {
        idx = xInds[i];

        modelFluxAxis[idx] = m_ContinuumFluxAxis[idx];
        spcFluxAxisNoContinuum[idx] = m_SpcFluxAxis[idx]-m_ContinuumFluxAxis[idx];
    }
}

/**
 * \brief Accumulates the squared differences between model and spectrum in the argument lambdaRange and returns the sum.
 **/
Float64 CLineModelElementList::getLeastSquareMerit(const TFloat64Range& lambdaRange)
{
    const CSpectrumSpectralAxis& spcSpectralAxis = m_SpectrumModel->GetSpectralAxis();
    const CSpectrumFluxAxis& spcFluxAxis = m_SpcFluxAxis;
    const CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();

    //Int32 numDevs = 0;
    Float64 fit = 0.0;
    const Float64* Ymodel = modelFluxAxis.GetSamples();
    const Float64* Yspc = spcFluxAxis.GetSamples();
    Float64 diff = 0.0;

    Int32 imin = spcSpectralAxis.GetIndexAtWaveLength(lambdaRange.GetBegin());
    Int32 imax = spcSpectralAxis.GetIndexAtWaveLength(lambdaRange.GetEnd());
    for(Int32 j=imin; j<imax; j++)
    {
        //numDevs++;
        diff = (Yspc[j] - Ymodel[j]);
        fit += (diff*diff) / (m_ErrorNoContinuum[j]*m_ErrorNoContinuum[j]);
        //        if ( 1E6 * diff < m_ErrorNoContinuum[j] )
        //        {
        //            Log.LogDebug( "Warning: noise is at least 6 orders greater than the residue!" );
        //            Log.LogDebug( "CLineModelElementList::getLeastSquareMerit diff = %f", diff );
        //            Log.LogDebug( "CLineModelElementList::getLeastSquareMerit m_ErrorNoContinuum[%d] = %f", j, m_ErrorNoContinuum[j] );
        //        }
    }
    Log.LogDebug( "CLineModelElementList::getLeastSquareMerit fit = %f", fit );
    if(std::isnan(fit))
    {
        Log.LogError( "CLineModelElementList::getLeastSquareMerit: NaN value found on the lambdarange = (%f, %f)", lambdaRange.GetBegin(), lambdaRange.GetEnd() );
        Log.LogError( "CLineModelElementList::getLeastSquareMerit: NaN value found on the true observed spectral axis lambdarange = (%f, %f)", spcSpectralAxis[imin], spcSpectralAxis[imax] );
        for(Int32 j=imin; j<imax; j++)
        {
            if(std::isnan(Yspc[j]))
            {
                Log.LogError( "CLineModelElementList::getLeastSquareMerit: NaN value found for the observed spectrum at lambda=%f", spcSpectralAxis[j] );
                break;
            }
            if(std::isnan(Ymodel[j]))
            {
                Log.LogError( "CLineModelElementList::getLeastSquareMerit: NaN value found for the model at lambda=%f", spcSpectralAxis[j] );
                break;
            }

            if(std::isnan(m_ErrorNoContinuum[j]))
            {
                Log.LogError( "CLineModelElementList::getLeastSquareMerit: NaN value found for the sqrt(variance) at lambda=%f", spcSpectralAxis[j] );
                break;
            }
            if(m_ErrorNoContinuum[j]==0.0)
            {
                Log.LogError( "CLineModelElementList::getLeastSquareMerit: 0 value found for the sqrt(variance) at lambda=%f", spcSpectralAxis[j] );
                break;
            }
        }
        throw runtime_error("CLineModelElementList::getLeastSquareMerit: NaN value found");
    }
    return fit;
}

Float64 CLineModelElementList::getLeastSquareContinuumMerit(const TFloat64Range& lambdaRange)
{
    const CSpectrumSpectralAxis& spcSpectralAxis = m_SpectrumModel->GetSpectralAxis();
    const CSpectrumFluxAxis& spcFluxAxis = m_SpcFluxAxis;

    Int32 numDevs = 0;
    Float64 fit = 0.0;
    const Float64* YCont = m_ContinuumFluxAxis.GetSamples();
    const Float64* Yspc = spcFluxAxis.GetSamples();
    Float64 diff = 0.0;

    Float64 imin = spcSpectralAxis.GetIndexAtWaveLength(lambdaRange.GetBegin());
    Float64 imax = spcSpectralAxis.GetIndexAtWaveLength(lambdaRange.GetEnd());
    for( UInt32 j=imin; j<imax; j++ )
    {
        numDevs++;
        diff = (Yspc[j] - YCont[j]);
        fit += (diff*diff) / (m_ErrorNoContinuum[j]*m_ErrorNoContinuum[j]);
    }
    return fit;
}


/**
 * \brief Get the squared difference by fast method proposed by D. Vibert
 **/
Float64 CLineModelElementList::getLeastSquareMeritFast(Int32 idxLine)
{
    Float64 fit;

    fit = getLeastSquareContinuumMeritFast();

    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        if(idxLine!=-1 && idxLine!=iElts)
        {
            continue;
        }
        Float64 dtm = m_Elements[iElts]->GetSumCross();
        Float64 mtm = m_Elements[iElts]->GetSumGauss();
        Float64 a = m_Elements[iElts]->GetFitAmplitude();
        Float64 term1 = a*a*mtm;
        Float64 term2 = - 2.*a*dtm;
        fit += term1 + term2;
    }

    Log.LogDebug( "CLineModelElementList::getLeastSquareMerit fit fast = %f", fit );
    return fit;
}

Float64 CLineModelElementList::getLeastSquareContinuumMeritFast()
{
    Float64 fit;

    if( m_ContinuumComponent=="tplfit" )
    {
        fit = m_dTransposeDRaw;

        Float64 term1 = m_fitContinuum_tplFitAmplitude*m_fitContinuum_tplFitAmplitude*m_fitContinuum_tplFitMtM;
        Float64 term2 = - 2.*m_fitContinuum_tplFitAmplitude*m_fitContinuum_tplFitDtM;
        fit += term1 + term2;
    }else{
        fit = m_dTransposeDNocontinuum;
    }

    return fit;
}



/**
 * \brief Get the scale marginalization correction
 *
 * WARNING: (todo-check) for lm-rules or lm-free, if hybrid method was used to fit, mtm and dtm are not estimated for now...
 **/
Float64 CLineModelElementList::getScaleMargCorrection(Int32 idxLine)
{
    Float64 corr=0.0;

    //scale marg for continuum
    corr += getContinuumScaleMargCorrection();



    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        if(idxLine!=-1 && idxLine!=iElts)
        {
            continue;
        }
        if(m_Elements[iElts]->IsOutsideLambdaRange() == true){
            continue;
        }
        if(m_Elements[iElts]->GetElementAmplitude()<=0.0){
            continue;
        }

        Float64 mtm = m_Elements[iElts]->GetSumGauss();
        if(mtm>0.0){
            corr += log(mtm);
        }
    }

    return corr;
}


Float64 CLineModelElementList::getContinuumScaleMargCorrection()
{
    Float64 corr=0.0;

    //scale marg for continuum
    if(m_ContinuumComponent == "tplfit") //the support has to be already computed when LoadFitContinuum() is called
    {
        corr += log(m_fitContinuum_tplFitMtM);
    }

    return corr;
}

/**
 * \brief Returns the number of spectral samples between lambdaRange.
 **/
Int32 CLineModelElementList::getSpcNSamples(const TFloat64Range& lambdaRange){
    const CSpectrumSpectralAxis& spcSpectralAxis = m_SpectrumModel->GetSpectralAxis();

    Int32 numDevs = 0;

    Float64 imin = spcSpectralAxis.GetIndexAtWaveLength(lambdaRange.GetBegin());
    Float64 imax = spcSpectralAxis.GetIndexAtWaveLength(lambdaRange.GetEnd());

    for( UInt32 j=imin; j<imax; j++ )
    {
        numDevs++;
    }

    return numDevs;
}

/**
 * \brief Accumulates the squared differences between model and spectrum and returns the sum.
 **/
Float64 CLineModelElementList::getLeastSquareMeritUnderElements()
{
    const CSpectrumFluxAxis& spcFluxAxis = m_SpcFluxAxis;
    const CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();

    Int32 numDevs = 0;
    Float64 fit = 0;
    const Float64* Ymodel = modelFluxAxis.GetSamples();
    const Float64* Yspc = spcFluxAxis.GetSamples();
    Float64 diff = 0.0;

    TInt32RangeList support;
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        if(m_Elements[iElts]->IsOutsideLambdaRange()){
            continue;
        }
        TInt32RangeList s = m_Elements[iElts]->getSupport();
        for( UInt32 iS=0; iS<s.size(); iS++ )
        {
            support.push_back(s[iS]);
        }
    }


    for( UInt32 iS=0; iS<support.size(); iS++ )
    {
        for( UInt32 j=support[iS].GetBegin(); j<support[iS].GetEnd(); j++ )
        {
            numDevs++;
            diff = (Yspc[j] - Ymodel[j]);
            fit += (diff*diff) / (m_ErrorNoContinuum[j]*m_ErrorNoContinuum[j]);
        }
    }
    return fit;
}

/**
 * \brief Returns the error of the support for subelements under the element with the argument eltId as index.
 * Accumulate "fit", the squared difference between model and spectrum, divied by the square of the m_ErrorNoContinuum value.
 * Accumulate "sumErr" 1 / square of the m_ErrorNoContinuum value.
 * return the square root of fit / sumErr.
 **/
Float64 CLineModelElementList::getModelErrorUnderElement( UInt32 eltId )
{
    if(eltId<0)
    {
        return -1.0;
    }
    const CSpectrumFluxAxis& spcFluxAxis = m_SpcFluxAxis;
    const CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();

    Int32 numDevs = 0;
    Float64 fit = 0.0;
    const Float64* Ymodel = modelFluxAxis.GetSamples();
    const Float64* Yspc = spcFluxAxis.GetSamples();
    Float64 diff = 0.0;

    Float64 sumErr=0.0;

    TInt32RangeList support;
    UInt32 iElts=eltId;
    {
        if(m_Elements[iElts]->IsOutsideLambdaRange()){
            return 0.0;
        }
        TInt32RangeList s = m_Elements[iElts]->getSupport();
        for( UInt32 iS=0; iS<s.size(); iS++ )
        {
            support.push_back(s[iS]);
        }
    }


    Float64 w=0.0;
    for( UInt32 iS=0; iS<support.size(); iS++ )
    {
        for( UInt32 j=support[iS].GetBegin(); j<support[iS].GetEnd(); j++ )
        {
            numDevs++;
            diff = (Yspc[j] - Ymodel[j]);
            w = 1.0 / (m_ErrorNoContinuum[j]*m_ErrorNoContinuum[j]);
            fit += (diff*diff) * w;
            sumErr += w;
        }
    }
    return sqrt(fit/sumErr);
}


/**
 * \brief Returns the Stronger Multiple Emission Lines Amplitude Coefficient (SMELAC)
 * 1. retrieve the lines amplitudes list for the Strong, and the Weak lines
 * 2. TODO: estimate the coefficient to be used as prior to penalise solutions with less Strong lines
 **/
Float64 CLineModelElementList::getStrongerMultipleELAmpCoeff()
{
    TFloat64List AmpsStrong;
    TFloat64List AmpsWeak;
    Float64 sumAmps = 0.0;

    //Retrieve all the lines amplitudes in two lists (1 Strong, 1 weak)
    std::vector<UInt32> validEltsIdx = GetModelValidElementsIndexes();
    for( UInt32 iValidElts=0; iValidElts<validEltsIdx.size(); iValidElts++ )
    {
        Int32 iElts = validEltsIdx[iValidElts];
        UInt32 nlines =  m_Elements[iElts]->GetRays().size();
        for(UInt32 lineIdx=0; lineIdx<nlines; lineIdx++)
        {
            if( !m_RestRayList[m_Elements[iElts]->m_LineCatalogIndexes[lineIdx]].GetIsEmission() ){
                continue;
            }

            Float64 amp = m_Elements[iElts]->GetFittedAmplitude(lineIdx);
            sumAmps += amp;
            if( m_RestRayList[m_Elements[iElts]->m_LineCatalogIndexes[lineIdx]].GetIsStrong() )
            {
                AmpsStrong.push_back(amp);
            }
            else
            {
                AmpsWeak.push_back(amp);
            }
        }
    }

    Float64 sumAmpsStrong = 0.0;
    for(UInt32 k=0; k<AmpsStrong.size(); k++)
    {
        sumAmpsStrong += AmpsStrong[k];
    }

    return sumAmpsStrong;
}

/**
 * \brief Returns the cumulative SNR under the Strong Emission Lines
 * 1. retrieve the lines support
 * 2. process each
 **/
Float64 CLineModelElementList::getCumulSNRStrongEL()
{

  //Retrieve all the liens supports in a list of range
    TInt32RangeList supportList;
    std::vector<UInt32> validEltsIdx = GetModelValidElementsIndexes();
    for( UInt32 iValidElts=0; iValidElts<validEltsIdx.size(); iValidElts++ )
    {
        Int32 iElts = validEltsIdx[iValidElts];
        UInt32 nlines =  m_Elements[iElts]->GetRays().size();
        for(UInt32 lineIdx=0; lineIdx<nlines; lineIdx++)
        {
            if( !m_RestRayList[m_Elements[iElts]->m_LineCatalogIndexes[lineIdx]].GetIsStrong() ){
                continue;
            }
            if( !m_RestRayList[m_Elements[iElts]->m_LineCatalogIndexes[lineIdx]].GetIsEmission() ){
                continue;
            }

            TInt32Range support = m_Elements[iElts]->getTheoreticalSupportSubElt(lineIdx);
            supportList.push_back(support);
        }
    }

    //merge overlapping ranges
    TInt32RangeList nonOverlappingSupportList;
    std::vector<Int32> processedSupport;
    for(UInt32 k=0; k<supportList.size(); k++)
    {
        //skip if already fitted
        bool alreadyProcessed=false;
        for(Int32 i=0; i<processedSupport.size(); i++)
        {
            if(k == processedSupport[i])
            {
                alreadyProcessed=true;
                break;
            }
        }
        if(alreadyProcessed)
        {
            continue;
        }


        processedSupport.push_back(k);
        TInt32Range support = TInt32Range( supportList[k].GetBegin(), supportList[k].GetEnd());

        for(UInt32 l=0; l<supportList.size(); l++)
        {
            //skip if already fitted
            bool alreadyProcessed=false;
            for(Int32 i=0; i<processedSupport.size(); i++)
            {
                if(l == processedSupport[i])
                {
                    alreadyProcessed=true;
                    break;
                }
            }
            if(alreadyProcessed)
            {
                continue;
            }

            //try if current range is bluer than l and overlaps ?
            Float64 xinf = support.GetBegin();
            Float64 xsup = support.GetEnd();
            Float64 yinf = supportList[l].GetBegin();
            Float64 ysup = supportList[l].GetEnd();
            Float64 max = std::max(xinf,yinf);
            Float64 min = std::min(xsup,ysup);
            if( max-min < 0 ){
                processedSupport.push_back(l);
                support.SetBegin(std::min(xinf,yinf));
                support.SetEnd(std::max(xsup,ysup));
            }
        }
        nonOverlappingSupportList.push_back(support);
    }

    TFloat64List snrList;
    Float64 sumSNR = 0.0;
    //process SNR on the non overlapping ranges
    for(UInt32 k=0; k<nonOverlappingSupportList.size(); k++)
    {
        snrList.push_back(getCumulSNROnRange(nonOverlappingSupportList[k]));
        sumSNR += snrList[k];
    }
    std::sort( snrList.rbegin(), snrList.rend() );

    //compute the snr metric
    std::vector<Int32> snrIsrelevantList;
    for(UInt32 k=0; k<snrList.size(); k++)
    {
        snrIsrelevantList.push_back(0);
    }
    Float64 thresRatio = 0.8;
    Float64 curRatio = 0.0;
    for(UInt32 k=0; k<snrList.size(); k++)
    {
        // relevant if snr>8.0
        if(snrList[k]>8.0)
        {
            snrIsrelevantList[k]=1;
        }

        // relevant if contributes to the 'thresRatio'*100 percent (ex. 80%) of the SumSNR value
        curRatio += snrList[k]/sumSNR;
        if(curRatio<=thresRatio)
        {
            snrIsrelevantList[k]=1;
        }
    }

    Float64 sumIsRelevant = 0.0;
    for(UInt32 k=0; k<snrIsrelevantList.size(); k++)
    {
        sumIsRelevant += snrIsrelevantList[k];
    }

    Float64 snrMetric = sumSNR*sumIsRelevant;
    return snrMetric;
}

/**
 * @brief CLineModelElementList::GetModelStrongLinePresent
 * @return 1 if there is 1 strong emission line present
 */
bool CLineModelElementList::GetModelStrongEmissionLinePresent()
{
    bool isStrongPresent = false;

    std::vector<UInt32> validEltsIdx = GetModelValidElementsIndexes();
    for( UInt32 iValidElts=0; iValidElts<validEltsIdx.size(); iValidElts++ )
    {
        Int32 iElts = validEltsIdx[iValidElts];
        UInt32 nlines =  m_Elements[iElts]->GetRays().size();
        for(UInt32 lineIdx=0; lineIdx<nlines; lineIdx++)
        {
            if( !m_RestRayList[m_Elements[iElts]->m_LineCatalogIndexes[lineIdx]].GetIsEmission() )
            {
                continue;
            }
            if( !m_RestRayList[m_Elements[iElts]->m_LineCatalogIndexes[lineIdx]].GetIsStrong() )
            {
                continue;
            }

            Float64 amp = m_Elements[iElts]->GetFittedAmplitude(lineIdx);
            if( amp>0.0)
            {
                isStrongPresent = true;
                Log.LogDebug("    model: GetModelStrongEmissionLinePresent - found Strong EL: %s", m_RestRayList[m_Elements[iElts]->m_LineCatalogIndexes[lineIdx]].GetName().c_str());
                break;
            }
        }
    }

    return isStrongPresent;
}

/**
 * \brief Returns the cumulative SNR on the idxRange
 **/
Float64 CLineModelElementList::getCumulSNROnRange( TInt32Range idxRange  )
{
    Int32 n = idxRange.GetEnd() - idxRange.GetBegin() +1;
    if(n<2){
        return -1;
    }

    const CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();
    const Float64* Ymodel = modelFluxAxis.GetSamples();

    Int32 idx = 0;
    Float64 maxAmp = 0.0;
    Float64 sumM = 0.0;
    for(Int32 i = 0; i < n; i++)
    {
        idx = i + idxRange.GetBegin();
        Float64 amp = Ymodel[idx]-m_ContinuumFluxAxis[idx]; //using only the no-continuum component to estimate SNR
        if(maxAmp<amp)
        {
            maxAmp = amp;
        }
        sumM +=  m_ErrorNoContinuum[idx];
    }
    Float64 Err = sumM/Float64(n);
    Float64 rangeSNR = maxAmp/Err;

    return rangeSNR;
}


Float64 CLineModelElementList::getContinuumMeanUnderElement(UInt32 eltId)
{
    Int32 n = 0;
    Float64 m=0.0;
    //Float64 sumErr=0.0;

    TInt32RangeList support;
    UInt32 iElts=eltId;
    {
        if(m_Elements[iElts]->IsOutsideLambdaRange()){
            return 0.0;
        }
        TInt32RangeList s = m_Elements[iElts]->getSupport();
        for( UInt32 iS=0; iS<s.size(); iS++ )
        {
            support.push_back(s[iS]);
        }
    }


    for( UInt32 iS=0; iS<support.size(); iS++ )
    {
        for( UInt32 j=support[iS].GetBegin(); j<support[iS].GetEnd(); j++ )
        {
            n++;
            //w = 1.0 / m_ErrorNoContinuum[j];
            //sumErr += w;
            //m += m_ContinuumFluxAxis[j] * w;
            m += m_ContinuumFluxAxis[j];
        }
    }


    return m/Float64(n);
}

/**
 * \brief Search the line catalog for lines whose name match the argument strTag.
 **/
std::vector<UInt32> CLineModelElementList::findLineIdxInCatalog(const CRayCatalog::TRayVector& restRayList, std::string strTag, Int32 type)
{
    std::vector<UInt32> indexes;
    for( UInt32 iRestRay=0; iRestRay<restRayList.size(); iRestRay++ )
    {
        if(restRayList[iRestRay].GetType() != type){
            continue;
        }
        std::string name = restRayList[iRestRay].GetName();
        std::size_t foundstra = name.find(strTag.c_str());
        if (foundstra!=std::string::npos){
            indexes.push_back(iRestRay);
        }
    }
    return indexes;
}

/**
 * \brief Adds an entry to m_Elements as a CMultiLine constructed from the arguments.
 **/
void CLineModelElementList::addDoubleLine(const CRay &r1, const CRay &r2, Int32 index1, Int32 index2, Float64 nominalWidth, Float64 a1, Float64 a2)
{
    std::vector<CRay> lines;
    lines.push_back(r1);
    lines.push_back(r2);
    std::vector<Float64> amps;
    amps.push_back(a1);
    amps.push_back(a2);
    std::vector<UInt32> a;
    a.push_back(index1);
    a.push_back(index2);
    m_Elements.push_back(boost::shared_ptr<CLineModelElement> (new CMultiLine(lines, m_LineWidthType, m_resolution, m_velocityEmission, m_velocityAbsorption, amps, nominalWidth, a)));
}

/**
 * /brief Calls the rules' methods depending on the JSON options.
 * If m_rulesoption is "no", do nothing.
 * If either "balmer" or "all" is in the rules string, call ApplyBalmerRuleLinSolve.
 * If "all" or "ratiorange" is in the rules string, call ApplyAmplitudeRatioRangeRule parameterized for OII.
 * If "all" or "strongweak" is in the rules string, call ApplyStrongHigherWeakRule for emission and then for absorption.
 **/
void CLineModelElementList::applyRules( bool enableLogs )
{
  if( m_rulesoption=="no" )
    {
      return;
    }

  m_Regulament-> EnableLogs(enableLogs);
  m_Regulament->Apply( *this );

}

TStringList CLineModelElementList::GetModelRulesLog()
{
  return m_Regulament->GetLogs();
}

/*
Reset all the model value to the previous solution found.
return 0 if every was ok; else -1
*/
Int32 CLineModelElementList::LoadModelSolution(const CLineModelSolution&  modelSolution){
  if(m_RestRayList.size() != modelSolution.Rays.size() ){
    Log.LogError("Unable to load the model solution because m_restRayList has not the same size as the number of ray in  modelSolution");
    return -1;
  }
  for( UInt32 iRestRay=0; iRestRay<m_RestRayList.size(); iRestRay++ )
  {
    Int32 eIdx = modelSolution.ElementId[iRestRay];
    Int32 subeIdx = m_Elements[eIdx]->FindElementIndex(iRestRay);

    if(m_RestRayList[iRestRay] !=modelSolution.Rays[iRestRay]){
      Log.LogError("ray in m_restRest doesn't correspond to the modelSolution ray for index %d", iRestRay);
    }
    if(eIdx != -1 and subeIdx==-1){
      m_Elements[eIdx]->SetFittedAmplitude(subeIdx,modelSolution.Amplitudes[iRestRay], -1.0);
    }else{
      m_Elements[eIdx]->SetFittedAmplitude(subeIdx,modelSolution.Amplitudes[iRestRay], modelSolution.Errors[iRestRay]);
    }
    //Log.LogInfo("LoadModelSolution eidx %d subeIdx %d, iRestRay %d, Amp %f", eIdx, subeIdx, iRestRay,modelSolution.Amplitudes[iRestRay]);
  }

  if(modelSolution.LyaWidthCoeff != -1.0 or modelSolution.LyaAlpha !=-1.0 or modelSolution.LyaDelta !=-1.0){
    linetags ltags;
    std::string lyaTag = ltags.lya_em;
    Int32 idxLyaE = FindElementIndex(lyaTag);
    if( idxLyaE>-1 )
    {
        m_Elements[idxLyaE]-> SetAsymfitWidthCoeff(modelSolution.LyaWidthCoeff );
        m_Elements[idxLyaE]-> SetAsymfitAlphaCoeff(modelSolution.LyaAlpha );
        m_Elements[idxLyaE]-> SetAsymfitDelta(modelSolution.LyaDelta);
    }
  }
  SetVelocityEmission(modelSolution.EmissionVelocity);
  SetVelocityAbsorption(modelSolution.AbsorptionVelocity);
  setRedshift(modelSolution.Redshift, false);
  return 0;
}

// return error: 1=can't find element index, 2=Abs_width not high enough compared to Em_width
Int32 CLineModelElementList::improveBalmerFit()
{
    linetags ltags;
    //Emission Balmer lines
    std::vector<std::string> linetagsE;
    linetagsE.push_back( ltags.halpha_em );
    linetagsE.push_back( ltags.hbeta_em );
    linetagsE.push_back( ltags.hgamma_em );
    linetagsE.push_back( ltags.hdelta_em );
    //Absorption Balmer lines
    std::vector<std::string> linetagsA;
    linetagsA.push_back( ltags.halpha_abs );
    linetagsA.push_back( ltags.hbeta_abs );
    linetagsA.push_back( ltags.hgamma_abs );
    linetagsA.push_back( ltags.hdelta_abs );
    //Additional lines to be fitted with the Balmer lines, WARNING: only EMISSION for now !!
    std::vector<std::string> linetagsNII;
    linetagsNII.push_back( ltags.niia_em );
    linetagsNII.push_back( ltags.niib_em );
    std::vector<std::string> linetagsVoid;
    std::vector<TStringList> linetagsMore;
    linetagsMore.push_back( linetagsNII );
    linetagsMore.push_back( linetagsVoid );
    linetagsMore.push_back( linetagsVoid );
    linetagsMore.push_back( linetagsVoid );

    if(linetagsE.size()!=linetagsA.size() || linetagsE.size()!=linetagsMore.size())
    {
        return -1;
    }

    for( Int32 itag=0; itag<linetagsE.size(); itag++)
    {
        std::string tagE = linetagsE[itag];
        std::string tagA = linetagsA[itag];

        Int32 ilineE = FindElementIndex( tagE, CRay::nType_Emission );
        Int32 ilineA = FindElementIndex( tagA, CRay::nType_Absorption );
        // Were the lines indexes found ?
        if(ilineE<0 || ilineA<0)
        {
            continue;
        }
        // for now only allow this process if Em and Abs line are single lines
        if(m_Elements[ilineE]->m_Rays.size()>1 || m_Elements[ilineA]->m_Rays.size()>1)
        {
            continue;
        }
        Int32 subeIdxE = 0;
        Int32 subeIdxA = 0;

        //find the linesMore unique elements indexes
        std::vector<UInt32> ilinesMore;
        for( UInt32 imore=0; imore<linetagsMore[itag].size(); imore++)
        {
            std::string tagMore = linetagsMore[itag][imore];
            Int32 ilineMore = FindElementIndex( tagMore, CRay::nType_Emission );
            if(ilineMore<0)
            {
                continue;
            }
            ilinesMore.push_back(ilineMore);
        }
        std::sort(ilinesMore.begin(), ilinesMore.end());
        ilinesMore.erase( std::unique( ilinesMore.begin(), ilinesMore.end() ), ilinesMore.end() );
        for( Int32 imore=0; imore<ilinesMore.size(); imore++)
        {
            Log.LogDebug("    model: balmerImprove more tags = %d", ilinesMore[imore]);
        }

        //try if the width is significantly different: abs > em
        Float64 AbsVSEmWidthCoeffThreshold = 2.0;
        Float64 sigmaE = m_Elements[ilineE]->GetWidth(subeIdxE, m_Redshift);
        Float64 sigmaA = m_Elements[ilineA]->GetWidth(subeIdxA, m_Redshift);
        if(sigmaA<AbsVSEmWidthCoeffThreshold*sigmaE)
        {
            continue;
        }

        //simulatneous fit with linsolve
        Float64 modelErr_init = getModelErrorUnderElement(ilineA);
        Float64 ampA = m_Elements[ilineA]->GetFittedAmplitude(0);
        Float64 amp_errorA = m_Elements[ilineA]->GetFittedAmplitudeErrorSigma(0);
        Float64 ampE = m_Elements[ilineE]->GetFittedAmplitude(0);
        Float64 amp_errorE = m_Elements[ilineE]->GetFittedAmplitudeErrorSigma(0);
        TFloat64List ampsMore;
        TFloat64List ampErrorsMore;
        for( Int32 imore=0; imore<ilinesMore.size(); imore++)
        {
            Float64 amp = m_Elements[ilinesMore[imore]]->GetFittedAmplitude(0);
            Float64 ampErr = m_Elements[ilinesMore[imore]]->GetFittedAmplitudeErrorSigma(0);
            ampsMore.push_back(amp);
            ampErrorsMore.push_back(ampErr);
        }

        std::vector<UInt32> eltsIdx;
        eltsIdx.push_back(ilineA);
        eltsIdx.push_back(ilineE);
        for( Int32 imore=0; imore<ilinesMore.size(); imore++)
        {
            eltsIdx.push_back(ilinesMore[imore]);
        }
        std::vector<Float64> ampsfitted;
        std::vector<Float64> errorsfitted;
        fitAmplitudesLinSolve(eltsIdx, m_SpectrumModel->GetSpectralAxis(), m_spcFluxAxisNoContinuum, m_ContinuumFluxAxis, ampsfitted, errorsfitted);

        //decide if the fit is better than previous amps
        std::vector<UInt32> elts;
        elts.push_back(ilineA);
        elts.push_back(ilineE);
        for( UInt32 imore=0; imore<ilinesMore.size(); imore++)
        {
            elts.push_back(ilinesMore[imore]);
        }
        refreshModelUnderElements(elts);
        Float64 modelErr_withfit = getModelErrorUnderElement(ilineA);
        if(modelErr_withfit>modelErr_init)
        {
            Float64 nominal_ampA = m_Elements[ilineA]->GetNominalAmplitude(0);
            Float64 nominal_ampE = m_Elements[ilineE]->GetNominalAmplitude(0);
            m_Elements[ilineA]->SetFittedAmplitude(ampA/nominal_ampA, amp_errorA/nominal_ampA);
            m_Elements[ilineE]->SetFittedAmplitude(ampE/nominal_ampE, amp_errorE/nominal_ampE);
            for( Int32 imore=0; imore<ilinesMore.size(); imore++)
            {
                Float64 nominal_ampMore = m_Elements[ilinesMore[imore]]->GetNominalAmplitude(0);
                m_Elements[ilinesMore[imore]]->SetFittedAmplitude(ampsMore[imore]/nominal_ampMore, ampErrorsMore[imore]/nominal_ampMore);
            }
        }
    }

    return 0;
}



/**
 * \brief Returns a SLineModelSolution object populated with the current solutions.
 **/
CLineModelSolution CLineModelElementList::GetModelSolution(Int32 opt_level)
{
    CLineModelSolution modelSolution;
    modelSolution.nDDL = GetModelNonZeroElementsNDdl();
    modelSolution.Rays = m_RestRayList;
    modelSolution.ElementId.resize(m_RestRayList.size());
    modelSolution.snrHa = -1.0;
    modelSolution.lfHa = -1.0;
    modelSolution.snrOII = -1.0;
    modelSolution.lfOII = -1.0;

    TInt32List eIdx_oii;
    TInt32List subeIdx_oii;


    for( UInt32 iRestRay=0; iRestRay<m_RestRayList.size(); iRestRay++ )
    {
        Int32 eIdx = FindElementIndex(iRestRay);
        Int32 subeIdx = m_Elements[eIdx]->FindElementIndex(iRestRay);
        modelSolution.ElementId[iRestRay] = eIdx;

        if(eIdx==-1 || subeIdx==-1 || m_Elements[eIdx]->IsOutsideLambdaRange(subeIdx))
        {
            modelSolution.Amplitudes.push_back(m_Elements[eIdx]->GetFittedAmplitude(subeIdx));
            modelSolution.Errors.push_back(-1.0);
            modelSolution.FittingError.push_back(-1.0);
            modelSolution.LambdaObs.push_back(-1.0);
            modelSolution.Offset.push_back(-1.0);
            modelSolution.Velocity.push_back(-1.0);
            modelSolution.CenterContinuumFlux.push_back(-1.0);
            modelSolution.ContinuumError.push_back(-1.0);
            modelSolution.Sigmas.push_back(-1.0);
            modelSolution.Fluxs.push_back(-1.0);
            modelSolution.FluxErrors.push_back(-1.0);
            modelSolution.FluxDirectIntegration.push_back(-1.0);
            modelSolution.OutsideLambdaRange.push_back(true);
            modelSolution.fittingGroupInfo.push_back("-1");
        }else{
            Float64 amp = m_Elements[eIdx]->GetFittedAmplitude(subeIdx);
            modelSolution.Amplitudes.push_back(amp);
            Float64 ampError = m_Elements[eIdx]->GetFittedAmplitudeErrorSigma(subeIdx);
            modelSolution.Errors.push_back(ampError);

            modelSolution.LambdaObs.push_back(m_Elements[eIdx]->GetObservedPosition(subeIdx, m_Redshift));
            modelSolution.Velocity.push_back(m_Elements[eIdx]->GetVelocity());
            modelSolution.Offset.push_back(m_Elements[eIdx]->m_Rays[subeIdx].GetOffset());

            if(opt_level!=0)// brief, to save processing time, do not estimate fluxes and high level line properties
            {
                modelSolution.FittingError.push_back(getModelErrorUnderElement(eIdx));
                Float64 cont = m_Elements[eIdx]->GetContinuumAtCenterProfile(subeIdx, m_SpectrumModel->GetSpectralAxis(), m_Redshift, m_ContinuumFluxAxis);
                modelSolution.CenterContinuumFlux.push_back(cont);
                modelSolution.ContinuumError.push_back(GetContinuumError(eIdx, subeIdx));
                Float64 sigma = m_Elements[eIdx]->GetWidth(subeIdx, m_Redshift);
                Float64 flux = -1;
                Float64 fluxError = -1;
                Float64 fluxDI = -1;
                Float64 snrDI = -1;
                TInt32List eIdx_line(1, eIdx);
                TInt32List subeIdx_line(1, subeIdx);
                Int32 retdi = GetFluxDirectIntegration(eIdx_line, subeIdx_line, fluxDI, snrDI);
                if(amp>=0.0)
                {
                    if(m_RestRayList[iRestRay].GetType()==CRay::nType_Emission)
                    {
                        flux = m_Elements[eIdx]->GetLineFlux(m_RestRayList[iRestRay].GetProfile(), sigma, amp);
                        fluxError = m_Elements[eIdx]->GetLineFlux(m_RestRayList[iRestRay].GetProfile(), sigma, ampError);
                    }else{
                        Float64 _amp = cont*amp;
                        flux = -m_Elements[eIdx]->GetLineFlux(m_RestRayList[iRestRay].GetProfile(), sigma, _amp);
                        fluxError = m_Elements[eIdx]->GetLineFlux(m_RestRayList[iRestRay].GetProfile(), sigma, ampError); //to be checked
                    }
                }
                modelSolution.Sigmas.push_back(sigma);
                modelSolution.Fluxs.push_back(flux);
                modelSolution.FluxErrors.push_back(fluxError);
                modelSolution.FluxDirectIntegration.push_back(fluxDI);

                //rough estimation of SNR_Ha, using the given model and fitting method
                //(warning: Ha flux and error could have been obtained by global fitting of the model, which leads to different results than fitted individually...)
                linetags ltags;
                if(m_RestRayList[iRestRay].GetName()==ltags.halpha_em && m_RestRayList[iRestRay].GetType()==CRay::nType_Emission)
                {
                    if(false)
                    {
                        if(fluxError>0.0)
                        {
                            modelSolution.snrHa = flux/fluxError;
                        }
                        if(flux>0.0)
                        {
                            modelSolution.lfHa = log10(flux);
                        }
                    }else{
                        modelSolution.snrHa = snrDI;
                        if(fluxDI>0.0)
                        {
                            modelSolution.lfHa = log10(fluxDI);
                        }
                    }
                }
                if((m_RestRayList[iRestRay].GetName()==ltags.oII3726_em || m_RestRayList[iRestRay].GetName()==ltags.oII3729_em)
                        && m_RestRayList[iRestRay].GetType()==CRay::nType_Emission)
                {
                    //here we only cover the fluxDI case.
                    eIdx_oii.push_back(eIdx);
                    subeIdx_oii.push_back(subeIdx);
                    if(true)
                    {
                        fluxDI = -1;
                        Float64 snrDI = -1;
                        Int32 retdi = GetFluxDirectIntegration(eIdx_oii, subeIdx_oii, fluxDI, snrDI);

                        modelSolution.snrOII = snrDI;
                        if(fluxDI>0.0)
                        {
                            modelSolution.lfOII = log10(fluxDI);
                        }
                    }
                }
            }
            else{
                modelSolution.FittingError.push_back(-1.0);
                modelSolution.CenterContinuumFlux.push_back(-1.0);
                modelSolution.ContinuumError.push_back(-1.0);
                modelSolution.Sigmas.push_back(-1.0);
                modelSolution.Fluxs.push_back(-1.0);
                modelSolution.FluxErrors.push_back(-1.0);
                modelSolution.FluxDirectIntegration.push_back(-1.0);
            }

            modelSolution.fittingGroupInfo.push_back(m_Elements[eIdx]->m_fittingGroupInfo);
            modelSolution.OutsideLambdaRange.push_back(m_Elements[eIdx]->IsOutsideLambdaRange(subeIdx));
        }
    }
    //retrieve Lya params if fitted
    modelSolution.LyaWidthCoeff = -1.0;
    modelSolution.LyaAlpha = -1.0;
    modelSolution.LyaDelta = -1.0;
    linetags ltags;
    std::string lyaTag = ltags.lya_em;
    Int32 idxLyaE = FindElementIndex(lyaTag);
    if( idxLyaE>-1 )
    {
        modelSolution.LyaWidthCoeff = m_Elements[idxLyaE]->GetAsymfitWidthCoeff();
        modelSolution.LyaAlpha = m_Elements[idxLyaE]->GetAsymfitAlphaCoeff();
        modelSolution.LyaDelta = m_Elements[idxLyaE]->GetAsymfitDelta();
    }
    modelSolution.EmissionVelocity = m_velocityEmission;
    modelSolution.AbsorptionVelocity = m_velocityAbsorption;
    modelSolution.Redshift = m_Redshift;
    return modelSolution;
}


CContinuumModelSolution CLineModelElementList::GetContinuumModelSolution()
{
    CContinuumModelSolution continuumModelSolution;

    continuumModelSolution.tplName = m_fitContinuum_tplName;
    continuumModelSolution.tplDustCoeff = m_fitContinuum_tplFitDustCoeff;
    continuumModelSolution.tplMeiksinIdx = m_fitContinuum_tplFitMeiksinIdx;
    continuumModelSolution.tplAmplitude = m_fitContinuum_tplFitAmplitude;
    continuumModelSolution.tplMerit = m_fitContinuum_tplFitMerit;
    continuumModelSolution.tplDtm = m_fitContinuum_tplFitDtM;
    continuumModelSolution.tplMtm = m_fitContinuum_tplFitMtM;
    continuumModelSolution.tplRedshift = m_fitContinuum_tplFitRedshift;
    continuumModelSolution.pCoeffs = m_fitContinuum_tplFitPolyCoeffs;

    return continuumModelSolution;
}

/**
 * \brief Returns the size of m_Elements.
 **/
Int32 CLineModelElementList::GetNElements()
{
    Int32 nddl = m_Elements.size();
    return nddl;
}

/**
 * \brief Returns the number of m_Elements that fail IsOutsideLambdaRange().
 **/
Int32 CLineModelElementList::GetModelValidElementsNDdl()
{
    Int32 nddl = 0;
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        if(m_Elements[iElts]->IsOutsideLambdaRange() == true){
            continue;
        }

        nddl++;
    }
    return nddl;
}

/**
 * \brief Returns the number of elements that have only subelements with non-positive amplitude.
 **/
Int32 CLineModelElementList::GetModelNonZeroElementsNDdl()
{
    Int32 nddl = 0;
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        if(m_Elements[iElts]->IsOutsideLambdaRange() == true){
            continue;
        }
        bool isAllZero=true;
        for(Int32 ie=0; ie<m_Elements[iElts]->GetSize(); ie++){
            if(m_Elements[iElts]->GetFittedAmplitude(ie) > 0.0){
                isAllZero=false;
            }
        }

        if(isAllZero==false){
            nddl++;
        }
    }
    return nddl;
}

/**
 * \brief Returns the list of indexes of elements that fail IsOutsideLambdaRange.
 **/
std::vector<UInt32> CLineModelElementList::GetModelValidElementsIndexes()
{
    std::vector<UInt32> nonZeroIndexes;
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        if(m_Elements[iElts]->IsOutsideLambdaRange() == true){
            continue;
        }
        if(IsElementIndexInDisabledList(iElts))
        {
            continue;
        }

        nonZeroIndexes.push_back(iElts);
    }
    return nonZeroIndexes;
}

/**
 * \brief Returns the list of groups, with each group being a set of line indexes with the velcocity to be jointly
 * TEMPORARY-DEV: return all the indexes individually as  agroup
**/
std::vector<std::vector<Int32>> CLineModelElementList::GetModelVelfitGroups( Int32 lineType )
{
    Bool verbose = false;
    if(verbose)
    {
        Log.LogInfo("    model: group tags for lineType=%d", lineType);
    }
    std::vector<std::string> tags;
    std::vector<UInt32> nonGroupedLines;

    std::vector<UInt32> nonZeroIndexes = GetModelValidElementsIndexes();
    for(Int32 i=0; i<nonZeroIndexes.size(); i++)
    {
        Int32 iElts = nonZeroIndexes[i];
        Int32 nRays = m_Elements[iElts]->GetSize();
        for(Int32 iSubElts=0; iSubElts<nRays; iSubElts++)
        {
            if(verbose)
            {
                Log.LogInfo("    model: group tags - lineType=%d", lineType);
                Log.LogInfo("    model: group tags - m_Elements[iElts]->m_Rays[iSubElts].GetType()=%d", m_Elements[iElts]->m_Rays[iSubElts].GetType());
            }
            if(lineType == m_Elements[iElts]->m_Rays[iSubElts].GetType())
            {
                std::string _tag = m_Elements[iElts]->m_Rays[iSubElts].GetVelGroupName();
                if(_tag != "-1"){
                    tags.push_back(_tag);
                }else{
                    nonGroupedLines.push_back(iElts);
                }
            }
        }
    }

    if(verbose)
    {
        Log.LogInfo("    model: group tags non unique found=%d", tags.size());
        for( Int32 itag = 0; itag<tags.size(); itag++){
            Log.LogInfo("    model: non unique tag %d/%d = %s", itag+1, tags.size(), tags[itag].c_str());
        }
    }

    // create the group tag set by removing duplicates
    std::sort( tags.begin(), tags.end() );
    tags.erase( std::unique( tags.begin(), tags.end() ), tags.end() );
    //*
    if(verbose)
    {
        //print the tags
        for( Int32 itag = 0; itag<tags.size(); itag++){
            Log.LogInfo("    model: velfit group Tag %d/%d = %s", itag+1, tags.size(), tags[itag].c_str());
        }
    }
    //*/

    //add the grouped lines
    std::vector<std::vector<Int32>> groups;
    std::vector<std::string> groupsTags;
    for( Int32 itag = 0; itag<tags.size(); itag++)
    {
        std::vector<Int32> _group;
        for(Int32 i=0; i<nonZeroIndexes.size(); i++)
        {
            Int32 iElts = nonZeroIndexes[i];
            Int32 nRays = m_Elements[iElts]->GetSize();
            for(Int32 iSubElts=0; iSubElts<nRays; iSubElts++)
            {
                if(lineType == m_Elements[iElts]->m_Rays[iSubElts].GetType())
                {
                    std::string _tag = m_Elements[iElts]->m_Rays[iSubElts].GetVelGroupName();
                    if(_tag == tags[itag]){
                        _group.push_back(iElts);
                    }
                }
            }
        }
        //add the grouped lines, no duplicates
        std::sort( _group.begin(), _group.end() );
        _group.erase( std::unique( _group.begin(), _group.end() ), _group.end() );

        groups.push_back(_group);
        groupsTags.push_back(tags[itag]);
    }
    //add the non grouped lines, no duplicates
    std::sort( nonGroupedLines.begin(), nonGroupedLines.end() );
    nonGroupedLines.erase( std::unique( nonGroupedLines.begin(), nonGroupedLines.end() ), nonGroupedLines.end() );
    for( Int32 i = 0; i<nonGroupedLines.size(); i++)
    {
        std::vector<Int32> _group;
        _group.push_back(nonGroupedLines[i]);
        groups.push_back(_group);
        groupsTags.push_back("-1");
    }

    if(true)
    {
        //print the groups
        for( Int32 igr = 0; igr<groups.size(); igr++){
            Log.LogInfo("    model: Group %d/%d: nlines=%d, tag=%s", igr+1, groups.size(), groups[igr].size(), groupsTags[igr].c_str());
            for(Int32 i=0; i<groups[igr].size(); i++)
            {
                Log.LogInfo("    model: \t%d: iElt=%d", i+1, groups[igr][i]);
            }
        }
    }

    //Override velGroups from Catalog: => Individual lines as groups
    /*
    std::vector<std::vector<Int32>> groups;
    std::vector<Int32> nonZeroIndexes = GetModelValidElementsIndexes();
    for(Int32 i=0; i<nonZeroIndexes.size(); i++)
    {
        if(lineType == m_Elements[nonZeroIndexes[i]]->m_Rays[0].GetType())
        {
            std::vector<Int32> gr;
            gr.push_back(nonZeroIndexes[i]);
            groups.push_back(gr);
            //Log.LogInfo("Group %d, idx=%d", groups.size(), groups[groups.size()-1][0]);
        }
    }
    //*/

    return groups;
}

bool CLineModelElementList::IsElementIndexInDisabledList(Int32 index)
{
    for( UInt32 i=0; i<m_elementsDisabledIndexes.size(); i++ )
    {
        if( m_elementsDisabledIndexes[i]== index){
            return true;
        }
    }
    return false;
}

/**
 * @brief CLineModelElementList::SetElementIndexesDisabledAuto
 * Disables all the elements that have all sub-elements (lines) amplitudes equal to zero
 */
void CLineModelElementList::SetElementIndexesDisabledAuto()
{
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        if(m_Elements[iElts]->IsOutsideLambdaRange() == true){
            continue;
        }
        bool isAllZero=true;
        for(Int32 ie=0; ie<m_Elements[iElts]->GetSize(); ie++){
            if(m_Elements[iElts]->GetFittedAmplitude(ie) > 0.0){
                isAllZero=false;
            }
        }

        if(isAllZero==true){
            m_elementsDisabledIndexes.push_back(iElts);
        }
    }
}

void CLineModelElementList::ResetElementIndexesDisabled()
{
    m_elementsDisabledIndexes.clear();
}


/**
 * \brief Returns the first index of m_Elements where calling the element's FindElementIndex method with LineCatalogIndex argument does not return -1.
 **/
Int32 CLineModelElementList::FindElementIndex(Int32 LineCatalogIndex)
{
    Int32 idx = -1;
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        if(m_Elements[iElts]->FindElementIndex(LineCatalogIndex) !=-1){
            idx = iElts;
            break;
        }
    }
    return idx;
}

/**
 * \brief Returns the first index of m_Elements where calling the element's FindElementIndex method with LineTagStr argument does not return -1.
 **/
Int32 CLineModelElementList::FindElementIndex(std::string LineTagStr, Int32 linetype, Int32& lineIdx )
{
    Int32 idx = -1;
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        lineIdx = m_Elements[iElts]->FindElementIndex(LineTagStr) ;
        if( lineIdx!=-1 ){
            if( linetype!=-1 )
            {
                if(m_RestRayList[m_Elements[iElts]->m_LineCatalogIndexes[lineIdx]].GetType() != linetype){
                    continue;
                }
            }
            idx = iElts;
            break;
        }
    }
    return idx;
}

/**
 * \brief If argument j is a valid index of m_Elements, updates the element in that index calling its SetFittedAmplitude with arguments a and snr.
 **/
void CLineModelElementList::SetElementAmplitude(Int32 j, Float64 a, Float64 snr)
{
    if(j>=0 && j<m_Elements.size())
    {
        m_Elements[j]->SetFittedAmplitude(a, snr);
    }
    return;
}

/**
 * \brief If j is a valid index of m_Elements, returns a call to that element's GetElementAmplitude. If not, returns -1.
 **/
Float64 CLineModelElementList::GetElementAmplitude(Int32 j)
{
    Float64 a=-1.0;
    if(j>=0 && j<m_Elements.size())
    {
        a = m_Elements[j]->GetElementAmplitude();
    }
    return a;
}


void CLineModelElementList::SetSourcesizeDispersion(Float64 sizeArcsec)
{
    for(Int32 j=0; j<m_Elements.size(); j++)
    {
        m_Elements[j]->SetSourcesizeDispersion(sizeArcsec);
    }
}

void CLineModelElementList::SetVelocityEmission(Float64 vel)
{
    m_velocityEmission = vel;
    for(Int32 j=0; j<m_Elements.size(); j++)
    {
        m_Elements[j]->SetVelocityEmission(vel);
    }
}

void CLineModelElementList::SetVelocityEmissionOneElement(Float64 vel, Int32 idxElt)
{
    m_velocityEmission = vel;
    if(idxElt<m_Elements.size())
    {
        m_Elements[idxElt]->SetVelocityEmission(vel);
    }
}

void CLineModelElementList::SetVelocityAbsorption(Float64 vel)
{
    m_velocityAbsorption = vel;
    for(Int32 j=0; j<m_Elements.size(); j++)
    {
        m_Elements[j]->SetVelocityAbsorption(vel);
    }
}

void CLineModelElementList::SetVelocityAbsorptionOneElement(Float64 vel, Int32 idxElt)
{
    m_velocityAbsorption = vel;
    if(idxElt<m_Elements.size())
    {
        m_Elements[idxElt]->SetVelocityAbsorption(vel);
    }
}

Float64 CLineModelElementList::GetVelocityEmission()
{
    return m_velocityEmission;
}
Float64 CLineModelElementList::GetVelocityAbsorption()
{
    return m_velocityAbsorption;
}

Float64 CLineModelElementList::GetVelocityInfFromInstrumentResolution()
{
    static Float64 c = 300000.0;
    static Float64 tolCoeff = 2.0;
    Float64 vel = c/m_resolution/tolCoeff;
    Float64 roundingVal = 10.0;
    vel = Float64(Int32(vel/roundingVal))*roundingVal;
//    if(vel<5.0){
//        vel=5.0;
//    }
    return vel;
}

Float64 CLineModelElementList::GetRedshift()
{
    return m_Redshift;
}

Int32 CLineModelElementList::ApplyVelocityBound(Float64 inf, Float64 sup)
{

    Int32 corrected=false;
    static Float64 velInfFromInstrument = inf;
    static Float64 velSupEmission = sup;
    static Float64 velSupAbsorption = sup;

    Float64 vel;
    vel = GetVelocityEmission();
    if(vel>velSupEmission || vel<velInfFromInstrument)
    {
        SetVelocityEmission(m_velocityEmissionInit);
        corrected = true;
        Log.LogInfo( "\nLineModel Infos: Reset Velocity Emission, to v = %.1f", m_velocityEmissionInit);
    }
    vel = GetVelocityAbsorption();
    if(vel>velSupAbsorption || vel<velInfFromInstrument)
    {
        SetVelocityAbsorption(m_velocityAbsorptionInit);
        corrected = true;
        Log.LogInfo( "\nLineModel Infos: Reset Velocity Absorption, to v = %.1f", m_velocityAbsorptionInit);
    }
    return corrected;
}


/**
 * \brief this function estimates the continuum after removal(interpolation) of the flux samples under the lines for a given redshift 
 * //todo: use lambdaRange in order to speed up the continuum estimation process. (rmq: use lambdarange with some margin in order to not detetriorate cont. estimation close to the borders)
 **/
void CLineModelElementList::EstimateSpectrumContinuum( Float64 opt_enhance_lines, const TFloat64Range& lambdaRange )
{
    std::vector<UInt32> validEltsIdx = GetModelValidElementsIndexes();
    //std::vector<UInt32> xInds = getSupportIndexes( validEltsIdx );
    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();

    //create new spectrum, which is corrected under the lines
    CSpectrum spcCorrectedUnderLines(*m_SpectrumModel);
    CSpectrumFluxAxis& fluxAxisNothingUnderLines = spcCorrectedUnderLines.GetFluxAxis();
    Float64* Y = fluxAxisNothingUnderLines.GetSamples();

    /*
    //1. interp from previous values
    for( Int32 t=0;t<spectralAxis.GetSamplesCount();t++)
    {
        Y[t] = m_SpcFluxAxis[t];
    }
    Int32 idx = 0;
    Float64 valf = m_SpcFluxAxis[0];
    Float64 corrRatio=0.8;
    for (int i = 0; i < xInds.size(); i++)
    {
        idx = xInds[i];
        if(idx>0){
            if ( std::find(xInds.begin(), xInds.end(), idx-1) == xInds.end() )
            {
                valf=m_SpcFluxAxis[idx-1];
            }
        }
        Y[idx]= corrRatio*valf + (1.0-corrRatio)*m_SpcFluxAxis[idx];
    }
    //*/

    //2. subtract lines from model
    //model for subtraction
    CSpectrum spcmodel4linefitting = GetModelSpectrum();
    for(UInt32 i=0; i<spcmodel4linefitting.GetFluxAxis().GetSamplesCount(); i++){
        spcmodel4linefitting.GetFluxAxis()[i] = spcmodel4linefitting.GetFluxAxis()[i]-m_ContinuumFluxAxis[i];
    }
    //optionnaly enhance the abs model component
    if(opt_enhance_lines>0.0){

        bool filterOnlyAbs = false;
        for( Int32 t=0;t<spectralAxis.GetSamplesCount();t++)
        {
            if(filterOnlyAbs)
            {
                if(spcmodel4linefitting.GetFluxAxis()[t]<0.0)
                {
                    spcmodel4linefitting.GetFluxAxis()[t] *= opt_enhance_lines;
                }
            }else{
                spcmodel4linefitting.GetFluxAxis()[t] *= opt_enhance_lines;
            }
        }
    }
    // subtract the lines component
    for( Int32 t=0;t<spectralAxis.GetSamplesCount();t++)
    {
        Y[t] = m_SpcFluxAxis[t]-spcmodel4linefitting.GetFluxAxis()[t];
    }

    /*
    // export for debug
    FILE* fspc = fopen( "continuum_correctedSpc_dbg.txt", "w+" );
    Float64 coeffSaveSpc = 1e16;
    for( Int32 t=0;t<spectralAxis.GetSamplesCount();t++)
    {
        fprintf( fspc, "%f %f %f\n", t, spectralAxis[t], (m_SpcFluxAxis[t])*coeffSaveSpc, (spcmodel4linefitting.GetFluxAxis()[t])*coeffSaveSpc);
    }
    fclose( fspc );
    //*/

    // TODO: use the continuum remover defined in the CSpectrum continuum member, with params defined in the smae place
    // Remove continuum
    CSpectrumFluxAxis fluxAxisWithoutContinuumCalc;
    if(1)
    {
        CContinuumIrregularSamplingMedian continuum;
        Float64 opt_medianKernelWidth = m_ContinuumWinsize;
        continuum.SetMedianKernelWidth(opt_medianKernelWidth);
        continuum.RemoveContinuum( spcCorrectedUnderLines, fluxAxisWithoutContinuumCalc );
    }else
    {
        Int64 nscales = 6;
        std::string dfBinPath="/home/aschmitt/gitlab/amazed/extern/df_linux/";
        CContinuumDF continuum(dfBinPath);
        spcCorrectedUnderLines.SetDecompScales(nscales);
        continuum.RemoveContinuum( spcCorrectedUnderLines, fluxAxisWithoutContinuumCalc );
    }



    CSpectrumFluxAxis fluxAxisNewContinuum;
    fluxAxisNewContinuum.SetSize( fluxAxisNothingUnderLines.GetSamplesCount() );

    for( Int32 t=0;t<spectralAxis.GetSamplesCount();t++)
    {
        fluxAxisNewContinuum[t] = fluxAxisNothingUnderLines[t];
    }
    fluxAxisNewContinuum.Subtract(fluxAxisWithoutContinuumCalc);

    /*
    // export for debug
    FILE* f = fopen( "continuum_estimated_dbg.txt", "w+" );
    Float64 coeffSave = 1e16;
    for( Int32 t=0;t<spectralAxis.GetSamplesCount();t++)
    {
        fprintf( f, "%f %f %f\n", t, spectralAxis[t], (m_SpcFluxAxis[t])*coeffSave, (fluxAxisNewContinuum[t])*coeffSave);
    }
    fclose( f );
    //*/

    /*
    // export for debug
    FILE* f = fopen( "continuumfree_estimated_dbg.txt", "w+" );
    Float64 coeffSave = 1e16;
    for( Int32 t=0;t<spectralAxis.GetSamplesCount();t++)
    {
        fprintf( f, "%f %f %f\n", t, spectralAxis[t], (m_SpcFluxAxis[t]-m_SpcContinuumFluxAxis[t])*coeffSave, (m_SpcFluxAxis[t]-fluxAxisNewContinuum[t])*coeffSave);
    }
    fclose( f );
    //*/

//    //modify m_SpcFluxAxis
//    CSpectrumFluxAxis& fluxAxisModified = m_SpcFluxAxis;
//    Float64* Y2 = fluxAxisModified.GetSamples();
//    for( Int32 t=0;t<spectralAxis.GetSamplesCount();t++)
//    {
//        Y2[t] = m_SpcFluxAxis[t]-fluxAxisNewContinuum[t];
//        //Y2[t] = m_SpcFluxAxis[t]-m_SpcContinuumFluxAxis[t];
//    }

    //modify m_ContinuumFluxAxis
    CSpectrumFluxAxis& fluxAxisModified = m_ContinuumFluxAxis;
    Float64* Y2 = fluxAxisModified.GetSamples();
    for( Int32 t=0;t<spectralAxis.GetSamplesCount();t++)
    {
        Y2[t] = fluxAxisNewContinuum[t];
        //Y2[t] = m_SpcContinuumFluxAxis[t];
    }
}


/**
 * \brief this function returns the dtd value withing the wavelength range for a given spcComponent
 *
 **/
Float64 CLineModelElementList::getDTransposeD(const TFloat64Range& lambdaRange, std::string spcComponent)
{
    if(! (m_dTransposeDLambdaRange.GetBegin()==lambdaRange.GetBegin() && m_dTransposeDLambdaRange.GetEnd()==lambdaRange.GetEnd() ) )
    {
        initDtd(lambdaRange);
    }

    if(spcComponent=="nocontinuum")
    {
        return m_dTransposeDNocontinuum;
    }else
    {
        return m_dTransposeDRaw;
    }
}

/**
 * \brief this function returns the dtd value withing the wavelength range for a given spcComponent
 *
 **/
Float64 CLineModelElementList::getLikelihood_cstLog(const TFloat64Range& lambdaRange)
{
    if(! (m_dTransposeDLambdaRange.GetBegin()==lambdaRange.GetBegin() && m_dTransposeDLambdaRange.GetEnd()==lambdaRange.GetEnd() ) )
    {
        initDtd(lambdaRange);
    }

    return m_likelihood_cstLog;
}

/**
 * \brief this function estimates the dtd value withing the wavelength range
 **/
Float64 CLineModelElementList::EstimateDTransposeD(const TFloat64Range& lambdaRange, std::string spcComponent)
{
    const CSpectrumSpectralAxis& spcSpectralAxis = m_SpectrumModel->GetSpectralAxis();
    const CSpectrumFluxAxis& spcFluxAxis = m_SpcFluxAxis;
    const CSpectrumFluxAxis& spcFluxAxisNoContinuum = m_spcFluxAxisNoContinuum;

    Int32 numDevs = 0;
    Float64 dtd = 0.0;
    const Float64* Yspc = spcFluxAxis.GetSamples();
    const Float64* YspcNoContinuum = spcFluxAxisNoContinuum.GetSamples();
    Float64 flux = 0.0;

    Float64 imin = spcSpectralAxis.GetIndexAtWaveLength(lambdaRange.GetBegin());
    Float64 imax = spcSpectralAxis.GetIndexAtWaveLength(lambdaRange.GetEnd());
    for( UInt32 j=imin; j<imax; j++ )
    {
        numDevs++;
        if(spcComponent=="nocontinuum")
        {
            flux = YspcNoContinuum[j];
        }else
        {
            flux = Yspc[j];
        }
        dtd += (flux*flux) / (m_ErrorNoContinuum[j]*m_ErrorNoContinuum[j]);
    }
    Log.LogDebug( "CLineModelElementList::EstimateDTransposeD val = %f", dtd );

    return dtd;
}

/**
 * \brief this function estimates the mtm value withing the wavelength range
 **/
Float64 CLineModelElementList::EstimateMTransposeM(const TFloat64Range& lambdaRange)
{
    const CSpectrumSpectralAxis& spcSpectralAxis = m_SpectrumModel->GetSpectralAxis();
    const CSpectrumFluxAxis& spcFluxAxis = m_SpectrumModel->GetFluxAxis();

    Int32 numDevs = 0;
    Float64 mtm = 0.0;
    const Float64* Yspc = spcFluxAxis.GetSamples();
    Float64 diff = 0.0;

    Float64 imin = spcSpectralAxis.GetIndexAtWaveLength(lambdaRange.GetBegin());
    Float64 imax = spcSpectralAxis.GetIndexAtWaveLength(lambdaRange.GetEnd());
    for( UInt32 j=imin; j<imax; j++ )
    {
        numDevs++;
        {
            diff = Yspc[j];
        }
        mtm += (diff*diff) / (m_ErrorNoContinuum[j]*m_ErrorNoContinuum[j]);
    }
    //Log.LogDebug( "CLineModelElementList::EstimateMTransposeM val = %f", mtm );


    return mtm;
}

/**
 * \brief this function estimates the mtm value withing the wavelength range
 **/
Int32 CLineModelElementList::getMTransposeMCumulative(const TFloat64Range& lambdaRange, std::vector<Float64> lbda, std::vector<Float64> mtmCumul)
{
    mtmCumul.clear();
    lbda.clear();
    const CSpectrumSpectralAxis& spcSpectralAxis = m_SpectrumModel->GetSpectralAxis();
    const CSpectrumFluxAxis& spcFluxAxis = m_SpectrumModel->GetFluxAxis();

    Int32 numDevs = 0;
    Float64 mtm = 0.0;
    const Float64* Yspc = spcFluxAxis.GetSamples();
    Float64 diff = 0.0;

    Float64 imin = spcSpectralAxis.GetIndexAtWaveLength(lambdaRange.GetBegin());
    Float64 imax = spcSpectralAxis.GetIndexAtWaveLength(lambdaRange.GetEnd());
    for( UInt32 j=imin; j<imax; j++ )
    {
        numDevs++;
        {
            diff = Yspc[j];
        }
        mtm += (diff*diff) / (m_ErrorNoContinuum[j]*m_ErrorNoContinuum[j]);
        lbda.push_back(spcSpectralAxis[j]);
        mtmCumul.push_back(mtm);
    }
    //Log.LogDebug( "CLineModelElementList::EstimateMTransposeM val = %f", mtm );


    return 0;
}


/**
 * \brief this function estimates the likelihood_cstLog term withing the wavelength range
 **/
Float64 CLineModelElementList::EstimateLikelihoodCstLog(const TFloat64Range& lambdaRange)
{
    const CSpectrumSpectralAxis& spcSpectralAxis = m_SpectrumModel->GetSpectralAxis();

    Int32 numDevs = 0;
    Float64 cstLog = 0.0;
    Float64 sumLogNoise = 0.0;

    Float64 imin = spcSpectralAxis.GetIndexAtWaveLength(lambdaRange.GetBegin());
    Float64 imax = spcSpectralAxis.GetIndexAtWaveLength(lambdaRange.GetEnd());
    for( UInt32 j=imin; j<imax; j++ )
    {
        numDevs++;
        sumLogNoise += log( m_ErrorNoContinuum[j] );
    }
    //Log.LogDebug( "CLineModelElementList::EstimateMTransposeM val = %f", mtm );

    cstLog = -numDevs*0.5*log(2*M_PI) - sumLogNoise;

    return cstLog;
}


