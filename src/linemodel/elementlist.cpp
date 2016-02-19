#include <epic/redshift/linemodel/elementlist.h>
#include <epic/redshift/linemodel/singleline.h>
#include <epic/redshift/linemodel/multiline.h>

#include <epic/redshift/gaussianfit/multigaussianfit.h>
#include <gsl/gsl_multifit.h>

#include <epic/redshift/spectrum/io/genericreader.h>
#include <epic/redshift/spectrum/template/template.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>


#include <epic/redshift/continuum/irregularsamplingmedian.h>
#include <epic/redshift/spectrum/io/fitswriter.h>

#include <epic/core/debug/assert.h>
#include <epic/core/log/log.h>

#include <math.h>

#include <boost/format.hpp>
#include <boost/chrono/thread_clock.hpp>

#include <algorithm>

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include "lmfitfunctions.c"



using namespace NSEpic;

CLineModelElementList::CLineModelElementList(const CSpectrum& spectrum, const CSpectrum &spectrumContinuum, const CRayCatalog::TRayVector& restRayList, const std::string& opt_fittingmethod, const std::string& opt_continuumcomponent, const std::string& widthType, const Float64 resolution, const Float64 velocityEmission, const Float64 velocityAbsorption, const std::string& opt_rules)
{
    m_ContinuumComponent = opt_continuumcomponent;
    m_LineWidthType = widthType;
    m_resolution = resolution;
    m_velocityEmission = velocityEmission;
    m_velocityAbsorption = velocityAbsorption;
    m_fittingmethod = opt_fittingmethod;
    m_rulesoption = opt_rules;

    // to be deleted: nominal width
    //m_nominalWidthDefaultEmission = 3.4;//3.4; //suited to PFS RJLcont simulations
    m_nominalWidthDefaultEmission = 1.15;// suited to new pfs simulations
    m_nominalWidthDefaultAbsorption = m_nominalWidthDefaultEmission;

    LoadCatalog(restRayList);
    //LoadCatalogSingleLines(restRayList);

    //LogCatalogInfos();
    m_RestRayList = restRayList;
    m_SpectrumModel = std::shared_ptr<CSpectrum>( new CSpectrum(spectrum) );

    Int32 spectrumSampleCount = spectrum.GetSampleCount();
    m_SpcFluxAxis.SetSize( spectrumSampleCount );
    m_SpcContinuumFluxAxis = spectrumContinuum.GetFluxAxis();
    m_ContinuumFluxAxis.SetSize( spectrumSampleCount );
    m_SpcFluxAxisModelDerivSigma.SetSize( spectrumSampleCount );

    CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();
    const CSpectrumFluxAxis& spectrumFluxAxis = spectrum.GetFluxAxis();

    m_spcFluxAxisNoContinuum.SetSize( spectrumSampleCount );
    const Float64* error = spectrumFluxAxis.GetError();
    m_ErrorNoContinuum = m_spcFluxAxisNoContinuum.GetError();
    Float64* errorSpc = m_SpcFluxAxis.GetError();
    Float64* errorSpcContinuum = m_SpcContinuumFluxAxis.GetError();
    for(UInt32 i=0; i<spectrumSampleCount; i++){
        m_ErrorNoContinuum[i] = error[i];
        errorSpc[i] = error[i];
        errorSpcContinuum[i] = error[i];
    }

    if(m_ContinuumComponent == "nocontinuum"){
        //the continuum is set to zero and the observed spectrum is the spectrum without continuum
        for(UInt32 i=0; i<modelFluxAxis.GetSamplesCount(); i++){
            modelFluxAxis[i] = 0.0;
            m_ContinuumFluxAxis[i] = 0.0;
            m_spcFluxAxisNoContinuum[i] = spectrumFluxAxis[i]-m_SpcContinuumFluxAxis[i];
            m_SpcFluxAxis[i] = m_spcFluxAxisNoContinuum[i];
        }
    }else if(m_ContinuumComponent == "fromspectrum"){
        //the continuum is set to the SpcContinuum and the observed spectrum is the raw spectrum
        CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();
        for(UInt32 i=0; i<modelFluxAxis.GetSamplesCount(); i++){

            m_ContinuumFluxAxis[i] = m_SpcContinuumFluxAxis[i];
            m_spcFluxAxisNoContinuum[i] = spectrumFluxAxis[i]-m_SpcContinuumFluxAxis[i];

            modelFluxAxis[i] = m_ContinuumFluxAxis[i];
            m_SpcFluxAxis[i] = spectrumFluxAxis[i];
        }
    }
    m_precomputedFineGridContinuumFlux = NULL;


    /*
    // export continuum for debug
    FILE* fspc = fopen( "lm_continuum_dbg.txt", "w+" );
    Float64 coeffSaveSpc = 1e16;
    for(UInt32 i=0; i<spectrumSampleCount; i++){
        fprintf( fspc, "%f %f %f\n", m_SpectrumModel->GetSpectralAxis()[i], (m_SpcFluxAxis[i])*coeffSaveSpc, (m_ContinuumFluxAxis[i])*coeffSaveSpc);
    }
    fclose( fspc );
    //*/
}

CLineModelElementList::~CLineModelElementList()
{
}

const CSpectrum& CLineModelElementList::GetModelSpectrum() const
{
    return *m_SpectrumModel;
}


Float64 CLineModelElementList::getModelFluxVal(Int32 idx) const
{
    CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();
    if(idx<modelFluxAxis.GetSamplesCount()){
        return modelFluxAxis[idx];
    }

    return -1.0;
}


Float64 CLineModelElementList::getModelFluxDerivEltVal(Int32 DerivEltIdx, Int32 idx) const
{

    CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    Int32 iElts=DerivEltIdx;
    if(idx<spectralAxis.GetSamplesCount())
    {
        Float64 derivateVal = m_Elements[iElts]->GetModelDerivAmplitudeAtLambda(spectralAxis[idx], m_Redshift);
        return derivateVal;
    }

    return -1.0;
}


Float64 CLineModelElementList::getModelFluxDerivSigmaVal(Int32 idx) const
{
    if(idx<m_SpcFluxAxisModelDerivSigma.GetSamplesCount()){
        return m_SpcFluxAxisModelDerivSigma[idx];
    }

    return -1.0;
}


void CLineModelElementList::LoadCatalog(const CRayCatalog::TRayVector& restRayList)
{
    CRayCatalog crctlg;
    std::vector<CRayCatalog::TRayVector> groupList = crctlg.ConvertToGroupList(restRayList);
    for(Int32 ig=0; ig<groupList.size(); ig++){
        std::vector<CRay> lines;
        std::vector<Float64> amps;
        std::vector<Int32> inds;
        for(Int32 i=0; i<groupList[ig].size(); i++){
            std::vector<Int32> idx = findLineIdxInCatalog( restRayList, groupList[ig][i].GetName(), groupList[ig][i].GetType());
            inds.push_back(idx[0]);
            amps.push_back(groupList[ig][i].GetNominalAmplitude());
            lines.push_back(groupList[ig][i]);
        }
        if(lines.size()>0){
            m_Elements.push_back(boost::shared_ptr<CLineModelElement> (new CMultiLine(lines, m_LineWidthType, m_resolution, m_velocityEmission, m_velocityEmission, amps, m_nominalWidthDefaultAbsorption, inds)));
        }
    }
}

void CLineModelElementList::LogCatalogInfos()
{
    Log.LogInfo( "\n");
    Log.LogInfo( "LineModel Infos: %d elements", m_Elements.size());
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        Int32 nRays = m_Elements[iElts]->GetSize();
        if(nRays<1)
        {
            Log.LogInfo( "LineModel ctlg: elt %d (%s): no lines", iElts, m_Elements[iElts]->GetElementTypeTag().c_str());

        }
        for(UInt32 j=0; j<nRays; j++){
            std::string nominalAmpStr = "";
            if(nRays>1){
                nominalAmpStr = boost::str(boost::format("(nominal amp = %.4f)") % m_Elements[iElts]->GetNominalAmplitude(j));
            }
            Log.LogInfo( "LineModel ctlg: elt %d (%s): line %d = %s %s", iElts, m_Elements[iElts]->GetElementTypeTag().c_str(), j, m_Elements[iElts]->GetRayName(j).c_str(), nominalAmpStr.c_str());
        }
    }
    Log.LogInfo( "\n");
}

void CLineModelElementList::LoadCatalogSingleLines(const CRayCatalog::TRayVector& restRayList)
{
    //Load the rest of the single lines
    for( UInt32 iRestRay=0; iRestRay<restRayList.size(); iRestRay++ )
    {
        if ( FindElementIndex(iRestRay)==-1 )
        {
            addSingleLine(restRayList[iRestRay], iRestRay, m_nominalWidthDefaultEmission);
        }
    }
}

void CLineModelElementList::LoadContinuum()
{
    //path p( templatePath );
    //path name = p.leaf();
    //std::string templatePath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/Templates/ExtendedGalaxyEL2/emission/NEW_Im_extended.dat";
    //std::string templatePath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/Templates/ExtendedGalaxyEL2/galaxy/EW_SB2extended.dat";
    //std::string templatePath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/Templates/ExtendedGalaxyEL2/galaxy/BulgedataExtensionData.dat";

    std::string templatePath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/Templates/linemodel/emission/NEW_Im_extended_blue_continuum.txt";
    //std::string templatePath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/Templates/linemodel/emission/NEW_Im_extended_continuum.txt";
    //std::string templatePath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/Templates/linemodel/galaxy/EW_SB2extended.txt";
    //std::string templatePath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/Templates/linemodel/galaxy/BulgedataExtensionData.txt";

    //std::string templatePath = "/home/aschmitt/data/pfs/pfs_testsimu_20151009/47002690000013_flam.txt";


    //CRef<CTemplate> tmpl = new CTemplate();
    CTemplate tpl;
    CSpectrumIOGenericReader asciiReader;
    if( !asciiReader.Read( templatePath.c_str(), tpl ) ) {
        Log.LogError("Fail to read template: %s", templatePath.c_str());
        return;
    }

//    //save as a spectrum
//    CSpectrumIOFitsWriter writer;
//    Bool retVal1 = writer.Write( "tplexported.fits",  tpl);


//    //create a continuum tpl
//    if(1)
//    {
//        // Remove continuum
//        CContinuumIrregularSamplingMedian continuum;
//        CSpectrumFluxAxis fluxAxisWithoutContinuumCalc;

//        Int32 retVal = continuum.RemoveContinuum( tpl, fluxAxisWithoutContinuumCalc );
//        CSpectrumFluxAxis fluxAxis = tpl.GetFluxAxis();
//        fluxAxis.Subtract(fluxAxisWithoutContinuumCalc);
//        CSpectrumSpectralAxis tplSpectralAxis = tpl.GetSpectralAxis();
//        //*//debug:
//        // save continuum tpl and xmad,  flux data
//        FILE* f = fopen( "continuum_tpl_dbg.txt", "w+" );
//        for( Int32 t=0;t<fluxAxisWithoutContinuumCalc.GetSamplesCount();t++)
//        {
//            fprintf( f, "%f %f\n", t, tplSpectralAxis[t], fluxAxis[t]);//*1e12);
//        }
//        fclose( f );
//        //*/
//        return;
//    }


    //*/
    // Precalculate a fine grid template to be used for the 'closest value' rebin method
    Int32 n = tpl.GetSampleCount();
    CSpectrumFluxAxis tplFluxAxis = tpl.GetFluxAxis();
    CSpectrumSpectralAxis tplSpectralAxis = tpl.GetSpectralAxis();
    //Float64 dLambdaTgt =  1.0 * ( spectrum.GetMeanResolution()*0.9 )/( 1+sortedRedshifts[sortedRedshifts.size()-1] );
    Float64 dLambdaTgt =  0.1; //should be sufficient for the continuum
    //Float64 lmin = tplSpectralAxis[0];
    Float64 lmin = 0;
    Float64 lmax = tplSpectralAxis[n-1];
    Int32 nTgt = (lmax-lmin)/dLambdaTgt + 2.0/dLambdaTgt;

    m_precomputedFineGridContinuumFlux = (Float64*)malloc(nTgt*sizeof(Float64));

    //inialise and allocate the gsl objects
    Float64* Ysrc = tplFluxAxis.GetSamples();
    Float64* Xsrc = tplSpectralAxis.GetSamples();

    //spline
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n);
    gsl_spline_init (spline, Xsrc, Ysrc, n);
    gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();

    Int32 k = 0;
    Float64 x = 0.0;
    for(k=0; k<nTgt; k++){
        x = lmin + k*dLambdaTgt;
        if(x < tplSpectralAxis[0] || x > tplSpectralAxis[n-1]){
            m_precomputedFineGridContinuumFlux[k] = 0.0; //todo, make sure this is never used in the next steps...
        }else{
            m_precomputedFineGridContinuumFlux[k] = gsl_spline_eval (spline, x, accelerator);
        }
    }
    //*/

}

//this function prepares the continuum for use in the fit with the line elements
//1. Rebin with PFG buffer
//2. Find and apply amplitude factor from cross-corr
void CLineModelElementList::PrepareContinuum(Float64 z)
{
    const CSpectrumSpectralAxis& targetSpectralAxis = m_SpectrumModel->GetSpectralAxis();
    const Float64* Xtgt = targetSpectralAxis.GetSamples();
    Float64* Yrebin = m_ContinuumFluxAxis.GetSamples();

    if(m_precomputedFineGridContinuumFlux == NULL){
        for ( Int32 i = 0; i<targetSpectralAxis.GetSamplesCount(); i++)
        {
            //Yrebin[i] = m_SpcNoContinuumFluxAxis[i]; //
            Yrebin[i] = m_SpcContinuumFluxAxis[i];
        }
        return;
    }

    TFloat64Range currentRange(Xtgt[0], Xtgt[targetSpectralAxis.GetSamplesCount()-1]);

    // Move cursors up to lambda range start
    Int32 j = 0;
    while( j<targetSpectralAxis.GetSamplesCount() && Xtgt[j] < currentRange.GetBegin() )
    {
        Yrebin[j] = 0.0;
        j++;
    }
    //* // Precomputed FINE GRID nearest sample,
    Int32 k = 0;
    Float64 dl = 0.1;
    Float64 Coeffk = 1.0/dl/(1+z);
    // For each sample in the target spectrum
    while( j<targetSpectralAxis.GetSamplesCount() && Xtgt[j] <= currentRange.GetEnd() )
    {
        k = (int)(Xtgt[j]*Coeffk+0.5);
        Yrebin[j] = m_precomputedFineGridContinuumFlux[k];
        //Yrebin[j] = 0.0;
        j++;

    }
    //*/
    while( j < targetSpectralAxis.GetSamplesCount() )
    {
        Yrebin[j] = 0.0;
        j++;
    }


    //fit the continuum
    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    //const Float64* flux = m_SpcFluxAxis.GetSamples();
    const Float64* flux = m_SpcContinuumFluxAxis.GetSamples();

    const Float64* spectral = spectralAxis.GetSamples();

    Float64 sumCross = 0.0;
    Float64 sumGauss = 0.0;
    Float64 err2 = 0.0;
    Int32 num = 0;

    Float64 x=0.0;
    Float64 y=0.0;
    Float64 yg=0.0;

    //A estimation
    for ( Int32 i = 0; i<targetSpectralAxis.GetSamplesCount(); i++)
    {
        y = flux[i];
        x = spectral[i];
        yg = Yrebin[i];

        num++;
        err2 = 1.0 / (m_ErrorNoContinuum[i] * m_ErrorNoContinuum[i]);
        sumCross += yg*y*err2;
        sumGauss += yg*yg*err2;
    }

    if ( num==0 || sumGauss==0 )
    {
        return;
    }

    Float64 A = std::max(0.0, sumCross / sumGauss);
    for ( Int32 i = 0; i<targetSpectralAxis.GetSamplesCount(); i++)
    {
        Yrebin[i] *=A;
    }
}

Float64 CLineModelElementList::fit(Float64 redshift, const TFloat64Range& lambdaRange, CLineModelResult::SLineModelSolution& modelSolution, Int32 contreest_iterations)
{
    m_Redshift = redshift;

    if(m_ContinuumComponent != "nocontinuum"){
        //prepare the continuum
        PrepareContinuum(redshift);
    }

    //initialize the model spectrum
    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();


    //prepare the elements
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        m_Elements[iElts]->prepareSupport(spectralAxis, redshift, lambdaRange);
    }

    //EstimateSpectrumContinuum();

    if(m_ContinuumComponent == "nocontinuum"){
//        for(UInt32 i=0; i<modelFluxAxis.GetSamplesCount(); i++){
//            modelFluxAxis[i] = 0.0;
//        }
    }else{
        for(UInt32 i=0; i<modelFluxAxis.GetSamplesCount(); i++){
            modelFluxAxis[i] = m_ContinuumFluxAxis[i];
            m_spcFluxAxisNoContinuum[i] = m_SpcFluxAxis[i]-m_ContinuumFluxAxis[i];
        }
    }

    //fit the amplitudes of each element independently
    if(m_fittingmethod=="individual")
    {
        //fit the model amplitudes individually
        for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
        {
            m_Elements[iElts]->fitAmplitude(spectralAxis, m_spcFluxAxisNoContinuum, redshift);
        }
    }
    //else{
    //    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    //    {
    //        m_Elements[iElts]->SetFittedAmplitude(1e-19);
    //    }
    //}

    //fit the amplitude of all elements together with iterative solver: Nelder Mead Simplex
    if(0){
        //fit the amplitudes together
        fitAmplitudesSimplex();
    }

    //fit the amplitude of all elements together (but Emission or Absorption separately) with iterative   solver: lmfit
    if(m_fittingmethod=="lmfit")
    {
        std::vector<Int32> validEltsIdx = GetModelValidElementsIndexes();
        std::vector<Float64> ampsfitted;
        Int32 lineType = CRay::nType_Absorption;
        fitAmplitudesLmfit(validEltsIdx, m_spcFluxAxisNoContinuum, ampsfitted, lineType);
        lineType = CRay::nType_Emission;
        fitAmplitudesLmfit(validEltsIdx, m_spcFluxAxisNoContinuum, ampsfitted, lineType);
    }

    //fit the amplitude of all elements together with linear solver: gsl_multifit_wlinear
    if(m_fittingmethod=="svd")
    {
        std::vector<Int32> validEltsIdx = GetModelValidElementsIndexes();
        std::vector<Float64> ampsfitted;
        fitAmplitudesLinSolve(validEltsIdx, spectralAxis, m_spcFluxAxisNoContinuum, ampsfitted);
    }

    //fit the amplitudes of each element independently, unless there is overlap
    if(m_fittingmethod=="hybrid")
    {
        fitAmplitudesHybrid(spectralAxis, m_spcFluxAxisNoContinuum, redshift);

        //apply a continuum iterative re-estimation with lines removed from the initial spectrum
        Int32 nIt = contreest_iterations;
        Int32 it=0;
        while(it<nIt){
            applyRules();

            //*
            //iterative continuum estimation :: RAW SLOW METHOD
            refreshModel();
            EstimateSpectrumContinuum();

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


            fitAmplitudesHybrid(spectralAxis, m_spcFluxAxisNoContinuum, redshift);
            it++;
        }
    }


    //Apply rules,
    applyRules();

    refreshModel();

    //create spectrum model
    modelSolution = GetModelSolution();


    /*
    //model
    CSpectrum spcmodel = GetModelSpectrum();
    CSpectrumIOFitsWriter writer;
    Bool retVal1 = writer.Write( "model.fits",  spcmodel);

    if(retVal1){
        CSpectrum s(spcmodel);
        s.GetFluxAxis() = m_SpcFluxAxis;
        Bool retVal2 = writer.Write( "spectrum.fits",  s);
    }
    //*/


    /*
    //model for linefitting
    CSpectrum spcmodel4linefitting = GetModelSpectrum();
    for(UInt32 i=0; i<spcmodel4linefitting.GetFluxAxis().GetSamplesCount(); i++){
        spcmodel4linefitting.GetFluxAxis()[i] = spcmodel4linefitting.GetFluxAxis()[i]-m_ContinuumFluxAxis[i];
    }
    CSpectrumIOFitsWriter writer2;
    Bool retVal2 = writer2.Write( "model4linefit.fits",  spcmodel4linefitting);

    if(retVal1){
        CSpectrum s4linefitting(spcmodel);
        s4linefitting.GetFluxAxis() = spcFluxAxisNoContinuum;
        Bool retVal2 = writer2.Write( "spectrum4linefit.fits",  s4linefitting);
    }
    //*/
    Float64 merit = getLeastSquareMerit(lambdaRange);

    if(m_ContinuumComponent == "nocontinuum"){
        reinitModel();
    }

    return merit;
}

void CLineModelElementList::SetFittingMethod(std::string fitMethod)
{
    m_fittingmethod = fitMethod;
}

/**
 * @brief CLineModelElementList::fit with model selection
 * @param redshift
 * @param lambdaRange
 * @param modelSolution
 */
// WARNING - DEPRECATED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
void CLineModelElementList::fitWithModelSelection(Float64 redshift, const TFloat64Range& lambdaRange, CLineModelResult::SLineModelSolution& modelSolution)
{
    m_Redshift = redshift;

    //prepare the continuum
    PrepareContinuum(redshift);

    //initialize the model spectrum
    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();
    CSpectrumFluxAxis spcFluxAxisNoContinuum(m_SpcFluxAxis);

    const Float64* error = m_SpcFluxAxis.GetError();
    Float64* errorNoContinuum = spcFluxAxisNoContinuum.GetError();

    spcFluxAxisNoContinuum.SetSize( modelFluxAxis.GetSamplesCount() );


    //prepare the elements
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        m_Elements[iElts]->prepareSupport(spectralAxis, redshift, lambdaRange);
    }

    for(UInt32 i=0; i<modelFluxAxis.GetSamplesCount(); i++){
        modelFluxAxis[i] = m_ContinuumFluxAxis[i];
        spcFluxAxisNoContinuum[i] = m_SpcFluxAxis[i]-m_ContinuumFluxAxis[i];
        errorNoContinuum[i] = error[i];
    }

    //fit the amplitudes of each element independently, unless there is overlap
    if(1)
    {
        fitAmplitudesHybrid(spectralAxis, spcFluxAxisNoContinuum, redshift);

        //apply a continuum iterative re-estimation with lines removed from the initial spectrum
        Int32 nIt = 1;
        Int32 it=0;
        while(it<nIt){
            applyRules();

            /*
            //iterative continuum estimation :: RAW SLOW METHOD
            refreshModel();
            EstimateSpectrumContinuum();

            for(UInt32 i=0; i<modelFluxAxis.GetSamplesCount(); i++){
                modelFluxAxis[i] = m_ContinuumFluxAxis[i];
                spcFluxAxisNoContinuum[i] = m_SpcFluxAxis[i]-m_ContinuumFluxAxis[i];
            }
            //*/

            //*
            std::vector<Int32> validEltsIdx = GetModelValidElementsIndexes();
            std::vector<Int32> refreshIdxs = ReestimateContinuumApprox(validEltsIdx);
            refreshModelAfterContReestimation(refreshIdxs, modelFluxAxis, spcFluxAxisNoContinuum);
            //*/


            fitAmplitudesHybrid(spectralAxis, spcFluxAxisNoContinuum, redshift);
            it++;
        }
    }


    //Apply rules,
    applyRules();

    reinitModel();

    //create spectrum model by adding only the lines that improve the BIC
    CSpectrumFluxAxis modelFluxAxisTmp = m_SpectrumModel->GetFluxAxis();
    Int32 nddl = 0;
    Float64 bic = 1e12;
    Int32 nsamples = modelFluxAxisTmp.GetSamplesCount();
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        m_Elements[iElts]->addToSpectrumModel(spectralAxis, modelFluxAxis, m_Redshift);
        Float64 fitTmp = getLeastSquareMerit(lambdaRange);//getLeastSquareMeritUnderElements();
        Float64 bicTmp = fitTmp + (nddl+1)*log(nsamples);
        //Float64 bicTmp = fitTmp + 2*(nddl+1); //AIC

        if(iElts==0 || bicTmp<bic){
            m_Elements[iElts]->addToSpectrumModel(spectralAxis, modelFluxAxisTmp, m_Redshift);
            nddl++;
            bic = bicTmp;
        }else{
            SetElementAmplitude(iElts, 0.0, 0.0);
            modelFluxAxis = modelFluxAxisTmp;
        }
    }

    //create spectrum model
    modelSolution = GetModelSolution();


    //*
    //model
    CSpectrum spcmodel = GetModelSpectrum();
    CSpectrumIOFitsWriter writer;
    Bool retVal1 = writer.Write( "model.fits",  spcmodel);

    if(retVal1){
        CSpectrum s(spcmodel);
        s.GetFluxAxis() = m_SpcFluxAxis;
        Bool retVal2 = writer.Write( "spectrum.fits",  s);
    }
    //*/


    /*
    //model for linefitting
    CSpectrum spcmodel4linefitting = GetModelSpectrum();
    for(UInt32 i=0; i<spcmodel4linefitting.GetFluxAxis().GetSamplesCount(); i++){
        spcmodel4linefitting.GetFluxAxis()[i] = spcmodel4linefitting.GetFluxAxis()[i]-m_ContinuumFluxAxis[i];
    }
    CSpectrumIOFitsWriter writer2;
    Bool retVal2 = writer2.Write( "model4linefit.fits",  spcmodel4linefitting);

    if(retVal1){
        CSpectrum s4linefitting(spcmodel);
        s4linefitting.GetFluxAxis() = spcFluxAxisNoContinuum;
        Bool retVal2 = writer2.Write( "spectrum4linefit.fits",  s4linefitting);
    }
    //*/

}

void CLineModelElementList::reinitModel()
{
    CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();
    //init spectrum model with continuum
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        m_Elements[iElts]->initSpectrumModel(modelFluxAxis, m_ContinuumFluxAxis);
    }
}

void CLineModelElementList::reinitModelUnderElements(std::vector<Int32>  filterEltsIdx)
{
    Int32 iElts;
    CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();
    //init spectrum model with continuum
    for( UInt32 i=0; i<filterEltsIdx.size(); i++ )
    {
        iElts = filterEltsIdx[i];
        m_Elements[iElts]->initSpectrumModel(modelFluxAxis, m_ContinuumFluxAxis);
    }
}

void CLineModelElementList::refreshModel()
{
    reinitModel();
    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();
    //create spectrum model
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        m_Elements[iElts]->addToSpectrumModel(spectralAxis, modelFluxAxis, m_Redshift);
    }
}

void CLineModelElementList::refreshModelUnderElements(std::vector<Int32> filterEltsIdx)
{
    reinitModelUnderElements(filterEltsIdx);
    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();
    //create spectrum model
    Int32 iElts;
    for( UInt32 i=0; i<filterEltsIdx.size(); i++ )
    {
        iElts = filterEltsIdx[i];
        m_Elements[iElts]->addToSpectrumModel(spectralAxis, modelFluxAxis, m_Redshift);
    }
}

void CLineModelElementList::refreshModelDerivSigmaUnderElements(std::vector<Int32> filterEltsIdx)
{
    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    std::vector<Int32> supportIdxes = getSupportIndexes( filterEltsIdx );
    for( UInt32 i=0; i<supportIdxes.size(); i++ )
    {
        m_SpcFluxAxisModelDerivSigma[supportIdxes[i]]=0.0;
    }

    //create spectrum model partial derivate vs sigma
    Int32 iElts;
    for( UInt32 i=0; i<filterEltsIdx.size(); i++ )
    {
        iElts = filterEltsIdx[i];
        m_Elements[iElts]->addToSpectrumModelDerivSigma(spectralAxis, m_SpcFluxAxisModelDerivSigma, m_Redshift);
    }
}

Int32 CLineModelElementList::fitAmplitudesHybrid( const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& spcFluxAxisNoContinuum, Float64 redshift)
{
    std::vector<Int32> validEltsIdx = GetModelValidElementsIndexes();

    std::vector<Int32> indexesFitted;
    for( UInt32 iValidElts=0; iValidElts<validEltsIdx.size(); iValidElts++ )
    {
        Int32 iElts = validEltsIdx[iValidElts];

        //skip if already fitted
        bool alreadyfitted=false;
        for(Int32 i=0; i<indexesFitted.size(); i++){
            if(iElts == indexesFitted[i]){
                alreadyfitted=true;
                continue;
            }
        }
        if(alreadyfitted){
            continue;
        }

        //do the fit on the ovelapping elements
        Float64 overlapThres = 0.33;
        std::vector<Int32> overlappingInds = getOverlappingElements(iElts, indexesFitted, overlapThres);

//        Log.LogInfo( "Redshift: %f", m_Redshift);
//        Log.LogInfo( "hybrid fit: idx=%d - overlappingIdx=%d", iValidElts, overlappingInds.size());
//        for(Int32 ifit=0; ifit<overlappingInds.size(); ifit++)
//        {
//            Log.LogInfo( "hybrid fit: i=%d - Id=%d", ifit, overlappingInds[ifit]);
//        }

        if(overlappingInds.size()<2){
            m_Elements[iElts]->fitAmplitude(spectralAxis, spcFluxAxisNoContinuum, redshift);
        }else{
            std::vector<Float64> ampsfitted;
            Int32 retVal = fitAmplitudesLinSolve(overlappingInds, spectralAxis, spcFluxAxisNoContinuum, ampsfitted);
            // todo: if all the amplitudes fitted don't have the same sign, do it separately
            std::vector<Int32> overlappingIndsSameSign;
            if(retVal!=1){
                for(Int32 ifit=0; ifit<overlappingInds.size(); ifit++)
                {
                    if(ampsfitted[ifit]>0){
                        overlappingIndsSameSign.push_back(overlappingInds[ifit]);
                        //m_Elements[overlappingInds[ifit]]->fitAmplitude(spectralAxis, spcFluxAxisNoContinuum, redshift);
                    }else{
                        SetElementAmplitude(overlappingInds[ifit], 0.0, 0.0);
                    }
                }
                //fit the rest of the overlapping elements (same sign) together
                if(overlappingIndsSameSign.size()==1){
                    m_Elements[overlappingIndsSameSign[0]]->fitAmplitude(spectralAxis, spcFluxAxisNoContinuum, redshift);
                }else if(overlappingIndsSameSign.size()>1){
//                    for(Int32 ifit=0; ifit<overlappingIndsSameSign.size(); ifit++)
//                    {
//                        SetElementAmplitude(overlappingIndsSameSign[ifit], 0.0, 0.0);
//                    }
                    Int32 retVal2 = fitAmplitudesLinSolve(overlappingIndsSameSign, spectralAxis, spcFluxAxisNoContinuum, ampsfitted);
                    if(retVal2!=1){
                        for(Int32 ifit=0; ifit<overlappingIndsSameSign.size(); ifit++)
                        {
                            if(ampsfitted[ifit]>0){
                                m_Elements[overlappingIndsSameSign[ifit]]->fitAmplitude(spectralAxis, spcFluxAxisNoContinuum, redshift);
                            }else{
                                SetElementAmplitude(overlappingIndsSameSign[ifit], 0.0, 0.0);
                            }
                        }
                    }
                }
            }
        }

        //update the already fitted list
        for(Int32 i=0; i<overlappingInds.size(); i++){
            indexesFitted.push_back(overlappingInds[i]);
        }

    }

    return 0;
}

void CLineModelElementList::fitAmplitudesSimplex()
{
    CMultiGaussianFit fitter;

    Int32 status = fitter.Compute( *this );
//    if(status!=NSEpic::CMultiGaussianFit::nStatus_Success){
//        continue;
//    }

//    Float64 gaussAmp;
//    Float64 gaussPos;
//    Float64 gaussWidth;
//    fitter.GetResults( gaussAmp, gaussPos, gaussWidth );
//    Float64 gaussAmpErr;
//    Float64 gaussPosErr;
//    Float64 gaussWidthErr;
//    fitter.GetResultsError( gaussAmpErr, gaussPosErr, gaussWidthErr );


}

//temporary stuff, dev of the lmfit method: int

void print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
  printf ("iter: %3u x = % 15.8f % 15.8f % 15.8f "
          "|f(x)| = %g\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2),
          gsl_blas_dnrm2 (s->f));
}


Int32 CLineModelElementList::fitAmplitudesLmfit(std::vector<Int32> EltsIdx, const CSpectrumFluxAxis& fluxAxis, std::vector<Float64>& ampsfitted, Int32 lineType)
{
    //http://www.gnu.org/software/gsl/manual/html_node/Example-programs-for-Nonlinear-Least_002dSquares-Fitting.html
    Bool verbose=false;

    std::vector<Int32> filteredEltsIdx;
    for (Int32 iElt = 0; iElt < EltsIdx.size(); iElt++)
    {
        if(m_RestRayList[m_Elements[iElt]->m_LineCatalogIndexes[0]].GetType() != lineType)
        {
            continue;
        }
        filteredEltsIdx.push_back(EltsIdx[iElt]);
    }

    Int32 nddl = filteredEltsIdx.size();
    if(nddl<1){
        return -1;
    }
    std::vector<Int32> xInds = getSupportIndexes( filteredEltsIdx );

    const Float64* flux = fluxAxis.GetSamples();

    //create gsl solver object
    const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
    gsl_multifit_fdfsolver *s;
    int status, info;
    size_t i;
    size_t n; //n samples on the support, /* number of data points to fit */
    size_t p = nddl+1; //DOF = n amplitudes to fit (1 for each element) + 1 (EL velocity)

    n = xInds.size();
    if(n<nddl){
        ampsfitted.resize(nddl);
        for (Int32 iddl = 0; iddl < nddl; iddl++)
        {
            ampsfitted[iddl] = 0.0;
        }
        return -1;
    }

    Int32 idx = 0;
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
        fprintf(stderr, "normFactor = '%.3e'\n", normFactor);
    }

    gsl_matrix *J = gsl_matrix_alloc(n, p);
    gsl_matrix *covar = gsl_matrix_alloc (p, p);
    double y[n], weights[n];
    struct lmfitdata d = {n,y,this, filteredEltsIdx, xInds, normFactor, lineType};
    gsl_multifit_function_fdf f;

    Float64* x_init = (Float64*) calloc( p, sizeof( Float64 ) );
    //TODO: initialize lmfit with individual/hybrid fit method
    for(Int32 kp=0; kp<nddl; kp++)
    {
        x_init[kp] = 0;
    }
    if(lineType==CRay::nType_Emission)
    {
        x_init[nddl] = GetVelocityEmission();
    }else
    {
        x_init[nddl] = GetVelocityAbsorption();

    }

    gsl_vector_view x = gsl_vector_view_array (x_init, p);
    gsl_vector_view w = gsl_vector_view_array(weights, n);
    const gsl_rng_type * type;
    gsl_rng * r;
    gsl_vector *res_f;
    double chi, chi0;

    const double xtol = 1e-8;
    const double gtol = 1e-8;
    const double ftol = 0.0;

    gsl_rng_env_setup();

    type = gsl_rng_default;
    r = gsl_rng_alloc (type);

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

    /* solve the system with a maximum of 20 iterations */
    status = gsl_multifit_fdfsolver_driver(s, 50, xtol, gtol, ftol, &info);

    gsl_multifit_fdfsolver_jac(s, J);
    gsl_multifit_covar (J, 0.0, covar);

    /* compute final residual norm */
    chi = gsl_blas_dnrm2(res_f);

    if(verbose)
    {
#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

        fprintf(stderr, "summary from method '%s'\n",
                gsl_multifit_fdfsolver_name(s));
        fprintf(stderr, "number of iterations: %zu\n",
                gsl_multifit_fdfsolver_niter(s));
        fprintf(stderr, "function evaluations: %zu\n", f.nevalf);
        fprintf(stderr, "Jacobian evaluations: %zu\n", f.nevaldf);
        fprintf(stderr, "reason for stopping: %s\n",
                (info == 1) ? "small step size" : "small gradient");
        fprintf(stderr, "initial |f(x)| = %g\n", chi0);
        fprintf(stderr, "final   |f(x)| = %g\n", chi);

        {
            double dof = n - p;
            double c = GSL_MAX_DBL(1, chi / sqrt(dof));

            fprintf(stderr, "chisq/dof = %g\n",  pow(chi, 2.0) / dof);

            for(Int32 k=0; k<p; k++)
            {
                if(FIT(k)<1e-3)
                {
                    fprintf (stderr, "A      = %.3e +/- %.8f\n", FIT(k), c*ERR(k));
                }else{
                    fprintf (stderr, "A      = %.5f +/- %.8f\n", FIT(k), c*ERR(k));
                }

            }
        }
        fprintf (stderr, "status = %s\n", gsl_strerror (status));
    }

    // finally populate the fitting results to the linemodel
    {
        for (Int32 iddl = 0; iddl < nddl; iddl++)
        {
            Float64 a = gsl_vector_get(s->x,iddl)/normFactor;
            Float64 dof = n - p;
            Float64 c = GSL_MAX_DBL(1, chi / sqrt(dof));
            Float64 sigma = c*sqrt(gsl_matrix_get(covar,iddl,iddl))/normFactor;
            if(a<0.0){
                a=0.0;
            }
            SetElementAmplitude(filteredEltsIdx[iddl], a, sigma);
        }
        Float64 velEmission = gsl_vector_get(s->x,nddl);
        SetVelocityEmission(velEmission);
    }

    gsl_multifit_fdfsolver_free (s);
    gsl_matrix_free (covar);
    gsl_matrix_free (J);
    gsl_rng_free (r);

    return 0;
}

std::vector<Int32> CLineModelElementList::getSupportIndexes( std::vector<Int32> EltsIdx)
{
    std::vector<Int32> indexes;

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


std::vector<Int32> CLineModelElementList::getOverlappingElementsBySupport(  Int32 ind, Float64 overlapThres)
{
    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    std::vector<Int32> indexes;

    if(m_Elements[ind]->IsOutsideLambdaRange()){
        indexes.push_back(ind);
        return indexes;
    }
    TInt32RangeList refsupport = m_Elements[ind]->getSupport();
    CRay ray = m_RestRayList[m_Elements[ind]->m_LineCatalogIndexes[0]];
    Int32 linetype = ray.GetType();
    Float64 mu = ray.GetPosition()*(1+m_Redshift);
    Float64 c = m_Elements[ind]->GetLineWidth(mu, m_Redshift, ray.GetIsEmission());
    std::string profile = ray.GetProfile();
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
                if( max-min < -overlapThresholdMin ){
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

std::vector<Int32> CLineModelElementList::getOverlappingElements(Int32 ind, std::vector<Int32> excludedInd, Float64 overlapThres)
{
    std::vector<Int32> indexes;

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
                Float64 cRef = m_Elements[ind]->GetLineWidth(muRef, m_Redshift, raysRef[iRayRef].GetIsEmission());
                std::string profileRef = raysRef[iRayRef].GetProfile();
                Float64 winsizeRef = m_Elements[ind]->GetNSigmaSupport(profileRef)*cRef;
                Float64 overlapSizeMin = winsizeRef*overlapThres;
                xinf = muRef-winsizeRef/2.0;
                xsup = muRef+winsizeRef/2.0;

                Float64 muElt = raysElt[iRayElt].GetPosition()*(1+m_Redshift);
                Float64 cElt = m_Elements[iElts]->GetLineWidth(muElt, m_Redshift, raysElt[iRayElt].GetIsEmission());
                std::string profileElt = raysElt[iRayElt].GetProfile();
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

Int32 CLineModelElementList::fitAmplitudesLinSolve( std::vector<Int32> EltsIdx, const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis, std::vector<Float64>& ampsfitted)
{
    boost::chrono::thread_clock::time_point start_prep = boost::chrono::thread_clock::now();

    Int32 nddl = EltsIdx.size();
    if(nddl<1){
        return -1;
    }
    std::vector<Int32> xInds = getSupportIndexes( EltsIdx );

    for (Int32 iddl = 0; iddl < nddl; iddl++)
    {
        SetElementAmplitude(EltsIdx[iddl], 1.0, 0.0);
    }

    const Float64* spectral = spectralAxis.GetSamples();
    const Float64* flux = fluxAxis.GetSamples();

    //Linear fit
    int i, n;
    Float64 fval;
    double xi, yi, ei, chisq;
    gsl_matrix *X, *cov;
    gsl_vector *y, *w, *c;

    n = xInds.size();
    if(n<nddl){
        ampsfitted.resize(nddl);
        for (Int32 iddl = 0; iddl < nddl; iddl++)
        {
            ampsfitted[iddl] = 0.0;
        }
        return -1;
    }

    X = gsl_matrix_alloc (n, nddl);
    y = gsl_vector_alloc (n);
    w = gsl_vector_alloc (n);

    c = gsl_vector_alloc (nddl);
    cov = gsl_matrix_alloc (nddl, nddl);

    //
    //todo: normalize, center...
    //
    Int32 idx = 0;
    for (i = 0; i < n; i++)
    {
        idx = xInds[i];
        xi = spectral[idx];
        yi = flux[idx];
        ei = m_ErrorNoContinuum[idx];

        for (Int32 iddl = 0; iddl < nddl; iddl++)
        {
            fval =  m_Elements[EltsIdx[iddl]]->getModelAtLambda(xi, m_Redshift);
            gsl_matrix_set (X, i, iddl, fval);
        }

        gsl_vector_set (y, i, yi);
        gsl_vector_set (w, i, 1.0/(ei*ei));
    }

    //
    boost::chrono::thread_clock::time_point stop_prep = boost::chrono::thread_clock::now();
    Float64 duration_prep = boost::chrono::duration_cast<boost::chrono::microseconds>(stop_prep - start_prep).count();
    boost::chrono::thread_clock::time_point start_fit = boost::chrono::thread_clock::now();

    {
      gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, nddl);
      gsl_multifit_wlinear (X, w, y, c, cov, &chisq, work);
      gsl_multifit_linear_free (work);
    }

    //
    boost::chrono::thread_clock::time_point stop_fit = boost::chrono::thread_clock::now();
    Float64 duration_fit = boost::chrono::duration_cast<boost::chrono::microseconds>(stop_fit - start_fit).count();
    //Log.LogInfo( "LineModel linear fit: prep = %.3f - fit = %.3f", duration_prep, duration_fit);

//#define C(i) (gsl_vector_get(c,(i)))
//#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))
//    if(0){
//        Log.LogInfo("# best fit: Y = %g X1 + %g X2 ...", C(0), C(1));
//        Log.LogInfo("# covariance matrix:\n");
//        Log.LogInfo("[ %+.5e, %+.5e \n", COV(0,0), COV(0,1));
//        Log.LogInfo("  %+.5e, %+.5e \n", COV(1,0), COV(1,1));

////        Log.LogInfo("[ %+.5e, %+.5e, %+.5e  \n", COV(0,0), COV(0,1), COV(0,2));
////        Log.LogInfo("  %+.5e, %+.5e, %+.5e  \n", COV(1,0), COV(1,1), COV(1,2));
////        Log.LogInfo("  %+.5e, %+.5e, %+.5e ]\n", COV(2,0), COV(2,1), COV(2,2));

//        Log.LogInfo("# chisq/n = %g", chisq/n);
//    }

    Int32 sameSign = 1;
    Float64 a0 = gsl_vector_get(c,0);
    for (Int32 iddl = 1; iddl < nddl; iddl++)
    {
        Float64 a = gsl_vector_get(c,iddl);
        Float64 product = a0*a;
        if(product<0){
            sameSign = 0;
        }
    }


    if(sameSign){
        for (Int32 iddl = 0; iddl < nddl; iddl++)
        {
            Float64 a = gsl_vector_get(c,iddl);
            Float64 cova = gsl_matrix_get(cov,iddl,iddl);
            Float64 sigma = sqrt(cova);
            SetElementAmplitude(EltsIdx[iddl], a, sigma);
        }
        //refreshModel();
    }else{
        ampsfitted.resize(nddl);
        for (Int32 iddl = 0; iddl < nddl; iddl++)
        {
            Float64 a = gsl_vector_get(c,iddl);
            ampsfitted[iddl] = (a);
        }
    }
    gsl_matrix_free (X);
    gsl_vector_free (y);
    gsl_vector_free (w);
    gsl_vector_free (c);
    gsl_matrix_free (cov);

    return sameSign;
}

/**
* @brief CLineModelElementList::ReestimateContinuumUnderLines
* For each line, reestimate the continuum using the original median routines on a sub segment of the original spectrum
* - the subsegment is a contiguous segment around the support of all the elements
* @param EltsIdx
* @return
*/
std::vector<Int32> CLineModelElementList::ReestimateContinuumUnderLines(std::vector<Int32> EltsIdx)
{
    if(EltsIdx.size()<1){
        std::vector<Int32> empty;
        return empty;
    }

    //smoothing factor in continuum median filter
    Float64 smoof = 150;

    //modify m_ContinuumFluxAxis
    CSpectrumFluxAxis& fluxAxisModified = m_ContinuumFluxAxis;
    Float64* Ycont = fluxAxisModified.GetSamples();
    CSpectrumSpectralAxis spectralAxis = m_SpectrumModel->GetSpectralAxis();


    std::vector<Int32> xInds = getSupportIndexes( EltsIdx );
    Int32 minInd = xInds[0];
    Int32 maxInd = xInds[xInds.size()-1];
    //
    Int32 iminMerge = spectralAxis.GetIndexAtWaveLength(spectralAxis[minInd]-2*smoof);
    if(iminMerge==0){
        iminMerge = 1;
    }
    Int32 imaxMerge = spectralAxis.GetIndexAtWaveLength(spectralAxis[maxInd]+2*smoof);
    //
    //
    Int32 iminMerge2 = spectralAxis.GetIndexAtWaveLength(spectralAxis[minInd]-smoof);
    if(iminMerge2==0){
        iminMerge2 = 1;
    }
    Int32 imaxMerge2 = spectralAxis.GetIndexAtWaveLength(spectralAxis[maxInd]+smoof);
    //
    Int32 imin = spectralAxis.GetIndexAtWaveLength(spectralAxis[minInd]-3*smoof);
    if(imin==0){
        imin = 1;
    }
    Int32 imax = spectralAxis.GetIndexAtWaveLength(spectralAxis[maxInd]+3*smoof);
    Int32 spcSize = imax-imin+1;


    //prepare the line model with no continuum;
    CSpectrumFluxAxis modelFluxAxisTmp = m_SpectrumModel->GetFluxAxis();
    for(Int32 i=0; i<spcSize; i++){
        modelFluxAxisTmp[i+imin]=0.0;
    }
    for(Int32 idx=0; idx<EltsIdx.size(); idx++){
        Int32 eltIdx = EltsIdx[idx];
        m_Elements[eltIdx]->addToSpectrumModel(spectralAxis, modelFluxAxisTmp, m_Redshift);

    }

    //gather the modified indexes
    std::vector<Int32> modifiedIdxs;

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
        fprintf( fspc, "%f %f %f\n", t, spcBuffer.GetSpectralAxis()[t], (m_SpcFluxAxis[t+imin])*coeffSaveSpc, (spcBuffer.GetFluxAxis()[t])*coeffSaveSpc);//*1e12);
    }
    fclose( fspc );
    //*/

    //apply continuum routine on this spcbuffer
    CContinuumIrregularSamplingMedian continuum;
    CSpectrumFluxAxis fluxAxisWithoutContinuumCalc;
    Int32 retVal = continuum.RemoveContinuum( spcBuffer, fluxAxisWithoutContinuumCalc );


    /*
    // export for debug
    FILE* f = fopen( "continuum_reestimated_underlines_dbg.txt", "w+" );
    Float64 coeffSave = 1e16;
    for( Int32 t=iminMerge;t<imaxMerge;t++)
    {
        fprintf( f, "%f %f %f\n", t, spectralAxis[t], (m_SpcContinuumFluxAxis[t])*coeffSave, (spcBuffer.GetFluxAxis()[t-imin] - fluxAxisWithoutContinuumCalc[t-imin])*coeffSave);//*1e12);
    }
    fclose( f );
    //*/

    Float64 modified = 0.0;
    Float64 coeff=0.0;
    //merge raw continuum free with the newly calculated cont. under the line, (todo: with cross-fade on the borders)
    for(Int32 i=iminMerge; i<imaxMerge; i++){
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


std::vector<Int32> CLineModelElementList::ReestimateContinuumApprox(std::vector<Int32> EltsIdx)
{
    //smoothing factor in continuum median filter
    Float64 smoof = 150;

    //
    std::vector<Int32> modifiedIdxs;


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
            //integratedA /= 5.0;
            //Float64 integratedA = A/8.0;

            Float64 term=0.0;
            for(Int32 i=imin; i<imax; i++){
                Float64 dx=spectralAxis[imin]-spectralAxis[imin-1];
                if(i>smin && i<smax){
                    term = coeffA/dx;
                }
                else if(i>imin && i<=smin){
                    //term = integratedA/(Float64)sSize/dx;
                    term = ((Float64(i-imin)/Float64(smin-imin))) * coeffA/dx;
                }else if(i>=smax && i<imax){
                    //term = integratedA/(Float64)sSize/dx;
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

void CLineModelElementList::refreshModelAfterContReestimation(std::vector<Int32> xInds, CSpectrumFluxAxis& modelFluxAxis, CSpectrumFluxAxis& spcFluxAxisNoContinuum)
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

Float64 CLineModelElementList::getLeastSquareMerit(const TFloat64Range& lambdaRange)
{
    const CSpectrumSpectralAxis& spcSpectralAxis = m_SpectrumModel->GetSpectralAxis();
    const CSpectrumFluxAxis& spcFluxAxis = m_SpcFluxAxis;
    const CSpectrumFluxAxis& modelFluxAxis = m_SpectrumModel->GetFluxAxis();

    Int32 numDevs = 0;
    Float64 fit = 0.0;
    const Float64* Ymodel = modelFluxAxis.GetSamples();
    const Float64* Yspc = spcFluxAxis.GetSamples();
    Float64 diff = 0.0;

    Float64 imin = spcSpectralAxis.GetIndexAtWaveLength(lambdaRange.GetBegin());
    Float64 imax = spcSpectralAxis.GetIndexAtWaveLength(lambdaRange.GetEnd());

    for( UInt32 j=imin; j<imax; j++ )
    {
        numDevs++;
        // fit
        diff = (Yspc[j] - Ymodel[j]);
        fit += (diff*diff) / (m_ErrorNoContinuum[j]*m_ErrorNoContinuum[j]);
        //fit += (diff*diff) / (1e-16*1e-16);
        //fit += (diff*diff)/ (error[0]*error[0]);
        //fit += pow( Yspc[j] - Ymodel[j] , 2.0 );
    }
    //fit /= numDevs;

    return fit;
}

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
            // fit
            diff = (Yspc[j] - Ymodel[j]);
            fit += (diff*diff) / (m_ErrorNoContinuum[j]*m_ErrorNoContinuum[j]);
            //fit += diff*diff;

        }
    }


    return fit;
}

Float64 CLineModelElementList::getModelErrorUnderElement(Int32 eltId)
{
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
            // fit
            diff = (Yspc[j] - Ymodel[j]);
            w = 1.0 / (m_ErrorNoContinuum[j]*m_ErrorNoContinuum[j]);
            fit += (diff*diff) * w;
            //fit += diff*diff;
            sumErr += w;
        }
    }


    return sqrt(fit/sumErr);
}

std::vector<int> CLineModelElementList::findLineIdxInCatalog(const CRayCatalog::TRayVector& restRayList, std::string strTag, Int32 type)
{
    std::vector<Int32> indexes;
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

void CLineModelElementList::addSingleLine(const CRay &r, Int32 index, Float64 nominalWidth)
{
    //CSingleLine line = CSingleLine(r, m_LineWidthType, nominalWidth);
    std::vector<Int32> a;
    a.push_back(index);
    //CSingleLine c(r, nominalWidth, a);
    m_Elements.push_back(boost::shared_ptr<CLineModelElement> (new CSingleLine(r, m_LineWidthType, m_resolution, m_velocityEmission, m_velocityAbsorption, nominalWidth, a)));
    //m_Elements.push_back(new CSingleLine(r, m_LineWidthType, nominalWidth, a));
}

void CLineModelElementList::addDoubleLine(const CRay &r1, const CRay &r2, Int32 index1, Int32 index2, Float64 nominalWidth, Float64 a1, Float64 a2)
{
    std::vector<CRay> lines;
    lines.push_back(r1);
    lines.push_back(r2);

    std::vector<Float64> amps;
    amps.push_back(a1);
    amps.push_back(a2);

    std::vector<Int32> a;
    a.push_back(index1);
    a.push_back(index2);
    //CSingleLine c(r, nominalWidth, a);
    m_Elements.push_back(boost::shared_ptr<CLineModelElement> (new CMultiLine(lines, m_LineWidthType, m_resolution, m_velocityEmission, m_velocityAbsorption, amps, nominalWidth, a)));
    //m_Elements.push_back(new CSingleLine(r, nominalWidth, a));
}

void CLineModelElementList::applyRules()
{
    if(m_rulesoption=="no"){
        return;
    }

    Bool enableBalmer = m_rulesoption.find("balmer") != std::string::npos;
    if(m_rulesoption=="all" || enableBalmer){

        //*
        ApplyBalmerRuleLinSolve();
        //*/

        //*
        Apply2SingleLinesAmplitudeRule(CRay::nType_Emission, "Halpha", "Hbeta", 1.0/2.86*1.1);
        Apply2SingleLinesAmplitudeRule(CRay::nType_Emission, "Hbeta", "Hgamma", 0.47*1.1);
        Apply2SingleLinesAmplitudeRule(CRay::nType_Emission, "Hgamma", "Hdelta", 1.1);
        //Apply2SingleLinesAmplitudeRule(CRay::nType_Emission, "Hdelta", "Hepsilon", 1.1);
        //Apply2SingleLinesAmplitudeRule(CRay::nType_Emission, "Hepsilon", "H8", 1.1);
        Apply2SingleLinesAmplitudeRule(CRay::nType_Emission, "Hdelta", "H8", 1.1);
        Apply2SingleLinesAmplitudeRule(CRay::nType_Emission, "H8", "H9", 1.1);
        Apply2SingleLinesAmplitudeRule(CRay::nType_Emission, "H9", "H10", 1.1);
        Apply2SingleLinesAmplitudeRule(CRay::nType_Emission, "H10", "H11", 1.1);
        //*/

        /*
    Apply2SingleLinesAmplitudeRule(CRay::nType_Absorption, "HbetaA", "HgammaA", 1.0);
    Apply2SingleLinesAmplitudeRule(CRay::nType_Absorption, "HgammaA", "HdeltaA", 1.0);
    Apply2SingleLinesAmplitudeRule(CRay::nType_Absorption, "HdeltaA", "HepsilonA", 1.0);
    Apply2SingleLinesAmplitudeRule(CRay::nType_Absorption, "HepsilonA", "H8A", 1.0);
    Apply2SingleLinesAmplitudeRule(CRay::nType_Absorption, "H8A", "H9A", 1.0);
    Apply2SingleLinesAmplitudeRule(CRay::nType_Absorption, "H9A", "H10A", 1.0);
    Apply2SingleLinesAmplitudeRule(CRay::nType_Absorption, "H10A", "H11A", 1.0);
    //*/
    }

    Bool enableOIIRatioRange= m_rulesoption.find("oiiratio") != std::string::npos;
    if(m_rulesoption=="all" || enableOIIRatioRange){

        //*
        ApplyAmplitudeRatioRangeRule(CRay::nType_Emission, "[OII]3726", "[OII]3729", 2.0);
        //*/
    }

    //*
    //add rule, if OIII present, then OII should be there too
    //*/


    Bool enableStrongWeak= m_rulesoption.find("strongweak") != std::string::npos;
    if(m_rulesoption=="all" || enableStrongWeak){

        //*
        ApplyStrongHigherWeakRule(CRay::nType_Emission);
        //*/

        //*
        ApplyStrongHigherWeakRule(CRay::nType_Absorption);
        //*/
    }
}

Void CLineModelElementList::ApplyStrongHigherWeakRule( Int32 linetype )
{
    Float64 coeff = 1.0;
    Float64 erStrong=-1.0;
    Float64 maxiStrong = FindHighestStrongLineAmp(linetype, erStrong);
    if(maxiStrong == -1){
        return;
    }

    for( UInt32 iRestRayWeak=0; iRestRayWeak<m_RestRayList.size(); iRestRayWeak++ ) //loop on the strong lines
    {
        if(m_RestRayList[iRestRayWeak].GetForce() != CRay::nForce_Weak){
            continue;
        }
        if(m_RestRayList[iRestRayWeak].GetType() != linetype){
            continue;
        }
        Int32 eIdxWeak = FindElementIndex(iRestRayWeak);
        Int32 subeIdxWeak = m_Elements[eIdxWeak]->FindElementIndex(iRestRayWeak);

        if(m_Elements[eIdxWeak]->IsOutsideLambdaRange(subeIdxWeak) == true){
            continue;
        }

        Float64 nSigma = 1.0;
        Float64 ampA = maxiStrong;
        Float64 erA = erStrong;

        Float64 ampB = m_Elements[eIdxWeak]->GetFittedAmplitude(subeIdxWeak);
        Float64 erB = m_Elements[eIdxWeak]->GetFittedAmplitudeErrorSigma(subeIdxWeak);

        Float64 maxB = (coeff*ampA) + coeff*(erA*nSigma);

        m_Elements[eIdxWeak]->LimitFittedAmplitude(subeIdxWeak, maxB);

    }

}

Float64 CLineModelElementList::FindHighestStrongLineAmp( Int32 linetype , Float64 &er)
{
    Float64 maxi = -1.0;
    for( UInt32 iRestRayStrong=0; iRestRayStrong<m_RestRayList.size(); iRestRayStrong++ ) //loop on the strong lines
    {
        if(m_RestRayList[iRestRayStrong].GetForce() != CRay::nForce_Strong){
            continue;
        }
        if(m_RestRayList[iRestRayStrong].GetType() != linetype){
            continue;
        }

        Int32 eIdxStrong = FindElementIndex(iRestRayStrong);
        Int32 subeIdxStrong = m_Elements[eIdxStrong]->FindElementIndex(iRestRayStrong);

        if(m_Elements[eIdxStrong]->IsOutsideLambdaRange(subeIdxStrong) == true){
            continue;
        }


        Float64 ampStrong = m_Elements[eIdxStrong]->GetFittedAmplitude(subeIdxStrong);
        if(maxi<ampStrong){
            maxi = ampStrong;
            er = m_Elements[eIdxStrong]->GetFittedAmplitudeErrorSigma(subeIdxStrong);
        }
    }
    return maxi;
}

Void CLineModelElementList::Apply2SingleLinesAmplitudeRule( Int32 linetype, std::string lineA, std::string lineB, Float64 coeff )
{
    Int32 iA = FindElementIndex(lineA, linetype);
    if(iA==-1){
        return;
    }
    if(m_Elements[iA]->GetSize()>1){
        iA=-1;
    }
    Int32 iB = FindElementIndex(lineB, linetype);
    if(iB==-1){
        return;
    }
    if(m_Elements[iB]->GetSize()>1){
        iB=-1;
    }
    if(iA==-1 || iB==-1 || iA==iB){
        return;
    }

    if(m_Elements[iA]->IsOutsideLambdaRange() == false){
        Float64 nSigma = 1.0;
        Float64 ampA = m_Elements[iA]->GetFittedAmplitude(0);
        Float64 erA = m_Elements[iA]->GetFittedAmplitudeErrorSigma(0);

        Float64 ampB = m_Elements[iB]->GetFittedAmplitude(0);
        Float64 erB = m_Elements[iB]->GetFittedAmplitudeErrorSigma(0);

        /*
        //Method 1, limit the weakest line's amplitude
        Float64 maxB = (coeff*ampA) + (erA*nSigma*coeff);
        m_Elements[iB]->LimitFittedAmplitude(0, maxB);
        //*/

        //*
        //Method 2, correct both lines depending on their sigmas
        if(ampB!=0.0 && (erA!=0 && erB!=0) && std::abs(ampB) > std::abs(ampA*coeff) ){
            Float64 R = 1.0/coeff;
            Float64 wA = 0.0;
            if(erA!=0.0){
                wA = 1.0/(erA*erA);
            }
            Float64 wB = 0.0;
            if(erB!=0.0){
                wB = 1.0/(erB*erB*R*R);
            }
            Float64 correctedA = (ampA*wA + ampB*wB*R)/(wA+wB) ;
            Float64 correctedB = correctedA/R;

            m_Elements[iA]->SetFittedAmplitude(correctedA, erA); //check: keep the original error sigma ?
            m_Elements[iB]->SetFittedAmplitude(correctedB, erB); //check: keep the original error sigma ?
        }else if(ampB!=0.0 && ampA==0.0){
            Float64 maxB = erA;//*nSigma*coeff;
            m_Elements[iB]->LimitFittedAmplitude(0, maxB);
        }

        //*/

    }
}

Void CLineModelElementList::ApplyAmplitudeRatioRangeRule( Int32 linetype, std::string lineA, std::string lineB, Float64 coeff)
{
    Int32 iA = FindElementIndex(lineA, linetype);
    if(iA==-1){
        return;
    }
    if(m_Elements[iA]->GetSize()>1){
        iA=-1;
    }
    Int32 iB = FindElementIndex(lineB, linetype);
    if(iB==-1){
        return;
    }
    if(m_Elements[iB]->GetSize()>1){
        iB=-1;
    }
    if(iA==-1 || iB==-1 || iA==iB){
        return;
    }

    if(m_Elements[iA]->IsOutsideLambdaRange() == false && m_Elements[iB]->IsOutsideLambdaRange() == false){
        Float64 ampA = m_Elements[iA]->GetFittedAmplitude(0);
        Float64 erA = m_Elements[iA]->GetFittedAmplitudeErrorSigma(0);

        Float64 ampB = m_Elements[iB]->GetFittedAmplitude(0);
        Float64 erB = m_Elements[iB]->GetFittedAmplitudeErrorSigma(0);

        Int32 i1 = iA;
        Int32 i2 = iB;
        Float64 amp1 = ampA;
        Float64 er1 = erA;
        Float64 amp2 = ampB;
        Float64 er2 = erB;
        if( std::abs(ampA) > std::abs(ampB*coeff) ){
            i1 = iA;
            i2 = iB;
            amp1 = ampA;
            er1 = erA;
            amp2 = ampB;
            er2 = erB;
        }else if( std::abs(ampB) > std::abs(ampA*coeff) ){
            i1 = iB;
            i2 = iA;
            amp1 = ampB;
            er1 = erB;
            amp2 = ampA;
            er2 = erA;
        }else{
            return;
        }

        Float64 R = coeff;
        Float64 w1 = 0.0;
        if(er1!=0.0){
            w1 = 1.0/(er1*er1);
        }
        Float64 w2 = 0.0;
        if(er2!=0.0){
            w2 = 1.0/(er2*er2*R*R);
        }
        Float64 corrected1 = (amp1*w1 + amp2*w2*R)/(w1+w2) ;
        Float64 corrected2 = corrected1/R;

        m_Elements[i1]->SetFittedAmplitude(corrected1, er1);
        m_Elements[i2]->SetFittedAmplitude(corrected2, er2);

    }
}

Int32 CLineModelElementList::ApplyBalmerRuleLinSolve()
{
    const CSpectrumSpectralAxis& spectralAxis = m_SpectrumModel->GetSpectralAxis();
    std::vector<Float64> lambdax;
    std::vector<Float64> continuumx;
    std::vector<Float64> datax;
    std::vector<Float64> errdatax;
    std::vector<Int32> iEltE;
    std::vector<Int32> iEltA;

    std::vector<std::string> linetags;
    linetags.push_back("Halpha");
    linetags.push_back("Hbeta");
    linetags.push_back("Hgamma");
    linetags.push_back("Hdelta");
    //linetags.push_back("Hepsilon");
    linetags.push_back("H8");
    linetags.push_back("H9");
    linetags.push_back("H10");
    linetags.push_back("H11");

    Int32 ilineE;
    Int32 ilineA;
    std::string tagE;
    std::string tagA;
    TFloat64List ampsE;
    TFloat64List ersE;
    TFloat64List ersFitE;
    TFloat64List ampsA;
    TFloat64List ersA;
    TFloat64List ersFitA;
    TBoolList ampsEwasfitted;
    Int32 nLines = 0;
    for(Int32 itag=0; itag<linetags.size(); itag++){
        tagE = linetags[itag];
        tagA = linetags[itag];
        tagA.append("A");
        ilineE = FindElementIndex(tagE, CRay::nType_Emission);
        ilineA = FindElementIndex(tagA, CRay::nType_Absorption);
        if(ilineE==-1 || ilineA==-1){
            continue;
        }
        Float64 ampE = m_Elements[ilineE]->GetFittedAmplitude(0);
        Float64 erE = m_Elements[ilineE]->GetFittedAmplitudeErrorSigma(0);
        Float64 erFitE = getModelErrorUnderElement(ilineE);
        erE = sqrt(erE*erE+erFitE*erFitE);
        Float64 ampA = m_Elements[ilineA]->GetFittedAmplitude(0);
        Float64 erA = m_Elements[ilineA]->GetFittedAmplitudeErrorSigma(0);
        Float64 erFitA = getModelErrorUnderElement(ilineA);
        erA = sqrt(erA*erA+erFitA*erFitA);
        Float64 amp;
        Float64 er;
        if(ampE>0.0 && ampE>ampA){
            amp = ampE;
            er = erE;
            ampsEwasfitted.push_back(true);
        }else if(ampA>0.0){
            amp = -ampA;
            er = erA;
            ampsEwasfitted.push_back(false);
        }else{
            continue;
        }
        nLines++;
        ampsE.push_back(ampE);
        ersE.push_back(erE);
        ersFitE.push_back(erFitE);
        ampsA.push_back(ampA);
        ersA.push_back(erA);
        ersFitA.push_back(erFitA);
        iEltE.push_back(ilineE);
        iEltA.push_back(ilineA);
        Float64 lambda = m_Elements[ilineE]->GetRays()[0].GetPosition();
        lambdax.push_back(lambda);
        Int32 Idx = spectralAxis.GetIndexAtWaveLength(lambda);
        continuumx.push_back(m_SpcContinuumFluxAxis[Idx]);
        datax.push_back(m_SpcContinuumFluxAxis[Idx] + amp);
        errdatax.push_back(er);
    }

    /*
    TFloat64List coeffs = BalmerModelLinSolve( lambdax, continuumx, datax, errdatax );
    // export for debug
    FILE* fspc = fopen( "BalmerLinSolve_dbg.txt", "w+" );
    Float64 coeffSaveSpc = 1.0;
    for(Int32 i=0; i<nLines; i++){
        Float64 ampRegE = coeffs[0]*lambdax[i] + coeffs[3] ;
        Float64 coeffRegA = (coeffs[2]*lambdax[i] + coeffs[1]);
        Float64 coeffA = (continuumx[i] - ampsA[i])/continuumx[i];
        fprintf( fspc, "%d %f %f %f %f %f %f\n", i, lambdax[i], errdatax[i], ampRegE, ampsE[i], coeffRegA, coeffA);
    }
    fclose( fspc );
    //*/

    //apply the absorption rule
    for(Int32 j=0; j<nLines; j++){
        Int32 i1=nLines-1-j; //go through the lines, small wavelengths first
        Int32 i2 = i1-1;
        if(i2<0){
            break;
        }
        Float64 coeffA1 = (continuumx[i1] - ampsA[i1])/continuumx[i1];
        Float64 coeffA2 = (continuumx[i2] - ampsA[i2])/continuumx[i2];
        Float64 erA1 = ersA[i1]; //todo, check if this error value can be calculated better !
        Float64 erA2 = ersA[i2]; //...

        if(coeffA1<coeffA2){ //A1 should absorb less than A2
            Float64 R = 1.0;
            Float64 wA1 = 0.0;
            if(erA1!=0.0){
                wA1 = 1.0/(erA1*erA1);
            }
            Float64 wA2 = 0.0;
            if(erA2!=0.0){
                wA2 = 1.0/(erA2*erA2*R*R);
            }
            Float64 correctedCoeffA1 = (coeffA1*wA1 + coeffA2*wA2*R)/(wA1+wA2) ;
            Float64 correctedCoeffA2 = correctedCoeffA1/R;

            {  //make AL and EL higher...
                Float64 correctedAmpA2 = (1.0 - correctedCoeffA2)*continuumx[i2];
//                if(correctedAmpA2 > ampsE[i2]){
//                    correctedAmpA2 = ampsE[i2];
//                }
                if(correctedAmpA2 > 0){
                    Float64 correction = correctedAmpA2-ampsA[i2];
                    Float64 widthRatioAE = 1.0/3.0;
                    ampsA[i2] = ampsA[i2]+correction;
                    if(ampsEwasfitted[i2]){
                        ampsE[i2] = ampsE[i2]+correction;
                        m_Elements[iEltE[i2]]->SetFittedAmplitude(ampsE[i2]+ampsA[i2]*widthRatioAE, ersE[i2]);
                    }
                    m_Elements[iEltA[i2]]->SetFittedAmplitude(ampsA[i2], ersA[i2]);


                    //than re-fit the lines together, because the correction applied so far is for R=1.0, it could be R>1.0
                    if(ampsEwasfitted[i2]){
                        Float64 AStepRatio = 0.15;
                        Int32 maxIterations = 6;
                        Int32 iterations = 0;
                        Float64 ACorrected=ampsA[i2];

                        Float64 ValidAcorrected = ampsA[i2];
                        Float64 ValidEcorrected = ampsE[i2];

                        Float64 overlap = 0.33;
                        std::vector<Int32> indexesFitted;
                        std::vector<Int32> EOverlapIdx = getOverlappingElements(iEltE[i2], indexesFitted, overlap);
                        std::vector<Int32> EOverlapIdxA = getOverlappingElements(iEltA[i2], indexesFitted, overlap);
                        for(Int32 io=0; io<EOverlapIdxA.size(); io++){
                            EOverlapIdx.push_back(EOverlapIdxA[io]);
                        }
                        refreshModelUnderElements(EOverlapIdx);
                        Float64 previousFitErr = getModelErrorUnderElement(iEltE[i2]);

                        while(iterations<maxIterations){
                            iterations++;
                            ACorrected=ACorrected*(1.0+AStepRatio);
                            Float64 correction = ACorrected-ampsA[i2];
                            m_Elements[iEltA[i2]]->SetFittedAmplitude(ACorrected, ersA[i2]);
                            m_Elements[iEltE[i2]]->SetFittedAmplitude(ampsE[i2]+correction+ACorrected*widthRatioAE, ersE[i2]);
                            refreshModelUnderElements(EOverlapIdx);
                            Float64 newFitErr = getModelErrorUnderElement(iEltE[i2]);

                            if(newFitErr >= previousFitErr){ //todo put a ratio threshold ?
                                iterations = maxIterations;
                            }else{
                                ValidAcorrected = ACorrected;
                                ValidEcorrected = ampsE[i2]+correction;
                            }
                        }
                        ampsA[i2] = ValidAcorrected;
                        ampsE[i2] = ValidEcorrected+ValidAcorrected*widthRatioAE;
                        m_Elements[iEltA[i2]]->SetFittedAmplitude(ValidAcorrected, ersA[i2]);
                        m_Elements[iEltE[i2]]->SetFittedAmplitude(ValidEcorrected+ValidAcorrected*widthRatioAE, ersE[i2]);
                    }
                }
            }
        }
    }

    /*
    TFloat64List coeffs = BalmerModelLinSolve( lambdax, continuumx, datax, errdatax );
    // export for debug
    FILE* fspc = fopen( "BalmerLinSolve_dbg.txt", "w+" );
    Float64 coeffSaveSpc = 1.0;
    for(Int32 i=0; i<nLines; i++){
        Float64 ampRegE = coeffs[0]*lambdax[i] + coeffs[3] ;
        Float64 coeffRegA = (coeffs[2]*lambdax[i] + coeffs[1]);
        Float64 coeffA = (continuumx[i] - ampsA[i])/continuumx[i];
        fprintf( fspc, "%d %f %f %f %f %f %f\n", i, lambdax[i], errdatax[i], ampRegE, ampsE[i], coeffRegA, coeffA);
    }
    fclose( fspc );
    //*/

//    for(Int32 i=0; i<nLines; i++){
//        Float64 ampRegE = coeffs[0]*lambdax[i] + coeffs[3] ;
//        Float64 sigma = errdatax[i];
//        SetElementAmplitude(iEltE[i], ampRegE, sigma);

//        Float64 aa = (1-(coeffs[2]*lambdax[i] + coeffs[1]))*continuumx[i];
//        Float64 ampA = m_Elements[iEltA[i]]->GetFittedAmplitude(0);
//        //Float64 aa = (1-(coeffs[1]))*continuumx[i];
//        SetElementAmplitude(iEltA[i], aa, sigma);
//    }

    return -1;
}


TFloat64List CLineModelElementList::BalmerModelLinSolve( std::vector<Float64> lambdax, std::vector<Float64> continuumx, std::vector<Float64> datax, std::vector<Float64> errdatax )
{
    //Linear fit
    int i, n;
    Float64 fval;
    double yi, ei, chisq;
    gsl_matrix *X, *cov;
    gsl_vector *y, *w, *c;

    n = lambdax.size();
    Int32 nddl = 4;
    if(n<nddl){
        TFloat64List empty;
        return empty;
    }

    X = gsl_matrix_alloc (n, nddl);
    y = gsl_vector_alloc (n);
    w = gsl_vector_alloc (n);

    c = gsl_vector_alloc (nddl);
    cov = gsl_matrix_alloc (nddl, nddl);

    //
    //
    for (i = 0; i < n; i++)
    {
        yi = datax[i];
        ei = errdatax[i];

        for (Int32 iddl = 0; iddl < nddl; iddl++)
        {
            if(iddl==0){
               fval = lambdax[i];
            }else if(iddl==1){
                fval = continuumx[i];
            }else if(iddl==2){
                fval = lambdax[i]*continuumx[i];
            }else if(iddl==3){
                fval = 1.0;
            }
            gsl_matrix_set (X, i, iddl, fval);
        }

        gsl_vector_set (y, i, yi);
        gsl_vector_set (w, i, 1.0/(ei*ei));
    }

    //
    {
      gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, nddl);
      gsl_multifit_wlinear (X, w, y, c, cov, &chisq, work);
      gsl_multifit_linear_free (work);
    }

    //

#define C(i) (gsl_vector_get(c,(i)))
#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))
    if(1){
        Log.LogInfo("# best fit: Y = %g X + %g CX + %g XCX + %g", C(0), C(1), C(2), C(3));
//        Log.LogInfo("# covariance matrix:\n");
//        Log.LogInfo("[ %+.5e, %+.5e \n", COV(0,0), COV(0,1));
//        Log.LogInfo("  %+.5e, %+.5e \n", COV(1,0), COV(1,1));

//        Log.LogInfo("[ %+.5e, %+.5e, %+.5e  \n", COV(0,0), COV(0,1), COV(0,2));
//        Log.LogInfo("  %+.5e, %+.5e, %+.5e  \n", COV(1,0), COV(1,1), COV(1,2));
//        Log.LogInfo("  %+.5e, %+.5e, %+.5e ]\n", COV(2,0), COV(2,1), COV(2,2));

        Log.LogInfo("# chisq/n = %g", chisq/n);
    }

    TFloat64List coeffs;
    coeffs.push_back(gsl_vector_get(c,(0)));
    coeffs.push_back(gsl_vector_get(c,(1)));
    coeffs.push_back(gsl_vector_get(c,(2)));
    coeffs.push_back(gsl_vector_get(c,(3)));

    gsl_matrix_free (X);
    gsl_vector_free (y);
    gsl_vector_free (w);
    gsl_vector_free (c);
    gsl_matrix_free (cov);

    return coeffs;
}

CLineModelResult::SLineModelSolution CLineModelElementList::GetModelSolution()
{
    CLineModelResult::SLineModelSolution modelSolution;
    modelSolution.nDDL = GetModelNonZeroElementsNDdl();

    for( UInt32 iRestRay=0; iRestRay<m_RestRayList.size(); iRestRay++ )
    {
        Int32 eIdx = FindElementIndex(iRestRay);
        Int32 subeIdx = m_Elements[eIdx]->FindElementIndex(iRestRay);
        modelSolution.Rays.push_back(m_RestRayList[iRestRay]);
        modelSolution.ElementId.push_back( eIdx );
        modelSolution.Amplitudes.push_back(m_Elements[eIdx]->GetFittedAmplitude(subeIdx));
        modelSolution.Errors.push_back(m_Elements[eIdx]->GetFittedAmplitudeErrorSigma(subeIdx));
        modelSolution.FittingError.push_back(getModelErrorUnderElement(eIdx));

        //modelSolution.Widths.push_back(-1.0);
        //modelSolution.OutsideLambdaRange.push_back(true);
    }

    return modelSolution;
}


Int32 CLineModelElementList::GetNElements()
{
    Int32 nddl = m_Elements.size();
    return nddl;
}

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

std::vector<Int32> CLineModelElementList::GetModelValidElementsIndexes()
{
    std::vector<Int32> nonZeroIndexes;
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

Int32 CLineModelElementList::FindElementIndex(std::string LineTagStr, Int32 linetype)
{
    Int32 idx = -1;
    for( UInt32 iElts=0; iElts<m_Elements.size(); iElts++ )
    {
        if(m_RestRayList[m_Elements[iElts]->m_LineCatalogIndexes[0]].GetType() != linetype){
            continue;
        }

        if(m_Elements[iElts]->FindElementIndex(LineTagStr) !=-1){
            idx = iElts;
            break;
        }
    }
    return idx;
}


void CLineModelElementList::SetElementAmplitude(Int32 j, Float64 a, Float64 snr)
{
    if(j>=0 && j<m_Elements.size())
    {
        m_Elements[j]->SetFittedAmplitude(a, snr);
    }
    return;
}

Float64 CLineModelElementList::GetElementAmplitude(Int32 j)
{
    Float64 a=-1.0;
    if(j>=0 && j<m_Elements.size())
    {
        a = m_Elements[j]->GetElementAmplitude();
    }
    return a;
}


void CLineModelElementList::SetVelocityEmission(Float64 vel)
{
    m_velocityEmission = vel;
    for(Int32 j=0; j<m_Elements.size(); j++)
    {
        m_Elements[j]->SetVelocityEmission(vel);
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

Float64 CLineModelElementList::GetVelocityEmission()
{
    return m_velocityEmission;
}
Float64 CLineModelElementList::GetVelocityAbsorption()
{
    return m_velocityAbsorption;
}


//this function estimates the continuum after removal(interpolation) of the flux samples under the lines for a given redshift
void CLineModelElementList::EstimateSpectrumContinuum()
{
    std::vector<Int32> validEltsIdx = GetModelValidElementsIndexes();
    std::vector<Int32> xInds = getSupportIndexes( validEltsIdx );
    //m_SpcContinuumFluxAxis = m_SpcFluxAxis;

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
        fprintf( fspc, "%f %f %f\n", t, spectralAxis[t], (m_SpcFluxAxis[t])*coeffSaveSpc, (spcmodel4linefitting.GetFluxAxis()[t])*coeffSaveSpc);//*1e12);
    }
    fclose( fspc );
    //*/

    // Remove continuum
    CContinuumIrregularSamplingMedian continuum;
    CSpectrumFluxAxis fluxAxisWithoutContinuumCalc;
    Int32 retVal = continuum.RemoveContinuum( spcCorrectedUnderLines, fluxAxisWithoutContinuumCalc );

    CSpectrumFluxAxis fluxAxisNewContinuum;
    fluxAxisNewContinuum.SetSize( fluxAxisNothingUnderLines.GetSamplesCount() );

    for( Int32 t=0;t<spectralAxis.GetSamplesCount();t++)
    {
        fluxAxisNewContinuum[t] = fluxAxisNothingUnderLines[t];
    }
    fluxAxisNewContinuum.Subtract(fluxAxisWithoutContinuumCalc);

//    for (int i = 0; i < fluxAxisWithoutContinuumCalc.GetSamplesCount(); i++)
//    {
//        m_SpcFluxAxis[i]=fluxAxisWithoutContinuumCalc[i];
//    }


    //return;

    /*
    // export for debug
    FILE* f = fopen( "continuum_estimated_dbg.txt", "w+" );
    Float64 coeffSave = 1e16;
    for( Int32 t=0;t<spectralAxis.GetSamplesCount();t++)
    {
        fprintf( f, "%f %f %f\n", t, spectralAxis[t], (m_SpcFluxAxis[t])*coeffSave, (fluxAxisNewContinuum[t])*coeffSave);//*1e12);
    }
    fclose( f );
    //*/

    /*
    // export for debug
    FILE* f = fopen( "continuumfree_estimated_dbg.txt", "w+" );
    Float64 coeffSave = 1e16;
    for( Int32 t=0;t<spectralAxis.GetSamplesCount();t++)
    {
        fprintf( f, "%f %f %f\n", t, spectralAxis[t], (m_SpcFluxAxis[t]-m_SpcContinuumFluxAxis[t])*coeffSave, (m_SpcFluxAxis[t]-fluxAxisNewContinuum[t])*coeffSave);//*1e12);
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
