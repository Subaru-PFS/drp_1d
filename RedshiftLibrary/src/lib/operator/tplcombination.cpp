#include <RedshiftLibrary/operator/tplcombination.h>
#include <RedshiftLibrary/operator/templatefittingresult.h>

#include <RedshiftLibrary/spectrum/axis.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/tools.h>
#include <RedshiftLibrary/common/mask.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/common/quicksort.h>
#include <RedshiftLibrary/common/indexing.h>

#include <RedshiftLibrary/processflow/datastore.h>
#include <RedshiftLibrary/log/log.h>

#include <boost/numeric/conversion/bounds.hpp>
#include <math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <algorithm>    // std::sort
#include <float.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include <sstream>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/chrono/thread_clock.hpp>
//#include <boost/progress.hpp>
#include <assert.h>

#define NOT_OVERLAP_VALUE NAN
#include <stdio.h>

namespace bfs = boost::filesystem;
using namespace NSEpic;
using namespace std;


void COperatorTplcombination::BasicFit_preallocateBuffers(const CSpectrum& spectrum, const std::vector<CTemplate>& tplList)
{

    // Pre-Allocate the rebined template and mask with regard to the spectrum size
    m_templatesRebined_bf.resize(tplList.size());
    m_masksRebined_bf.resize(tplList.size());
    m_spcSpectralAxis_restframe.SetSize(spectrum.GetSampleCount());

    for(Int32 ktpl=0; ktpl<tplList.size(); ktpl++)
    {
        m_templatesRebined_bf[ktpl].GetSpectralAxis().SetSize(spectrum.GetSampleCount());
        m_templatesRebined_bf[ktpl].GetFluxAxis().SetSize(spectrum.GetSampleCount());
        m_masksRebined_bf[ktpl].SetSize(spectrum.GetSampleCount());

        m_templatesRebined_bf[ktpl].m_ismCorrectionCalzetti = tplList[ktpl].m_ismCorrectionCalzetti;
        m_templatesRebined_bf[ktpl].m_igmCorrectionMeiksin = tplList[ktpl].m_igmCorrectionMeiksin;
    }
}

// @BasicFit@ in its current content and the commented code corresponds to @::computeModelSpectrum@ since there is no
// loops over igm/ism corrections

void COperatorTplcombination::BasicFit(const CSpectrum& spectrum,
                                       const std::vector<CTemplate>& tplList,
                                       const TFloat64Range& lambdaRange,
                                       Float64 redshift,
                                       Float64 overlapThreshold,
                                       STplcombination_basicfitresult& fittingResults,
                                       std::string opt_interp,
                                       Float64 forcedAmplitude,
                                       Int32 opt_extinction,
                                       Int32 opt_dustFitting,
                                       CMask spcMaskAdditional)
{
    bool verbose = false;
    if(verbose)
    {
        Log.LogDebug("  Operator-tplcombination: BasicFit - for z=%f", redshift);
    }
    boost::chrono::thread_clock::time_point start_prep = boost::chrono::thread_clock::now();


    fittingResults.chisquare = boost::numeric::bounds<float>::highest();
    bool status_chisquareSetAtLeastOnce = false;
    //why is it hardcoded?
    Int32 nISM = 1; //no ISM fitting for now
    Int32 nIGM = 1; //no IGM fitting for now
    for(Int32 kism=0; kism<nISM; kism++)
    {

        TFloat64List _chi2List(nIGM, DBL_MAX);
        fittingResults.ChiSquareInterm.push_back(_chi2List);

    }

    fittingResults.fittingAmplitudes = std::vector<Float64>(tplList.size(), -1.);
    fittingResults.fittingErrors = std::vector<Float64>(tplList.size(), -1.);
    fittingResults.overlapRate = 0.0;
    fittingResults.status = COperator::nStatus_DataError;

    const CSpectrumSpectralAxis& spcSpectralAxis = spectrum.GetSpectralAxis();
    const CSpectrumFluxAxis& spcFluxAxis = spectrum.GetFluxAxis();
    const CSpectrumNoiseAxis& spcError = spcFluxAxis.GetError();

    if(spcMaskAdditional.GetMasksCount()!=spcFluxAxis.GetSamplesCount())
    {
        Log.LogError("  Operator-Tplcombination: spcMaskAdditional does not have the same size as the spectrum flux vector... (%d vs %d), aborting", spcMaskAdditional.GetMasksCount(), spcFluxAxis.GetSamplesCount());
        throw std::runtime_error("  Operator-Tplcombination: spcMaskAdditional does not have the same size as the spectrum flux vector. aborting");       
    }

    TFloat64Range currentRange;   
    Int32 ret = RebinTemplate(spectrum, tplList, 
                                redshift, lambdaRange, opt_interp,
                                currentRange,
                                fittingResults.overlapRate,
                                overlapThreshold); 
    if( ret == -1 ){
        fittingResults.status = COperator::nStatus_NoOverlap; 
        return;
    }
    /*if( ret == -2 ){//this is not coded
        fittingResults.status = COperator::nStatus_DataError;
        return;
    }*/

    // Linear fit
    Int32 imin_lbda = m_spcSpectralAxis_restframe.GetIndexAtWaveLength(currentRange.GetBegin());
    Int32 imax_lbda = m_spcSpectralAxis_restframe.GetIndexAtWaveLength(currentRange.GetEnd());

    Int32 n = imax_lbda - imin_lbda + 1;
    Log.LogDebug("  Operator-Tplcombination: prep. linear fitting with n=%d samples in the clamped lambdarange spectrum (imin=%d, lbda_min=%.3f - imax=%d, lbda_max=%.3f)", n, imin_lbda, spcSpectralAxis[imin_lbda], imax_lbda, spcSpectralAxis[imax_lbda]);
    Int32 nddl=tplList.size();
    gsl_matrix *X, *cov;
    gsl_vector *y, *w, *c;

    X = gsl_matrix_alloc (n, nddl);
    y = gsl_vector_alloc (n);
    w = gsl_vector_alloc (n);
    c = gsl_vector_alloc (nddl);
    cov = gsl_matrix_alloc (nddl, nddl);

    // Normalizing factor
    Float64 normFactor;
    Float64 maxabsval = DBL_MIN;
    for (Int32 k = 0; k < n; k++)
    {
        if(maxabsval<std::abs(spcFluxAxis[k+imin_lbda]))
        {
            maxabsval=std::abs(spcFluxAxis[k+imin_lbda]);
        }
    }
    normFactor = maxabsval;
    if(verbose)
    {
        Log.LogDetail("  Operator-Tplcombination: Linear fitting, found normalization Factor=%e", normFactor);
    }

    // Prepare the fit data
    Float64 yi, ei, chisq;
    for(Int32 i = 0; i < n; i++)
    {
        yi = spcFluxAxis[i+imin_lbda]/normFactor;
        ei = spcError[i+imin_lbda]/normFactor;

        for (Int32 iddl = 0; iddl < nddl; iddl++)
        {
            Float64 fval =  m_templatesRebined_bf[iddl].GetFluxAxis()[i+imin_lbda];
            gsl_matrix_set (X, i, iddl, fval);

            if(0 && verbose)
            {
                Log.LogDebug("  Operator-Tplcombination: Linear fitting, component[%d]_fval = %.3e", iddl, fval);
            }
        }

        gsl_vector_set (y, i, yi);
        gsl_vector_set (w, i, 1.0/(ei*ei));
    }

    // Now fitting
    boost::chrono::thread_clock::time_point stop_prep = boost::chrono::thread_clock::now();
    Float64 duration_prep = boost::chrono::duration_cast<boost::chrono::microseconds>(stop_prep - start_prep).count();
    Log.LogDebug( "  Operator-Tplcombination: Linear fitting, preparation time = %.3f microsec", duration_prep);
    boost::chrono::thread_clock::time_point start_fit = boost::chrono::thread_clock::now();
    {
      gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, nddl);
      gsl_multifit_wlinear (X, w, y, c, cov, &chisq, work);
      gsl_multifit_linear_free (work);
    }
    //
    boost::chrono::thread_clock::time_point stop_fit = boost::chrono::thread_clock::now();
    Float64 duration_fit = boost::chrono::duration_cast<boost::chrono::microseconds>(stop_fit - start_fit).count();
    Log.LogDebug( "  Operator-Tplcombination: Linear fitting, fit = %.3f microsec", duration_fit);
    boost::chrono::thread_clock::time_point start_postprocess = boost::chrono::thread_clock::now();

#define C(i) (gsl_vector_get(c,(i)))
#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))
    if(verbose)
    {
        if(1){
            Log.LogInfo("# best fit: Y = %g X1 + %g X2 ...", C(0), C(1));
            Log.LogInfo("# covariance matrix:");
            Log.LogInfo("[");
            Log.LogInfo("  %+.5e, %+.5e", COV(0,0), COV(0,1));
            Log.LogInfo("  %+.5e, %+.5e", COV(1,0), COV(1,1));
            Log.LogInfo("]");
            Log.LogInfo("# chisq/n = %g", chisq/n);
            }

        for (Int32 iddl = 0; iddl < nddl; iddl++)
        {
            Float64 a = gsl_vector_get(c,iddl)*normFactor;
            Log.LogInfo("# Found amplitude %d: %+.5e +- %.5e", iddl, a, COV(iddl,iddl)*normFactor);
        }
    }

    //save the fitted amps and fitErrors, etc...
    Float64 a, err;
    for (Int32 iddl = 0; iddl < nddl; iddl++)
    {
        a = gsl_vector_get(c,iddl)*normFactor;
        fittingResults.fittingAmplitudes[iddl]=a;
        err = COV(iddl,iddl)*normFactor;
        fittingResults.fittingErrors[iddl] = err;
    }
    if(fittingResults.fittingAmplitudes.size()!=nddl)
    {
        Log.LogDebug("  Operator-Tplcombination: Found nfittedamps(=%d) different than nddl(=%d)", fittingResults.fittingAmplitudes.size(), nddl);
    }

    //build the model
    TFloat64List spc_extract = TFloat64List(spcSpectralAxis.GetSamplesVector().begin() + imin_lbda, 
                                            spcSpectralAxis.GetSamplesVector().begin() + imax_lbda);
    TFloat64List modelFlux(n, 0.0);
    for (Int32 iddl = 0; iddl < nddl; iddl++){
        TFloat64List tmp = TFloat64List( m_templatesRebined_bf[iddl].GetFluxAxis().GetSamplesVector().begin() + imin_lbda, 
                                         m_templatesRebined_bf[iddl].GetFluxAxis().GetSamplesVector().begin() + imax_lbda);
        Float32 a = fittingResults.fittingAmplitudes[iddl];
        for(Int32 k = 0; k<n; k++){
            modelFlux[k]+= a*tmp[k];
        }
    }
    Int32 bestkDust = -1, bestidxMeiksin = -1; //these should be the result of iterating over ism/igm values and selecting based on the best chi2 value
    /*//TODO: uncomment as part of #6329
    CSpectrumSpectralAxis xAxis(spc_extract, n);
    CSpectrumFluxAxis yAxis(modelFlux, n);
    CTemplate modelSpec("combination", "", xAxis, yAxis);//m_templatesRebined_bf cant be used, thus we create a new CTemplate to apply corrections

    if ((bestkDust>-1) || (bestidxMeiksin>-1))
        modelSpec.InitIsmIgmConfig();

    modelSpec.SetIsmIgmLambdaRange(currentRange);


    //applyDustCoeff
    if (bestkDust>-1){
        modelSpec.ApplyDustCoeff(bestkDust);
    }
    //applyMeiksin
    if(bestidxMeiksin > -1){
        modelSpec.ApplyMeiksinCoeff(meiksinIdx, redshift);
    }  
    fittingResults.modelSpectrum = modelSpec;
    */
    fittingResults.modelSpectrum = CSpectrum(CSpectrumSpectralAxis(spc_extract), CSpectrumFluxAxis(modelFlux));
    fittingResults.igmIdx = bestidxMeiksin;
    fittingResults.ebmvCoeff = bestkDust;
    fittingResults.snr = -1.0;
       
    //estimate the lst-square brute force
    fittingResults.chisquare = .0;
    Float64 diff, err2;
    for(Int32 k=0; k<n; k++)
    {
        diff = modelFlux[k]-spcFluxAxis[k+imin_lbda];
        err2 = spcError[k+imin_lbda]*spcError[k+imin_lbda];
        fittingResults.chisquare += diff*diff/err2;
    }

    //save the interm chisquares: for now, ism and igm deactivated so that interm chi2=global chi2
    for(Int32 kism=0; kism<fittingResults.ChiSquareInterm.size(); kism++)
    {
        for(Int32 kigm=0; kigm<fittingResults.ChiSquareInterm[kism].size(); kigm++)
        {
            fittingResults.ChiSquareInterm[kism][kigm] = fittingResults.chisquare;
        }
    }

    boost::chrono::thread_clock::time_point stop_postprocess = boost::chrono::thread_clock::now();
    Float64 duration_postprocess = boost::chrono::duration_cast<boost::chrono::microseconds>(stop_postprocess - start_postprocess).count();
    Log.LogDebug( "  Operator-Tplcombination: Linear fitting, postprocess = %.3f microsec", duration_postprocess);

    gsl_matrix_free (X);
    gsl_vector_free (y);
    gsl_vector_free (w);
    gsl_vector_free (c);
    gsl_matrix_free (cov);


    if(status_chisquareSetAtLeastOnce)
    {
        fittingResults.status = COperator::nStatus_OK;
    }else{
        fittingResults.status = COperator::nStatus_LoopError;
    }
}


Int32  COperatorTplcombination::RebinTemplate( const CSpectrum& spectrum,
                                                const std::vector<CTemplate>& tplList,
                                                Float64 redshift,
                                                const TFloat64Range& lambdaRange,
                                                std::string opt_interp,
                                                TFloat64Range& currentRange,
                                                Float64& overlapRate,
                                                Float64 overlapThreshold)
{
    Float64 onePlusRedshift = 1.0 + redshift;

    //shift lambdaRange backward to be in restframe
    TFloat64Range spcLambdaRange_restframe;
    TFloat64Range lambdaRange_restframe( lambdaRange.GetBegin() / onePlusRedshift,
                                         lambdaRange.GetEnd() / onePlusRedshift );

    //redshift in restframe the tgtSpectralAxis,
    m_spcSpectralAxis_restframe.ShiftByWaveLength(spectrum.GetSpectralAxis(), onePlusRedshift, CSpectrumSpectralAxis::nShiftBackward);
    m_spcSpectralAxis_restframe.ClampLambdaRange( lambdaRange_restframe, spcLambdaRange_restframe );

    TFloat64Range intersectedAllLambdaRange(spcLambdaRange_restframe);

    // Now interpolating all the templates
    Log.LogDebug("  Operator-tplcombination: BasicFit - interpolating");
    for(Int32 ktpl=0; ktpl<tplList.size(); ktpl++)
    {
        const CSpectrumSpectralAxis& tplSpectralAxis = tplList[ktpl].GetSpectralAxis();
        const CSpectrumFluxAxis& tplFluxAxis = tplList[ktpl].GetFluxAxis();

        // Compute clamped lambda range over template
        TFloat64Range tplLambdaRange;
        tplSpectralAxis.ClampLambdaRange( lambdaRange_restframe, tplLambdaRange );

        // if there is any intersection between the lambda range of the spectrum and the lambda range of the template
        // Compute the intersected range
        TFloat64Range intersectedLambdaRange( 0.0, 0.0 );
        TFloat64Range::Intersect( tplLambdaRange, spcLambdaRange_restframe, intersectedLambdaRange );

        // find lambda range intersection common to all templates
        intersectedAllLambdaRange.IntersectWith(intersectedLambdaRange);

        CSpectrum & itplTplSpectrum = m_templatesRebined_bf[ktpl];
        CMask & itplMask = m_masksRebined_bf[ktpl];

        tplList[ktpl].Rebin( intersectedLambdaRange, m_spcSpectralAxis_restframe, itplTplSpectrum, itplMask, opt_interp );
        
        const CSpectrumSpectralAxis& itplTplSpectralAxis = itplTplSpectrum.GetSpectralAxis();
        Log.LogDebug("  Operator-Tplcombination: Rebinned template #%d has n=%d samples in lambdarange: %.2f - %.2f", 
                        ktpl, itplTplSpectralAxis.GetSamplesCount(), itplTplSpectralAxis[0], 
                        itplTplSpectralAxis[itplTplSpectralAxis.GetSamplesCount()-1]);

        overlapRate = m_spcSpectralAxis_restframe.IntersectMaskAndComputeOverlapRate( lambdaRange_restframe, itplMask );

        // Check for overlap rate
        if( overlapRate < overlapThreshold || overlapRate<=0.0 )
        {
            return  -1;
        }
    }

    currentRange = intersectedAllLambdaRange;

    Int32 imin_lbda = m_spcSpectralAxis_restframe.GetIndexAtWaveLength(currentRange.GetBegin());
    if(m_spcSpectralAxis_restframe[imin_lbda]<spcLambdaRange_restframe.GetBegin() && imin_lbda+1<m_spcSpectralAxis_restframe.GetSamplesCount())
    {
        imin_lbda += 1;
    }
    Int32 imax_lbda = m_spcSpectralAxis_restframe.GetIndexAtWaveLength(currentRange.GetEnd());
    if(m_spcSpectralAxis_restframe[imax_lbda]>spcLambdaRange_restframe.GetEnd() && imax_lbda-1>=imin_lbda)
    {
        imax_lbda -= 1;
    }
    currentRange.Set(m_spcSpectralAxis_restframe[imin_lbda], m_spcSpectralAxis_restframe[imax_lbda]);
    
    return 0;
}

/**
 * \brief
 *
 * input: if additional_spcMasks size is 0, no additional mask will be used, otherwise its size should match the redshifts list size
 * @lambdaRange is not clamped
 **/
std::shared_ptr<COperatorResult> COperatorTplcombination::Compute(const CSpectrum& spectrum,
                                                                  const std::vector<CTemplate> & tplList,
                                                                  const TFloat64Range& lambdaRange,
                                                                  const TFloat64List& redshifts,
                                                                  Float64 overlapThreshold,
                                                                  std::vector<CMask> additional_spcMasks,
                                                                  std::string opt_interp,
                                                                  Int32 opt_extinction,
                                                                  Int32 opt_dustFitting)
{
    Log.LogInfo("  Operator-tplcombination: starting computation with N-template = %d", tplList.size());

    if( spectrum.GetSpectralAxis().IsInLinearScale() == false )
    {
        Log.LogError("  Operator-tplcombination: input spectrum is not in log scale");
        throw std::runtime_error("  Operator-tplcombination: input spectrum is not in log scale");
    }

    for(Int32 ktpl=0; ktpl<tplList.size(); ktpl++)
    {
        if( tplList[ktpl].GetSpectralAxis().IsInLinearScale() == false )
        {
            Log.LogError("  Operator-tplcombination: input template k=%d are not in log scale", ktpl);
            throw std::runtime_error("  Operator-tplcombination: input template k=%d are not in log scale");
        }
        /*//temporarily commented the time ISM/IGM fitting is coded
        if( opt_dustFitting && tplList[ktpl].m_ismCorrectionCalzetti->calzettiInitFailed)
        {
            Log.LogError("  Operator-tplcombination: no calzetti calib. file loaded... aborting");
            throw std::runtime_error("  Operator-tplcombination: no calzetti calib. file loaded... aborting");
        }
        if( opt_extinction && tplList[ktpl].m_igmCorrectionMeiksin->meiksinInitFailed)
        {
            Log.LogError("  Operator-tplcombination: no meiksin calib. file loaded... aborting");
            throw std::runtime_error("  Operator-tplcombination: no meiksin calib. file loaded... aborting");
        }
        */

    }

    Log.LogDebug("  Operator-tplcombination: allocating memory for buffers (N = %d)", tplList.size());

    BasicFit_preallocateBuffers(spectrum, tplList);

    //sort the redshift and keep track of the indexes
    TFloat64List sortedRedshifts;
    TFloat64List sortedIndexes;
    // This is a vector of {value,index} pairs
    vector<pair<Float64,Int32> > vp;
    vp.reserve(redshifts.size());
    for (Int32 i = 0 ; i < redshifts.size() ; i++) {
        vp.push_back(make_pair(redshifts[i], i));
    }
    std::sort(vp.begin(), vp.end());
    for (Int32 i = 0 ; i < vp.size() ; i++) {
        sortedRedshifts.push_back(vp[i].first);
        sortedIndexes.push_back(vp[i].second);
    }

    Log.LogDebug("  Operator-tplcombination: prepare the results");
    std::shared_ptr<CTemplateFittingResult> result = std::shared_ptr<CTemplateFittingResult>( new CTemplateFittingResult() );
    Int32 nEbmvCoeffs=1;
    /*if(opt_dustFitting)
    {
        nEbmvCoeffs = tplList.front().m_ismCorrectionCalzetti->GetNPrecomputedDustCoeffs();
    }*/
    Log.LogDebug("  Operator-tplcombination: prepare N ism coeffs = %d", nEbmvCoeffs);
    Int32 nIGMCoeffs=1;
    /*if(opt_extinction)
    {
        nIGMCoeffs = tplList.front().m_igmCorrectionMeiksin->GetIdxCount();
    }*/
    Log.LogDebug("  Operator-tplcombination: prepare N igm coeffs = %d", nIGMCoeffs);

    result->Init(sortedRedshifts.size(), nEbmvCoeffs, nIGMCoeffs);
    result->Redshifts = sortedRedshifts;

    CMask additional_spcMask(spectrum.GetSampleCount());
    CMask default_spcMask(spectrum.GetSampleCount());
    //default mask
    for(Int32 km=0; km<default_spcMask.GetMasksCount(); km++)
    {
        default_spcMask[km] = 1.0;
    }
    bool useDefaultMask = 0;
    if(additional_spcMasks.size()!=sortedRedshifts.size())
    {
        useDefaultMask=true;
    }
    if(additional_spcMasks.size()!=sortedRedshifts.size() && additional_spcMasks.size()!=0)
    {
        Log.LogError("  Operator-Tplcombination: using default mask, masks-list size (%d) didn't match the input redshift-list (%d) !)", additional_spcMasks.size(), sortedRedshifts.size());
        throw std::runtime_error("  Operator-Tplcombination: using default mask, masks-list size didn't match the input redshift-list size");

    }

    TFloat64Range clampedlambdaRange;
    spectrum.GetSpectralAxis().ClampLambdaRange(lambdaRange, clampedlambdaRange );

    for (Int32 i=0;i<sortedRedshifts.size();i++)
    {
        //default mask
        if(useDefaultMask)
        {
            additional_spcMask = default_spcMask;
        }else{
            //masks from the input masks list
            additional_spcMask = additional_spcMasks[sortedIndexes[i]];
        }

        Float64 redshift = result->Redshifts[i];

        STplcombination_basicfitresult fittingResults;
        
        BasicFit( spectrum,
                  tplList,
                  clampedlambdaRange,
                  redshift,
                  overlapThreshold,
                  fittingResults,
                  opt_interp,
                  -1,
                  opt_extinction,
                  opt_dustFitting,
                  additional_spcMask);

        if(result->Status[i]==COperator::nStatus_InvalidProductsError)
        {
            Log.LogError("  Operator-Tplcombination: found invalid tplcombination products for z=%f. Now breaking z loop.", redshift);
            throw std::runtime_error("  Operator-Tplcombination: found invalid tplcombination products for z");      
        }

        result->ChiSquare[i]=fittingResults.chisquare;
        result->Overlap[i]=fittingResults.overlapRate;
        result->ChiSquareIntermediate[i]=fittingResults.ChiSquareInterm;

    }

    //overlap warning
    Float64 overlapValidInfZ = -1;
    for (Int32 i=0;i<sortedRedshifts.size();i++)
    {
        if(result->Overlap[i]>=overlapThreshold && overlapValidInfZ==-1)
        {
            overlapValidInfZ=sortedRedshifts[i];
            break;
        }
    }
    Float64 overlapValidSupZ = -1;
    for (Int32 i=sortedRedshifts.size()-1;i>=0;i--)
    {
        if(result->Overlap[i]>=overlapThreshold && overlapValidSupZ==-1)
        {
            overlapValidSupZ=sortedRedshifts[i];
            break;
        }
    }
    if(overlapValidInfZ!=sortedRedshifts[0] || overlapValidSupZ!=sortedRedshifts[sortedRedshifts.size()-1])
    {
        Log.LogInfo("  Operator-Tplcombination: overlap warning for: minz=%.3f, maxz=%.3f", overlapValidInfZ, overlapValidSupZ);
    }

    //only bad status warning
    Int32 oneValidStatusFoundIndex = -1;
    for (Int32 i=0;i<sortedRedshifts.size();i++)
    {
        if(result->Status[i]==COperator::nStatus_OK)
        {
            oneValidStatusFoundIndex=i;
            Log.LogDebug("  Operator-Tplcombination: STATUS VALID found at least at index=%d", i);
            break;
        }
    }if(oneValidStatusFoundIndex==-1)
    {
        Log.LogWarning("  Operator-Tplcombination: STATUS WARNING: Not even one single valid fit/merit value found");
    }


    //loop error status warning
    Int32 loopErrorStatusFoundIndex = -1;
    for (Int32 i=0;i<sortedRedshifts.size();i++)
    {
        if(result->Status[i]==COperator::nStatus_LoopError)
        {
            loopErrorStatusFoundIndex=i;
            Log.LogDebug("  Operator-Tplcombination: STATUS Loop Error found at least at index=%d", i);
            break;
        }
    }if(loopErrorStatusFoundIndex!=-1)
    {
        Log.LogWarning("    Tplcombination-operator: Loop Error - lst-square values not set even once");
    }

    //estimate CstLog for PDF estimation
    result->CstLog = EstimateLikelihoodCstLog(spectrum, clampedlambdaRange);
    
    //store spectrum results
/*    Int32 nMaxExtremaSpectraSave = 1, n = min(int(result->Extrema.size()), nMaxExtremaSpectraSave);
    TFloat64Index indexing;
    for( Int32 ie=0; ie<n; ie++ )
    {
        Int32 i=-1;
        i = indexing.getIndex(result->Redshifts, result->Extrema[ie]);
        if(i==-1){
            Log.LogError("  Operator-Tplcombination: Unable to find the z index for spectrum result save");
            throw std::runtime_error("  Operator-Tplcombination: Unable to find the z index for spectrum result save");
        }
        //default mask
        if(useDefaultMask)
        {
            additional_spcMask = default_spcMask;
        }else{
            //masks from the input masks list
            additional_spcMask = additional_spcMasks[sortedIndexes[i]];
        }

        STplcombination_basicfitresult fittingResults;

        BasicFit( spectrum,
                    tplList,
                    clampedlambdaRange,
                    result->Extrema[ie],
                    overlapThreshold,
                    fittingResults,
                    opt_interp,
                    -1,
                    opt_extinction,
                    opt_dustFitting,
                    additional_spcMask);

        if(result->Status[i]==COperator::nStatus_InvalidProductsError)
        {
            Log.LogError("  Operator-Tplcombination: found invalid tplcombination products for z=%f. Now breaking z loop.", result->Extrema[ie]);
            throw std::runtime_error("  Operator-Tplcombination: found invalid tplcombination products for z");
        }

        Log.LogDetail("  Operator-Tplcombination: creating spectrum/model tplcombination result for z=%f", result->Extrema[ie]);
        for(Int32 itpl=0; itpl<fittingResults.fittingAmplitudes.size(); itpl++)
        {
            Log.LogDetail("  Operator-Tplcombination: creating spectrum/model with tpl-%d amp = %.4e", itpl, fittingResults.fittingAmplitudes[itpl]);
        }
        std::shared_ptr<CModelSpectrumResult>  resultspcmodel = std::shared_ptr<CModelSpectrumResult>( new CModelSpectrumResult(fittingResults.modelSpectrum) );
        m_savedModelSpectrumResults.push_back(resultspcmodel);      
        m_savedModelContinuumFittingResults.push_back(std::make_shared<CModelContinuumFittingResult>(result->Extrema[ie],
                                                                                                    fittingResults.tplNames, 
                                                                                                    fittingResults.chisquare,
                                                                                                    fittingResults.fittingAmplitudes,
                                                                                                    fittingResults.fittingErrors,
                                                                                                    fittingResults.ebmvCoeff,
                                                                                                    fittingResults.igmIdx,
                                                                                                    fittingResults.snr));  

    }*/


    // Deallocate the rebined template and mask buffers
    m_templatesRebined_bf.clear();
    m_masksRebined_bf.clear();

    return result;

}

/**
 * \brief this function estimates the likelihood_cstLog term withing the wavelength range
 **/
Float64 COperatorTplcombination::EstimateLikelihoodCstLog(const CSpectrum& spectrum, const TFloat64Range& lambdaRange)
{
    const CSpectrumSpectralAxis& spcSpectralAxis = spectrum.GetSpectralAxis();
    const TFloat64List& error = spectrum.GetFluxAxis().GetError().GetSamplesVector();

    Int32 numDevs = 0;
    Float64 cstLog = 0.0;
    Float64 sumLogNoise = 0.0;

    Float64 imin = spcSpectralAxis.GetIndexAtWaveLength(lambdaRange.GetBegin());
    Float64 imax = spcSpectralAxis.GetIndexAtWaveLength(lambdaRange.GetEnd());
    for( UInt32 j=imin; j<imax; j++ )
    {
        numDevs++;
        sumLogNoise += log( error[j] );
    }
    //Log.LogDebug( "CLineModelElementList::EstimateMTransposeM val = %f", mtm );

    cstLog = -numDevs*0.5*log(2*M_PI) - sumLogNoise;

    return cstLog;
}

void COperatorTplcombination::SaveSpectrumResults(std::shared_ptr<COperatorResultStore> resultStore)
{
    Log.LogDetail("  Operator-Tplcombination: now saving spectrum/model tplcombination results n=%d", m_savedModelSpectrumResults.size());
    for(Int32 ie=0; ie<m_savedModelSpectrumResults.size(); ie++)
    {
        std::string fname_spc = (boost::format("tplcombinationmodel_spc_extrema_%1%") % ie).str();
       
        resultStore->StoreGlobalResult("tplcombinationsolve",fname_spc.c_str(), m_savedModelSpectrumResults[ie] );
    }
}
