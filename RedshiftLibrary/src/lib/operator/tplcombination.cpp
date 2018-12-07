#include <RedshiftLibrary/operator/tplcombination.h>

#include <RedshiftLibrary/spectrum/axis.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/tools.h>
#include <RedshiftLibrary/common/mask.h>
#include <RedshiftLibrary/operator/chisquareresult.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/common/quicksort.h>

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
#include <boost/progress.hpp>

#include <assert.h>

#define NOT_OVERLAP_VALUE NAN
#include <stdio.h>

namespace bfs = boost::filesystem;
using namespace NSEpic;
using namespace std;


COperatorTplcombination::COperatorTplcombination( std::string calibrationPath )
{

    //ISM
    m_ismCorrectionCalzetti = new CSpectrumFluxCorrectionCalzetti();
    m_ismCorrectionCalzetti->Init(calibrationPath, 0.0, 0.1, 10);
    //Allocate buffer for Ytpl reinit during Dust-fit loop
    m_YtplRawBufferMaxBufferSize = 10*1e6; //allows array from 0A to 100000A with dl=0.01
    m_YtplRawBuffer = new Float64[(int)m_YtplRawBufferMaxBufferSize]();

    //IGM
    m_igmCorrectionMeiksin = new CSpectrumFluxCorrectionMeiksin();
    m_igmCorrectionMeiksin->Init(calibrationPath);

}

COperatorTplcombination::~COperatorTplcombination()
{
    delete[] m_YtplRawBuffer;
    delete m_ismCorrectionCalzetti;
    delete m_igmCorrectionMeiksin;
}


void COperatorTplcombination::BasicFit(const CSpectrum& spectrum,
                                       std::vector<CTemplate>& tplList,
                                       Float64* pfgTplBuffer,
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
    const TFloat64List& spcError = spcFluxAxis.GetError();

    if(spcMaskAdditional.GetMasksCount()!=spcFluxAxis.GetSamplesCount())
    {
        Log.LogInfo("  Operator-Tplcombination: spcMaskAdditional does not have the same size as the spectrum flux vector... (%d vs %d), aborting!", spcMaskAdditional.GetMasksCount(), spcFluxAxis.GetSamplesCount());
        fittingResults.status = COperator::nStatus_DataError;
        return ;
    }
    // Compute clamped lambda range over spectrum
    TFloat64Range spcLambdaRange;
    spcSpectralAxis.ClampLambdaRange( lambdaRange, spcLambdaRange );
    Log.LogDebug("  Operator-Tplcombination: spectrum has n=%d samples in lambdarange: %.2f - %.2f", spcSpectralAxis.GetSamplesCount(), spcSpectralAxis[0], spcSpectralAxis[spcSpectralAxis.GetSamplesCount()-1]);


    // Now interpolating all the templates
    Log.LogDebug("  Operator-tplcombination: BasicFit - interpolating");
    for(Int32 ktpl=0; ktpl<tplList.size(); ktpl++)
    {
        const CSpectrumSpectralAxis& tplSpectralAxis = tplList[ktpl].GetSpectralAxis();
        const CSpectrumFluxAxis& tplFluxAxis = tplList[ktpl].GetFluxAxis();

        // Compute shifted template
        Float64 onePlusRedshift = 1.0 + redshift;
        m_shiftedTemplatesSpectralAxis_bf[ktpl]->ShiftByWaveLength( tplSpectralAxis, onePlusRedshift, CSpectrumSpectralAxis::nShiftForward );
        TFloat64Range intersectedLambdaRange( 0.0, 0.0 );

        // Compute clamped lambda range over template
        TFloat64Range tplLambdaRange;
        m_shiftedTemplatesSpectralAxis_bf[ktpl]->ClampLambdaRange( lambdaRange, tplLambdaRange );

        // if there is any intersection between the lambda range of the spectrum and the lambda range of the template
        // Compute the intersected range
        TFloat64Range::Intersect( tplLambdaRange, spcLambdaRange, intersectedLambdaRange );

        //UInt32 tgtn = spcSpectralAxis.GetSamplesCount() ;
        CSpectrumFluxAxis& itplTplFluxAxis = m_templatesRebined_bf[ktpl]->GetFluxAxis();
        CSpectrumSpectralAxis& itplTplSpectralAxis = m_templatesRebined_bf[ktpl]->GetSpectralAxis();
        CMask& itplMask = *m_masksRebined_bf[ktpl];

        //CSpectrumFluxAxis::Rebin( intersectedLambdaRange, tplFluxAxis, shiftedTplSpectralAxis, spcSpectralAxis, itplTplFluxAxis, itplTplSpectralAxis, itplMask );
        CSpectrumFluxAxis::Rebin2( intersectedLambdaRange, tplFluxAxis, pfgTplBuffer, redshift, *m_shiftedTemplatesSpectralAxis_bf[ktpl], spcSpectralAxis, itplTplFluxAxis, itplTplSpectralAxis, itplMask, opt_interp );

        Log.LogDebug("  Operator-Tplcombination: Rebinned template #%d has n=%d samples in lambdarange: %.2f - %.2f", ktpl, itplTplSpectralAxis.GetSamplesCount(), itplTplSpectralAxis[0], itplTplSpectralAxis[itplTplSpectralAxis.GetSamplesCount()-1]);

        /*//overlapRate, Method 1
        CMask mask;
        spcSpectralAxis.GetMask( lambdaRange, mask );
        itplMask &= mask;
        overlapRate = mask.CompouteOverlapRate( itplMask );
        //*/

        //overlapRate, Method 2
        //CMask mask;
        //spcSpectralAxis.GetMask( lambdaRange, mask );
        //overlapRate = mask.IntersectAndComputeOverlapRate( itplMask );

        //overlapRate, Method 3
        fittingResults.overlapRate = spcSpectralAxis.IntersectMaskAndComputeOverlapRate( lambdaRange, itplMask );

        // Check for overlap rate
        if( fittingResults.overlapRate < overlapThreshold || fittingResults.overlapRate<=0.0 )
        {
            fittingResults.status = COperator::nStatus_NoOverlap;
            return;
        }
    }

    // Linear fit
    Int32 imin_lbda = spcSpectralAxis.GetIndexAtWaveLength(spcLambdaRange.GetBegin());
    if(spcSpectralAxis[imin_lbda]<spcLambdaRange.GetBegin() && imin_lbda+1<spcSpectralAxis.GetSamplesCount())
    {
        imin_lbda += 1;
    }
    Int32 imax_lbda = spcSpectralAxis.GetIndexAtWaveLength(spcLambdaRange.GetEnd());
    if(spcSpectralAxis[imax_lbda]>spcLambdaRange.GetEnd() && imax_lbda-1>=imin_lbda)
    {
        imax_lbda -= 1;
    }
    Int32 n=imax_lbda-imin_lbda+1;
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
            Float64 fval =  m_templatesRebined_bf[iddl]->GetFluxAxis()[i+imin_lbda];
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
            Float64 a = gsl_vector_get(c,iddl)*normFactor;
            if(verbose)
            {
                Log.LogInfo("# Found amplitude %d: %+.5e +- %.5e", iddl, a, COV(iddl,iddl)*normFactor);
            }
        }
    }

    //save the fitted amps and fitErrors, etc...
    for (Int32 iddl = 0; iddl < nddl; iddl++)
    {
        Float64 a = gsl_vector_get(c,iddl)*normFactor;
        fittingResults.fittingAmplitudes[iddl]=a;
        Float64 err = COV(iddl,iddl)*normFactor;
        fittingResults.fittingErrors[iddl] = err;
    }
    if(fittingResults.fittingAmplitudes.size()!=nddl)
    {
        Log.LogDebug("  Operator-Tplcombination: Found nfittedamps(=%d) different than nddl(=%d)", fittingResults.fittingAmplitudes.size(), nddl);
    }

    //build the model
    CSpectrum modelSpectrum;
    modelSpectrum.GetSpectralAxis().SetSize(n);
    modelSpectrum.GetFluxAxis().SetSize(n);
    for(Int32 k=0; k<modelSpectrum.GetSampleCount(); k++)
    {
        modelSpectrum.GetSpectralAxis()[k]=spcSpectralAxis[k+imin_lbda];
        modelSpectrum.GetFluxAxis()[k]=.0;
        for (Int32 iddl = 0; iddl < nddl; iddl++)
        {
            Float32 a = fittingResults.fittingAmplitudes[iddl];
            if(k==0){
                Log.LogDebug("  Operator-Tplcombination: Building model with tpl/component #%d with amplitude a=%+.5e", iddl, a);
            }
            modelSpectrum.GetFluxAxis()[k] += a*m_templatesRebined_bf[iddl]->GetFluxAxis()[k+imin_lbda];
        }
    }
    fittingResults.modelSpectrum = modelSpectrum;

    //estimate the lst-square brute force
    fittingResults.chisquare = .0;
    Float64 diff, err2;
    for(Int32 k=0; k<n; k++)
    {
        diff = modelSpectrum.GetFluxAxis()[k]-spcFluxAxis[k+imin_lbda];
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

/**
 * \brief
 *
 * input: if additional_spcMasks size is 0, no additional mask will be used, otherwise its size should match the redshifts list size
 **/
std::shared_ptr<COperatorResult> COperatorTplcombination::Compute(const CSpectrum& spectrum,
                                                                  std::vector<CTemplate> tplList,
                                                                  const TFloat64Range& lambdaRange,
                                                                  const TFloat64List& redshifts,
                                                                  Float64 overlapThreshold,
                                                                  std::vector<CMask> additional_spcMasks,
                                                                  std::string opt_interp,
                                                                  Int32 opt_extinction,
                                                                  Int32 opt_dustFitting)
{
    Log.LogInfo("  Operator-tplcombination: starting computation with N-template = %d", tplList.size());


    if( opt_dustFitting && m_ismCorrectionCalzetti->calzettiInitFailed)
    {
        Log.LogError("  Operator-tplcombination: no calzetti calib. file loaded... aborting!");
        return NULL;
    }
    if( opt_extinction && m_igmCorrectionMeiksin->meiksinInitFailed)
    {
        Log.LogError("  Operator-tplcombination: no meiksin calib. file loaded... aborting!");
        return NULL;
    }
    if( spectrum.GetSpectralAxis().IsInLinearScale() == false )
    {
        Log.LogError("  Operator-tplcombination: input spectrum is not in log scale (ignored)");
        //return NULL;
    }


    for(Int32 ktpl=0; ktpl<tplList.size(); ktpl++)
    {
        if( tplList[ktpl].GetSpectralAxis().IsInLinearScale() == false )
        {
            Log.LogError("  Operator-tplcombination: input template k=%d are not in log scale (ignored)", ktpl);
            //return NULL;
        }
    }

    Log.LogDebug("  Operator-tplcombination: allocating memory for buffers (N = %d)", tplList.size());
    Float64* precomputedFineGridTplFlux = NULL; //Todo: define as a list for each tpl
    for(Int32 ktpl=0; ktpl<tplList.size(); ktpl++)
    {
        // Pre-Allocate the rebined template with regard to the spectrum size
        std::shared_ptr<CTemplate> templateRebined_bf = std::shared_ptr<CTemplate>( new CTemplate( tplList[ktpl].GetName(), tplList[ktpl].GetCategory() ) );
        templateRebined_bf->GetSpectralAxis().SetSize(spectrum.GetSampleCount());
        templateRebined_bf->GetFluxAxis().SetSize(spectrum.GetSampleCount());
        m_templatesRebined_bf.push_back(templateRebined_bf);

        // Pre-Allocate the rebined mask with regard to the spectrum size
        std::shared_ptr<CMask> mskRebined_bf = std::shared_ptr<CMask>( new CMask() );
        mskRebined_bf->SetSize(spectrum.GetSampleCount());
        m_masksRebined_bf.push_back(mskRebined_bf);

        //
        std::shared_ptr<CSpectrumSpectralAxis> shiftedTplSpectralAxis_bf = std::shared_ptr<CSpectrumSpectralAxis>( new CSpectrumSpectralAxis(  ));
        shiftedTplSpectralAxis_bf->SetSize(tplList[ktpl].GetSampleCount());
        m_shiftedTemplatesSpectralAxis_bf.push_back(shiftedTplSpectralAxis_bf);

        if(opt_interp=="precomputedfinegrid"){
            /*/
            // Precalculate a fine grid template to be used for the 'closest value' rebin method
            Int32 n = tpl.GetSampleCount();
            CSpectrumFluxAxis tplFluxAxis = tpl.GetFluxAxis();
            CSpectrumSpectralAxis tplSpectralAxis = tpl.GetSpectralAxis();
            //Float64 dLambdaTgt =  1.0 * ( spectrum.GetMeanResolution()*0.9 )/( 1+sortedRedshifts[sortedRedshifts.size()-1] );
            Float64 dLambdaTgt =  0.1;
            //Float64 lmin = tplSpectralAxis[0];
            Float64 lmin = 0;
            Float64 lmax = tplSpectralAxis[n-1];
            Int32 nTgt = (lmax-lmin)/dLambdaTgt + 2.0/dLambdaTgt;

            // pfg with std::vector
            //CTemplate       templateFine;
            //templateFine.GetSpectralAxis().SetSize(nTgt);
            //templateFine.GetFluxAxis().SetSize(nTgt);
            //Float64* precomputedFineGridTplFlux = templateFine.GetFluxAxis().GetSamples();
            // pfg with malloc
            precomputedFineGridTplFlux = (Float64*)malloc(nTgt*sizeof(Float64));
            if(precomputedFineGridTplFlux == NULL)
            {
                Log.LogError("  Operator-Chisquare2: unable to allocate the precomputed fine grid buffer... aborting!");
                return NULL;
            }
            // pfg with static array => doesn't work
            //nTgt = 999999;
            //Float64 precomputedFineGridTplFlux[999999];
            //Log.LogInfo( "nTgt: %d samples", nTgt);

            //inialise and allocate the gsl objects
            Float64* Ysrc = tplFluxAxis.GetSamples();
            Float64* Xsrc = tplSpectralAxis.GetSamples();
            // linear
            //gsl_interp *interpolation = gsl_interp_alloc (gsl_interp_linear,n);
            //gsl_interp_init(interpolation, Xsrc, Ysrc, n);
            //gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();

            //spline
            gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n);
            gsl_spline_init (spline, Xsrc, Ysrc, n);
            gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();

            Int32 k = 0;
            Float64 x = 0.0;
            for(k=0; k<nTgt; k++){
                x = lmin + k*dLambdaTgt;
                if(x < tplSpectralAxis[0] || x > tplSpectralAxis[n-1]){
                    precomputedFineGridTplFlux[k] = 0.0;
                }else{
                    //precomputedFineGridTplFlux[k] = gsl_interp_eval(interpolation, Xsrc, Ysrc, x, accelerator);
                    precomputedFineGridTplFlux[k] = gsl_spline_eval (spline, x, accelerator);
                }
            }

            gsl_spline_free (spline);
            gsl_interp_accel_free (accelerator);
            //*/
        }
        /*//debug:
        // save templateFine
        FILE* f = fopen( "template_precomputedFineGrid.txt", "w+" );
        for(Int32 m=0; m<nTgt; m++){
            Float64 x = lmin + k*dLambdaTgt;
            if( Yfine[m] < 0.0001 ){
                fprintf( f, "%e %e\n", x, precomputedFineGridTplFlux[m]);
            }else{
                fprintf( f, "%f %f\n", x, precomputedFineGridTplFlux[m]);
            }
        }
        fclose( f );
        //*/
    }

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
    std::shared_ptr<CTplcombinationResult> result = std::shared_ptr<CTplcombinationResult>( new CTplcombinationResult() );
    Int32 nDustCoeffs=1;
    if(opt_dustFitting)
    {
        nDustCoeffs = m_ismCorrectionCalzetti->GetNPrecomputedDustCoeffs();
    }
    Log.LogDebug("  Operator-tplcombination: prepare N ism coeffs = %d", nDustCoeffs);
    Int32 nIGMCoeffs=1;
    if(opt_extinction)
    {
        nIGMCoeffs = m_igmCorrectionMeiksin->GetIdxCount();
    }
    Log.LogDebug("  Operator-tplcombination: prepare N igm coeffs = %d", nIGMCoeffs);

    result->Init(sortedRedshifts.size(), nDustCoeffs, nIGMCoeffs);
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
    }


    boost::progress_display show_progress(Float64(sortedRedshifts.size())) ;
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
                  precomputedFineGridTplFlux,
                  lambdaRange,
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
            break;
        }

        result->ChiSquare[i]=fittingResults.chisquare;
        result->Overlap[i]=fittingResults.overlapRate;
        result->ChiSquareIntermediate[i]=fittingResults.ChiSquareInterm;

        ++show_progress;
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
    result->CstLog = EstimateLikelihoodCstLog(spectrum, lambdaRange);

    // extrema
    Int32 extremumCount = 10;
    if(result->Redshifts.size()>extremumCount)
    {
        TPointList extremumList;
        TFloat64Range redshiftsRange(result->Redshifts[0], result->Redshifts[result->Redshifts.size()-1]);
        CExtremum extremum( redshiftsRange, extremumCount, true);
        extremum.Find( result->Redshifts, result->ChiSquare, extremumList );
        // Refine Extremum with a second maximum search around the z candidates:
        // This corresponds to the finer xcorrelation in EZ Pandora (in standard_DP fctn in SolveKernel.py)
        Float64 radius = 0.001;
        for( Int32 i=0; i<extremumList.size(); i++ )
        {
            Float64 x = extremumList[i].X;
            Float64 left_border = max(redshiftsRange.GetBegin(), x-radius);
            Float64 right_border=min(redshiftsRange.GetEnd(), x+radius);

            TPointList extremumListFine;
            TFloat64Range rangeFine = TFloat64Range( left_border, right_border );
            CExtremum extremumFine( rangeFine , 1, true);
            extremumFine.Find( result->Redshifts, result->ChiSquare, extremumListFine );
            if(extremumListFine.size()>0){
                extremumList[i] = extremumListFine[0];
            }
        }
        // store extrema results
        result->Extrema.resize( extremumCount );
        for( Int32 i=0; i<extremumList.size(); i++ )
        {
            result->Extrema[i] = extremumList[i].X;
        }

        Log.LogDebug("  Operator-Tplcombination: EXTREMA found n=%d", extremumList.size());
    }else
    {
        // store extrema results
        result->Extrema.resize( result->Redshifts.size() );
        TFloat64List tmpX;
        TFloat64List tmpY;
        for( Int32 i=0; i<result->Redshifts.size(); i++ )
        {
            tmpX.push_back(result->Redshifts[i]);
            tmpY.push_back(result->ChiSquare[i]);
        }
        // sort the results by merit
        CQuickSort<Float64> sort;
        vector<Int32> sortedIndexes( result->Redshifts.size() );
        sort.SortIndexes( tmpY.data(), sortedIndexes.data(), sortedIndexes.size() );
        for( Int32 i=0; i<result->Redshifts.size(); i++ )
        {
            result->Extrema[i] = tmpX[sortedIndexes[i]];
        }
        Log.LogDebug("  Operator-Tplcombination: EXTREMA forced n=%d", result->Extrema.size());
    }

    //store spectrum results
    Int32 nMaxExtremaSpectraSave = 1;
    for( Int32 ie=0; ie<result->Extrema.size(); ie++ )
    {
        if(ie>=nMaxExtremaSpectraSave)
        {
            break;
        }

        Int32 i=-1;
        for( Int32 iz=0; iz<result->Redshifts.size(); iz++ )
        {
            if(result->Redshifts[iz]==result->Extrema[ie])
            {
                i=iz;
                break;
            }
        }

        if(i!=-1)
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
                      precomputedFineGridTplFlux,
                      lambdaRange,
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
                break;
            }

            Log.LogDetail("  Operator-Tplcombination: creating spectrum/model tplcombination result for z=%f", redshift);
            for(Int32 itpl=0; itpl<fittingResults.fittingAmplitudes.size(); itpl++)
            {
                Log.LogDetail("  Operator-Tplcombination: creating spectrum/model with tpl-%d amp = %.4e", itpl, fittingResults.fittingAmplitudes[itpl]);
            }
            std::shared_ptr<CModelSpectrumResult>  resultspcmodel = std::shared_ptr<CModelSpectrumResult>( new CModelSpectrumResult(fittingResults.modelSpectrum) );
            m_savedModelSpectrumResults.push_back(resultspcmodel);

        }else{
            Log.LogError("  Operator-Tplcombination: Unable to find the z index for spectrum result save");
        }
    }

    if(opt_interp=="precomputedfinegrid"){
        free(precomputedFineGridTplFlux);
    }
    return result;

}


/* @brief COperatorTplcombination::getDustCoeff: get the dust coeff at a fixed resolution of 1A
* @param dustCoeff
* @param maxLambda
* @return
*/
const Float64*  COperatorTplcombination::getDustCoeff(Float64 dustCoeff, Float64 maxLambda)
{
    //find kDust
    Int32 idxDust = -1;
    for(Int32 kDust=0; kDust<m_ismCorrectionCalzetti->GetNPrecomputedDustCoeffs(); kDust++)
    {
        Float64 coeffEBMV = m_ismCorrectionCalzetti->GetEbmvValue(kDust);
        if(dustCoeff==coeffEBMV)
        {
            idxDust = kDust;
            break;
        }
    }

    if(idxDust<0)
    {
        return 0;
    }

    Int32 nSamples = maxLambda+1; //+1 for security
    Float64* dustCoeffs = new Float64 [(int)nSamples]();


    for(Int32 kl=0; kl<nSamples; kl++)
    {
        Float64 restLambda = kl;
        Float64 coeffDust = m_ismCorrectionCalzetti->getDustCoeff( idxDust, restLambda);
        dustCoeffs[kl] = coeffDust;
    }
    return dustCoeffs;
}

/**
 * @brief COperatorTplcombination::getMeiksinCoeff: get the IGM Meiksin coeff at a fixed resolution of 1A
 * @param dustCoeff
 * @param maxLambda
 * @return
 */
const Float64*  COperatorTplcombination::getMeiksinCoeff(Int32 meiksinIdx, Float64 redshift, Float64 maxLambda)
{
    if(meiksinIdx<0 || meiksinIdx>m_igmCorrectionMeiksin->GetIdxCount()-1)
    {
        return 0;
    }

    //find redshiftIdx from redshift value
   Int32 redshiftIdx = m_igmCorrectionMeiksin->GetRedshiftIndex(redshift);

    Int32 nSamples = maxLambda+1; //+1 for security
    Float64* meiksinCoeffs = new Float64 [(int)nSamples]();


    for(Int32 kl=0; kl<nSamples; kl++)
    {
        Float64 restLambda = kl;
        Float64 coeffIGM = 1.0;
        if(restLambda <= m_igmCorrectionMeiksin->GetLambdaMax())
        {
            Int32 kLbdaMeiksin = 0;
            if(restLambda >= m_igmCorrectionMeiksin->GetLambdaMin())
            {
                kLbdaMeiksin = Int32(restLambda-m_igmCorrectionMeiksin->GetLambdaMin());
            }else //if lambda lower than min meiksin value, use lower meiksin value
            {
                kLbdaMeiksin = 0;
            }

            coeffIGM = m_igmCorrectionMeiksin->m_corrections[redshiftIdx].fluxcorr[meiksinIdx][kLbdaMeiksin];

        }
        meiksinCoeffs[kl] = coeffIGM;
    }
    return meiksinCoeffs;
}

/**
 * \brief this function estimates the likelihood_cstLog term withing the wavelength range
 **/
Float64 COperatorTplcombination::EstimateLikelihoodCstLog(const CSpectrum& spectrum, const TFloat64Range& lambdaRange)
{
    const CSpectrumSpectralAxis& spcSpectralAxis = spectrum.GetSpectralAxis();
    const TFloat64List& error = spectrum.GetFluxAxis().GetError();

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

void COperatorTplcombination::SaveSpectrumResults(CDataStore &dataStore)
{
    Log.LogDetail("  Operator-Tplcombination: now saving spectrum/model tplcombination results n=%d", m_savedModelSpectrumResults.size());
    for(Int32 ie=0; ie<m_savedModelSpectrumResults.size(); ie++)
    {
        std::string fname_spc = (boost::format("tplcombinationmodel_spc_extrema_%1%") % ie).str();
        dataStore.StoreScopedGlobalResult( fname_spc.c_str(), m_savedModelSpectrumResults[ie] );
    }
}
