#include  "RedshiftLibrary/operator/tplcombination.h"
#include  "RedshiftLibrary/operator/tplcombinationresult.h"
#include  "RedshiftLibrary/spectrum/axis.h"
#include  "RedshiftLibrary/spectrum/spectrum.h"
#include  "RedshiftLibrary/spectrum/template/template.h"
#include  "RedshiftLibrary/spectrum/tools.h"
#include  "RedshiftLibrary/common/mask.h"
#include  "RedshiftLibrary/extremum/extremum.h"
#include  "RedshiftLibrary/common/quicksort.h"
#include  "RedshiftLibrary/common/indexing.h"
#include  "RedshiftLibrary/processflow/datastore.h"
#include  "RedshiftLibrary/log/log.h"

#include <boost/numeric/conversion/bounds.hpp>
#include <math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <algorithm>    // std::sort
#include <float.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit_nlin.h>
#include <numeric>
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


void COperatorTplcombination::BasicFit_preallocateBuffers(const CSpectrum& spectrum, const TTemplateConstRefList& tplList)
{
    Int32 componentCount = tplList.size();
    // Pre-Allocate the rebined template and mask with regard to the spectrum size
    m_templatesRebined_bf.resize(componentCount);
    m_masksRebined_bf.resize(componentCount);
    m_spcSpectralAxis_restframe.SetSize(spectrum.GetSampleCount());

    for(Int32 ktpl=0; ktpl<componentCount; ktpl++)
    {
        m_templatesRebined_bf[ktpl].m_ismCorrectionCalzetti = tplList[ktpl]->m_ismCorrectionCalzetti;
        m_templatesRebined_bf[ktpl].m_igmCorrectionMeiksin = tplList[ktpl]->m_igmCorrectionMeiksin;
    }
}

// @BasicFit@ in its current content and the commented code corresponds to @::computeModelSpectrum@ since there is no
// loops over igm/ism corrections

void COperatorTplcombination::BasicFit(const CSpectrum& spectrum,
                                       const TTemplateConstRefList& tplList,
                                       const TFloat64Range& lambdaRange,
                                       Float64 redshift,
                                       Float64 overlapThreshold,
                                       STplcombination_basicfitresult& fittingResults,
                                       std::string opt_interp,
                                       Float64 forcedAmplitude,
                                       Int32 opt_extinction,
                                       Int32 opt_dustFitting,
                                       CMask spcMaskAdditional,
                                       CPriorHelper::TPriorEList logpriore,
                                       bool keepigmism,
                                       const TInt32List& MeiksinList)
{
    bool verbose = false;
    if(verbose)
    {
        Log.LogDebug("  Operator-tplcombination: BasicFit - for z=%f", redshift);
    }
    boost::chrono::thread_clock::time_point start_prep = boost::chrono::thread_clock::now();

    bool status_chisquareSetAtLeastOnce = false;

    const CSpectrumSpectralAxis& spcSpectralAxis = spectrum.GetSpectralAxis();
    const CSpectrumFluxAxis& spcFluxAxis = spectrum.GetFluxAxis();
    const CSpectrumNoiseAxis& spcError = spcFluxAxis.GetError();

    if(spcMaskAdditional.GetMasksCount()!=spcFluxAxis.GetSamplesCount())
    {
        Log.LogError("  Operator-Tplcombination: spcMaskAdditional does not have the same size as the spectrum flux vector... (%d vs %d), aborting", spcMaskAdditional.GetMasksCount(), spcFluxAxis.GetSamplesCount());
        throw std::runtime_error("  Operator-Tplcombination: spcMaskAdditional does not have the same size as the spectrum flux vector. aborting");       
    }

    TFloat64Range currentRange;   
    Int32 ret = RebinTemplate(spectrum, 
                              tplList, 
                              redshift, 
                              lambdaRange, 
                              opt_interp,
                              currentRange,
                              fittingResults.overlapRate,
                              overlapThreshold); 
    
    Int32 kStart = -1, kEnd = -1, kIgmEnd = -1;
    //I consider here that all templates share the same spectralAxis
    currentRange.getClosedIntervalIndices(m_templatesRebined_bf[0].GetSpectralAxis().GetSamplesVector(), kStart, kEnd);
    Int32 kStart_model = kStart; //mainly used at high redshifts, when desextincting spectrum is happenning with null coeffs
    Int32 nddl = tplList.size();
    if (opt_extinction || opt_dustFitting){
        for (Int32 iddl = 0; iddl <nddl; iddl++)
            m_templatesRebined_bf[iddl].InitIsmIgmConfig(kStart, kEnd, redshift);
    }
    if(opt_extinction)
        kIgmEnd = m_templatesRebined_bf[0].GetIgmEndIndex();


    if( ret == -1 ){
        fittingResults.status = COperator::nStatus_NoOverlap; 
        return;
    }
    /*if( ret == -2 ){//this is not coded
        fittingResults.status = COperator::nStatus_DataError;
        return;
    }*/

    //determine min and max value of ebmv coeff
    Int32 nISM = fittingResults.IsmCalzettiCoeffInterm.size(); 
    Int32 nIGM = fittingResults.IgmMeiksinIdxInterm[0].size(); 

    Int32 iEbmvCoeffMin = 0;
    Int32 iEbmvCoeffMax = iEbmvCoeffMin+nISM;
    if(keepigmism){
        iEbmvCoeffMin = m_templatesRebined_bf[0].m_ismCorrectionCalzetti->GetEbmvIndex(fittingResults.EbmvCoeff);
        iEbmvCoeffMax = iEbmvCoeffMin;
    }
    // Linear fit
    Int32 n = kEnd-kStart +1;
    Log.LogDebug("  Operator-Tplcombination: prep. linear fitting with n=%d samples in the clamped lambdarange spectrum (imin=%d, lbda_min=%.3f - imax=%d, lbda_max=%.3f)", n, kStart, spcSpectralAxis[kStart], kEnd, spcSpectralAxis[kEnd]);
    
    gsl_matrix *X, *cov;
    gsl_vector *y, *w, *c;

    X = gsl_matrix_alloc (n, nddl);
    y = gsl_vector_alloc (n);
    w = gsl_vector_alloc (n);
    c = gsl_vector_alloc (nddl);
    cov = gsl_matrix_alloc (nddl, nddl);

    // Normalizing factor
    Float64 normFactor = GetNormFactor(spcFluxAxis, kStart, n);

    if(verbose)
    {
        Log.LogDetail("  Operator-Tplcombination: Linear fitting, found normalization Factor=%e", normFactor);
    }

    Bool option_igmFastProcessing = (MeiksinList.size()==1 ? false : true); //TODO
    Bool igmLoopUseless_WavelengthRange = false;
    fittingResults.chisquare = boost::numeric::bounds<float>::highest();//final best Xi2 value
    Float64 chisq;
    Float64 dtd_complement = 0.;//value mainly relevant with DisextinctData method
    /**
     * there are two ways to apply extinction:
     * 1. apply the inverse of extinction on both spcFlux and error on flux
     * 2. apply extinction on all templates, using original spcFlux an error
     * 
     * TODO: keep the first method once we completely validate the equivalence of both method
    */
    Bool DisextinctData = false; // true = option 1; false = option 2
    //create a template with cte flux = 1, to be used only when disextincting data and noise
    CTemplate identityTemplate("identity", "idle", m_templatesRebined_bf[0].GetSpectralAxis(), std::move(CSpectrumFluxAxis(m_templatesRebined_bf[0].GetSampleCount(), 1)));
    // Prepare the fit data, once for all
    if(!DisextinctData){
        Float64 yi, ei;
        for(Int32 i = 0; i < n; i++)
        {
            yi = spcFluxAxis[i+kStart]/normFactor;
            ei = spcError[i+kStart]/normFactor;

            gsl_vector_set (y, i, yi); //y[i] = yi
            gsl_vector_set (w, i, 1.0/(ei*ei));//w[i] = 1/(ei*ei)
        }
    } 
    for(Int32 kigm = 0; kigm<nIGM; kigm++)
    {
        if(igmLoopUseless_WavelengthRange)
        {
            //Now copy from the already calculated k>0 igm values
            for(Int32 kism=0; kism<fittingResults.ChiSquareInterm.size(); kism++)
            {
                for(Int32 kigm=1; kigm<fittingResults.ChiSquareInterm[kism].size(); kigm++)
                {
                    fittingResults.ChiSquareInterm[kism][kigm] = fittingResults.ChiSquareInterm[kism][0];
                    fittingResults.IsmCalzettiCoeffInterm[kism][kigm] = fittingResults.IsmCalzettiCoeffInterm[kism][0];
                    fittingResults.IgmMeiksinIdxInterm[kism][kigm] = fittingResults.IgmMeiksinIdxInterm[kism][0];
                }
            }
            break;
        }
        Int32 meiksinIdx = MeiksinList[kigm];

        bool igmCorrectionAppliedOnce = false;
        //applyMeiksin on all templates
        if(opt_extinction == 1 && !DisextinctData){
            igmCorrectionAppliedOnce = true;//since we are multiplying the bools, rather init to true
            for (Int32 iddl = 0; iddl <nddl; iddl++)
            {
                igmCorrectionAppliedOnce &= m_templatesRebined_bf[iddl].ApplyMeiksinCoeff(meiksinIdx);
            } 
            if(!igmCorrectionAppliedOnce)
                igmLoopUseless_WavelengthRange = true; 
        }

        for(Int32 kEbmv = iEbmvCoeffMin; kEbmv<nISM; kEbmv++)
        {
            Int32 kEbmv_ = kEbmv - iEbmvCoeffMin; //index used to fill some arrays
            Float64 coeffEBMV = -1.; //no ism by default (ie DustCoeff=1.)
            //apply ism on all templates, once for all
            if (opt_dustFitting == 1 && !DisextinctData){
                coeffEBMV = m_templatesRebined_bf[0].m_ismCorrectionCalzetti->GetEbmvValue(kEbmv);
                for (Int32 iddl = 0; iddl <nddl; iddl++)
                    m_templatesRebined_bf[iddl].ApplyDustCoeff(kEbmv);
            } 

            if((opt_extinction || opt_dustFitting) && DisextinctData)
            {
                identityTemplate.InitIsmIgmConfig(kStart, kEnd, redshift, tplList[0]->m_ismCorrectionCalzetti, tplList[0]->m_igmCorrectionMeiksin);
                if(opt_extinction){ 
                    igmCorrectionAppliedOnce = identityTemplate.ApplyMeiksinCoeff(meiksinIdx);
                    if(!igmCorrectionAppliedOnce)
                        igmLoopUseless_WavelengthRange = true;
                }
                if(opt_dustFitting){
                    coeffEBMV = m_templatesRebined_bf[0].m_ismCorrectionCalzetti->GetEbmvValue(kEbmv);
                    identityTemplate.ApplyDustCoeff(kEbmv);
                }
                const CSpectrumFluxAxis& extinction = identityTemplate.GetFluxAxis();
                Float64 yi, ei;
                //Important: hereinafter we assume that kStart_up can only be equal or higher than kStart when moving throw the igm curves
                //this is valid only iif igm curves are loaded in the increasing ordre,i.e., from the least extinction curve to the highest extinction curve
                UInt32 kStart_up = kStart;
                for(Int32 i = 0; i < n; i++)
                {
                    if(!extinction[i]){
                        kStart_up++; //find the first non-null value and update kStart accordingly
                        continue;
                    }
                    if(kStart_up>kStart)
                    {
                        //compute the dtd part than wont be computed using gsl
                        dtd_complement += ComputeDtD(spcFluxAxis, TInt32Range(kStart_model, kStart_up - 1));
                        n = kEnd - kStart_up + 1; //update the length value
                        kStart = kStart_up;
                        //free previously allocated buffers
                        gsl_matrix_free (X);
                        gsl_vector_free (y);
                        gsl_vector_free (w);
                        //reallocate buffers with the newly calculated value
                        X = gsl_matrix_alloc (n, nddl);
                        y = gsl_vector_alloc (n);
                        w = gsl_vector_alloc (n);
                    }
                    yi = spcFluxAxis[i+kStart]/extinction[i]/normFactor;
                    ei = spcError[i+kStart]/extinction[i]/normFactor;

                    gsl_vector_set (y, i, yi); //y[i] = yi
                    gsl_vector_set (w, i, 1.0/(ei*ei));//w[i] = 1/(ei*ei)
                }
            }     

            //preparing the template matrix for computing the Xi2 value
            //the fastigm has its effect mainly on the number of spectrum to consider for computing the amplitudes 
            //and then the fit
            for(Int32 i = 0; i < n; i++)
            {
                for (Int32 iddl = 0; iddl < nddl; iddl++)
                {
                    Float64 fval =  m_templatesRebined_bf[iddl].GetFluxAxis()[i+kStart];
                    gsl_matrix_set (X, i, iddl, fval);//i.e., X[i,iddl]=fval -> X is an extract of m_templatesRebinned_bf 
                }
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
            chisq+= dtd_complement;//complement with the dtd value, only in the case of DisextinctData.  
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
            Float64 a, err2, err, sA = 0, sE = 0;
            for (Int32 iddl = 0; iddl < nddl; iddl++)
            {
                a = gsl_vector_get(c,iddl);
                fittingResults.fittingAmplitudes[iddl] = a*normFactor;
                err2 = COV(iddl,iddl);
                err = std::sqrt(err2);
                fittingResults.fittingAmplitudeErrors[iddl] = err*normFactor;
                fittingResults.fittingAmplitudeSigmas[iddl] = a/err;
                sA += a*a;
                sE += err2;
            }
            fittingResults.SNR = std::sqrt(sA/sE);
            fittingResults.fittingAmplitudesInterm[kEbmv_][kigm] = fittingResults.fittingAmplitudes;//saving 

            //save covariance matrix into MtM
            for (Int32 iddl = 0; iddl < nddl; iddl++)
            {
               for (Int32 iddc = 0; iddc < nddl; iddc++)
               {
                   fittingResults.COV[iddl][iddc]=COV(iddl,iddc);
               } 
            }
            if(fittingResults.fittingAmplitudes.size()!=nddl)
            {
                Log.LogDebug("  Operator-Tplcombination: Found nfittedamps(=%d) different than nddl(=%d)", fittingResults.fittingAmplitudes.size(), nddl);
            }
            UInt32 modelSize = kEnd - kStart_model + 1;
            //build the combined template model: spcmodel should have the same size as input spc
            TFloat64List spc_extract = TFloat64List(spcSpectralAxis.GetSamplesVector().begin() + kStart_model, 
                                                    spcSpectralAxis.GetSamplesVector().begin() + kEnd + 1);
            TFloat64List modelFlux(modelSize, 0.0);
            for (Int32 iddl = 0; iddl < nddl; iddl++){
                const CSpectrumFluxAxis & tmp = m_templatesRebined_bf[iddl].GetFluxAxis(); 
                const CSpectrumFluxAxis & ext = identityTemplate.GetFluxAxis(); //works for both cases
                
                Float32 a = fittingResults.fittingAmplitudes[iddl];
                for(Int32 k = 0; k<n; k++){
                    modelFlux[k]+= a*tmp[k+kStart]*ext[k+kStart];
                }
            }
            //TODO: optimize the below
            CSpectrumSpectralAxis xAxis(std::move(spc_extract), n);
            CSpectrumFluxAxis yAxis(std::move(modelFlux));
            CTemplate modelSpec("combination", "", std::move(xAxis), std::move(yAxis));//m_templatesRebined_bf cant be used, thus we create a new CTemplate to apply corrections
            fittingResults.modelSpectrum = modelSpec;

            if(chisq < fittingResults.chisquare)
            {
                fittingResults.chisquare = chisq;
                fittingResults.IGMIdx = igmCorrectionAppliedOnce?meiksinIdx:-1;
                fittingResults.EbmvCoeff = coeffEBMV;
                status_chisquareSetAtLeastOnce = true;
            }

            //save the interm chisquares in the intermediate vector
            fittingResults.ChiSquareInterm[kEbmv_][kigm] = fittingResults.chisquare;
            fittingResults.IsmCalzettiCoeffInterm[kEbmv_][kigm] = coeffEBMV;
            fittingResults.IgmMeiksinIdxInterm[kEbmv_][kigm] = igmCorrectionAppliedOnce?meiksinIdx:-1;

            boost::chrono::thread_clock::time_point stop_postprocess = boost::chrono::thread_clock::now();
            Float64 duration_postprocess = boost::chrono::duration_cast<boost::chrono::microseconds>(stop_postprocess - start_postprocess).count();
            Log.LogDebug( "  Operator-Tplcombination: Linear fitting, postprocess = %.3f microsec", duration_postprocess);

        }//end iterating over ISM
    }//end iterating over IGM

    //fittingResults.modelSpectrum = CSpectrum(CSpectrumSpectralAxis(std::move(spc_extract)), CSpectrumFluxAxis(std::move(modelFlux)));

    fittingResults.SNR = -1.0;

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

//estimate the lst-square brute force
Float64 COperatorTplcombination::ComputeXi2_bruteForce(const CSpectrumFluxAxis& correctedFlux, 
                                                       const CSpectrumFluxAxis& spcFluxAxis,
                                                       const Int32 kStart)
{
    const CSpectrumNoiseAxis& spcError = spcFluxAxis.GetError();

    Float64 diff, err2;
    Float64 chi2Value = .0;
    for(Int32 k=0; k<correctedFlux.GetSamplesCount(); k++)
    {
        diff = correctedFlux[k]-spcFluxAxis[k+kStart];//indeces should be verified
        err2 = spcError[k+kStart]*spcError[k+kStart];
        chi2Value += diff*diff/err2;
    }
    return chi2Value;
}

Int32  COperatorTplcombination::RebinTemplate( const CSpectrum& spectrum,
                                                const TTemplateConstRefList& tplList,
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
        const CSpectrumSpectralAxis& tplSpectralAxis = tplList[ktpl]->GetSpectralAxis();
        const CSpectrumFluxAxis& tplFluxAxis = tplList[ktpl]->GetFluxAxis();

        // Compute clamped lambda range over template
        TFloat64Range tplLambdaRange;
        tplSpectralAxis.ClampLambdaRange( lambdaRange_restframe, tplLambdaRange );

        // if there is any intersection between the lambda range of the spectrum and the lambda range of the template
        // Compute the intersected range
        TFloat64Range intersectedLambdaRange( 0.0, 0.0 );
        TFloat64Range::Intersect( tplLambdaRange, spcLambdaRange_restframe, intersectedLambdaRange );

        // find lambda range intersection common to all templates
        intersectedAllLambdaRange.IntersectWith(intersectedLambdaRange);

        CTemplate & itplTplSpectrum = m_templatesRebined_bf[ktpl];
        CMask & itplMask = m_masksRebined_bf[ktpl];

        tplList[ktpl]->Rebin( intersectedLambdaRange, m_spcSpectralAxis_restframe, itplTplSpectrum, itplMask, opt_interp );
        
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
    return 0;
}

/**
 * \brief
 *
 * input: if additional_spcMasks size is 0, no additional mask will be used, otherwise its size should match the redshifts list size
 * @lambdaRange is not clamped
 **/
std::shared_ptr<COperatorResult> COperatorTplcombination::Compute(const CSpectrum& spectrum,
                                                                  const TTemplateConstRefList& tplList,
                                                                  const TFloat64Range& lambdaRange,
                                                                  const TFloat64List& redshifts,
                                                                  Float64 overlapThreshold,
                                                                  std::vector<CMask> additional_spcMasks,
                                                                  std::string opt_interp,
                                                                  Int32 opt_extinction,
                                                                  Int32 opt_dustFitting,
                                                                  CPriorHelper::TPriorZEList logpriorze,
                                                                  Bool keepigmism,
                                                                  Float64 FitEbmvCoeff,
                                                                  Float64 FitMeiksinIdx)
{
    Int32 componentCount = tplList.size();
    Log.LogInfo("  Operator-tplcombination: starting computation with N-template = %d", componentCount);

    if( spectrum.GetSpectralAxis().IsInLinearScale() == false )
    {
        Log.LogError("  Operator-tplcombination: input spectrum is not in log scale");
        throw std::runtime_error("  Operator-tplcombination: input spectrum is not in log scale");
    }

    for(Int32 ktpl=0; ktpl<componentCount; ktpl++)
    {
        if( tplList[ktpl]->GetSpectralAxis().IsInLinearScale() == false )
        {
            Log.LogError("  Operator-tplcombination: input template k=%d are not in log scale", ktpl);
            throw std::runtime_error("  Operator-tplcombination: input template k=%d are not in log scale");
        }
        if( opt_dustFitting && tplList[ktpl]->m_ismCorrectionCalzetti->calzettiInitFailed)
        {
            Log.LogError("  Operator-tplcombination: no calzetti calib. file loaded... aborting");
            throw std::runtime_error("  Operator-tplcombination: no calzetti calib. file loaded... aborting");
        }
        if( opt_extinction && tplList[ktpl]->m_igmCorrectionMeiksin->meiksinInitFailed)
        {
            Log.LogError("  Operator-tplcombination: no meiksin calib. file loaded... aborting");
            throw std::runtime_error("  Operator-tplcombination: no meiksin calib. file loaded... aborting");
        }
    }

    Log.LogDebug("  Operator-tplcombination: allocating memory for buffers (N = %d)", componentCount);

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
    std::shared_ptr<CTplCombinationResult> result = std::shared_ptr<CTplCombinationResult>( new CTplCombinationResult() );
    Int32 nEbmvCoeffs = 1;
    if(opt_dustFitting && !keepigmism)
    {
        nEbmvCoeffs = tplList.front()->m_ismCorrectionCalzetti->GetNPrecomputedEbmvCoeffs();
    }
    Log.LogDebug("  Operator-tplcombination: prepare N ism coeffs = %d", nEbmvCoeffs);
    
    Int32 nIGMCoeffs = 1;
    if(opt_extinction && !keepigmism)
    {
        nIGMCoeffs = tplList.front()->m_igmCorrectionMeiksin->GetIdxCount();
    }
    Log.LogDebug("  Operator-tplcombination: prepare N igm coeffs = %d", nIGMCoeffs);
    
    //create meikinList
    TInt32List MeiksinList(nIGMCoeffs);
    if(opt_extinction)
    {
        if(keepigmism)
            MeiksinList[0] = FitMeiksinIdx;
        else
            std::iota(MeiksinList.begin(), MeiksinList.end(), 0);
    }else{//at least have one element
        MeiksinList[0] = -1;
    }

    result->Init(sortedRedshifts.size(), nEbmvCoeffs, nIGMCoeffs, componentCount);
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

    TFloat64List _chi2List(nIGMCoeffs, DBL_MAX);
    TFloat64List _ismList(nIGMCoeffs, NAN);
    TInt32List   _igmList(nIGMCoeffs, -1);
    TFloat64List _ampList(componentCount, NAN);
    TFloat64List _covList(componentCount, NAN);
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
        CPriorHelper::TPriorEList logp;
        if(logpriorze.size()>0 && logpriorze.size()==sortedRedshifts.size())
        {
            logp = logpriorze[i];
        }
        Float64 redshift = result->Redshifts[i];

        //initializing fittingResults: could be moved to another function
        STplcombination_basicfitresult fittingResults;
        if(keepigmism && opt_dustFitting && opt_extinction)
        {
            fittingResults.IGMIdx = FitMeiksinIdx;
            fittingResults.EbmvCoeff = FitEbmvCoeff;
        }
        //init fittingResult intermediate values before passing to ::BasicFit
        fittingResults.fittingAmplitudesInterm.resize(nEbmvCoeffs);
        for(Int32 kism=0; kism<nEbmvCoeffs; kism++)
        {   
            fittingResults.ChiSquareInterm.push_back(_chi2List);
            fittingResults.IsmCalzettiCoeffInterm.push_back(_ismList);
            fittingResults.IgmMeiksinIdxInterm.push_back(_igmList);

            fittingResults.fittingAmplitudesInterm[kism].resize(nIGMCoeffs);
            for(Int32 kigm = 0; kigm<nIGMCoeffs; kigm++)
                fittingResults.fittingAmplitudesInterm[kism][kigm]=_ampList;
        }
        fittingResults.chisquare = boost::numeric::bounds<float>::highest();
        fittingResults.fittingAmplitudes = TFloat64List(componentCount, NAN);
        fittingResults.fittingAmplitudeErrors = TFloat64List(componentCount, NAN);
        fittingResults.fittingAmplitudeSigmas = TFloat64List(componentCount, NAN);
        fittingResults.overlapRate = 0.0;
        fittingResults.status = COperator::nStatus_DataError;
        fittingResults.COV.resize(componentCount);//MtM = COV-1
        for(Int32 iddl = 0; iddl<componentCount; iddl++){
            fittingResults.COV[iddl] = _covList;
        }

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
                  additional_spcMask,
                  logp,
                  keepigmism,
                  MeiksinList);

        if(result->Status[i]==COperator::nStatus_InvalidProductsError)
        {
            Log.LogError("  Operator-Tplcombination: found invalid tplcombination products for z=%f. Now breaking z loop.", redshift);
            throw std::runtime_error("  Operator-Tplcombination: found invalid tplcombination products for z");      
        }

        result->ChiSquare[i] = fittingResults.chisquare;
        result->Overlap[i] = fittingResults.overlapRate;
        result->ChiSquareIntermediate[i] = fittingResults.ChiSquareInterm;
        result->FitAmplitude[i] = fittingResults.fittingAmplitudes;
        result->FitAmplitudeSigma[i] =fittingResults.fittingAmplitudeSigmas;
        result->FitAmplitudeError[i] = fittingResults.fittingAmplitudeErrors;
        result->SNR[i] = fittingResults.SNR;
        result->FitCOV[i] = fittingResults.COV;
        //result->LogPrior[i]=NAN: //not yet calculated
        result->FitEbmvCoeff[i] = fittingResults.EbmvCoeff; 
        result->FitMeiksinIdx[i] = fittingResults.IGMIdx;
        result->Status[i] = fittingResults.status;
        result->ChiSquareIntermediate[i] = fittingResults.ChiSquareInterm;
        result->IsmEbmvCoeffIntermediate[i] = fittingResults.IsmCalzettiCoeffInterm;
        result->IgmMeiksinIdxIntermediate[i] = fittingResults.IgmMeiksinIdxInterm;

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

    // Deallocate the rebined template and mask buffers
    m_templatesRebined_bf.clear();
    m_masksRebined_bf.clear();

    return result;

}

Int32   COperatorTplcombination::ComputeSpectrumModel( const CSpectrum& spectrum,
                                                        const TTemplateConstRefList& tplList,
                                                        Float64 redshift,
                                                        Float64 EbmvCoeff,
                                                        Int32 meiksinIdx,
                                                        const TFloat64List& amplitudes,
                                                        std::string opt_interp,
                                                        const TFloat64Range& lambdaRange,
                                                        Float64 overlapThreshold,
                                                        std::shared_ptr<CModelSpectrumResult> & spcPtr)
{
    Log.LogDetail("  Operator-COperatorTplCombination: building spectrum model tptCombination for candidate Zcand=%f", redshift);
    
    BasicFit_preallocateBuffers(spectrum, tplList);
    Int32 nddl = tplList.size();

    //Estatus status; 
    Float64 overlapRate = 0.0;
    TFloat64Range currentRange;   
    Int32 ret = RebinTemplate(spectrum, 
                              tplList, 
                              redshift, 
                              lambdaRange, 
                              opt_interp,
                              currentRange,
                              overlapRate,
                              overlapThreshold); 
    if( ret == -1 ){
        //status = nStatus_NoOverlap; 
        return -1;
    }
    if( ret == -2 ){
        //status = nStatus_DataError;
        return -1;
    }
    Int32 kStart = -1, kEnd = -1, kIgmEnd = -1;

    currentRange.getClosedIntervalIndices(m_templatesRebined_bf[0].GetSpectralAxis().GetSamplesVector(), kStart, kEnd);

    //create identityTemplate on which we apply meiksin and ism, once for all tpllist
    CTemplate identityTemplate("identity", "idle", m_templatesRebined_bf[0].GetSpectralAxis(), CSpectrumFluxAxis(m_templatesRebined_bf[0].GetSampleCount(), 1));
    
    if ((EbmvCoeff>0.) || (meiksinIdx>-1))
    {
        identityTemplate.InitIsmIgmConfig(kStart, kEnd, redshift, tplList[0]->m_ismCorrectionCalzetti, tplList[0]->m_igmCorrectionMeiksin);
    }

    if(EbmvCoeff>0.)
    {
        Int32 idxEbmv = -1;
        idxEbmv = identityTemplate.m_ismCorrectionCalzetti->GetEbmvIndex(EbmvCoeff); 

        if (idxEbmv!=-1)
            identityTemplate.ApplyDustCoeff(idxEbmv);
    }

    if(meiksinIdx>-1)
    { 
        identityTemplate.ApplyMeiksinCoeff(meiksinIdx);
    }

    const CSpectrumFluxAxis& extinction = identityTemplate.GetFluxAxis();

    UInt32 modelSize = kEnd - kStart + 1;
    //build the combined template model: spcmodel should have the same size as input spc
    CSpectrumSpectralAxis modelSpcAxis; 
    modelSpcAxis.extractFrom(spectrum.GetSpectralAxis().GetSamplesVector(), kStart, kEnd);
    CSpectrumFluxAxis modelFlux(modelSize, 0.0);
    for (Int32 iddl = 0; iddl < nddl; iddl++)
    {
        const CSpectrumFluxAxis & tmp = m_templatesRebined_bf[iddl].GetFluxAxis();       
        for(Int32 k = 0; k<modelSize; k++){
            modelFlux[k]+= amplitudes[iddl]*tmp[k+kStart]*extinction[k+kStart];
        }
    }
    spcPtr = std::make_shared<CModelSpectrumResult>(CSpectrum(std::move(modelSpcAxis), std::move(modelFlux)));

    // Deallocate the rebined template and mask buffers
    m_templatesRebined_bf.clear();
    m_masksRebined_bf.clear();   
    return 0;
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

Float64 COperatorTplcombination::GetNormFactor(const CSpectrumFluxAxis spcFluxAxis, UInt32 kStart, UInt32 n)
{
    Float64 maxabsval = DBL_MIN;
    for (Int32 k = 0; k < n; k++)
    {
        if(maxabsval<std::abs(spcFluxAxis[k+kStart]))
        {
            maxabsval=std::abs(spcFluxAxis[k+kStart]);
        }
    }
    return maxabsval;
}
//compute dtd on a specific range
Float64 COperatorTplcombination::ComputeDtD(const CSpectrumFluxAxis& spcFluxAxis,
                                            const TInt32Range& range)
{
    const CSpectrumNoiseAxis& spcError = spcFluxAxis.GetError();

    Float64 err2;
    Float64 dtd = .0;
    for(Int32 k=range.GetBegin(); k<=range.GetEnd(); k++)
    {
        err2 = spcError[k]*spcError[k];
        dtd += spcFluxAxis[k]*spcFluxAxis[k]/err2;
    }
    return dtd;
}