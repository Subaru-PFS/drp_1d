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
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>
//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_blas.h>
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


void COperatorTplcombination::BasicFit_preallocateBuffers(const CSpectrum& spectrum, const std::vector<CTemplate>& tplList)
{

    // Pre-Allocate the rebined template and mask with regard to the spectrum size
    m_templatesRebined_bf.resize(tplList.size());
    m_masksRebined_bf.resize(tplList.size());
    m_spcSpectralAxis_restframe.SetSize(spectrum.GetSampleCount());

    for(Int32 ktpl=0; ktpl<tplList.size(); ktpl++)
    {
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
        iEbmvCoeffMin = m_templatesRebined_bf[0].m_ismCorrectionCalzetti->GetEbmvIndex(fittingResults.ebmvCoeff);
        iEbmvCoeffMax = iEbmvCoeffMin;
    }
    // Linear fit
    //imin/max_lbda : valid spectrum samples within currentRange 
    //Int32 imin_lbda = m_spcSpectralAxis_restframe.GetIndexAtWaveLength(currentRange.GetBegin());
    //Int32 imax_lbda = m_spcSpectralAxis_restframe.GetIndexAtWaveLength(currentRange.GetEnd());

    Int32 n = kEnd-kStart +1;
    Log.LogDebug("  Operator-Tplcombination: prep. linear fitting with n=%d samples in the clamped lambdarange spectrum (imin=%d, lbda_min=%.3f - imax=%d, lbda_max=%.3f)", n, kStart, spcSpectralAxis[kStart], kEnd, spcSpectralAxis[kEnd]);
    
    gsl_matrix *X, *cov, *E;
    gsl_vector *y, *w, *c;

    X = gsl_matrix_alloc (n, nddl);
    y = gsl_vector_alloc (n);
    w = gsl_vector_alloc (n);
    c = gsl_vector_alloc (nddl);
    cov = gsl_matrix_alloc (nddl, nddl);
    E = gsl_matrix_alloc (n, n);//Extinction diagonal matrix

    // Normalizing factor
    Float64 normFactor;
    Float64 maxabsval = DBL_MIN;
    for (Int32 k = 0; k < n; k++)
    {
        if(maxabsval<std::abs(spcFluxAxis[k+kStart]))
        {
            maxabsval=std::abs(spcFluxAxis[k+kStart]);
        }
    }
    normFactor = maxabsval;
    if(verbose)
    {
        Log.LogDetail("  Operator-Tplcombination: Linear fitting, found normalization Factor=%e", normFactor);
    }

    // Prepare the fit data, once for all
    Float64 yi, ei, chisq;
    for(Int32 i = 0; i < n; i++)
    {
        yi = spcFluxAxis[i+kStart]/normFactor;
        ei = spcError[i+kStart]/normFactor;

        gsl_vector_set (y, i, yi); //y[i] = yi
        gsl_vector_set (w, i, 1.0/(ei*ei));//w[i] = 1/(ei*ei)
    }

    
    Bool option_igmFastProcessing = (MeiksinList.size()==1 ? false : true);
    TFloat64List chi2_outsideIGM(nIGM, 0.0);//save fit values to save computation
    //iterate over ism/igm; apply corrections as part of vector preparations for gsl; then call GSL
    fittingResults.chisquare = 0.;//final best Xi2 value
       
    //apply ism/igm corrections
    Bool igmLoopUseless_WavelengthRange = false;
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
        if(opt_extinction == 1){
            for (Int32 iddl = 0; iddl <nddl; iddl++)
            {
                m_templatesRebined_bf[iddl].ApplyMeiksinCoeff(meiksinIdx);
            } 
            if(!igmCorrectionAppliedOnce)
                igmLoopUseless_WavelengthRange = true; 
        }

        for(Int32 kEbmv = iEbmvCoeffMin; kEbmv<nISM; kEbmv++)
        {
            Int32 kEbmv_ = kEbmv - iEbmvCoeffMin; //index used to fill some arrays
            Float64 coeffEBMV = -1.; //no ism by default (ie DustCoeff=1.)
            //apply ism on all templates, once for all
            if (opt_dustFitting == 1){
                coeffEBMV = m_templatesRebined_bf[0].m_ismCorrectionCalzetti->GetEbmvValue(kEbmv);
                for (Int32 iddl = 0; iddl <nddl; iddl++)
                    m_templatesRebined_bf[iddl].ApplyDustCoeff(kEbmv);
            } 

            if(opt_extinction || opt_dustFitting)
            {
                const TFloat64List & ismCoeffList = m_templatesRebined_bf[0].GetcomputedDustCoeffs();
                const TFloat64List & igmCoeffList = m_templatesRebined_bf[0].GetcomputedMeiksinCoeffs();
                
                //init to zero and then fill the diag by the multi of ism/igmcoeff
                gsl_matrix_set_zero(E);
                for(Int32 i = 0; i<ismCoeffList.size(); i++)
                {
                    gsl_matrix_set (E, i, i, ismCoeffList[i]*igmCoeffList[i]);
                }

                // inverse the extinction diag matrix
                gsl_matrix *invE = InvertMatrix(E, n);
#define INVEXT(i,j) (gsl_matrix_get(invE,(i),(j)))
                //
                for(Int32 i = 0; i < n; i++)
                {
                    yi = spcFluxAxis[i+kStart]/normFactor*INVEXT(i,i);
                    ei = spcError[i+kStart]/normFactor*INVEXT(i,i);

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
                fittingResults.fittingAmplitudes[iddl] = a;
                err = COV(iddl,iddl)*normFactor;
                fittingResults.fittingErrors[iddl] = err;
            }
            fittingResults.fittingAmplitudesInterm[kEbmv_][kigm] = fittingResults.fittingAmplitudes;//saving 

            if(fittingResults.fittingAmplitudes.size()!=nddl)
            {
                Log.LogDebug("  Operator-Tplcombination: Found nfittedamps(=%d) different than nddl(=%d)", fittingResults.fittingAmplitudes.size(), nddl);
            }

            //build the combined template model: can be optimized
            TFloat64List spc_extract = TFloat64List(spcSpectralAxis.GetSamplesVector().begin() + kStart, 
                                                    spcSpectralAxis.GetSamplesVector().begin() + kEnd + 1);
            TFloat64List modelFlux(n, 0.0);
            for (Int32 iddl = 0; iddl < nddl; iddl++){
                const CSpectrumFluxAxis & tmp = m_templatesRebined_bf[iddl].GetFluxAxis();
                Float32 a = fittingResults.fittingAmplitudes[iddl];
                for(Int32 k = 0; k<n; k++){
                    modelFlux[k]+= a*tmp[k+kStart];
                }
            }
            //TODO: optimize the below
            CSpectrumSpectralAxis xAxis(std::move(spc_extract), n);
            CSpectrumFluxAxis yAxis(std::move(modelFlux));
            CTemplate modelSpec("combination", "", std::move(xAxis), std::move(yAxis));//m_templatesRebined_bf cant be used, thus we create a new CTemplate to apply corrections
            fittingResults.modelSpectrum = modelSpec;

            const CSpectrumFluxAxis& correctedFlux = modelSpec.GetFluxAxis();

            Bool optimizeComputation = true;
            Float64 optimized_Chi2 = NAN;
            if(optimizeComputation){
                const TInt32Range range(kStart, kEnd);
                    optimized_Chi2 = ComputeChi2_invCovBased(cov, 
                                                            fittingResults.fittingAmplitudes, 
                                                            spcFluxAxis, 
                                                            nddl, 
                                                            range, normFactor);
            }

            //estimate the lst-square brute force
            //now if we are in the configuration of fastigm
            Int32 chi2IgmSaved = 0;//indicator that chi2_igm is saved
            Float64 chi2_IGM = 0.;
            Int32 kEndloop = correctedFlux.GetSamplesCount()-1; //kEnd;
            if (option_igmFastProcessing && kigm>0) kEndloop = kIgmEnd - kStart;//i have a doubt about kStart and kStart/kEnd

            Float64 diff, err2;
            Float64 chi2Value = 0.;
            for(Int32 j=0; j<=kEndloop; j++)
            {
                if(option_igmFastProcessing && chi2IgmSaved == 0 && j>(kIgmEnd-kStart)) //since modelFkux is an extract of tpl
                {
                    chi2_IGM = chi2Value;//backup value for future use
                    chi2IgmSaved = 1;
                }

                diff = correctedFlux[j]-spcFluxAxis[j+kStart];//indeces should be verified
                err2 = spcError[j+kStart]*spcError[j+kStart];
                chi2Value += diff*diff/err2;

                if(option_igmFastProcessing && kigm==0)
                {
                    chi2_outsideIGM[kEbmv_] = chi2Value - chi2_IGM;
                }
                if(option_igmFastProcessing && kigm>0)
                {
                    chi2Value += chi2_outsideIGM[kEbmv_];
                }
            }
            /*
             bool xi2BrutForce = true;
             Float64 chi2Value = NAN;
             if(xi2BrutForce)
                chi2Value = ComputeXi2_bruteForce(correctedFlux, spcFluxAxis, kStart);
            */
            //temporary control
            if((optimized_Chi2 - chi2Value)>1E-6)
            {
                Log.LogError("  Operator-Tplcombination: OptimizedChi2: %f,  vs chi2Value: %f", optimized_Chi2, chi2Value);
                throw std::runtime_error("  Operator-Tplcombination: problem validating the Chi2 values" );
            }
            if(chi2Value < fittingResults.chisquare)
            {
                fittingResults.chisquare = chi2Value;
                fittingResults.igmIdx = meiksinIdx;
                fittingResults.ebmvCoeff = coeffEBMV;
            }

            //save the interm chisquares in the intermediate vector
            fittingResults.ChiSquareInterm[kEbmv_][kigm] = fittingResults.chisquare;
            fittingResults.IsmCalzettiCoeffInterm[kEbmv_][kigm] = coeffEBMV;
            fittingResults.IgmMeiksinIdxInterm[kEbmv_][kigm] = meiksinIdx;

            boost::chrono::thread_clock::time_point stop_postprocess = boost::chrono::thread_clock::now();
            Float64 duration_postprocess = boost::chrono::duration_cast<boost::chrono::microseconds>(stop_postprocess - start_postprocess).count();
            Log.LogDebug( "  Operator-Tplcombination: Linear fitting, postprocess = %.3f microsec", duration_postprocess);

        }//end iterating over ISM
    }//end iterating over IGM

    //fittingResults.modelSpectrum = CSpectrum(CSpectrumSpectralAxis(std::move(spc_extract)), CSpectrumFluxAxis(std::move(modelFlux)));

    fittingResults.snr = -1.0;

    gsl_matrix_free (X);
    gsl_vector_free (y);
    gsl_vector_free (w);
    gsl_vector_free (c);
    gsl_matrix_free (cov);
    gsl_matrix_free (E);

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

        CTemplate & itplTplSpectrum = m_templatesRebined_bf[ktpl];
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

    Int32 kStart = m_spcSpectralAxis_restframe.GetIndexAtWaveLength(currentRange.GetBegin());
    if(m_spcSpectralAxis_restframe[kStart]<spcLambdaRange_restframe.GetBegin() && kStart+1<m_spcSpectralAxis_restframe.GetSamplesCount())
    {
        kStart += 1;
    }
    Int32 kEnd = m_spcSpectralAxis_restframe.GetIndexAtWaveLength(currentRange.GetEnd());
    if(m_spcSpectralAxis_restframe[kEnd]>spcLambdaRange_restframe.GetEnd() && kEnd-1>=kStart)
    {
        kEnd -= 1;
    }
    currentRange.Set(m_spcSpectralAxis_restframe[kStart], m_spcSpectralAxis_restframe[kEnd]);
    
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
                                                                  Int32 opt_dustFitting,
                                                                  CPriorHelper::TPriorZEList logpriorze,
                                                                  Bool keepigmism,
                                                                  Float64 FitEbmvCoeff,
                                                                  Float64 FitMeiksinIdx)
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
        //temporarily commented the time ISM/IGM fitting is coded
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
    Int32 nEbmvCoeffs = 1;
    if(opt_dustFitting && !keepigmism)
    {
        nEbmvCoeffs = tplList.front().m_ismCorrectionCalzetti->GetNPrecomputedEbmvCoeffs();
    }
    Log.LogDebug("  Operator-tplcombination: prepare N ism coeffs = %d", nEbmvCoeffs);
    
    Int32 nIGMCoeffs = 1;
    if(opt_extinction && !keepigmism)
    {
        nIGMCoeffs = tplList.front().m_igmCorrectionMeiksin->GetIdxCount();
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

    TFloat64List _chi2List(nIGMCoeffs, DBL_MAX);
    TFloat64List _ismList(nIGMCoeffs, -1.0);
    TInt32List   _igmList(nIGMCoeffs, -1);
    TFloat64List _ampList(tplList.size(), NAN);
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
        //
        STplcombination_basicfitresult fittingResults;
        if(keepigmism && opt_dustFitting && opt_extinction)
        {
            fittingResults.igmIdx = FitMeiksinIdx;
            fittingResults.ebmvCoeff = FitEbmvCoeff;
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
        fittingResults.fittingAmplitudes = TFloat64List(tplList.size(), NAN);
        fittingResults.fittingErrors = TFloat64List(tplList.size(), NAN);
        fittingResults.overlapRate = 0.0;
        fittingResults.status = COperator::nStatus_DataError;
        //

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

gsl_matrix * COperatorTplcombination::InvertMatrix(gsl_matrix* m, UInt32 dim)
{
    int s;
    // Compute the  inverse of the LU decomposition
    gsl_permutation * p = gsl_permutation_alloc (dim);

    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp (m, p, &s);  

    gsl_matrix *inv = gsl_matrix_alloc(dim, dim);
    gsl_linalg_LU_invert(m, p, inv);

    gsl_permutation_free (p);
    return inv;
}
/**
 * Xi2 = dt.N-1.d -at.C-1.a (where all elements are vect or matrices)
 * xi2Value (z,ism,igm) = SumD -SumAC
 *       SumD = Sum(di*di*ni); //di and ni = 1/err2 correspond to data and error on data
 *       SumAC= Sum(Sum((ai*aj*bij))); //ai corresponds to gsl-computed amplitude per template and bij corresponds to elements of the inverse of covariance matrix
 **/
Float64 COperatorTplcombination::ComputeChi2_invCovBased(gsl_matrix* cov, 
                                                         const TFloat64List& amplitudes, 
                                                         const CSpectrumFluxAxis& spcFlux, 
                                                         UInt32 nTpl,
                                                         const TInt32Range& spcRange, Float64 normFactor)
{
    //compute the inverse of covariance matrix:
    gsl_matrix *invCov = InvertMatrix(cov, nTpl);
#define INVCOV(i,j) (gsl_matrix_get(invCov,(i),(j)))
#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))
    const CSpectrumNoiseAxis& spcError = spcFlux.GetError();

    UInt32 kStart = spcRange.GetBegin();
    UInt32 kEnd = spcRange.GetEnd();

    Float64 sumD = 0.;
    for(Int32 i = kStart; i<= kEnd; i++)
    {
        sumD+= spcFlux[i]*spcFlux[i]/spcError[i]/spcError[i];
    }

    Float64 sumAC = 0.;
    for(Int32 i = 0; i<amplitudes.size(); i++)
    {
        for(Int32 j = 0; j<amplitudes.size(); j++)
        {
            sumAC+= amplitudes[i]*amplitudes[j]/normFactor/normFactor*INVCOV(i,j);
        }
    }

    Float64 fit = sumD - sumAC;

    gsl_matrix_free (invCov);
    return fit;
}