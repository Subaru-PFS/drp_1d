#include <RedshiftLibrary/operator/templatefitting.h>

#include <RedshiftLibrary/spectrum/axis.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/tools.h>
#include <RedshiftLibrary/common/mask.h>
#include <RedshiftLibrary/operator/templatefittingresult.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/common/quicksort.h>
#include <RedshiftLibrary/log/log.h>

#include <boost/numeric/conversion/bounds.hpp>

#include <math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <algorithm>    // std::sort
#include <float.h>

#include <sstream>
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <assert.h>
#include <RedshiftLibrary/processflow/datastore.h>
#define NOT_OVERLAP_VALUE NAN
#include <stdio.h>

namespace bfs = boost::filesystem;
using namespace NSEpic;
using namespace std;

/**
 * @brief COperatorTemplateFitting::BasicFit
 * @param spectrum
 * @param tpl
 * @param pfgTplBuffer
 * @param lambdaRange
 * @param redshift
 * @param overlapThreshold
 * @param overlapRate
 * @param chiSquare
 * @param fittingAmplitude
 * @param fittingAmplitudeError
 * @param fittingAmplitudeNegative
 * @param fittingDtM
 * @param fittingMtM
 * @param fittingDustCoeff
 * @param fittingMeiksinIdx
 * @param status
 * @param ChiSquareInterm
 * @param IsmCalzettiCoeffInterm
 * @param IgmMeiksinIdxInterm
 * @param opt_interp
 * @param forcedAmplitude
 * @param opt_extinction
 * @param opt_dustFitting : -1 = disabled, -10 = fit over all available indexes, positive integer 0, 1 or ... will be used as ism-calzetti index as initialized in constructor.
 * @param spcMaskAdditional
 * @param priorjoint_pISM_tpl_z : vector size = nISM, joint prior p(ISM, TPL, Z)
 */
void COperatorTemplateFitting::BasicFit(const CSpectrum& spectrum,
                                   const CTemplate& tpl,
                                   const TFloat64Range& lambdaRange,
                                   Float64 redshift,
                                   Float64 overlapThreshold,
                                   Float64& overlapRate,
                                   Float64& chiSquare,
                                   Float64& fittingAmplitude,
                                   Float64& fittingAmplitudeError,
                                   Bool& fittingAmplitudeNegative,
                                   Float64& fittingDtM,
                                   Float64& fittingMtM,
                                   Float64& fittingLogprior,
                                   Float64 &fittingDustCoeff,
                                   Float64 &fittingMeiksinIdx,
                                   EStatus& status,
                                   std::vector<TFloat64List>& ChiSquareInterm,
                                   std::vector<TFloat64List>& IsmCalzettiCoeffInterm,
                                   std::vector<TInt32List>& IgmMeiksinIdxInterm,
                                   std::string opt_interp,
                                   Float64 forcedAmplitude,
                                   Int32 opt_extinction,
                                   Int32 opt_dustFitting,
                                   CMask spcMaskAdditional,
                                   CPriorHelper::TPriorEList logpriore,
                                   bool keepigmism)
{
    bool verbose = false;
    bool amplForcePositive=true;
    chiSquare = boost::numeric::bounds<float>::highest();
    bool status_chisquareSetAtLeastOnce = false;

    for(Int32 kism=0; kism<ChiSquareInterm.size(); kism++)
    {
        for(Int32 kigm=0; kigm<ChiSquareInterm[kism].size(); kigm++)
        {
            ChiSquareInterm[kism][kigm] = boost::numeric::bounds<float>::highest();
            IsmCalzettiCoeffInterm[kism][kigm] = -1.0;
            IgmMeiksinIdxInterm[kism][kigm] = -1;
        }
    }
    fittingAmplitude = -1.0;
    fittingAmplitudeError = -1.0;
    fittingAmplitudeNegative = 0;
    overlapRate = 0.0;
    status = nStatus_DataError;

    const CSpectrumFluxAxis& spcFluxAxis = spectrum.GetFluxAxis();
    const TAxisSampleList & Yspc = spcFluxAxis.GetSamplesVector();

    if(spcMaskAdditional.GetMasksCount()!=spcFluxAxis.GetSamplesCount())
    {
        Log.LogInfo("  Operator-TemplateFitting: spcMaskAdditional does not have the same size as the spectrum flux vector... (%d vs %d), aborting!", spcMaskAdditional.GetMasksCount(),  spectrum.GetFluxAxis().GetSamplesCount());
        status = nStatus_DataError;
        return ;
    }

    TFloat64Range currentRange;

    Int32 ret = RebinTemplate(  spectrum, 
                                tpl,
                                redshift, 
                                lambdaRange,
                                opt_interp, 
                                currentRange,
                                overlapRate,
                                overlapThreshold);
    
    const TAxisSampleList & Xspc = m_spcSpectralAxis_restframe.GetSamplesVector();
    bool apply_ism = ( (opt_dustFitting==-10 || opt_dustFitting>0) ? true : false);
    
    if (apply_ism || opt_extinction){             
        m_templateRebined_bf.InitIsmIgmConfig(tpl.m_ismCorrectionCalzetti, tpl.m_igmCorrectionMeiksin);
    }

    if( ret == -1 ){
        status = nStatus_NoOverlap; 
        return;
    }
    if( ret == -2 ){
        status = nStatus_DataError;
        return;
    }

    // Optionally Apply some Calzetti Extinction for DUST
    Int32 nDustCoeffs=1;
    Int32 iDustCoeffMin = 0;
    Int32 iDustCoeffMax = iDustCoeffMin+1;
    if(apply_ism)
    {
        if(opt_dustFitting>0)
        {
            nDustCoeffs = 1;
            iDustCoeffMin = opt_dustFitting;
            iDustCoeffMax = iDustCoeffMin;
        }else if(opt_dustFitting==-10)
        {
            nDustCoeffs = m_templateRebined_bf.m_ismCorrectionCalzetti->GetNPrecomputedDustCoeffs();
            iDustCoeffMin = 0;
            iDustCoeffMax = iDustCoeffMin+nDustCoeffs-1;
        }
    }else{
        nDustCoeffs = 1;
        iDustCoeffMin = -1;
        iDustCoeffMax = iDustCoeffMin;
    }

    //Optionally apply some IGM absorption
    Int32 nIGMCoeffs=1;
    TInt32List MeiksinList;
    if(opt_extinction)
    {
        if(keepigmism){
            nIGMCoeffs = 1;
            MeiksinList.push_back(fittingMeiksinIdx);//fill it with only the index passed as argument   
        }
        else{
            nIGMCoeffs = m_templateRebined_bf.m_igmCorrectionMeiksin->GetIdxCount();
            for(Int32 mk = 0; mk<nIGMCoeffs; mk++){
                MeiksinList.push_back(mk);
            }
        }
    }else{//at least have one element
        MeiksinList.push_back(-1);
    }

    Bool option_igmFastProcessing = true; //todo: find a way to unit-test this acceleration
    //Prepare the wavelengthRange Limits
    Log.LogDebug( "  Operator-TemplateFitting: currentRange_lbda_min=%f, currentRange_lbda_max=%f", currentRange.GetBegin(), currentRange.GetEnd());
    std::vector<Float64> sumCross_outsideIGM(nDustCoeffs, 0.0);
    std::vector<Float64>  sumT_outsideIGM(nDustCoeffs, 0.0);
    std::vector<Float64>  sumS_outsideIGM(nDustCoeffs, 0.0);
    
    Float64 lbda_max, lbdaMax_IGM;
    if (opt_extinction)
    {
        lbdaMax_IGM = m_templateRebined_bf.m_igmCorrectionMeiksin->GetLambdaMax()*(1+redshift);
    }
    //Loop on the meiksin Idx
    Bool igmLoopUseless_WavelengthRange = false;
    for(Int32 kM=0; kM<nIGMCoeffs; kM++)
    {
        if(igmLoopUseless_WavelengthRange)
        {
            //Now copy from the already calculated k>0 igm values
            for(Int32 kism=0; kism<ChiSquareInterm.size(); kism++)
            {
                for(Int32 kigm=1; kigm<ChiSquareInterm[kism].size(); kigm++)
                {
                    ChiSquareInterm[kism][kigm] = ChiSquareInterm[kism][0];
                    IsmCalzettiCoeffInterm[kism][kigm] = IsmCalzettiCoeffInterm[kism][0];
                    IgmMeiksinIdxInterm[kism][kigm] = IgmMeiksinIdxInterm[kism][0];
                }
            }
            break;
        }
        Int32 meiksinIdx = MeiksinList[kM]; //index for the Meiksin curve (0-6; 3 being the median extinction value)
        lbda_max = currentRange.GetEnd();  
        if(option_igmFastProcessing && meiksinIdx>0){
            lbda_max = std::min(lbdaMax_IGM, lbda_max);
        }

        TFloat64Range range(currentRange.GetBegin(), lbda_max);
        Int32 kStart = -1, kEnd = -1;
        m_templateRebined_bf.SetIsmIgmLambdaRange( range );
        m_templateRebined_bf.GetIsmIgmRangeIndex(kStart, kEnd);

        if(kStart==-1 || kEnd==-1)
        {
            Log.LogDebug( "  Operator-TemplateFitting: kStart=%d, kEnd=%d ! Aborting.", kStart, kEnd);
            break;
        }

        bool igmCorrectionAppliedOnce = false;
        //Meiksin IGM extinction
        if(opt_extinction)
        {
            igmCorrectionAppliedOnce = m_templateRebined_bf.ApplyMeiksinCoeff(meiksinIdx, redshift);
            if(!igmCorrectionAppliedOnce)
                igmLoopUseless_WavelengthRange = true; 
        }
        //Loop on the EBMV dust coeff
        for(Int32 kDust=iDustCoeffMin; kDust<=iDustCoeffMax; kDust++)
        {
            Int32 kDust_ = kDust - iDustCoeffMin; //index used to fill some arrays
            Float64 coeffEBMV = -1.; //no ism by default (ie DustCoeff=1.)

            if (apply_ism)
            {
                coeffEBMV = m_templateRebined_bf.m_ismCorrectionCalzetti->GetEbmvValue(kDust);
                //check that we got the same coeff:
                /*if(keepigmism && (coeffEBMV - fittingDustCoeff)< DBL_EPSILON){//comparing floats
                    Log.LogDebug("Keepigmism: coeffEBMW corresponds to passed param: %f vs %f", coeffEBMV, fittingDustCoeff);
                }*/
                m_templateRebined_bf.ApplyDustCoeff(kDust);
            }
            //*/
            /*//debug:
            // save final ISM/IGM template model
            if(redshift>=2.4 && redshift<2.4001 && meiksinIdx==6){
                FILE* f = fopen( "templateFitting_template_finalModel.txt", "w+" );
                for(Int32 m=0; m<itplTplSpectralAxis.GetSamplesCount(); m++){
                    fprintf( f, "%e %e\n", Xtpl[m], Ytpl[m]);
                }
                fclose( f );
            }
            //*/

            /*/
            // Optionally mask pixels far from the breaks
            bool opt_onlyBreaks = false;
            if(opt_onlyBreaks)
            {
                //add the ranges to be processed
                TFloat64RangeList restLambdaRanges_A;
                TFloat64RangeList restLambdaRanges_B;

                restLambdaRanges_A.push_back(TFloat64Range( 1043.0, 1174.0 )); //Lya, A
                restLambdaRanges_B.push_back(TFloat64Range( 1304.0, 1369.0 )); //Lya, B

                restLambdaRanges_A.push_back(TFloat64Range( 3200.0, 3600.0 )); //OII, A
                restLambdaRanges_B.push_back(TFloat64Range( 4000.0, 4200.0 )); //OII, B

                restLambdaRanges_A.push_back(TFloat64Range( 4290.0, 4830.0 )); //OIII, A
                restLambdaRanges_B.push_back(TFloat64Range( 5365.0, 5635.0 )); //OIII, B

                restLambdaRanges_A.push_back(TFloat64Range( 5632.0, 6341.0 )); //Halpha, A
                restLambdaRanges_B.push_back(TFloat64Range( 7043.0, 7397.6 )); //Halpha, B

                Float64 z = redshift;
                for(Int32 k=0; k<itplTplSpectralAxis.GetSamplesCount(); k++)
                {
                    if(Xtpl[k] < lbda_min){
                        continue;
                    }
                    if(Xtpl[k] > lbda_max){
                        continue;
                    }

                    Float64 restLambda = Xtpl[k]/(1.0+z);
                    bool inBreakRange=false;
                    for(Int32 kRanges=0; kRanges<restLambdaRanges_A.size(); kRanges++)
                    {
                        if(restLambda>=restLambdaRanges_A[kRanges].GetBegin() && restLambda<=restLambdaRanges_B[kRanges].GetEnd())
                        {
                            inBreakRange = true;
                            break;
                        }
                    }
                    if(!inBreakRange)
                    {
                        spcMaskAdditional[k]=false;
                    }

                }
            }
            //*/
            //const CTemplate& tmp = m_templateRebined_bf;//obligatory, otherwise non-const GetFluxAxis is called
            const CSpectrumFluxAxis& tplIsmIgm = m_templateRebined_bf.GetFluxAxis();

            const TAxisSampleList & Ytpl = tplIsmIgm.GetSamplesVector();
            Float64 sumCross = 0.0;
            Float64 sumT = 0.0;
            Float64 sumS = 0.0;

            Float64 sumCross_IGM = 0.0;
            Float64 sumT_IGM = 0.0;
            Float64 sumS_IGM = 0.0;
            Int32 sumsIgmSaved = 0;

            Float64 err2 = 0.0;
            Float64 fit = 0;
            Int32 numDevs = 0;
            Int32 numDevsFull = 0;
            const TFloat64List& error = spcFluxAxis.GetError();

            for(Int32 j=kStart; j<=kEnd; j++)
            {
                numDevsFull++;
                if(spcMaskAdditional[j]){

                    /*
                    //check for invalid data
                    if(error[j]!=error[j]){
                        Log.LogDebug("  Operator-TemplateFitting: noise invalid=%e for i=%d", error[j], j);
                        Log.LogDebug("  Operator-TemplateFitting: noise invalid, w=%e for i=%d", spcSpectralAxis_restframe[j], j);
                        Log.LogDebug("  Operator-TemplateFitting: noise invalid,samples count=%d", spcSpectralAxis_restframe.GetSamplesCount(), j);
                    }
                    //*/

                    numDevs++;
                    err2 = 1.0 / (error[j] * error[j]);

                    // Tonry&Davis formulation
                    sumCross+=Yspc[j]*Ytpl[j]*err2;
                    sumT+=Ytpl[j]*Ytpl[j]*err2;
                    //sumCross+=Yspc[j]*Ytpl[j];
                    //sumT+=Ytpl[j]*Ytpl[j];

                    sumS+= Yspc[j]*Yspc[j]*err2;

                    if( std::isinf(err2) || std::isnan(err2) ){
                        Log.LogError("  Operator-TemplateFitting: found invalid inverse variance : err2=%e, for index=%d at restframe wl=%f", err2, j, m_spcSpectralAxis_restframe[j]);
                        status = nStatus_InvalidProductsError;
                        return ;
                    }

                    if( std::isinf(sumS) || std::isnan(sumS) || sumS!=sumS ){
                        Log.LogError("  Operator-TemplateFitting: found invalid dtd : dtd=%e, for index=%d at restframe wl=%f", sumS, j, m_spcSpectralAxis_restframe[j]);
                        Log.LogError("  Operator-TemplateFitting: found invalid dtd : Yspc=%e, for index=%d at restframe wl=%f", Yspc[j], j, m_spcSpectralAxis_restframe[j]);
                        Log.LogError("  Operator-TemplateFitting: found invalid dtd : err2=%e, for index=%d at restframe wl=%f", err2, j, m_spcSpectralAxis_restframe[j]);
                        Log.LogError("  Operator-TemplateFitting: found invalid dtd : error=%e, for index=%d at restframe wl=%f", error[j], j, m_spcSpectralAxis_restframe[j]);
                        status = nStatus_InvalidProductsError;
                        return ;
                    }

                    if( std::isinf(sumT) || std::isnan(sumT) ){
                        Log.LogError("  Operator-TemplateFitting: found invalid mtm : mtm=%e, for index=%d at restframe wl=%f", sumT, j, m_spcSpectralAxis_restframe[j]);
                        status = nStatus_InvalidProductsError;
                        return ;
                    }

                    if(option_igmFastProcessing && meiksinIdx==0)
                    {
                        //store intermediate sums for IGM range
                        if(sumsIgmSaved==0)
                        {
                            if(Xspc[j]>lbdaMax_IGM)
                            {
                                sumCross_IGM = sumCross;
                                sumT_IGM = sumT;
                                sumS_IGM = sumS;
                                sumsIgmSaved = 1;
                            }
                        }
                    }
                }
            }
            if(option_igmFastProcessing && meiksinIdx==0 && sumsIgmSaved==1)
            {
                sumCross_outsideIGM[kDust_] = sumCross-sumCross_IGM;
                sumT_outsideIGM[kDust_] = sumT-sumT_IGM;
                sumS_outsideIGM[kDust_] = sumS-sumS_IGM;
            }
            if(option_igmFastProcessing && meiksinIdx>0)
            {
                sumCross += sumCross_outsideIGM[kDust_];
                sumT += sumT_outsideIGM[kDust_];
                sumS += sumS_outsideIGM[kDust_];
            }

            if ( numDevs==0 )
            {
                status = nStatus_DataError;
                return;
            }

            Float64 ampl = 0.0;
            Float64 ampl_err = 0.0;
            Bool ampl_neg = 0;
            bool apply_priore = false;
            if(!logpriore.empty() && !m_templateRebined_bf.CalzettiInitFailed())
            { 
                if(logpriore.size()==m_templateRebined_bf.m_ismCorrectionCalzetti->GetNPrecomputedDustCoeffs())
                {
                    if(logpriore[kDust].A_sigma>0.0 && logpriore[kDust].betaA>0.0)
                        apply_priore = true;
                }
            }

            if( sumT==0 )
            {
                ampl = 0.0;
                ampl_err = 0.0;
                fit = sumS;
                ampl_neg = 0;
                //status = nStatus_DataError;
                //return;
            }else{
                
                ampl = sumCross/sumT;
                ampl_err = sqrt(1./sumT);
                
                if (apply_priore)
                {
                    Float64 bss2 = logpriore[kDust].betaA/(logpriore[kDust].A_sigma*logpriore[kDust].A_sigma);
                    ampl = (sumCross+logpriore[kDust].A_mean*bss2)/(sumT+bss2);
                    ampl_err = sqrt(sumT)/(sumT+bss2); 
                    /*
                    Log.LogInfo("  Operator-TemplateFitting: Constrained amplitude (betaA=%e):  s2b=%e, mtm=%e", logpriore[kDust].betaA, s2b, sumT);
                    //*/
                }

                if(forcedAmplitude !=-1){
                    ampl = forcedAmplitude;
                    ampl_err = 0.;
                }

                if (ampl < -3*ampl_err)
                {
                    ampl_neg = 1;
                }

                if(amplForcePositive)
                {
                    ampl = max(0.0, ampl);
                }

                //Generalized method (ampl can take any value now) for chi2 estimate
                fit = sumS + sumT*ampl*ampl - 2.*ampl*sumCross;
            }

            //*
            if(verbose)
            {
                Log.LogDebug("");
                Log.LogDebug("  Operator-TemplateFitting: z=%f", redshift);
                Log.LogDebug("  Operator-TemplateFitting: fit=%e", fit);
                Log.LogDebug("  Operator-TemplateFitting: sumT=%e", sumT);
                Log.LogDebug("  Operator-TemplateFitting: sumS=%e", sumS);
                Log.LogDebug("  Operator-TemplateFitting: sumCross=%e", sumCross);
                Log.LogDebug("  Operator-TemplateFitting: coeffEBMV=%.2f", coeffEBMV);
                Log.LogDebug("  Operator-TemplateFitting: meiksinIdx=%d", meiksinIdx);

            }
            //*/

            Float64 logprior = 0.;
            if (apply_priore)
            {
                logprior += -2.*logpriore[kDust].betaTE*logpriore[kDust].logprior_precompTE;
                logprior += -2.*logpriore[kDust].betaA*logpriore[kDust].logprior_precompA;
                logprior += -2.*logpriore[kDust].betaZ*logpriore[kDust].logprior_precompZ;
                if(logpriore[kDust].A_sigma>0.0)
                {
                    Float64 logPa = logpriore[kDust].betaA*(ampl-logpriore[kDust].A_mean)*(ampl-logpriore[kDust].A_mean)/(logpriore[kDust].A_sigma*logpriore[kDust].A_sigma);
                    logprior += logPa;
                }
                fit += logprior;
            }

            /*
            if(redshift==0.0)
            {
                Log.LogInfo( "templateFitting operator: for z=%f dtd = %e", redshift, sumS);
                Log.LogInfo( "templateFitting operator: for z=%f dtm = %e", redshift, sumCross);
                Log.LogInfo( "templateFitting operator: for z=%f mtm = %e", redshift, sumT);
            }
            //*/

            /*
            //mask correction coefficient for the masked samples
            Float64 maskedSamplesCorrection = (Float64)numDevsFull/(Float64)numDevs;
            fit *= maskedSamplesCorrection;
            //*/

            /*
            Float64 overlapCorrection = 1.0/overlapRate;
            fit *= overlapCorrection;
            //*/
            //Note: to avoid coresegmentation errors
            ChiSquareInterm[kDust_][kM] = fit;
            IsmCalzettiCoeffInterm[kDust_][kM] = coeffEBMV;
            IgmMeiksinIdxInterm[kDust_][kM] = meiksinIdx;

            if(fit<chiSquare)
            {
                chiSquare = fit;
                fittingDustCoeff = coeffEBMV;
                fittingMeiksinIdx = meiksinIdx;
                fittingAmplitude = ampl;
                fittingAmplitudeError = ampl_err;
                fittingAmplitudeNegative = ampl_neg;
                fittingDtM = sumCross;
                fittingMtM = sumT;
                fittingLogprior = logprior;

                status_chisquareSetAtLeastOnce = true;
            }
        }
    }

    if(status_chisquareSetAtLeastOnce)
    {
        status = nStatus_OK;
    }else{
        status = nStatus_LoopError;
    }
}

Int32  COperatorTemplateFitting::RebinTemplate( const CSpectrum& spectrum,
                                const CTemplate& tpl, 
                                Float64 redshift,
                                const TFloat64Range& lambdaRange,
                                std::string opt_interp,
                                //return variables
                                TFloat64Range& currentRange,
                                Float64& overlapRate,
                                Float64 overlapThreshold)// const
{
    Float64 onePlusRedshift = 1.0 + redshift;

    //shift lambdaRange backward to be in restframe
    TFloat64Range spcLambdaRange_restframe;
    TFloat64Range lambdaRange_restframe( lambdaRange.GetBegin() / onePlusRedshift,
                                         lambdaRange.GetEnd() / onePlusRedshift );

    //redshift in restframe the tgtSpectralAxis, i.e., division by (1+Z)
    m_spcSpectralAxis_restframe.ShiftByWaveLength(spectrum.GetSpectralAxis(), onePlusRedshift, CSpectrumSpectralAxis::nShiftBackward);
    m_spcSpectralAxis_restframe.ClampLambdaRange( lambdaRange_restframe, spcLambdaRange_restframe );
                                         
    // Compute clamped lambda range over template in restframe
    TFloat64Range tplLambdaRange;
    const CSpectrumSpectralAxis& tplSpectralAxis = tpl.GetSpectralAxis();
    tplSpectralAxis.ClampLambdaRange( lambdaRange_restframe, tplLambdaRange );
    // Compute the intersected range
    TFloat64Range intersectedLambdaRange( 0.0, 0.0 );
    TFloat64Range::Intersect( tplLambdaRange, spcLambdaRange_restframe, intersectedLambdaRange );

    tpl.Rebin( intersectedLambdaRange, m_spcSpectralAxis_restframe, m_templateRebined_bf, m_mskRebined_bf, opt_interp);   

    //overlapRate
    overlapRate = m_spcSpectralAxis_restframe.IntersectMaskAndComputeOverlapRate( lambdaRange_restframe, m_mskRebined_bf );

    // Check for overlap rate
    if( overlapRate < overlapThreshold || overlapRate<=0.0 )
    {
        //status = nStatus_NoOverlap; 
        return -1 ;
    }

    TFloat64Range logIntersectedLambdaRange( log( intersectedLambdaRange.GetBegin() ), log( intersectedLambdaRange.GetEnd() ) );
    //the spectral axis should be in the same scale
    currentRange = logIntersectedLambdaRange;
    if( m_spcSpectralAxis_restframe.IsInLinearScale() != tplSpectralAxis.IsInLinearScale() )
    {
        Log.LogError( "    chisquare operator: data and model not in the same scale (lin/log) ! Aborting.");
        //status = nStatus_DataError;
        return -2;
    }
    if(m_spcSpectralAxis_restframe.IsInLinearScale()){
        currentRange = intersectedLambdaRange;
    }
    return 0;
}
/**
 * \brief
 *
 * input: if additional_spcMasks size is 0, no additional mask will be used, otherwise its size should match the redshifts list size
 *
 * opt_dustFitting: -1 = disabled, -10 = fit over all available indexes, positive integer 0, 1 or ... will be used as ism-calzetti index as initialized in constructor.
 **/
std::shared_ptr<COperatorResult> COperatorTemplateFitting::Compute(const CSpectrum& spectrum,
                                                              const CTemplate& tpl,
                                                              const TFloat64Range& lambdaRange,
                                                              const TFloat64List& redshifts,
                                                              Float64 overlapThreshold,
                                                              std::vector<CMask> additional_spcMasks,
                                                              std::string opt_interp,
                                                              Int32 opt_extinction,
                                                              Int32 opt_dustFitting,
                                                              CPriorHelper::TPriorZEList logpriorze,
                                                              Bool keepigmism,
                                                              Float64 FitDustCoeff,
                                                              Float64 FitMeiksinIdx)
{
    Log.LogDetail("  Operator-TemplateFitting: starting computation for template: %s", tpl.GetName().c_str());

    if(0)
    {
        //CSpectrumFluxAxis tmp_tplFluxAxis = tpl.GetFluxAxis();
        CSpectrumSpectralAxis tmp_tplSpectralAxis = tpl.GetSpectralAxis();
        for(UInt32 k=0; k<std::min(Int32(tmp_tplSpectralAxis.GetSamplesCount()), Int32(10)); k++)
        {
            Log.LogDebug("  Operator-TemplateFitting: tpl_SpectralAxis[%d] = %f",
                         k,
                         tmp_tplSpectralAxis[k]);
        }
    }

    if( (opt_dustFitting==-10 || opt_dustFitting>-1) && tpl.CalzettiInitFailed())
    {
        Log.LogError("  Operator-TemplateFitting: no calzetti calib. file loaded in template... aborting!");
        throw std::runtime_error("  Operator-TemplateFitting: no calzetti calib. file in template");
    }
    if( opt_dustFitting>-1 && opt_dustFitting>tpl.m_ismCorrectionCalzetti->GetNPrecomputedDustCoeffs()-1)
    {
        Log.LogError("  Operator-TemplateFitting: calzetti index overflow (opt=%d, while NPrecomputedDustCoeffs=%d)... aborting!",
                     opt_dustFitting,
                     tpl.m_ismCorrectionCalzetti->GetNPrecomputedDustCoeffs());
        throw std::runtime_error("  Operator-TemplateFitting: calzetti index overflow");
    }

    if( opt_extinction && tpl.MeiksinInitFailed())
    {
        Log.LogError("  Operator-TemplateFitting: no meiksin calib. file loaded in template... aborting!");
        throw std::runtime_error("  Operator-TemplateFitting: no meiksin calib. file in template");
    }

    if( spectrum.GetSpectralAxis().IsInLinearScale() == false || tpl.GetSpectralAxis().IsInLinearScale() == false )
    {
        Log.LogError("  Operator-TemplateFitting: input spectrum or template are not in log scale (ignored)");
        //return NULL;
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

    std::shared_ptr<CTemplateFittingResult> result = std::shared_ptr<CTemplateFittingResult>( new CTemplateFittingResult() );
    Int32 nDustCoeffs=1;
    if(opt_dustFitting==-10)
    {
        nDustCoeffs = tpl.m_ismCorrectionCalzetti->GetNPrecomputedDustCoeffs();
    }
    Int32 nIGMCoeffs=1;
    if(opt_extinction && !keepigmism)
    { 
        nIGMCoeffs = tpl.m_igmCorrectionMeiksin->GetIdxCount();
    }

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
        Log.LogError("  Operator-TemplateFitting: using default mask, masks-list size (%d) didn't match the input redshift-list (%d) !)", additional_spcMasks.size(), sortedRedshifts.size());
    }
    if(logpriorze.size()>0 && logpriorze.size()!=sortedRedshifts.size())
    {
        Log.LogError("  Operator-TemplateFitting: prior list size (%d) didn't match the input redshift-list (%d) !)", logpriorze.size(), sortedRedshifts.size());
        throw std::runtime_error("  Operator-TemplateFitting: prior list size didn't match the input redshift-list size");
    }

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
        if( keepigmism /*&& opt_extinction && opt_dustFitting>-1 */){
            result->FitDustCoeff[i] = FitDustCoeff;
            result->FitMeiksinIdx[i] = FitMeiksinIdx;
        }
        //TODO: reinitialize m_fluxAxisIsmIgm prior to calling BasicFit
        BasicFit( spectrum,
                  tpl,
                  lambdaRange,
                  redshift,
                  overlapThreshold,
                  result->Overlap[i],
                  result->ChiSquare[i],
                  result->FitAmplitude[i],
                  result->FitAmplitudeError[i],
                  result->FitAmplitudeNegative[i],
                  result->FitDtM[i],
                  result->FitMtM[i],
                  result->LogPrior[i],
                  result->FitDustCoeff[i],
                  result->FitMeiksinIdx[i],
                  result->Status[i],
                  result->ChiSquareIntermediate[i],
                  result->IsmDustCoeffIntermediate[i],
                  result->IgmMeiksinIdxIntermediate[i],
                  opt_interp,
                  -1,
                  opt_extinction,
                  opt_dustFitting,
                  additional_spcMask,
                  logp,
                  keepigmism);

        if(result->Status[i]==nStatus_InvalidProductsError)
        {
            Log.LogError("  Operator-TemplateFitting: found invalid chisquare products for z=%f. Now breaking z loop.", redshift);
            break;
        }

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
        Log.LogInfo("  Operator-TemplateFitting: overlap warning for %s: minz=%.3f, maxz=%.3f", tpl.GetName().c_str(), overlapValidInfZ, overlapValidSupZ);
    }

    //only bad status warning
    Int32 oneValidStatusFoundIndex = -1;
    for (Int32 i=0;i<sortedRedshifts.size();i++)
    {
        if(result->Status[i]==nStatus_OK)
        {
            oneValidStatusFoundIndex=i;
            Log.LogDebug("  Operator-TemplateFitting: STATUS VALID found for %s: at least at index=%d", tpl.GetName().c_str(), i);
            break;
        }
    }if(oneValidStatusFoundIndex==-1)
    {
        Log.LogWarning("  Operator-TemplateFitting: STATUS WARNING for %s: Not even one single valid fit/merit value found", tpl.GetName().c_str());
    }


    //loop error status warning
    Int32 loopErrorStatusFoundIndex = -1;
    for (Int32 i=0;i<sortedRedshifts.size();i++)
    {
        if(result->Status[i]==nStatus_LoopError)
        {
            loopErrorStatusFoundIndex=i;
            Log.LogDebug("  Operator-TemplateFitting: STATUS Loop Error found for %s: at least at index=%d", tpl.GetName().c_str(), i);
            break;
        }
    }if(loopErrorStatusFoundIndex!=-1)
    {
        Log.LogWarning("    templateFitting-operator: Loop Error - chisquare values not set even once");
    }

    //estimate CstLog for PDF estimation
    result->CstLog = EstimateLikelihoodCstLog(spectrum, lambdaRange);

    return result;

}




const COperatorResult* COperatorTemplateFitting::ExportChi2versusAZ(const CSpectrum& spectrum, const CTemplate& tpl,
                          const TFloat64Range& lambdaRange, const TFloat64List& redshifts,
                          Float64 overlapThreshold )
{
    if( spectrum.GetSpectralAxis().IsInLinearScale() == false || tpl.GetSpectralAxis().IsInLinearScale() == false )
    {
        Log.LogError("  Operator-TemplateFitting: input spectrum or template are not in log scale (ignored)");
        //return NULL;
    }

    /*//debug:
    // save templateFine
    FILE* f = fopen( "template_fine.txt", "w+" );
    for(Int32 m=0; m<nTgt; m++){
        if( Yfine[m] < 0.0001 ){
            fprintf( f, "%e %e\n", Xfine[m], Yfine[m]);
        }else{
            fprintf( f, "%f %f\n", Xfine[m], Yfine[m]);
        }
    }
    fclose( f );
    //*/

    TFloat64List sortedRedshifts;// = redshifts;
    // for sc_020086471_F02P016_vmM1_red_107_1_atm_clean : zref = 1.3455, aref=4.6772031621836956e-17, zcalc = 2.1634
    Float64 zcenter = 1.3455;
    Float64 acenter = 4.6772031621836956e-17;
    Float64 arange = acenter/2.0/2.0; //
    //Float64 zcenter = 2.1634;
    //Float64 acenter = 3.3228859582551038e-17;
    //Float64 arange = acenter/4.0/2.0; //
    // for
    //Float64 zcenter = 1.134;
    //Float64 acenter = 4.6772031621836956e-17;
    //Float64 arange = acenter/2.0/2.0; //

    // fill the redshift grid
    //sortedRedshifts.push_back(zcenter);
    Float64 zmin = zcenter - 0.025;
    Float64 zmax = zcenter + 0.025;
    Float64 zstep = 0.0001;
    Int32 nz = (Int32)((zmax-zmin)/zstep);
    for (Int32 i=0;i<nz;i++)
    {
        Float64 z = zmin + zstep*i;
        sortedRedshifts.push_back(z);
    }
    std::sort(sortedRedshifts.begin(), sortedRedshifts.end());

    // fill the amplitude grid
    TFloat64List sortedAmplitudes;
    //sortedAmplitudes.push_back(-1); //auto amplitude
    Float64 amin = acenter - arange;
    Float64 amax = acenter + arange;
    Int32 na = 200;
    Float64 astep = ((amax-amin)/na);
    for (Int32 i=0;i<na;i++)
    {
        Float64 a = amin + astep*i;
        sortedAmplitudes.push_back(a);
    }
    std::sort(sortedAmplitudes.begin(), sortedAmplitudes.end());

    CTemplateFittingResult* result = new CTemplateFittingResult();
    result->ChiSquare.resize( sortedRedshifts.size() );
    result->FitAmplitude.resize( sortedRedshifts.size() );
    result->FitAmplitudeError.resize( sortedRedshifts.size() );
    result->FitAmplitudeNegative.resize( sortedRedshifts.size() );
    result->FitDtM.resize( sortedRedshifts.size() );
    result->FitMtM.resize( sortedRedshifts.size() );
    result->FitDustCoeff.resize( sortedRedshifts.size() );
    result->FitMeiksinIdx.resize( sortedRedshifts.size() );
    result->Redshifts.resize( sortedRedshifts.size() );
    result->Overlap.resize( sortedRedshifts.size() );
    result->Status.resize( sortedRedshifts.size() );

    result->Redshifts = sortedRedshifts;

    FILE* fa = fopen( "chi2_versus_ampl_z_axes_ampl.txt", "w+" );
    for (Int32 j=0;j<sortedAmplitudes.size();j++)
    {
        Float64 ampl = sortedAmplitudes[j];
        fprintf( fa, "%.15e\n", ampl);
    }
    fclose( fa );

    FILE* fz = fopen( "chi2_versus_ampl_z_axes_z.txt", "w+" );
    for (Int32 i=0;i<sortedRedshifts.size();i++)
    {
        fprintf( fz, "%.15e\n", sortedRedshifts[i]);
    }
    fclose( fz );

    FILE* f = fopen( "chi2_versus_ampl_z.txt", "w+" );
    for (Int32 i=0;i<sortedRedshifts.size();i++)
    {
        for (Int32 j=0;j<sortedAmplitudes.size();j++)
        {
            Float64 ampl = sortedAmplitudes[j];
            BasicFit( spectrum, tpl, lambdaRange, result->Redshifts[i], overlapThreshold,
                      result->Overlap[i],
                      result->ChiSquare[i],
                      result->FitAmplitude[i],
                      result->FitAmplitudeError[i],
                      result->FitAmplitudeNegative[i],
                      result->FitDtM[i],
                      result->FitMtM[i],
                      result->LogPrior[i],
                      result->FitDustCoeff[i],
                      result->FitMeiksinIdx[i],
                      result->Status[i],
                      result->ChiSquareIntermediate[i],
                      result->IsmDustCoeffIntermediate[i],
                      result->IgmMeiksinIdxIntermediate[i],
                      "lin",
                      ampl );

            fprintf( f, "%.15e", result->ChiSquare[i]);
            if(j<sortedAmplitudes.size()-1){
                fprintf( f, "\t");
            }
        }
        fprintf( f, "\n");
    }
    fclose( f );

    return result;

}

/**
 * \brief this function estimates the likelihood_cstLog term withing the wavelength range
 **/
Float64 COperatorTemplateFitting::EstimateLikelihoodCstLog(const CSpectrum& spectrum, const TFloat64Range& lambdaRange)
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


Int32   COperatorTemplateFitting::ComputeSpectrumModel(const CSpectrum& spectrum,
                                           const CTemplate& tpl,
                                           Float64 redshift,
                                           Float64 DustCoeff,
                                           Int32 meiksinIdx,
                                           Float64 amplitude,
                                           std::string opt_interp,
					                       std::string opt_extinction,
                                           const TFloat64Range& lambdaRange,
                                           Float64 overlapThreshold,
                                           std::shared_ptr<CModelSpectrumResult> & spcPtr)
{
    Log.LogDetail("  Operator-COperatorTemplateFitting: building spectrum model templateFitting for candidate Zcand=%f", redshift);
    EStatus status;
    TFloat64Range currentRange;
    Float64 overlapRate = 0.0;

    Int32 ret = RebinTemplate(  spectrum, 
                                tpl,
                                redshift, 
                                lambdaRange,
                                opt_interp,
                                currentRange,
                                overlapRate,
                                overlapThreshold);
    if( ret == -1 ){
        status = nStatus_NoOverlap; 
        return -1;
    }
    if( ret == -2 ){
        status = nStatus_DataError;
        return -1;
    }
    const TAxisSampleList & Xspc = m_spcSpectralAxis_restframe.GetSamplesVector();
    m_templateRebined_bf.SetIsmIgmLambdaRange(currentRange);

    if ((DustCoeff>0.) || (meiksinIdx>0)){
        m_templateRebined_bf.InitIsmIgmConfig(tpl.m_ismCorrectionCalzetti, tpl.m_igmCorrectionMeiksin);
    }

    if (DustCoeff>0.)
    {
        if (m_templateRebined_bf.CalzettiInitFailed())
        {
            Log.LogError("  Operator-TemplateFitting: asked model with Dust extinction with no calzetti calib. file loaded in template" );
            return -1;
        }
        Int32 idxDust = -1;
        for(Int32 kDust=0; m_templateRebined_bf.m_ismCorrectionCalzetti->GetNPrecomputedDustCoeffs(); kDust++)
        {
            Float64 kDustCoeff= m_templateRebined_bf.m_ismCorrectionCalzetti->GetEbmvValue(kDust);
            if(DustCoeff==kDustCoeff)
            {
                idxDust = kDust;
                break;
            }
        }

        if (idxDust!=-1)
            m_templateRebined_bf.ApplyDustCoeff(idxDust);
    }

    if(opt_extinction == "yes")
    {
        if (m_templateRebined_bf.MeiksinInitFailed())
        {
            Log.LogError("  Operator-TemplateFitting: asked model with IGM extinction with no Meikin calib. file loaded in template" );
            return -1;
        }
        Bool igmCorrectionAppliedOnce = m_templateRebined_bf.ApplyMeiksinCoeff(meiksinIdx, redshift);
    } 
    m_templateRebined_bf.ScaleFluxAxis(amplitude);
    spcPtr = std::make_shared<CModelSpectrumResult>(m_templateRebined_bf);
    return 0;
}



