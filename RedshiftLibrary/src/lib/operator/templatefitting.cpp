// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#include "RedshiftLibrary/operator/templatefitting.h"

#include "RedshiftLibrary/spectrum/axis.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "RedshiftLibrary/spectrum/tools.h"
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/operator/templatefittingresult.h"
#include "RedshiftLibrary/extremum/extremum.h"
#include "RedshiftLibrary/common/quicksort.h"
#include "RedshiftLibrary/log/log.h"

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
#include "RedshiftLibrary/processflow/datastore.h"
#define NOT_OVERLAP_VALUE NAN
#include <stdio.h>
#include <numeric>
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
 * @param fittingAmplitudeSigma
 * @param fittingDtM
 * @param fittingMtM
 * @param fittingEbmvCoeff
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
                                   Float64& fittingAmplitudeSigma,
                                   Float64& fittingDtM,
                                   Float64& fittingMtM,
                                   Float64& fittingLogprior,
                                   Float64 &fittingEbmvCoeff,
                                   Int32 &fittingMeiksinIdx,
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
                                   const TInt32List& MeiksinList,
                                   const TInt32List& EbmvList)
{
    bool verbose = false;
    bool amplForcePositive=true;
    chiSquare = INFINITY;
    bool status_chisquareSetAtLeastOnce = false;

    fittingAmplitude = -1.0;
    fittingAmplitudeError = -1.0;
    fittingAmplitudeSigma = 0.;
    overlapRate = 0.0;
    status = nStatus_DataError;

    const CSpectrumFluxAxis& spcFluxAxis = spectrum.GetFluxAxis();
    const TAxisSampleList & Yspc = spcFluxAxis.GetSamplesVector();

    if(spcMaskAdditional.GetMasksCount()!=spcFluxAxis.GetSamplesCount())
    {
        Log.LogInfo("  Operator-TemplateFitting: spcMaskAdditional does not have the same size as the spectrum flux vector... (%d vs %d), aborting", spcMaskAdditional.GetMasksCount(),  spectrum.GetFluxAxis().GetSamplesCount());
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
    
    Int32 kStart = -1, kEnd = -1, kIgmEnd = -1;
    bool kStartEnd_ok = currentRange.getClosedIntervalIndices(m_templateRebined_bf.GetSpectralAxis().GetSamplesVector(), kStart, kEnd);
    if (!kStartEnd_ok){
        Log.LogError("COperatorTemplateFitting::BasicFit: impossible to get valid kstart or kend");
        throw std::runtime_error("COperatorTemplateFitting::BasicFit: impossible to get valid kstart or kend");
    }
    if (apply_ism || opt_extinction){          
        m_templateRebined_bf.InitIsmIgmConfig(kStart, kEnd, redshift, tpl.m_ismCorrectionCalzetti, tpl.m_igmCorrectionMeiksin);
    }
    if (opt_extinction)
        kIgmEnd = m_templateRebined_bf.GetIgmEndIndex();
        
    if( ret == -1 ){
        status = nStatus_NoOverlap; 
        return;
    }
    if( ret == -2 ){
        status = nStatus_DataError;
        return;
    }

    Int32 EbmvListSize = ChiSquareInterm.size();
    Int32 MeiksinListSize = ChiSquareInterm[0].size();

    Int32 iEbmvCoeffMin = EbmvList[0];
    Int32 iEbmvCoeffMax = EbmvList[EbmvListSize-1];

    Bool option_igmFastProcessing = (MeiksinList.size()==1 ? false : true);
    TFloat64List  sumCross_outsideIGM(EbmvListSize, 0.0);
    TFloat64List  sumT_outsideIGM(EbmvListSize, 0.0);
    TFloat64List  sumS_outsideIGM(EbmvListSize, 0.0);
    //Loop on the meiksin Idx
    Bool igmLoopUseless_WavelengthRange = false;
    for(Int32 kM=0; kM<MeiksinListSize; kM++)
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

        if(kStart==-1 || kEnd==-1)
        {
            Log.LogDebug( "  Operator-TemplateFitting: kStart=%d, kEnd=%d ! Aborting.", kStart, kEnd);
            break;
        }

        bool igmCorrectionAppliedOnce = false;
        //Meiksin IGM extinction
        if(opt_extinction)
        {
            //check if lya belongs to current range. If not, no
            //zigm_min = GetIGMStartingRedshiftValue(currentRange.GetBegin())
            if(currentRange.GetBegin()/(1+redshift) > 1216)
                igmCorrectionAppliedOnce = false;
            else
                igmCorrectionAppliedOnce = m_templateRebined_bf.ApplyMeiksinCoeff(meiksinIdx);
            if(!igmCorrectionAppliedOnce)
                igmLoopUseless_WavelengthRange = true; 
        }
        //Loop on the EBMV dust coeff
        for(Int32 kEbmv=iEbmvCoeffMin; kEbmv<=iEbmvCoeffMax; kEbmv++)
        {
            Int32 kEbmv_ = kEbmv - iEbmvCoeffMin; //index used to fill some arrays
            Float64 coeffEBMV = -1.; //no ism by default (ie DustCoeff=1.)

            if (apply_ism)
            {
                coeffEBMV = m_templateRebined_bf.m_ismCorrectionCalzetti->GetEbmvValue(kEbmv);
                m_templateRebined_bf.ApplyDustCoeff(kEbmv);
            }
 
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
            const CSpectrumNoiseAxis& error = spcFluxAxis.GetError();

            Int32 kEndloop = kEnd;
            if (option_igmFastProcessing && kM>0) kEndloop = kIgmEnd;
            for(Int32 j=kStart; j<=kEndloop; j++)
            {
                if(option_igmFastProcessing && sumsIgmSaved==0 && j>kIgmEnd )
                {
                    //store intermediate sums for IGM range
                    sumCross_IGM = sumCross;
                    sumT_IGM = sumT;
                    sumS_IGM = sumS;
                    sumsIgmSaved = 1;
                }
                
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
                }
            }
            if(option_igmFastProcessing && kM==0)
            {
                sumCross_outsideIGM[kEbmv_] = sumCross-sumCross_IGM;
                sumT_outsideIGM[kEbmv_] = sumT-sumT_IGM;
                sumS_outsideIGM[kEbmv_] = sumS-sumS_IGM;
            }
            if(option_igmFastProcessing && kM>0)
            {
                sumCross += sumCross_outsideIGM[kEbmv_];
                sumT += sumT_outsideIGM[kEbmv_];
                sumS += sumS_outsideIGM[kEbmv_];
            }

            if ( numDevs==0 )
            {
                status = nStatus_DataError;
                return;
            }

            Float64 ampl = 0.0;
            Float64 ampl_err = 0.0;
            Float64 ampl_sigma = 0.0;
            bool apply_priore = false;
            if(!logpriore.empty() && !m_templateRebined_bf.CalzettiInitFailed())
            { 
                if(logpriore.size()==m_templateRebined_bf.m_ismCorrectionCalzetti->GetNPrecomputedEbmvCoeffs())
                {
                    if(logpriore[kEbmv].A_sigma>0.0 && logpriore[kEbmv].betaA>0.0)
                        apply_priore = true;
                }
            }

            if( sumT==0 )
            {
                ampl = 0.0;
                ampl_err = 0.0;
                fit = sumS;
                ampl_sigma = 0.0;
                //status = nStatus_DataError;
                //return;
            }else{
                
                ampl = sumCross/sumT;
                ampl_err = sqrt(1./sumT);
                
                if (apply_priore)
                {
                    Float64 bss2 = logpriore[kEbmv].betaA/(logpriore[kEbmv].A_sigma*logpriore[kEbmv].A_sigma);
                    ampl = (sumCross+logpriore[kEbmv].A_mean*bss2)/(sumT+bss2);
                    ampl_err = sqrt(sumT)/(sumT+bss2); 
                    /*
                    Log.LogInfo("  Operator-TemplateFitting: Constrained amplitude (betaA=%e):  s2b=%e, mtm=%e", logpriore[kEbmv].betaA, s2b, sumT);
                    //*/
                }

                if(forcedAmplitude !=-1){
                    ampl = forcedAmplitude;
                    ampl_err = 0.;
                }

                ampl_sigma = ampl/ampl_err;

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
                logprior += -2.*logpriore[kEbmv].betaTE*logpriore[kEbmv].logprior_precompTE;
                logprior += -2.*logpriore[kEbmv].betaA*logpriore[kEbmv].logprior_precompA;
                logprior += -2.*logpriore[kEbmv].betaZ*logpriore[kEbmv].logprior_precompZ;
                if(logpriore[kEbmv].A_sigma>0.0)
                {
                    Float64 logPa = logpriore[kEbmv].betaA*(ampl-logpriore[kEbmv].A_mean)*(ampl-logpriore[kEbmv].A_mean)/(logpriore[kEbmv].A_sigma*logpriore[kEbmv].A_sigma);
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
            ChiSquareInterm[kEbmv_][kM] = fit;
            IsmCalzettiCoeffInterm[kEbmv_][kM] = coeffEBMV;
            IgmMeiksinIdxInterm[kEbmv_][kM] = meiksinIdx;

            if(fit<chiSquare)
            {
                chiSquare = fit;
                fittingEbmvCoeff = coeffEBMV;
                fittingMeiksinIdx = igmCorrectionAppliedOnce==true?meiksinIdx:-1;
                fittingAmplitude = ampl;
                fittingAmplitudeError = ampl_err;
                fittingAmplitudeSigma = ampl_sigma;
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

/**
 * \brief
 *
 * input: if additional_spcMasks size is 0, no additional mask will be used, otherwise its size should match the redshifts list size
 *
 * opt_dustFitting: -1 = disabled, -10 = fit over all available indexes, positive integer 0, 1 or ... will be used as ism-calzetti index as initialized in constructor.
 * @lambdaRange is not clamped
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
                                                              Float64 FitEbmvCoeff,
                                                              Int32 FitMeiksinIdx)
{
    Log.LogDetail("  Operator-TemplateFitting: starting computation for template: %s", tpl.GetName().c_str());

    if(0)
    {
        //CSpectrumFluxAxis tmp_tplFluxAxis = tpl.GetFluxAxis();
        const CSpectrumSpectralAxis & tmp_tplSpectralAxis = tpl.GetSpectralAxis();
        for(UInt32 k=0; k<std::min(Int32(tmp_tplSpectralAxis.GetSamplesCount()), Int32(10)); k++)
        {
            Log.LogDebug("  Operator-TemplateFitting: tpl_SpectralAxis[%d] = %f",
                         k,
                         tmp_tplSpectralAxis[k]);
        }
    }

    if( (opt_dustFitting==-10 || opt_dustFitting>-1) && tpl.CalzettiInitFailed())
    {
        Log.LogError("  Operator-TemplateFitting: no calzetti calib. file loaded in template... aborting");
        throw std::runtime_error("  Operator-TemplateFitting: no calzetti calib. file in template");
    }
    if( opt_dustFitting>-1 && opt_dustFitting>tpl.m_ismCorrectionCalzetti->GetNPrecomputedEbmvCoeffs()-1)
    {
        Log.LogError("  Operator-TemplateFitting: calzetti index overflow (opt=%d, while NPrecomputedEbmvCoeffs=%d)... aborting",
                     opt_dustFitting,
                     tpl.m_ismCorrectionCalzetti->GetNPrecomputedEbmvCoeffs());
        throw std::runtime_error("  Operator-TemplateFitting: calzetti index overflow");
    }

    if( opt_extinction && tpl.MeiksinInitFailed())
    {
        Log.LogError("  Operator-TemplateFitting: no meiksin calib. file loaded in template... aborting");
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
    TInt32List MeiksinList;
    TInt32List EbmvList;
    tpl.GetIsmIgmIdxList( opt_extinction, opt_dustFitting, MeiksinList, EbmvList, keepigmism, FitEbmvCoeff, FitMeiksinIdx);
    Int32 MeiksinListSize = MeiksinList.size();
    Int32 EbmvListSize = EbmvList.size();

    result->Init(sortedRedshifts.size(), EbmvListSize, MeiksinListSize);
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
        CPriorHelper::TPriorEList logp;
        if(logpriorze.size()>0 && logpriorze.size()==sortedRedshifts.size())
        {
            logp = logpriorze[i];
        }

        Float64 redshift = result->Redshifts[i];

        result->ChiSquareIntermediate[i] = std::vector<TFloat64List>(EbmvListSize, TFloat64List(MeiksinListSize, DBL_MAX));
        result->IsmEbmvCoeffIntermediate[i] = std::vector<TFloat64List>(EbmvListSize, TFloat64List(MeiksinListSize, NAN));
        result->IgmMeiksinIdxIntermediate[i] = std::vector<TInt32List>(EbmvListSize, TInt32List(MeiksinListSize, -1));

        BasicFit( spectrum,
                  tpl,
                  clampedlambdaRange,
                  redshift,
                  overlapThreshold,
                  result->Overlap[i],
                  result->ChiSquare[i],
                  result->FitAmplitude[i],
                  result->FitAmplitudeError[i],
                  result->FitAmplitudeSigma[i],
                  result->FitDtM[i],
                  result->FitMtM[i],
                  result->LogPrior[i],
                  result->FitEbmvCoeff[i],
                  result->FitMeiksinIdx[i],
                  result->Status[i],
                  result->ChiSquareIntermediate[i],
                  result->IsmEbmvCoeffIntermediate[i],
                  result->IgmMeiksinIdxIntermediate[i],
                  opt_interp,
                  -1,
                  opt_extinction,
                  opt_dustFitting,
                  additional_spcMask,
                  logp,
                  MeiksinList,
                  EbmvList);

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
        Log.LogWarning("    Operator-TemplateFitting: Loop Error - chisquare values not set even once");
    }

    //estimate CstLog for PDF estimation
    result->CstLog = EstimateLikelihoodCstLog(spectrum, clampedlambdaRange);

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
    result->FitAmplitudeSigma.resize( sortedRedshifts.size() );
    result->FitDtM.resize( sortedRedshifts.size() );
    result->FitMtM.resize( sortedRedshifts.size() );
    result->FitEbmvCoeff.resize( sortedRedshifts.size() );
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
                      result->FitAmplitudeSigma[i],
                      result->FitDtM[i],
                      result->FitMtM[i],
                      result->LogPrior[i],
                      result->FitEbmvCoeff[i],
                      result->FitMeiksinIdx[i],
                      result->Status[i],
                      result->ChiSquareIntermediate[i],
                      result->IsmEbmvCoeffIntermediate[i],
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



