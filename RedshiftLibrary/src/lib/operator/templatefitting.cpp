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
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/common/flag.h"
#include "RedshiftLibrary/common/formatter.h"

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

#define NOT_OVERLAP_VALUE NAN

using namespace NSEpic;
using namespace std;



/**
 * @brief COperatorTemplateFitting::BasicFit
 * @param tpl
 * @param redshift
 * @param overlapThreshold
 * @param opt_interp
 * @param forcedAmplitude
 * @param opt_extinction
 * @param opt_dustFitting : -1 = disabled, -10 = fit over all available indexes, positive integer 0, 1 or ... will be used as ism-calzetti index as initialized in constructor.
 * @param spcMaskAdditional
 * @param priorjoint_pISM_tpl_z : vector size = nISM, joint prior p(ISM, TPL, Z)
 */
TFittingIsmIgmResult
COperatorTemplateFitting::BasicFit(const std::shared_ptr<const CTemplate>& tpl,
                                   Float64 redshift,
                                   Float64 overlapThreshold,
                                   std::string opt_interp,
                                   Float64 forcedAmplitude,
                                   Int32 opt_extinction,
                                   Int32 opt_dustFitting,
                                   CMask spcMaskAdditional,
                                   const CPriorHelper::TPriorEList & logpriore,
                                   const TInt32List& MeiksinList,
                                   const TInt32List& EbmvList)                                   
{
    bool verbose = false;
    bool amplForcePositive=true;
    bool status_chisquareSetAtLeastOnce = false;

    Int32 EbmvListSize = EbmvList.size();
    Int32 MeiksinListSize = MeiksinList.size();

    TFittingIsmIgmResult result(EbmvListSize, MeiksinListSize);

    if(spcMaskAdditional.GetMasksCount()!= m_spectrum.GetSampleCount())
    {
        Log.LogInfo("COperatorTemplateFitting::BasicFit: spcMaskAdditional does not have the same size as the spectrum flux vector... (%d vs %d), aborting", spcMaskAdditional.GetMasksCount(),  m_spectrum.GetFluxAxis().GetSamplesCount());
        result.status = nStatus_DataError;
        return result ;
    }

    TFloat64Range currentRange; // restframe
    RebinTemplate(  tpl,
                    redshift, 
                    opt_interp,
                    currentRange,
                    result.overlapRate,
                    overlapThreshold);
    
    bool apply_ism = ( (opt_dustFitting==-10 || opt_dustFitting>-1) ? true : false);
    
    bool kStartEnd_ok = currentRange.getClosedIntervalIndices(m_templateRebined_bf.GetSpectralAxis().GetSamplesVector(), m_kStart, m_kEnd);
    if (!kStartEnd_ok)
      throw GlobalException(INTERNAL_ERROR,"COperatorTemplateFitting::BasicFit: impossible to get valid kstart or kend");
    if(m_kStart==-1 || m_kEnd==-1)
      throw GlobalException(INTERNAL_ERROR,Formatter()<<"COperatorTemplateFitting::BasicFit:: kStart="<<m_kStart<<", kEnd="<<m_kEnd);

    if (apply_ism || opt_extinction){          
        InitIsmIgmConfig(redshift, tpl->m_ismCorrectionCalzetti, tpl->m_igmCorrectionMeiksin, EbmvListSize);
    }
    
    Int32 iEbmvCoeffMin = EbmvList[0];
    Int32 iEbmvCoeffMax = EbmvList[EbmvListSize-1];

    m_option_igmFastProcessing = (MeiksinList.size()>1 ? true : false);

    m_forcedAmplitude = forcedAmplitude;

    bool apply_priore = !logpriore.empty() && !tpl->CalzettiInitFailed() 
                        && (logpriore.size()==tpl->m_ismCorrectionCalzetti->GetNPrecomputedEbmvCoeffs());
    
    CPriorHelper::SPriorTZE logpriorTZEempty = {};

    //Loop on the meiksin Idx
    bool igmLoopUseless_WavelengthRange = false;
    for(Int32 kM=0; kM<MeiksinListSize; kM++)
    {
        if(igmLoopUseless_WavelengthRange)
        {
            //Now copy from the already calculated k>0 igm values
            for(Int32 kism=0; kism<result.ChiSquareInterm.size(); kism++)
            {
                for(Int32 kigm=1; kigm<result.ChiSquareInterm[kism].size(); kigm++)
                {
                    result.ChiSquareInterm[kism][kigm] = result.ChiSquareInterm[kism][0];
                    result.IsmCalzettiCoeffInterm[kism][kigm] = result.IsmCalzettiCoeffInterm[kism][0];
                    result.IgmMeiksinIdxInterm[kism][kigm] = result.IgmMeiksinIdxInterm[kism][0];
                }
            }
            break;
        }
        Int32 meiksinIdx = MeiksinList[kM]; //index for the Meiksin curve (0-6; 3 being the median extinction value)

        bool igmCorrectionAppliedOnce = false;
        //Meiksin IGM extinction
        if(opt_extinction)
        {
            //check if lya belongs to current range. 
            if (CheckLyaIsInCurrentRange(currentRange))
                igmCorrectionAppliedOnce = false;
            else
                igmCorrectionAppliedOnce = ApplyMeiksinCoeff(meiksinIdx);
            if(!igmCorrectionAppliedOnce)
                igmLoopUseless_WavelengthRange = true; 
        }
        //Loop on the EBMV dust coeff
        for (Int32 kEbmv_=0; kEbmv_<EbmvListSize; kEbmv_++)
        {
            Int32 kEbmv = EbmvList[kEbmv_];
            Float64 coeffEBMV = -1.; //no ism by default (ie DustCoeff=1.)

            if (apply_ism)
            {
                coeffEBMV = tpl->m_ismCorrectionCalzetti->GetEbmvValue(kEbmv);
                ApplyDustCoeff(kEbmv);
            }

            const CPriorHelper::SPriorTZE & logpriorTZE = apply_priore ?  logpriore[kEbmv] : logpriorTZEempty;
            TFittingResult fitRes =  ComputeLeastSquare(kM, kEbmv_, logpriorTZE, spcMaskAdditional);

            if(verbose)
            {
                Log.LogDebug("");
                Log.LogDebug("  Operator-TemplateFitting: z=%f", redshift);
                Log.LogDebug("  Operator-TemplateFitting: fit=%e", fitRes.chiSquare);
                Log.LogDebug("  Operator-TemplateFitting: sumT=%e", fitRes.sumT);
                Log.LogDebug("  Operator-TemplateFitting: sumS=%e", fitRes.sumS);
                Log.LogDebug("  Operator-TemplateFitting: sumCross=%e", fitRes.sumCross);
                Log.LogDebug("  Operator-TemplateFitting: coeffEBMV=%.2f", coeffEBMV);
                Log.LogDebug("  Operator-TemplateFitting: meiksinIdx=%d", meiksinIdx);

            }

            result.ChiSquareInterm[kEbmv_][kM] = fitRes.chiSquare;
            result.IsmCalzettiCoeffInterm[kEbmv_][kM] = coeffEBMV;
            result.IgmMeiksinIdxInterm[kEbmv_][kM] = meiksinIdx;

            if(fitRes.chiSquare < result.chiSquare)
            {
                TFittingResult & result_base = result;// upcasting (slicing)
                result_base = fitRes; // slicing, preserving specific TFittingIsmIGmResult members

                result.EbmvCoeff = coeffEBMV;
                result.MeiksinIdx = igmCorrectionAppliedOnce==true ? meiksinIdx : -1;
                status_chisquareSetAtLeastOnce = true;
            }
        }
    }

    if(status_chisquareSetAtLeastOnce)
    {
        result.status = nStatus_OK;
    }else{
        result.status = nStatus_LoopError;
    }

    return result;
}

void COperatorTemplateFitting::InitIsmIgmConfig(
        Float64 redshift,
        const std::shared_ptr<CSpectrumFluxCorrectionCalzetti>& ismCorrectionCalzetti,
        const std::shared_ptr<CSpectrumFluxCorrectionMeiksin>& igmCorrectionMeiksin,
        Int32 EbmvListSize)
{
    m_templateRebined_bf.InitIsmIgmConfig(m_kStart, m_kEnd, redshift, ismCorrectionCalzetti, igmCorrectionMeiksin);
    
    m_sumCross_outsideIGM = TFloat64List(EbmvListSize, 0.0);
    m_sumT_outsideIGM = TFloat64List(EbmvListSize, 0.0);
    m_sumS_outsideIGM = TFloat64List(EbmvListSize, 0.0);

}

TFittingResult 
COperatorTemplateFitting::ComputeLeastSquare(   Int32 kM,
                                                Int32 kEbmv_,
                                                const CPriorHelper::SPriorTZE & logpriorTZE,
                                                const CMask & spcMaskAdditional)
{
    TFittingResult fitResult = ComputeCrossProducts(kM, kEbmv_, spcMaskAdditional);

    ComputeAmplitudeAndChi2(fitResult, logpriorTZE);

    return fitResult;
}

TFittingResult 
COperatorTemplateFitting::ComputeCrossProducts( Int32 kM,
                                                Int32 kEbmv_,
                                                const CMask & spcMaskAdditional)
{
    const CSpectrumFluxAxis & spcFluxAxis = m_spectrum.GetFluxAxis();
    const TAxisSampleList & Yspc = spcFluxAxis.GetSamplesVector();
    const TAxisSampleList & Ytpl = m_templateRebined_bf.GetFluxAxis().GetSamplesVector();
    const TAxisSampleList & Xtpl = m_templateRebined_bf.GetSpectralAxis().GetSamplesVector();

    TFittingResult fitResult;
    
    Float64 & sumCross = fitResult.sumCross;
    Float64 & sumT = fitResult.sumT;
    Float64 & sumS = fitResult.sumS;

    Float64 sumCross_IGM = 0.0;
    Float64 sumT_IGM = 0.0;
    Float64 sumS_IGM = 0.0;
    Int32 sumsIgmSaved = 0;

    Float64 err2 = 0.0;
    Int32 numDevs = 0;
    Int32 numDevsFull = 0;
    const CSpectrumNoiseAxis& error = spcFluxAxis.GetError();

    Int32 kIgmEnd = m_option_igmFastProcessing ? m_templateRebined_bf.GetIgmEndIndex() :-1;
    Int32 kEndloop = m_option_igmFastProcessing && kM>0  ? kIgmEnd : m_kEnd;
    for(Int32 j=m_kStart; j<=kEndloop; j++)
    {
        if(m_option_igmFastProcessing && sumsIgmSaved==0 && j>kIgmEnd )
        {
            //store intermediate sums for IGM range
            sumCross_IGM = sumCross;
            sumT_IGM = sumT;
            sumS_IGM = sumS;
            sumsIgmSaved = 1;
        }
        
        numDevsFull++;

        if(spcMaskAdditional[j]){

            numDevs++;
            err2 = 1.0 / (error[j] * error[j]);

            // Tonry&Davis formulation
            sumCross+=Yspc[j]*Ytpl[j]*err2;
            sumT+=Ytpl[j]*Ytpl[j]*err2;
            sumS+= Yspc[j]*Yspc[j]*err2;

            if( std::isinf(err2) || std::isnan(err2) ){
                throw GlobalException(INTERNAL_ERROR, Formatter() << "COperatorTemplateFitting::ComputeCrossProducts: found invalid inverse variance : err2="<<err2<<", for index="<< j<<" at restframe wl="<<Xtpl[j]);
            }

            if( std::isinf(sumS) || std::isnan(sumS) || sumS!=sumS ){
                Log.LogError("COperatorTemplateFitting::ComputeCrossProducts: found invalid dtd : dtd=%e, for index=%d at restframe wl=%f", sumS, j, Xtpl[j]);
                Log.LogError("COperatorTemplateFitting::ComputeCrossProducts: found invalid dtd : Yspc=%e, for index=%d at restframe wl=%f", Yspc[j], j, Xtpl[j]);
                Log.LogError("COperatorTemplateFitting::ComputeCrossProducts: found invalid dtd : err2=%e, for index=%d at restframe wl=%f", err2, j, Xtpl[j]);
                Log.LogError("COperatorTemplateFitting::ComputeCrossProducts: found invalid dtd : error=%e, for index=%d at restframe wl=%f", error[j], j, Xtpl[j]);
                throw GlobalException(INTERNAL_ERROR, Formatter() << "COperatorTemplateFitting::ComputeCrossProducts: found invalid dtd");
            }

            if( std::isinf(sumT) || std::isnan(sumT) ){
                throw GlobalException(INTERNAL_ERROR, Formatter() << "COperatorTemplateFitting::ComputeCrossProducts: found invalid mtm : mtm=" << sumT
                    << " for index=" << j << " at restframe wl=" << Xtpl[j]);
            }
        }
    }

    if(m_option_igmFastProcessing)
    {
        if (kM==0){ 
            m_sumCross_outsideIGM[kEbmv_] = sumCross-sumCross_IGM;
            m_sumT_outsideIGM[kEbmv_] = sumT-sumT_IGM;
            m_sumS_outsideIGM[kEbmv_] = sumS-sumS_IGM;
        } else {
            sumCross += m_sumCross_outsideIGM[kEbmv_];
            sumT += m_sumT_outsideIGM[kEbmv_];
            sumS += m_sumS_outsideIGM[kEbmv_];
        }

    }

    if ( numDevs==0 )
    {
        throw GlobalException(INTERNAL_ERROR,  "COperatorTemplateFitting::ComputeCrossProducts: empty leastsquare sum");
    }

    return fitResult;
}

void 
COperatorTemplateFitting::ComputeAmplitudeAndChi2( TFittingResult & fitResult,
                                            const CPriorHelper::SPriorTZE & logpriorTZE) const
{ 
    const Float64 & sumCross = fitResult.sumCross;
    const Float64 & sumT = fitResult.sumT;
    const Float64 & sumS = fitResult.sumS;

    Float64 & ampl = fitResult.ampl;
    Float64 & ampl_err = fitResult.ampl_err;
    Float64 & ampl_sigma = fitResult.ampl_sigma;
    Float64 & fit = fitResult.chiSquare;

    bool apply_priore = (logpriorTZE.A_sigma>0.0 && logpriorTZE.betaA>0.0) ? true : false;

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
            Float64 bss2 = logpriorTZE.betaA/(logpriorTZE.A_sigma*logpriorTZE.A_sigma);
            ampl = (sumCross+logpriorTZE.A_mean*bss2)/(sumT+bss2);
            ampl_err = sqrt(sumT)/(sumT+bss2); 
        }

        if(m_forcedAmplitude !=-1){
            ampl = m_forcedAmplitude;
            ampl_err = 0.;
        }

        ampl_sigma = ampl/ampl_err;

        if(m_amplForcePositive)
        {
            ampl = max(0.0, ampl);
        }

        //Generalized method (ampl can take any value now) for chi2 estimate
        fit = sumS + sumT*ampl*ampl - 2.*ampl*sumCross;
    }

    Float64 & logprior = fitResult.logprior;
    if (apply_priore)
    {
        logprior += -2.*logpriorTZE.betaTE*logpriorTZE.logprior_precompTE;
        logprior += -2.*logpriorTZE.betaA*logpriorTZE.logprior_precompA;
        logprior += -2.*logpriorTZE.betaZ*logpriorTZE.logprior_precompZ;
        if(logpriorTZE.A_sigma>0.0)
        {
            Float64 logPa = logpriorTZE.betaA*(ampl-logpriorTZE.A_mean)*(ampl-logpriorTZE.A_mean)/(logpriorTZE.A_sigma*logpriorTZE.A_sigma);
            logprior += logPa;
        }
        fit += logprior;
    }
}

/**
 * \brief
 *
 * input: if additional_spcMasks size is 0, no additional mask will be used, otherwise its size should match the redshifts list size
 *
 * opt_dustFitting: -1 = disabled, -10 = fit over all available indexes, positive integer 0, 1 or ... will be used as ism-calzetti index as initialized in constructor.
 *
 **/
std::shared_ptr<COperatorResult> COperatorTemplateFitting::Compute(const std::shared_ptr<const CTemplate> & tpl,
                                                              Float64 overlapThreshold,
                                                              const std::vector<CMask> &additional_spcMasks,
                                                              std::string opt_interp,
                                                              Int32 opt_extinction,
                                                              Int32 opt_dustFitting,
                                                              const CPriorHelper::TPriorZEList &logpriorze,
                                                              bool keepigmism,
                                                              Float64 FitEbmvCoeff,
                                                              Int32 FitMeiksinIdx)
{
    Log.LogDetail("  Operator-TemplateFitting: starting computation for template: %s", tpl->GetName().c_str());

    if( (opt_dustFitting==-10 || opt_dustFitting>-1) && tpl->CalzettiInitFailed())
    {
        throw GlobalException(INTERNAL_ERROR,"  Operator-TemplateFitting: no calzetti calib. file in template");
    }
    if( opt_dustFitting>-1 && opt_dustFitting>tpl->m_ismCorrectionCalzetti->GetNPrecomputedEbmvCoeffs()-1)
    {
      throw GlobalException(INTERNAL_ERROR,Formatter()<<"Operator-TemplateFitting: calzetti index overflow (dustfitting="<<opt_dustFitting<<",while NPrecomputedEbmvCoeffs="<<tpl->m_ismCorrectionCalzetti->GetNPrecomputedEbmvCoeffs()<<")");
    }

    if( opt_extinction && tpl->MeiksinInitFailed())
    {
        throw GlobalException(INTERNAL_ERROR,"  Operator-TemplateFitting: no meiksin calib. file in template");
    }

    if( m_spectrum.GetSpectralAxis().IsInLinearScale() == false || tpl->GetSpectralAxis().IsInLinearScale() == false )
    {
        Log.LogError("  Operator-TemplateFitting: input spectrum or template are not in linear scale (ignored)");
    }

    //sort the redshift and keep track of the indexes
    TFloat64List sortedRedshifts;
    TFloat64List sortedIndexes;
    // This is a vector of {value,index} pairs
    vector<pair<Float64,Int32> > vp;
    vp.reserve(m_redshifts.size());
    for (Int32 i = 0 ; i < m_redshifts.size() ; i++) {
        vp.push_back(make_pair(m_redshifts[i], i));
    }
    std::sort(vp.begin(), vp.end());
    for (Int32 i = 0 ; i < vp.size() ; i++) {
        sortedRedshifts.push_back(vp[i].first);
        sortedIndexes.push_back(vp[i].second);
    }

    std::shared_ptr<CTemplateFittingResult> result = std::make_shared<CTemplateFittingResult>();
    TInt32List MeiksinList;
    TInt32List EbmvList;
    tpl->GetIsmIgmIdxList( opt_extinction, opt_dustFitting, MeiksinList, EbmvList, keepigmism, FitEbmvCoeff, FitMeiksinIdx);
    Int32 MeiksinListSize = MeiksinList.size();
    Int32 EbmvListSize = EbmvList.size();

    result->Init(sortedRedshifts.size());
    result->Redshifts = sortedRedshifts;

    //default mask
    bool useDefaultMask = additional_spcMasks.size()!=sortedRedshifts.size() ? true : false;
    CMask default_spcMask(m_spectrum.GetSampleCount());
    if (useDefaultMask)
        for(Int32 km=0; km<default_spcMask.GetMasksCount(); km++)
            default_spcMask[km] = 1.0;

    if(additional_spcMasks.size()!=sortedRedshifts.size() && additional_spcMasks.size()!=0)
        throw GlobalException(INTERNAL_ERROR,Formatter()<<"Operator-TemplateFitting: masks-list size ("<<additional_spcMasks.size()<<") didn't match the input redshift-list ("<<sortedRedshifts.size()<<") !)");
    
    if(logpriorze.size()>0 && logpriorze.size()!=sortedRedshifts.size())
        throw GlobalException(INTERNAL_ERROR,Formatter()<<"Operator-TemplateFitting: prior list size("<<logpriorze.size()<<") didn't match the input redshift-list size :"<< sortedRedshifts.size());

    for (Int32 i=0;i<sortedRedshifts.size();i++)
    {
        const CPriorHelper::TPriorEList &logp = logpriorze.size()>0 && logpriorze.size()==sortedRedshifts.size() ? logpriorze[i] : CPriorHelper::TPriorEList();

        Float64 redshift = result->Redshifts[i];

        const CMask &additional_spcMask = useDefaultMask ? default_spcMask : additional_spcMasks[sortedIndexes[i]];

        TFittingIsmIgmResult result_z = BasicFit(   tpl,
                                                    redshift,
                                                    overlapThreshold,
                                                    opt_interp,
                                                    -1,
                                                    opt_extinction,
                                                    opt_dustFitting,
                                                    additional_spcMask,
                                                    logp,
                                                    MeiksinList,
                                                    EbmvList);
        
        result->set_at_redshift(i, std::move(result_z));

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
        Log.LogInfo("  Operator-TemplateFitting: overlap warning for %s: minz=%.3f, maxz=%.3f", tpl->GetName().c_str(), overlapValidInfZ, overlapValidSupZ);
    }

    //only bad status warning
    Int32 oneValidStatusFoundIndex = -1;
    for (Int32 i=0;i<sortedRedshifts.size();i++)
    {
        if(result->Status[i]==nStatus_OK)
        {
            oneValidStatusFoundIndex=i;
            Log.LogDebug("  Operator-TemplateFitting: STATUS VALID found for %s: at least at index=%d", tpl->GetName().c_str(), i);
            break;
        }
    }if(oneValidStatusFoundIndex==-1)
    {
        Flag.warning(Flag.INVALID_MERIT_VALUES, Formatter() << "  COperatorTemplateFitting::"<<__func__<<": STATUS WARNING for "<<tpl->GetName().c_str()<<": Not even one single valid fit/merit value found");
    }


    //loop error status warning
    Int32 loopErrorStatusFoundIndex = -1;
    for (Int32 i=0;i<sortedRedshifts.size();i++)
    {
        if(result->Status[i]==nStatus_LoopError)
        {
            loopErrorStatusFoundIndex=i;
            Log.LogDebug("  Operator-TemplateFitting: STATUS Loop Error found for %s: at least at index=%d", tpl->GetName().c_str(), i);
            break;
        }
    }if(loopErrorStatusFoundIndex!=-1)
    {
        Flag.warning(Flag.INVALID_MERIT_VALUES, Formatter() <<"    COperatorTemplateFitting::"<<__func__<<": Loop Error - chisquare values not set even once");
    }

    //estimate CstLog for PDF estimation
    result->CstLog = EstimateLikelihoodCstLog(m_spectrum, m_lambdaRange);

    return result;

}
