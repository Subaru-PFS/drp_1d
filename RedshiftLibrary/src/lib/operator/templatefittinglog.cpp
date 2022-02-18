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
#include "RedshiftLibrary/operator/templatefittinglog.h"

#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/common/indexing.h"
#include "RedshiftLibrary/extremum/extremum.h"
#include "RedshiftLibrary/operator/templatefittingresult.h"
#include "RedshiftLibrary/spectrum/axis.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "RedshiftLibrary/spectrum/tools.h"

#include "RedshiftLibrary/log/log.h"

#include <boost/numeric/conversion/bounds.hpp>

#include <algorithm> // std::sort
#include <float.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <math.h>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <sstream>
#include <numeric>
#include <assert.h>

#define NOT_OVERLAP_VALUE NAN
#include <stdio.h>

namespace bfs = boost::filesystem;
using namespace NSEpic;
using namespace std;

COperatorTemplateFittingLog::COperatorTemplateFittingLog(const CSpectrum & spectrum,
                                                         const CSpectrum & logSampledSpectrum,
                                                         const TFloat64Range& lambdaRange, 
                                                         const TFloat64List& redshifts):
    COperatorTemplateFittingBase(spectrum, lambdaRange, redshifts),
    m_logSampledSpectrum(logSampledSpectrum)
{
    if (!m_logSampledSpectrum.GetSpectralAxis().IsLogSampled())
    {
        throw GlobalException(INTERNAL_ERROR,"  Operator-TemplateFittingLog: input spectrum is not log sampled");
    }

    m_logSampledSpectrum.GetSpectralAxis().ClampLambdaRange(lambdaRange, m_logSampledLambdaRange );

    CheckRedshifts();
};

void COperatorTemplateFittingLog::SetRedshifts(TFloat64List redshifts)
{
    COperatorTemplateFittingBase::SetRedshifts(std::move(redshifts));
    CheckRedshifts();
}

void COperatorTemplateFittingLog::CheckRedshifts()
{
    if(m_redshifts.size()<2)
    {
      throw GlobalException(INTERNAL_ERROR,Formatter()<<"Operator-TemplateFittingLog: Cannot compute on a redshift array " <<m_redshifts.size()<<" <2" );
    }

    //check if spectrum sampling is a multiple of redshift step
    m_logstep = log((m_redshifts[1]+1)/(m_redshifts[0]+1));
    Float64 modulo;
    m_ssRatio = m_logSampledSpectrum.GetSpectralAxis().GetLogSamplingIntegerRatio(m_logstep, modulo);
    if(std::abs(modulo)>1E-12)
        throw GlobalException(INTERNAL_ERROR,"TFLogOperator: spc and tpl do not have a lambdastep multiple of redshift step");

    //subsample spectrum if necessary (coarse redshift grid)
    if(m_ssRatio==1){//no required subsampling
        m_ssSpectrum = m_logSampledSpectrum;
    }else{
        TFloat64List mask_spc = m_logSampledSpectrum.GetSpectralAxis().GetSubSamplingMask(m_ssRatio, m_logSampledLambdaRange);
        m_ssSpectrum = CSpectrum(m_logSampledSpectrum, mask_spc);
        // scale the variance by ssratio
        CSpectrumNoiseAxis scaledNoise = m_ssSpectrum.GetErrorAxis();
        scaledNoise  *= 1./sqrt(m_ssRatio);
        m_ssSpectrum.SetErrorAxis(std::move(scaledNoise));

        //double make sure that subsampled spectrum are is well sampled
        if (!m_ssSpectrum.GetSpectralAxis().IsLogSampled(m_logstep))
        {
            throw GlobalException(INTERNAL_ERROR,"  Operator-TemplateFittingLog: subsampled spectrum is not log sampled the redshift step");
        }
    }

    // set the lambda range to the clamped subsampled spectrum
    m_ssSpectrum.GetSpectralAxis().ClampLambdaRange(m_logSampledLambdaRange, m_ssLambdaRange );
}

COperatorTemplateFittingLog::~COperatorTemplateFittingLog()
{
    freeFFTPlans();
}

Int32 COperatorTemplateFittingLog::EstimateXtYSlow(const TFloat64List& X,
                                                   const TFloat64List& Y,
                                                   UInt32 nShifts,
                                                   TFloat64List &XtY)
{
    XtY.resize(nShifts);

    UInt32 nX = X.size();
    Float64 xty = 0.0;
    for (Int32 k = 0; k < nShifts; k++)
    {
        xty = 0.0;
        for (Int32 j = 0; j < nX; j++)
        {
            xty += X[j] * Y[j + k];
        }
        XtY[k] = xty;
    }
    return 0;
}

// only works for mtm, Y=model^2, X=1.
Int32 COperatorTemplateFittingLog::EstimateMtMFast(const TFloat64List &X,
                                                   const TFloat64List &Y,
                                                   UInt32 nShifts,
                                                   TFloat64List &XtY)
{
    XtY.resize(nShifts);

    UInt32 nX = X.size();
    Float64 xty = 0.0;
    for (Int32 j = 0; j < nX; j++)
    {
        xty += X[j] * Y[j];
    }
    XtY[0] = xty;

    for (Int32 k = 1; k < nShifts; k++)
    {
        xty = XtY[k - 1];
        xty -= X[0] * Y[0 + k - 1];
        xty += X[nX - 1] * Y[nX - 1 + k - 1];

        XtY[k] = xty;
    }
    return 0;
}

Int32 COperatorTemplateFittingLog::EstimateXtY(const TFloat64List &X,
                                               const TFloat64List &Y,
                                               UInt32 nshifts,
                                               TFloat64List &XtY,
                                               Int32 precomputedFFT)
{
    // Processing the FFT
    Int32 nSpc = X.size();
    Int32 nTpl = Y.size();
    Int32 nPadded = m_nPaddedSamples;

    Int32 nPadBeforeSpc = nPadded - nSpc; //(Int32)nPadded/2.0;
    Int32 nPadBeforeTpl = 0;

    if (verboseLogXtYFFT)
    {
        Log.LogDebug("  Operator-TemplateFittingLog: FitAllz: Processing spc-fft with n=%d, padded to n=%d", nSpc, nPadded);
    }

    //    for(Int32 k=0; k<nSpc; k++)
    //    {
    //        inSpc[k][0] = X[k];
    //        inSpc[k][1] = 0.0;
    //    }
    //    for(Int32 k=nSpc; k<nPadded; k++)
    //    {
    //        inSpc[k][0] = 0.0;
    //        inSpc[k][1] = 0.0;
    //    }

    bool computeSpcFFT = true;
    if (precomputedFFT == 0 && precomputedFFT_spcFluxOverErr2 != 0)
    {
        computeSpcFFT = false;
    }
    if (precomputedFFT == 1 && precomputedFFT_spcOneOverErr2 != 0)
    {
        computeSpcFFT = false;
    }

    if (computeSpcFFT)
    {

        if (verboseLogXtYFFT)
        {
            Log.LogDebug("  Operator-TemplateFittingLog: FitAllz: Processing spc-fft with nPadBeforeSpc=%d", nPadBeforeSpc);
        }
        for (Int32 k = 0; k < nPadBeforeSpc; k++)
        {
            inSpc[k] = 0.0;
        }
        for (Int32 k = nPadBeforeSpc; k < nPadBeforeSpc + nSpc; k++)
        {
            inSpc[k] = X[nSpc - 1 - (k - nPadBeforeSpc)]; // X[k-nPadBeforeSpc];
        }
        for (Int32 k = nPadBeforeSpc + nSpc; k < nPadded; k++)
        {
            inSpc[k] = 0.0;
        }
        if (verboseExportXtYFFT)
        {
            // save spc-input data
            FILE *f_fftinput = fopen("loglbda_fitallz_xfftInput_dbg.txt", "w+");
            for (Int32 t = 0; t < nPadded; t++)
            {
                fprintf(f_fftinput, "%f\t%e\n", (Float64)t, inSpc[t]);
            }
            fclose(f_fftinput);
        }
 
        fftw_execute(pSpc);
        if (verboseLogXtYFFT)
        {
            Log.LogDebug("  Operator-TemplateFittingLog: FitAllz: spc-fft done");
        }
        if (verboseExportXtYFFT)
        {
            // save spc-fft data
            FILE *f_fftoutput = fopen("loglbda_fitallz_xfftOutput_dbg.txt", "w+");
            for (Int32 t = 0; t < nPadded; t++)
            {
                fprintf(f_fftoutput, "%f\t%e\n", (Float64)t, outSpc[t][0] * outSpc[t][0] + outSpc[t][1] * outSpc[t][1]);
            }
            fclose(f_fftoutput);
        }

        // save computed FFT into precomputed buffer
        if (precomputedFFT == 0)
        {
            precomputedFFT_spcFluxOverErr2 =
                (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nPadded);
            if (precomputedFFT_spcFluxOverErr2 == 0)
            {
                Log.LogError("  Operator-TemplateFittingLog: InitFFT: Unable to allocate precomputedFFT_spcFluxOverErr2");
                return -1;
            }
            for (Int32 k = 0; k < nPadded; k++)
            {
                precomputedFFT_spcFluxOverErr2[k][0] = outSpc[k][0];
                precomputedFFT_spcFluxOverErr2[k][1] = outSpc[k][1];
            }
        }
        if (precomputedFFT == 1)
        {

            precomputedFFT_spcOneOverErr2 =
                (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nPadded);
            if (precomputedFFT_spcOneOverErr2 == 0)
            {
                Log.LogError("  Operator-TemplateFittingLog: InitFFT: Unable to allocate precomputedFFT_spcOneOverErr2");
                return -1;
            }
            for (Int32 k = 0; k < nPadded; k++)
            {
                precomputedFFT_spcOneOverErr2[k][0] = outSpc[k][0];
                precomputedFFT_spcOneOverErr2[k][1] = outSpc[k][1];
            }
        }

    } else
    {
        if (precomputedFFT == 0)
        {
            for (Int32 k = 0; k < nPadded; k++)
            {
                outSpc[k][0] = precomputedFFT_spcFluxOverErr2[k][0];
                outSpc[k][1] = precomputedFFT_spcFluxOverErr2[k][1];
            }
        }
        if (precomputedFFT == 1)
        {
            for (Int32 k = 0; k < nPadded; k++)
            {
                outSpc[k][0] = precomputedFFT_spcOneOverErr2[k][0];
                outSpc[k][1] = precomputedFFT_spcOneOverErr2[k][1];
            }
        }
    }

    if (verboseLogXtYFFT)
    {
        Log.LogDebug("  Operator-TemplateFittingLog: FitAllz: Processing tpl-fft with n=%d, padded to n=%d", nTpl, nPadded);
    }
    //    for(Int32 k=0; k<nTpl; k++)
    //    {
    //        inTpl_padded[k][0] = Y[k];//Y[nTpl-1-k];//
    //        //inTpl_padded[k][1] = 0.0;
    //    }
    //    for(Int32 k=nTpl; k<nPadded; k++)
    //    {
    //        inTpl_padded[k][0] = 0.0;
    //        //inTpl_padded[k][1] = 0.0;
    //    }
    for (Int32 k = 0; k < nPadBeforeTpl; k++)
    {
        inTpl_padded[k] = 0.0;
    }
    for (Int32 k = nPadBeforeTpl; k < nPadBeforeTpl + nTpl; k++)
    {
        inTpl_padded[k] = Y[k - nPadBeforeTpl];
    }
    for (Int32 k = nPadBeforeTpl + nTpl; k < nPadded; k++)
    {
        inTpl_padded[k] = 0.0;
    }

    for (Int32 k = 0; k < nPadded; k++)
    {
        inTpl[k] = inTpl_padded[k]; // inTpl_padded[nPadded-1-k][0];//
    }

    if (verboseExportXtYFFT)
    {
        // save tpl input data
        FILE *f_fftinput = fopen("loglbda_fitallz_yfftInput_dbg.txt", "w+");
        for (Int32 t = 0; t < nPadded; t++)
        {
            fprintf(f_fftinput, "%f\t%e\n", (Float64)t, inTpl[t]);
        }
        fclose(f_fftinput);
    }
    fftw_execute(pTpl);
    if (verboseLogXtYFFT)
    {
        Log.LogInfo("  Operator-TemplateFittingLog: FitAllz: tpl-fft done");
    }
    if (verboseExportXtYFFT)
    {
        // save tpl-fft data
        FILE *f_fftoutput = fopen("loglbda_fitallz_yfftOutput_dbg.txt", "w+");
        for (Int32 t = 0; t < nPadded; t++)
        {
            fprintf(f_fftoutput, "%f\t%e\n", (Float64)t,
                    outTpl[t][0] * outTpl[t][0] + outTpl[t][1] * outTpl[t][1]);
        }
        fclose(f_fftoutput);
    }

    // Multiplying the FFT outputs
    for (Int32 k = 0; k < nPadded; k++)
    {
        outCombined[k][0] = (outTpl[k][0] * outSpc[k][0] - outTpl[k][1] * outSpc[k][1]);
        outCombined[k][1] = (outTpl[k][0] * outSpc[k][1] + outTpl[k][1] * outSpc[k][0]);
        // Y conjugate
        // outCombined[k][0] =
        // (outTpl[k][0]*outSpc[k][0]+outTpl[k][1]*outSpc[k][1]);
        // outCombined[k][1] =
        // (outTpl[k][0]*outSpc[k][1]-outTpl[k][1]*outSpc[k][0]);
    }

    fftw_execute(pBackward);
    if (verboseLogXtYFFT)
    {
        Log.LogDebug("  Operator-TemplateFittingLog: FitAllz: backward-fft done");
    }

    XtY.resize(nshifts);
    Int32 offsetSamples = 0; // nPadded/2.0;//nSpc;//(nTpl+nSpc);
    for (Int32 k = 0; k < nshifts; k++)
    {
        // XtY[k] = inCombined[k+(Int32)(nPadded/2.0)][0]/nPadded;
        Int32 IndexCyclic = k + offsetSamples;
        if (IndexCyclic < 0)
        {
            IndexCyclic += nPadded;
        } else if (IndexCyclic >= nPadded)
        {
            IndexCyclic -= nPadded;
        }

        XtY[k] = inCombined[IndexCyclic] / (Float64)nPadded;
    }

    if (verboseExportXtYFFT)
    {
        // save rebinned data
        FILE *f_fftoutput = fopen("loglbda_fitallz_xtyfft-output_dbg.txt", "w+");
        for (Int32 t = 0; t < nshifts; t++)
        {
            fprintf(f_fftoutput, "%f\t%e\n", (Float64)t, XtY[t]);
        }
        fclose(f_fftoutput);
    }

    return 0;
}

Int32 COperatorTemplateFittingLog::InitFFT(Int32 nPadded)
{
    freeFFTPlans();

    inSpc = (Float64 *)fftw_malloc(sizeof(Float64) * nPadded);
    outSpc = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nPadded);
    pSpc = fftw_plan_dft_r2c_1d(nPadded, inSpc, outSpc, FFTW_ESTIMATE);
    if (inSpc == 0)
    {
        throw GlobalException(INTERNAL_ERROR,"  Operator-TemplateFittingLog: InitFFT: Unable to allocate inSpc");
    }
    if (outSpc == 0)
    {
        throw GlobalException(INTERNAL_ERROR,"  Operator-TemplateFittingLog: InitFFT: Unable to allocate outSpc");
    }

    inTpl_padded = (Float64 *)fftw_malloc(sizeof(Float64) * nPadded);
    inTpl = (Float64 *)fftw_malloc(sizeof(Float64) * nPadded);
    outTpl = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nPadded);
    pTpl = fftw_plan_dft_r2c_1d(nPadded, inTpl, outTpl, FFTW_ESTIMATE);
    if (inTpl_padded == 0)
    {
        throw GlobalException(INTERNAL_ERROR,"  Operator-TemplateFittingLog: InitFFT: Unable to allocate inTpl_padded");
    }
    if (inTpl == 0)
    {
        throw GlobalException(INTERNAL_ERROR,"  Operator-TemplateFittingLog: InitFFT: Unable to allocate inTpl");
    }
    if (outTpl == 0)
    {
        throw GlobalException(INTERNAL_ERROR,"  Operator-TemplateFittingLog: InitFFT: Unable to allocate outTpl");
    }

    outCombined = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nPadded);
    inCombined = (Float64 *)fftw_malloc(sizeof(Float64) * nPadded);
    pBackward = fftw_plan_dft_c2r_1d(nPadded, outCombined, inCombined, FFTW_ESTIMATE);
    if (outCombined == 0)
    {
        throw GlobalException(INTERNAL_ERROR,"  Operator-TemplateFittingLog: InitFFT: Unable to allocate outCombined");
    }
    if (inCombined == 0)
    {
        throw GlobalException(INTERNAL_ERROR,"  Operator-TemplateFittingLog: InitFFT: Unable to allocate inCombined");
    }

    // reinit fft precomputed buffers
    freeFFTPrecomputedBuffers();

    return 0;
}

void COperatorTemplateFittingLog::freeFFTPrecomputedBuffers()
{
    if (precomputedFFT_spcFluxOverErr2 != NULL)
    {
        fftw_free(precomputedFFT_spcFluxOverErr2);
        precomputedFFT_spcFluxOverErr2 = 0;
    }
    if (precomputedFFT_spcOneOverErr2 != NULL)
    {
        fftw_free(precomputedFFT_spcOneOverErr2);
        precomputedFFT_spcOneOverErr2 = 0;
    }
}

void COperatorTemplateFittingLog::freeFFTPlans()
{
    if (pSpc)
    {
        fftw_destroy_plan(pSpc);
        pSpc = 0;
    }
    if (inSpc)
    {
        fftw_free(inSpc);
        inSpc = 0;
    }
    if (outSpc)
    {
        fftw_free(outSpc);
        outSpc = 0;
    }
    if (pTpl)
    {
        fftw_destroy_plan(pTpl);
        pTpl = 0;
    }
    if (inTpl_padded)
    {
        fftw_free(inTpl_padded);
        inTpl_padded = 0;
    }
    if (inTpl)
    {
        fftw_free(inTpl);
        inTpl = 0;
    }
    if (outTpl)
    {
        fftw_free(outTpl);
        outTpl = 0;
    }
    if (pBackward)
    {
        fftw_destroy_plan(pBackward);
        pBackward = 0;
    }
    if (inCombined)
    {
        fftw_free(inCombined);
        inCombined = 0;
    }
    if (outCombined)
    {
        fftw_free(outCombined);
        outCombined = 0;
    }

    freeFFTPrecomputedBuffers();
}

/**
 * @brief COperatorTemplateFittingLog::FitAllz
 *
 * @param lambdaRange
 * @param result
 * @param opt_extinction
 * @param opt_dustFitting
 * @param spcMaskAdditional
 * @param logpriorze: if size=0, prior is deactivated
 * @return
 */
Int32 COperatorTemplateFittingLog::FitAllz(std::shared_ptr<CTemplateFittingResult> result,
                                           const TInt32List &MeiksinList,
                                           const TInt32List &EbmvList,
                                           const CPriorHelper::TPriorZEList &logpriorze)
{
    bool verboseLogFitAllz = true;

    // prepare list of redshifts that need a full lbdarange lst-square
    // calculation
    TInt32RangeList izrangelist;
    TInt32List zindexesFullLstSquare;
    if (m_enableIGM && result->Redshifts.size() > 1)
    { 
        Float64 zmax_firstRange = GetIGMStartingRedshiftValue(m_ssSpectrum.GetSpectralAxis()[0]);
        zindexesFullLstSquare.push_back(0); // first index is always a mandatory full Lstsq Calculation case
        TFloat64List zlistsegments = m_templateRebined_bf.m_igmCorrectionMeiksin->GetSegmentsStartRedshiftList();
        for (Int32 k = 0; k < zlistsegments.size(); k++)
        {
            Int32 i_min = -1;
            bool b;
            if(zlistsegments[k]<zmax_firstRange) 
                b = CIndexing<Float64>::getClosestLowerIndex(result->Redshifts,zmax_firstRange, i_min);
            else 
                b = CIndexing<Float64>::getClosestLowerIndex(result->Redshifts,zlistsegments[k], i_min);
            if(b)
                zindexesFullLstSquare.push_back(i_min);
        }
        zindexesFullLstSquare.erase(std::unique(zindexesFullLstSquare.begin(),
                                                zindexesFullLstSquare.end()),
                                                zindexesFullLstSquare.end());

        if (verboseLogFitAllz)
        {
            Log.LogDebug("  Operator-TemplateFittingLog: FitAllz: indexes for full LstSquare calculation, count = %d", zindexesFullLstSquare.size());
            for (Int32 k = 0; k < zindexesFullLstSquare.size(); k++)
            {
                Log.LogDebug("  Operator-TemplateFittingLog: FitAllz: indexes ranges: for i=%d, zindexesFullLstSquare=%d", k, zindexesFullLstSquare[k]);
            }
        }

        UInt32 izmin = zindexesFullLstSquare[0];
        for (Int32 k = 1; k < zindexesFullLstSquare.size(); k++)
        {
            UInt32 izmax =  zindexesFullLstSquare[k];
            izrangelist.push_back(TInt32Range(izmin, izmax));
            izmin = izmax +1; //setting min for next range (one value overlapping)        
        }
        if (izrangelist.size() == 0)
            izrangelist.push_back(TInt32Range(izmin, result->Redshifts.size() - 1));
    } else
    {
        UInt32 izmin = 0;
        UInt32 izmax = result->Redshifts.size() - 1;
        izrangelist.push_back(TInt32Range(izmin, izmax));
    }
    UInt32 nzranges = izrangelist.size();

    if (verboseLogFitAllz)
    {
        Log.LogDebug("  Operator-TemplateFittingLog: FitAllz: indexes - izrangelist calculation, count = %d",
                     izrangelist.size());
        for (Int32 k = 0; k < nzranges; k++)
        {
            Log.LogDebug("  Operator-TemplateFittingLog: FitAllz: indexes ranges: for i=%d, zmin=%f, zmax=%f",
                         k, result->Redshifts[izrangelist[k].GetBegin()],
                         result->Redshifts[izrangelist[k].GetEnd()]);
        }
    }    
    //since dtd is cte, better compute it here
    const TAxisSampleList & error =  m_ssSpectrum.GetFluxAxis().GetError().GetSamplesVector();;
    const TAxisSampleList & spectrumRebinedFluxRaw = m_ssSpectrum.GetFluxAxis().GetSamplesVector();
    Float64 dtd = 0.0;
    TFloat64List inv_err2(error.size()); 
    for (Int32 j = 0; j < error.size(); j++)
    {
        inv_err2[j] = 1.0 / (error[j] * error[j]);
        dtd += spectrumRebinedFluxRaw[j] * spectrumRebinedFluxRaw[j] * inv_err2[j];
    }

    for (Int32 k = 0; k < nzranges; k++)
    {
        // prepare the zrange-result container
        std::shared_ptr<CTemplateFittingResult> subresult = std::shared_ptr<CTemplateFittingResult>(new CTemplateFittingResult());
        TFloat64Range zrange = TFloat64Range(result->Redshifts[izrangelist[k].GetBegin()],
                                             result->Redshifts[izrangelist[k].GetEnd()]);
        Float64 redshiftStep;   
        if (m_enableIGM && result->Redshifts.size() > 1)
        {
            TFloat64List::const_iterator first = result->Redshifts.begin() + izrangelist[k].GetBegin(),
                                         last = result->Redshifts.begin() + izrangelist[k].GetEnd()+1;
            TFloat64List subRedshifts(first, last);
            subresult->Init(subRedshifts.size(), EbmvList.size(), MeiksinList.size());
            subresult->Redshifts = subRedshifts;
            // slice the template
            redshiftStep = log((subRedshifts[1]+1.)/(subRedshifts[0]+1.));          
        } else{       
            redshiftStep = log((result->Redshifts[1]+1.)/(result->Redshifts[0]+1.));
            subresult->Init(result->Redshifts.size(), EbmvList.size(), MeiksinList.size());
            subresult->Redshifts = result->Redshifts;
        }
        TInt32Range ilbda = FindTplSpectralIndex(zrange);
        if (verboseLogFitAllz)
        {
            Log.LogDebug("  Operator-TemplateFittingLog: FitAllz: zrange min=%f, max=%f",
                         zrange.GetBegin(), zrange.GetEnd());
            Log.LogDebug("  Operator-TemplateFittingLog: FitAllz: full zmin=%f, full zmax=%f",
                         result->Redshifts[0],
                         result->Redshifts[result->Redshifts.size() - 1]);
            Log.LogDebug("  Operator-TemplateFittingLog: FitAllz: indexes tpl crop: "
                         "lbda min=%d, max=%d",
                         ilbda.GetBegin(), ilbda.GetEnd());
            Log.LogDebug("  Operator-TemplateFittingLog: FitAllz: indexes tpl full: "
                         "lbda min=%d, max=%d",
                         0, m_templateRebined_bf.GetSampleCount()-1);
            /*Log.LogDebug("  Operator-TemplateFittingLog: FitAllz: tpl lbda "
                         "min*zmax=%f, max*zmin=%f",
                         tplRebinedLambdaGlobal[ilbda.GetBegin()] * (1.0 + zrange.GetEnd()),
                         tplRebinedLambdaGlobal[ilbda.GetEnd()] * (1.0 + zrange.GetBegin()));*/
        }

        FitRangez(inv_err2, ilbda, subresult, MeiksinList, EbmvList, dtd);

        // copy subresults into global results
        for (UInt32 isubz = 0; isubz < subresult->Redshifts.size(); isubz++)
        {
            UInt32 fullResultIdx = isubz + izrangelist[k].GetBegin();
            result->ChiSquare[fullResultIdx] = subresult->ChiSquare[isubz];
            result->FitAmplitude[fullResultIdx] = subresult->FitAmplitude[isubz];
            result->FitAmplitudeError[fullResultIdx] = subresult->FitAmplitudeError[isubz];
            result->FitAmplitudeSigma[fullResultIdx] = subresult->FitAmplitudeSigma[isubz];
            result->FitDtM[fullResultIdx] = subresult->FitDtM[isubz];
            result->FitMtM[fullResultIdx] = subresult->FitMtM[isubz];


            Float64 logprior = 0.;
            if(logpriorze.size()>0)
            {
                bool verbose_priorA = false;
                Int32 kism_best = -1;
                if(subresult->FitEbmvCoeff[isubz]==-1)
                {
                    kism_best=0;
                }else{
                    kism_best = m_templateRebined_bf.m_ismCorrectionCalzetti->GetEbmvIndex(subresult->FitEbmvCoeff[isubz]);
                }

                const CPriorHelper::SPriorTZE &pTZE = logpriorze[fullResultIdx][kism_best];
                logprior += -2.0*pTZE.betaTE*pTZE.logprior_precompTE;
                logprior += -2.0*pTZE.betaA*pTZE.logprior_precompA;
                logprior += -2.0*pTZE.betaZ*pTZE.logprior_precompZ;


                if(pTZE.A_sigma>0.0 && pTZE.A_mean>0.0)
                {
                    //now update the amplitude if there is any constraints from the priors
                    Float64 ampl = result->FitAmplitude[fullResultIdx];
                    Float64 ampl_err = result->FitAmplitudeError[fullResultIdx];
                    Float64 ampl_sigma = result->FitAmplitudeSigma[fullResultIdx];
                    if(pTZE.betaA>0.0)
                    {
                        Float64 bss2 = pTZE.betaA/(pTZE.A_sigma*pTZE.A_sigma);
                        ampl = (result->FitDtM[fullResultIdx]+pTZE.A_mean*bss2)/(result->FitMtM[fullResultIdx]+bss2);
                        ampl_err = sqrt(result->FitMtM[fullResultIdx])/(result->FitMtM[fullResultIdx]+bss2); 

                    }else{
                        ampl = result->FitDtM[fullResultIdx]/result->FitMtM[fullResultIdx];
                        ampl_err = sqrt(1./result->FitMtM[fullResultIdx]);
                    }

                    if(verbose_priorA)
                    {
                        Log.LogDebug("  Operator-TemplateFittingLog: update the amplitude (a_mean=%e, a_sigma=%e)", pTZE.A_mean, pTZE.A_sigma);
                        Log.LogDebug("  Operator-TemplateFittingLog: update the amplitude (ampl was = %e, updated to %e)",
                                        result->FitAmplitude[fullResultIdx], ampl);
                    }

                    // check negative amplitude
                    ampl_sigma = ampl/ampl_err; 
                    
                    // force positivity
                    ampl = max(0., ampl);

                    result->FitAmplitude[fullResultIdx] = ampl;
                    result->FitAmplitudeError[fullResultIdx] = ampl_err;
                    result->FitAmplitudeSigma[fullResultIdx] = ampl_sigma;
                    result->ChiSquare[fullResultIdx] = dtd + result->FitMtM[fullResultIdx]*ampl*ampl - 2.*ampl*result->FitDtM[fullResultIdx];

                    Float64 logPa = pTZE.betaA*(ampl-pTZE.A_mean)*(ampl-pTZE.A_mean)/(pTZE.A_sigma*pTZE.A_sigma);
                    if(std::isnan(logPa) || logPa!=logPa || std::isinf(logPa))
                    {
		      throw GlobalException(INTERNAL_ERROR,Formatter()<<"Operator-TemplateFittingLog: logPa is NAN (a_mean="<<pTZE.A_mean<<", a_sigma="<< pTZE.A_sigma);
                    }
                    logprior += logPa;
                }else{
                    if(verbose_priorA)
                    {
                        Log.LogDebug("  Operator-TemplateFittingLog: NOT updating the amplitude (a_mean=%e, a_sigma=%e)",
                                     pTZE.A_mean,
                                     pTZE.A_sigma);
                    }
                }
                if(std::isnan(logprior) || logprior!=logprior || std::isinf(logprior))
                {
		  throw GlobalException(INTERNAL_ERROR,Formatter()<<"Operator-TemplateFittingLog: logPa is NAN (a_mean="<<pTZE.A_mean<<", a_sigma="<<pTZE.A_sigma<<", precompA="<<pTZE.logprior_precompA);
                }
                result->ChiSquare[fullResultIdx] += logprior;
            }
            result->Overlap[fullResultIdx] = subresult->Overlap[isubz];
            result->LogPrior[fullResultIdx] = logprior;
            result->FitEbmvCoeff[fullResultIdx] = subresult->FitEbmvCoeff[isubz];
            result->FitMeiksinIdx[fullResultIdx] = subresult->FitMeiksinIdx[isubz];
            result->Status[fullResultIdx] = subresult->Status[isubz];

            for (Int32 kism = 0; kism < result->ChiSquareIntermediate[fullResultIdx].size(); kism++)
            {
                for (Int32 kigm = 0;kigm < result->ChiSquareIntermediate[fullResultIdx][kism].size(); kigm++)
                {
                    result->ChiSquareIntermediate[fullResultIdx][kism][kigm] =
                            subresult->ChiSquareIntermediate[isubz][kism][kigm];
                    if(logpriorze.size()>0)
                    {
                        Float64 logprior = 0.; //not implemented -> not a problem for fullmodel, but will be necessary for tplmodel method for example
                        result->ChiSquareIntermediate[fullResultIdx][kism][kigm] += logprior;
                    }
                }
            }
        }
    }

    return 0;
}

/**
  // TODO : many vectors allocated in this function. Check if the allocation
 time is significant, and eventually use preallocated member buffers...
 * @brief COperatorTemplateFittingLog::FitRangez
 * @param spectrumRebinedLambda
 * @param spectrumRebinedFluxRaw
 * @param error
 * @param tplRebinedLambda
 * @param tplRebinedFluxRaw
 * @param nSpc
 * @param nTpl
 * @param result
 * @param MeiksinList
 * @param EbmvList
 * @return
 */
Int32 COperatorTemplateFittingLog::FitRangez(const TFloat64List &inv_err2,
                                             const TInt32Range &currentRange,
                                             const std::shared_ptr<CTemplateFittingResult> &result,
                                             const TInt32List &MeiksinList,
                                             const TInt32List &EbmvList,
                                             const Float64& dtd)
{
    const TAxisSampleList & spectrumRebinedLambda = m_ssSpectrum.GetSpectralAxis().GetSamplesVector();
    const TAxisSampleList & spectrumRebinedFluxRaw = m_ssSpectrum.GetFluxAxis().GetSamplesVector();
    UInt32 nSpc = spectrumRebinedLambda.size();

    const TAxisSampleList & tplRebinedLambdaGlobal = m_templateRebined_bf.GetSpectralAxis().GetSamplesVector();

    Int32 kstart, kend;
    kstart = currentRange.GetBegin(); kend = currentRange.GetEnd();
    UInt32 nTpl = kend - kstart + 1;

    TAxisSampleList spcRebinedFluxOverErr2(nSpc);
    for (Int32 j = 0; j < nSpc; j++)
    {
        spcRebinedFluxOverErr2[j] = spectrumRebinedFluxRaw[j] * inv_err2[j];
    }

   Float64 redshiftValueMeiksin = result->Redshifts[0];

   if (verboseLogFitFitRangez)
    {
        Log.LogDebug("  Operator-TemplateFittingLog: FitRangez: redshiftValueMeiksin = %f", redshiftValueMeiksin);
        Log.LogDebug("  Operator-TemplateFittingLog: FitRangez: spc[0] = %f", spectrumRebinedLambda[0]);
        Log.LogDebug("  Operator-TemplateFittingLog: FitRangez: spc[max] = %f", spectrumRebinedLambda[nSpc - 1]);
        Log.LogDebug("  Operator-TemplateFittingLog: FitRangez: tpl[0]*zmax = %f", tplRebinedLambdaGlobal[kstart] * (1.0 + result->Redshifts[result->Redshifts.size() - 1]));
        Log.LogDebug("  Operator-TemplateFittingLog: FitRangez: tpl[max]*zmin = %f", tplRebinedLambdaGlobal[kstart+nTpl - 1] * (1 + result->Redshifts[0]));
    }

    Int32 nshifts = nTpl - nSpc + 1;
    m_nPaddedSamples = nTpl;

    Log.LogDetail("  Operator-TemplateFittingLog: Now fitting using the FFT on "
                  "nshifts=%d values, for Meiksin redshift=%f",
                  nshifts, redshiftValueMeiksin);

    if (verboseLogFitFitRangez)
    {
        Log.LogDebug("  Operator-TemplateFittingLog: FitRangez: initializing FFT with n = %d points",
                     m_nPaddedSamples);
    }
    InitFFT(m_nPaddedSamples);

    TFloat64List z_vect = result->Redshifts;
    std::reverse(z_vect.begin(), z_vect.end());
    // prepare z array
    TFloat64List z_vect_verif(nshifts, 0.0);
    Float64 relative_zgrid_error_max = 0.;
    for (Int32 t = 0; t < nshifts; t++)
    {
        z_vect_verif[t] = (spectrumRebinedLambda[0] - tplRebinedLambdaGlobal[t + kstart]) /
                    tplRebinedLambdaGlobal[t + kstart];
        //compare with z_vect
        Float64 relative_zgrid_error = std::abs(z_vect[t] - z_vect_verif[t])/(1+z_vect_verif[t]);
        if (relative_zgrid_error > relative_zgrid_error_max) relative_zgrid_error_max = relative_zgrid_error;
    }
    if (verboseLogFitFitRangez)
        Log.LogDebug("  Operator-TemplateFittingLog: FitRangez: max diff in zgrid=%e",relative_zgrid_error_max);
    if (relative_zgrid_error_max>5E-7){
            throw GlobalException(INTERNAL_ERROR,"z_vect and z_vect_verification do not correspond.");
        }
    
    //check borders
    if(z_vect.size()!=z_vect_verif.size())
        throw GlobalException(INTERNAL_ERROR,"z_vect size and z_vect_verification size do not match.");
    Int32 nISM = EbmvList.size();
    Int32 nIGM = MeiksinList.size();

    // disable IGM if the redshift range and lambda range do not make the IGM
    // wavelength appear
    Int32 enableIGM = m_enableIGM;
    Int32 overrideNIGMTobesaved = -1;
    if (tplRebinedLambdaGlobal[kstart] > 1216. && nIGM > 1)
    {
        overrideNIGMTobesaved = nIGM;
        nIGM = 1;
        enableIGM = 0; 
        if (verboseLogFitFitRangez)
        {
            Log.LogDebug("  Operator-TemplateFittingLog: FitRangez: IGM disabled, min-tpl-lbda=%f", tplRebinedLambdaGlobal[kstart]);
        }
    }

    // prepare best fit data buffer
    TFloat64List bestChi2(nshifts, DBL_MAX);
    TFloat64List bestFitAmp(nshifts, -1.0);
    TFloat64List bestFitAmpErr(nshifts, -1.0);
    TFloat64List bestFitAmpSigma(nshifts);
    TFloat64List bestFitDtm(nshifts, -1.0);
    TFloat64List bestFitMtm(nshifts, -1.0);
    TFloat64List bestISMCoeff(nshifts, -1.0);
    TInt32List   bestIGMIdx(nshifts, -1);

    // prepare intermediate fit data buffer
    std::vector<std::vector<TFloat64List>> intermediateChi2;
    Int32 nIGMFinal = nIGM;
    if (overrideNIGMTobesaved > nIGM)
    {
        nIGMFinal = overrideNIGMTobesaved;
    }
    for (Int32 k = 0; k < nshifts; k++)
    {
        std::vector<TFloat64List> _ChiSquareISMList;
        for (Int32 kism = 0; kism < nISM; kism++)
        {
            TFloat64List _chi2List(nIGMFinal, DBL_MAX);
            _ChiSquareISMList.push_back(_chi2List);
        }
        intermediateChi2.push_back(_ChiSquareISMList);
    }

   //note that there is no need to copy the ism/igm cause they already exist in the rebinned template
    if(m_enableIGM || m_enableISM){
        m_templateRebined_bf.InitIsmIgmConfig(kstart, kend, redshiftValueMeiksin);
    }
    for (Int32 kIGM = 0; kIGM < nIGM; kIGM++)
    {
        if (verboseLogFitFitRangez && enableIGM)
        {
            Log.LogDebug("  Operator-TemplateFittingLog: FitRangez: IGM index=%d", kIGM);
        }

        if (enableIGM)
        {
            Int32 meiksinIdx = MeiksinList[kIGM];
            m_templateRebined_bf.ApplyMeiksinCoeff(meiksinIdx);
        }

        for (Int32 kISM = 0; kISM < nISM; kISM++)
        {
            if (verboseLogFitFitRangez && m_enableISM)
            {
                Log.LogDebug("  Operator-TemplateFittingLog: FitRangez: ISM index =%d", kISM);
            }

            if (m_enableISM)
            {
                Int32 kDust = EbmvList[kISM];
                m_templateRebined_bf.ApplyDustCoeff(kDust);
            }

            const TAxisSampleList & tplRebinedFluxcorr = m_templateRebined_bf.GetFluxAxis().GetSamplesVector();
            //extract only the relevant part
            TFloat64List::const_iterator first = tplRebinedFluxcorr.begin() + kstart,
                                last = tplRebinedFluxcorr.begin() + kend+1;
            const TAxisSampleList tplRebinedFluxcorr_cropped(first, last);
            TAxisSampleList tpl2RebinedFlux(nTpl);

            if(tplRebinedFluxcorr_cropped.size()!=nTpl){
                throw GlobalException(INTERNAL_ERROR,"prob with retrieved vector size");
            }
            //compute the square of the corrected flux
            for (Int32 j = 0; j < nTpl; j++)
            {
                tpl2RebinedFlux[j] = tplRebinedFluxcorr_cropped[j] * tplRebinedFluxcorr_cropped[j];
            }

            if (verboseExportFitRangez_model)
            {
                if ((m_enableISM && exportISMIdx == EbmvList[kISM]) ||(!m_enableISM))
                {
                    if ((enableIGM && exportIGMIdx == MeiksinList[kIGM]) || (!enableIGM))
                    {
                        // save chi2 data
                        FILE *f = fopen("loglbda_model_dbg.txt", "w+");
                        for (Int32 j = 0; j < nTpl; j++)
                        {
                            fprintf(f, "%f\t%e\n", tplRebinedLambdaGlobal[kstart+j], tplRebinedFluxcorr_cropped[j]);
                        }
                        fclose(f);
                    }
                }
            }

            // Estimate DtM: sumCross
            TFloat64List dtm_vec;
            EstimateXtY(spcRebinedFluxOverErr2, tplRebinedFluxcorr_cropped, nshifts, dtm_vec, 0);
            
	        Int32 dtm_vec_size = dtm_vec.size();
            if (verboseExportFitRangez)
            {
                // save chi2 data
                FILE *f = fopen("loglbda_dtm_dbg.txt", "w+");
                for (Int32 t = 0; t < dtm_vec_size; t++)
                {
                    fprintf(f, "%f\t%e\n", (Float64)t, dtm_vec[t]);
                }
                fclose(f);
            }

            // Estimate MtM: sumT
            TFloat64List mtm_vec;
            EstimateXtY(inv_err2, tpl2RebinedFlux, nshifts, mtm_vec, 1);

            if (verboseExportFitRangez)
            {
                // save chi2 data
                FILE *f = fopen("loglbda_mtm_dbg.txt", "w+");
                for (Int32 t = 0; t < mtm_vec.size(); t++)
                {
                    fprintf(f, "%f\t%e\n", (Float64)t, mtm_vec[t]);
                }
                fclose(f);
            }

            if (verboseExportFitRangez)
            {
                Log.LogDebug("  Operator-TemplateFittingLog: FitRangez: dtd = %e", dtd);
            }

            // Estimate Chi2
            if (mtm_vec.size() != dtm_vec_size)
            {
                freeFFTPlans();
                throw GlobalException(INTERNAL_ERROR,Formatter()<<"Operator-TemplateFittingLog: FitRangez: xty vectors sizes don't match: dtm size = "<< dtm_vec_size <<", mtm size ="<< mtm_vec.size());
            }
            //Log.LogDetail("  Operator-TemplateFittingLog: FitRangez: kISM = %d, kIGM = %d", kISM, kIGM);
            TFloat64List chi2(dtm_vec_size, DBL_MAX);
            TFloat64List amp(dtm_vec_size, DBL_MAX);
            TFloat64List amp_sigma(dtm_vec_size);
            TFloat64List amp_err(dtm_vec_size, DBL_MAX);
            for (Int32 k = 0; k < dtm_vec_size; k++)
            {
                if (mtm_vec[k] == 0.0)
                {
                    amp[k] = 0.0;
                    amp_err[k] = 0.0;
                    amp_sigma[k] = 0.0;
                    chi2[k] = dtd;
                } else
                {
                    amp[k] = dtm_vec[k] / mtm_vec[k];
                    amp_err[k] = sqrt(1./mtm_vec[k]);
                    amp_sigma[k] = amp[k]/amp_err[k];
                    amp[k] = max(0.0, amp[k]);

                    chi2[k] = dtd - 2 * dtm_vec[k] * amp[k] + mtm_vec[k] * amp[k] * amp[k];
                }
                //Log.LogDetail("  Operator-TemplateFittingLog: FitRangez: chi2[%d] = %f", k, chi2[k]);
            }

            for (Int32 k = 0; k < dtm_vec_size; k++)
            {
                intermediateChi2[k][kISM][kIGM] = chi2[k];
                // in the case of 1215A is not in the range, no need to
                // recompute with varying IGM coeff.
                if (overrideNIGMTobesaved > 1 && kIGM == 0)
                {
                    for (Int32 koigm = 1; koigm < overrideNIGMTobesaved; koigm++)
                    {
                        intermediateChi2[k][kISM][koigm] = chi2[k];
                    }
                }

                if (bestChi2[k] > chi2[k])
                {
                    bestChi2[k] = chi2[k];
                    bestFitAmp[k] = amp[k];
                    bestFitAmpErr[k] = amp_err[k];
                    bestFitAmpSigma[k] = amp_sigma[k];
                    bestFitDtm[k] = dtm_vec[k];
                    bestFitMtm[k] = mtm_vec[k];
                    bestISMCoeff[k] = m_enableISM ? m_templateRebined_bf.m_ismCorrectionCalzetti->GetEbmvValue(EbmvList[kISM]) : -1;
                    bestIGMIdx[k] = enableIGM ? MeiksinList[kIGM] : -1;
                }
            }

            if (verboseExportFitRangez)
            {
                Log.LogDebug("  Operator-TemplateFittingLog: FitRangez: spc lbda 0 =%f", spectrumRebinedLambda[0]);
                Log.LogDebug("  Operator-TemplateFittingLog: FitRangez: tpl lbda 0 =%f", tplRebinedLambdaGlobal[kstart]);
                Float64 z_O = (spectrumRebinedLambda[0] - tplRebinedLambdaGlobal[kstart]) / tplRebinedLambdaGlobal[kstart];
                Log.LogDebug("  Operator-TemplateFittingLog: FitRangez: z 0 =%f", z_O);

                // save chi2 data
                FILE *f_chi2 = fopen("loglbda_chi2output_dbg.txt", "w+");
                for (Int32 t = 0; t < chi2.size(); t++)
                {
                    fprintf(f_chi2, "%f\t%e\n", z_vect[t], chi2[t]);
                }
                fclose(f_chi2);
            }
        }
    }

    // reversing all vectors
    std::reverse(z_vect.begin(), z_vect.end());
    std::reverse(bestChi2.begin(), bestChi2.end());
    std::reverse(bestFitAmp.begin(), bestFitAmp.end());
    std::reverse(bestFitAmpErr.begin(), bestFitAmpErr.end());
    std::reverse(bestFitAmpSigma.begin(), bestFitAmpSigma.end());
    std::reverse(bestFitDtm.begin(), bestFitDtm.end());
    std::reverse(bestFitMtm.begin(), bestFitMtm.end());
    std::reverse(bestISMCoeff.begin(), bestISMCoeff.end());
    std::reverse(bestIGMIdx.begin(), bestIGMIdx.end());
    TFloat64List intermChi2BufferReversed_array(intermediateChi2.size());
    for (Int32 kism = 0; kism < nISM; kism++){
        for (Int32 kigm = 0; kigm < nIGMFinal; kigm++){
            for (Int32 t = 0; t < nshifts; t++)
                intermChi2BufferReversed_array[t] = intermediateChi2[nshifts - 1 - t][kism][kigm];
            for (Int32 t = 0; t < nshifts; t++)
                result->ChiSquareIntermediate[t][kism][kigm] = intermChi2BufferReversed_array[t];
        }
    }
    for (Int32 k = 0; k < result->Redshifts.size(); k++)
    {
        result->Overlap[k] = 1.0;
        result->Status[k] = nStatus_OK;
    }
    result->ChiSquare = bestChi2;
    result->FitAmplitude = bestFitAmp;
    result->FitAmplitudeError = bestFitAmpErr;
    result->FitAmplitudeSigma  = bestFitAmpSigma;
    result->FitDtM = bestFitDtm;
    result->FitMtM = bestFitMtm;
    result->FitEbmvCoeff = bestISMCoeff;
    result->FitMeiksinIdx = bestIGMIdx;

    freeFFTPlans();
    return 0;
}

Int32 COperatorTemplateFittingLog::InterpolateResult(const TFloat64List& in,
                                                     TFloat64List& inGrid,
                                                     const TFloat64List& tgtGrid,
                                                     TFloat64List &out,
                                                     Float64 defaultValue)
{
    bool debug = 0;
    Int32 n = inGrid.size();
    Int32 tgtn = tgtGrid.size();
    out.resize(tgtn);
    for (Int32 j = 0; j < tgtn; j++)
    {
        out[j] = defaultValue;
    }

    //* // GSL method spline
    // initialise and allocate the gsl objects
    // lin
    gsl_interp *interpolation = gsl_interp_alloc(gsl_interp_linear, n);
    gsl_interp_accel *accelerator = gsl_interp_accel_alloc();
    Int32 status;
    // spline
    // gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n);
    // gsl_spline_init (spline, inGrid, in, n);
    // gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();

    gsl_error_handler_t *gsl_error_handler_old = gsl_set_error_handler_off();

    // Deal with exp/log accuracy
    Float64 epsilon = 1E-8;

    if (tgtGrid[0] < inGrid[0])
    {
        if ( (inGrid[0] - tgtGrid[0]) < epsilon)
        {
            inGrid[0] = tgtGrid[0] - epsilon;
        }
        else
        {
	  throw GlobalException(INTERNAL_ERROR,Formatter()<<"Error while interpolating loglambda chi2 result : xrebin("<< tgtGrid[0] <<") < x[0]("<< inGrid[0] <<")");
        }
    }
    if (tgtGrid[tgtn-1] > inGrid[n-1])
    {
        if ((tgtGrid[tgtn-1] - inGrid[n-1]) < epsilon)
        {
            inGrid[n-1] = tgtGrid[tgtn-1] + epsilon;
        }
        else
        {
	  throw GlobalException(INTERNAL_ERROR,Formatter()<<"Error while interpolating loglambda chi2 result : xrebin("<< tgtGrid[tgtn-1] <<") > x[n-1]("<< inGrid[n-1] << ")");
        }
    }

    gsl_interp_init(interpolation, &inGrid.at(0), &in.at(0), n);

    for (Int32 j = 0; j < tgtn; j++)
    {
        Float64 Xrebin = tgtGrid[j];
        if(debug)
        {
            Log.LogDebug("  Operator-TemplateFittingLog: InterpolateResult, j=%d, xrebin=%f",
                     j,
                     Xrebin);
        }

        status = gsl_interp_eval_e(interpolation,
                                   &inGrid.at(0),
                                   &in.at(0),
                                   Xrebin,
                                   accelerator,
                                   &out[j]); // lin

        if (status != GSL_SUCCESS) {
	  throw GlobalException(EXTERNAL_LIB_ERROR,Formatter()<<"Operator-TemplateFittingLog: InterpolateError: GSL Error : "<<" gsl: "<< __FILENAME__<<":"<<__LINE__<<": ERROR:"<< status<<" (Errtype: "<<gsl_strerror(status)<<")");
        }
        // out[j] = gsl_spline_eval (spline, Xrebin, accelerator); //spline
        // Log.LogInfo("  Operator-TemplateFittingLog: FitAllz: interpolating
        // gsl-spline z result, , ztgt=%f, rebinY=%f", tgtGrid[j], out[j]);
    }

    //TODO why should there be an old gsl_error_handler ?
    gsl_set_error_handler(gsl_error_handler_old);

    gsl_interp_free(interpolation);
    // gsl_spline_free (spline);
    gsl_interp_accel_free(accelerator);

    return 0;
}
//find indexes in templateSpectra for which Z falls into the redshift range
/**
 * Method logic: 
 * 1. count the number of steps (Integer) between start/end borders of redshiftRange and 
 * the min/max Z values corresponding to the templateSpectralAxis min/max borders
 * 2. this number of zsteps correspond exactly to the right/left offsets on the spectralAxis in log
*/
TInt32Range COperatorTemplateFittingLog::FindTplSpectralIndex( const TFloat64Range & redshiftrange) const
{
    return FindTplSpectralIndex(m_ssSpectrum.GetSpectralAxis(), m_templateRebined_bf.GetSpectralAxis(), redshiftrange);
}

TInt32Range COperatorTemplateFittingLog::FindTplSpectralIndex( const CSpectrumSpectralAxis& spcSpectralAxis,
                                                               const CSpectrumSpectralAxis& tplSpectralAxis,
                                                               const TFloat64Range & redshiftrange) const
{

    const TFloat64Range spcRange = spcSpectralAxis.GetLambdaRange();//this is correct only if spcSpectralAxis is intersected with the input lambdaRange
    const TFloat64Range tplRange = tplSpectralAxis.GetLambdaRange();
    const UInt32 tplsize = tplSpectralAxis.GetSamplesCount();
    const Float64 logstep = tplSpectralAxis.GetlogGridStep();

    Float64 zmax = (spcRange.GetBegin() - tplRange.GetBegin()) / tplRange.GetBegin(); //get maximum z reachable with given template & spectra
    Float64 offset_left =  log((zmax+1)/(redshiftrange.GetEnd()+1))/logstep;  //  should be integer at numerical precision before round
    Int32 ilbdamin = round(offset_left); // deduce min lambda from max reachable z and max z in current range.
    Float64 zmin = (spcRange.GetEnd() - tplRange.GetEnd()) / tplRange.GetEnd(); // get minimum reachable z with given template & spectra

    Float64 offset_right = log((redshiftrange.GetBegin()+1)/(zmin+1))/logstep;
    Int32 ilbdamax = tplsize- 1 - round(offset_right); // deduce max lambda from min reachable z and min min z in current range.
 
    if( ilbdamax>=tplsize || ilbdamin < 0)
        throw GlobalException(INTERNAL_ERROR,"FindTplSpectralIndex failed to find indexes");
    
    if (ilbdamin > ilbdamax)
    {
      throw GlobalException(INTERNAL_ERROR,Formatter()<<"  Operator-TemplateFittingLog: Problem with tpl indexes for zranges, found ilbdamin="<< ilbdamin<<" > ilbdamax="<< ilbdamax);
    }

    return TInt32Range(ilbdamin, ilbdamax);
}
/**
 * \brief COperatorTemplateFittingLog::Compute
 *
 * This method computes the log_likelihood for the input spc and the tpl on a
 *given redshift range (should be a regular grid): 0. checks :
 *      - is the redshift list a regular grid ? (assumed in the rest of the
 *method)
 *      - is the overlap always >100% in the given redshift range ?
 * 1. resample the input spectrum/tpl on a loglambda regular grid (always
 *resampling as of 2017-06-13, option todo: use the already log-sampled input
 *spectrum grid) 2. input: if additional_spcMasks size is 0, no additional mask
 *will be used, otherwise its size should match the redshifts list size
 *
 * opt_dustFitting: -1 = disabled, -10 = fit over all available indexes, positive integer 0, 1 or ... will be used as ism-calzetti index as initialized in constructor
 * lambdaRange is not clamped
 **/
std::shared_ptr<COperatorResult> COperatorTemplateFittingLog::Compute(const std::shared_ptr<const CTemplate> &logSampledTpl,
                                                                    Float64 overlapThreshold,
                                                                    const std::vector<CMask> & additional_spcMasks,
                                                                    std::string opt_interp,
                                                                    Int32 opt_extinction,
                                                                    Int32 opt_dustFitting,
                                                                    const CPriorHelper::TPriorZEList &logpriorze,
                                                                    bool keepigmism,
                                                                    Float64 FitEbmvCoeff,
                                                                    Int32 FitMeiksinIdx)
{
    Log.LogDetail("  Operator-TemplateFittingLog: starting computation for template: %s", logSampledTpl->GetName().c_str());

    if ((opt_dustFitting==-10 || opt_dustFitting>-1) && logSampledTpl->CalzettiInitFailed())
    {
        throw GlobalException(INTERNAL_ERROR,"  Operator-TemplateFittingLog: no calzetti calib. file loaded... aborting");
    }
    if( opt_dustFitting>-1 && opt_dustFitting>logSampledTpl->m_ismCorrectionCalzetti->GetNPrecomputedEbmvCoeffs()-1)
    {
        Log.LogError("  Operator-TemplateFittingLog: calzetti index overflow (opt=%d, while NPrecomputedDustCoeffs=%d)... aborting",
                     opt_dustFitting,
                     logSampledTpl->m_ismCorrectionCalzetti->GetNPrecomputedEbmvCoeffs());
	      throw GlobalException(INTERNAL_ERROR,Formatter()<<"Operator-TemplateFitting: calzetti index overflow (dustfitting="<<opt_dustFitting<<",while NPrecomputedEbmvCoeffs="<<logSampledTpl->m_ismCorrectionCalzetti->GetNPrecomputedEbmvCoeffs()<<")");
    }

    if (opt_extinction && logSampledTpl->MeiksinInitFailed())
    {
        throw GlobalException(INTERNAL_ERROR,"  Operator-TemplateFittingLog: no meiksin calib. file loaded... aborting");
    }
    
    if (!logSampledTpl->GetSpectralAxis().IsLogSampled())
    {
        throw GlobalException(INTERNAL_ERROR,"  Operator-TemplateFittingLog: template is not log sampled");
    }

    //check if spc and tpl have same step
    const Float64 epsilon = 1E-8;
    if(std::abs(m_logSampledSpectrum.GetSpectralAxis().GetlogGridStep() - logSampledTpl->GetSpectralAxis().GetlogGridStep())>epsilon)
        throw GlobalException(INTERNAL_ERROR,"TFLogOperator: tpl and spc are not sampled the same step");

    //subsample template if necessary
    if(m_ssRatio==1){//no required subsampling
        m_templateRebined_bf = *logSampledTpl;
    }else{
        TInt32Range ilbda = FindTplSpectralIndex(m_ssSpectrum.GetSpectralAxis(), logSampledTpl->GetSpectralAxis(), TFloat64Range(m_redshifts));
        TFloat64List mask_tpl = logSampledTpl->GetSpectralAxis().GetSubSamplingMask(m_ssRatio, ilbda);
        
        m_templateRebined_bf = CTemplate(*logSampledTpl, mask_tpl);
        //double make sure that subsampled spectrum is well sampled
        if (!m_templateRebined_bf.GetSpectralAxis().IsLogSampled(m_logstep))
        {
            throw GlobalException(INTERNAL_ERROR,"  Operator-TemplateFittingLog: subsampled template is not log sampled with the redshift step");
        }
    }

    //**************** Fitting at all redshifts ****************//
    //Note: below corresponds to ::BasicFit code except that redshift loop belongs to ::compute 
    // Optionally apply some IGM absorption
    TInt32List MeiksinList;
    TInt32List EbmvList;
    m_templateRebined_bf.GetIsmIgmIdxList( opt_extinction, opt_dustFitting, MeiksinList, EbmvList, keepigmism, FitEbmvCoeff, FitMeiksinIdx);
    Int32 nIGMCoeffs = MeiksinList.size();
    Int32 nISMCoeffs = EbmvList.size();

    m_enableIGM = (MeiksinList[0] == -1)?0:1;
    m_enableISM = (EbmvList[0] == -1)?0:1;
    std::shared_ptr<CTemplateFittingResult> result = std::make_shared<CTemplateFittingResult>();
    result->Init(m_redshifts.size(), nISMCoeffs, nIGMCoeffs);
    result->Redshifts = m_redshifts;

    // WARNING: no additional masks coded for use as of 2017-06-13
    if (additional_spcMasks.size() != 0)
    {
        Log.LogError("  Operator-TemplateFittingLog: No additional masks used. "
                     "Feature not coded for this log-lambda operator!)");
        /*throw std::runtime_error("  Operator-TemplateFittingLog: No additional masks used. "
                     "Feature not coded for this log-lambda operator!)");*/
    }

    if(logpriorze.size()>0 && logpriorze.size()!=m_redshifts.size())
    {
      throw GlobalException(INTERNAL_ERROR,Formatter()<<"Operator-TemplateFittingLog: prior list size("<<logpriorze.size()<<") didn't match the input redshift-list size :"<< m_redshifts.size());
    }

    Int32 retFit = FitAllz(result, MeiksinList, EbmvList, logpriorze);
    
    if (retFit != 0)
    {
      throw GlobalException(INTERNAL_ERROR,Formatter()<<"Operator-TemplateFittingLog: FitAllz failed with error "<<retFit);
    }

    //**************** End Fitting at all redshifts ****************//

    // overlap warning
    Float64 overlapValidInfZ = -1;
    for (Int32 i = 0; i < m_redshifts.size(); i++)
    {
        if (result->Overlap[i] >= overlapThreshold && overlapValidInfZ == -1)
        {
            overlapValidInfZ = m_redshifts[i];
            break;
        }
    }
    Float64 overlapValidSupZ = -1;
    for (Int32 i = m_redshifts.size() - 1; i >= 0; i--)
    {
        if (result->Overlap[i] >= overlapThreshold && overlapValidSupZ == -1)
        {
            overlapValidSupZ = m_redshifts[i];
            break;
        }
    }
    if (overlapValidInfZ != m_redshifts.front() ||
        overlapValidSupZ != m_redshifts.back())
    {
        Log.LogInfo("  Operator-TemplateFittingLog: overlap warning for %s: minz=%.3f, maxz=%.3f",
                    logSampledTpl->GetName().c_str(), overlapValidInfZ, overlapValidSupZ);
    }

    // estimate CstLog for PDF estimation
    result->CstLog = EstimateLikelihoodCstLog(m_ssSpectrum, m_ssLambdaRange); // 0.0;//Todo: check how to estimate
                                            // that value for loglambda//

    return result;
}
