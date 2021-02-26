#include <RedshiftLibrary/operator/templatefittinglog.h>

#include <RedshiftLibrary/common/mask.h>
#include <RedshiftLibrary/common/quicksort.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/operator/templatefittingresult.h>
#include <RedshiftLibrary/spectrum/axis.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/tools.h>

#include <RedshiftLibrary/log/log.h>

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

#include <assert.h>

#define NOT_OVERLAP_VALUE NAN
#include <stdio.h>

namespace bfs = boost::filesystem;
using namespace NSEpic;
using namespace std;

COperatorTemplateFittingLog::COperatorTemplateFittingLog(
    std::string calibrationPath)
{
    m_opt_spcrebin = true;

    // ISM
    m_ismCorrectionCalzetti = std::unique_ptr<CSpectrumFluxCorrectionCalzetti> (new CSpectrumFluxCorrectionCalzetti());
    m_ismCorrectionCalzetti->Init(calibrationPath, 0.0, 0.1, 10);
    //m_ismCorrectionCalzetti->Init(calibrationPath, -0.6, 0.1, 16);

    // IGM
    m_igmCorrectionMeiksin = std::unique_ptr<CSpectrumFluxCorrectionMeiksin> (new CSpectrumFluxCorrectionMeiksin());
    m_igmCorrectionMeiksin->Init(calibrationPath);

    m_nPaddedSamples = 0;
    inSpc = 0;
    outSpc = 0;
    pSpc = 0;
    inTpl = 0;
    inTpl_padded = 0;
    pTpl = 0;
    outTpl = 0;
    outCombined = 0;
    inCombined = 0;
    pBackward = 0;
    precomputedFFT_spcFluxOverErr2 = 0;
    precomputedFFT_spcOneOverErr2 = 0;
}

COperatorTemplateFittingLog::~COperatorTemplateFittingLog()
{
    freeFFTPlans();
}

Int32 COperatorTemplateFittingLog::EstimateXtYSlow(const std::vector<Float64>& X,
                                                   const std::vector<Float64>& Y,
                                                   UInt32 nShifts,
                                                   std::vector<Float64> &XtY)
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
Int32 COperatorTemplateFittingLog::EstimateMtMFast(const std::vector<Float64> &X,
                                                   const std::vector<Float64> &Y,
                                                   UInt32 nShifts,
                                                   std::vector<Float64> &XtY)
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

Int32 COperatorTemplateFittingLog::EstimateXtY(const std::vector<Float64> &X,
                                               const std::vector<Float64> &Y,
                                               UInt32 nshifts,
                                               std::vector<Float64> &XtY,
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
        Log.LogError("  Operator-TemplateFittingLog: InitFFT: Unable to allocate inSpc");
        return -1;
    }
    if (outSpc == 0)
    {
        Log.LogError("  Operator-TemplateFittingLog: InitFFT: Unable to allocate outSpc");
        return -1;
    }

    inTpl_padded = (Float64 *)fftw_malloc(sizeof(Float64) * nPadded);
    inTpl = (Float64 *)fftw_malloc(sizeof(Float64) * nPadded);
    outTpl = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nPadded);
    pTpl = fftw_plan_dft_r2c_1d(nPadded, inTpl, outTpl, FFTW_ESTIMATE);
    if (inTpl_padded == 0)
    {
        Log.LogError("  Operator-TemplateFittingLog: InitFFT: Unable to allocate inTpl_padded");
        return -1;
    }
    if (inTpl == 0)
    {
        Log.LogError("  Operator-TemplateFittingLog: InitFFT: Unable to allocate inTpl");
        return -1;
    }
    if (outTpl == 0)
    {
        Log.LogError("  Operator-TemplateFittingLog: InitFFT: Unable to allocate outTpl");
        return -1;
    }

    outCombined = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nPadded);
    inCombined = (Float64 *)fftw_malloc(sizeof(Float64) * nPadded);
    pBackward = fftw_plan_dft_c2r_1d(nPadded, outCombined, inCombined, FFTW_ESTIMATE);
    if (outCombined == 0)
    {
        Log.LogError("  Operator-TemplateFittingLog: InitFFT: Unable to allocate outCombined");
        return -1;
    }
    if (inCombined == 0)
    {
        Log.LogError("  Operator-TemplateFittingLog: InitFFT: Unable to allocate inCombined");
        return -1;
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
Int32 COperatorTemplateFittingLog::FitAllz(const TFloat64Range &lambdaRange,
                                           std::shared_ptr<CTemplateFittingResult> result,
                                           std::vector<Int32> igmMeiksinCoeffs,
                                           std::vector<Int32> ismEbmvCoeffs,
                                           CMask spcMaskAdditional,
                                           CPriorHelper::TPriorZEList logpriorze)
{
    bool verboseLogFitAllz = true;

    const CSpectrumFluxAxis &spectrumRebinedFluxAxis =
        m_spectrumRebinedLog.GetFluxAxis();
    const TAxisSampleList & error = spectrumRebinedFluxAxis.GetError().GetSamplesVector();
    const CSpectrumSpectralAxis &spectrumRebinedSpectralAxis =
        m_spectrumRebinedLog.GetSpectralAxis();
    const CSpectrumFluxAxis &tplRebinedFluxAxis = m_templateRebinedLog.GetFluxAxis();
    const CSpectrumSpectralAxis &tplRebinedSpectralAxis =
        m_templateRebinedLog.GetSpectralAxis();

    bool enableIGM = true;
    if (igmMeiksinCoeffs.size() == 0)
    {
        enableIGM = false;
    }

    // prepare list of redshifts that need a full lbdarange lst-square
    // calculation
    TInt32RangeList izrangelist;
    std::vector<Int32> zindexesFullLstSquare;
    if (enableIGM && result->Redshifts.size() > 1)
    {
        zindexesFullLstSquare.push_back(
            0); // first index is always a mandatory full Lstsq Calculation case
        std::vector<Float64> zlistsegments =
            m_igmCorrectionMeiksin->GetSegmentsStartRedshiftList();
        for (Int32 k = 0; k < zlistsegments.size(); k++)
        {
            for (Int32 i = 0; i < result->Redshifts.size() - 1; i++)
            {
                if (zlistsegments[k] >= result->Redshifts[i] &&
                    zlistsegments[k] < result->Redshifts[i + 1])
                {
                    zindexesFullLstSquare.push_back(i);
                    // Log.LogInfo("  Operator-TemplateFittingLog:
                    // zindexesFullLstSquare: index found for zlistsegments_k =
                    // %f, redshift_i = %f, i = %d", zlistsegments[k], i);
                    break; // index found, go to next zlistsegments item
                }
            }
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
            UInt32 izmax = zindexesFullLstSquare[k] + 1;
            izrangelist.push_back(TInt32Range(izmin, izmax));
            if (k < zindexesFullLstSquare.size() - 1)
            {
                izmin = zindexesFullLstSquare[k] + 1;
            }
        }
        // add the last range until result->Redshifts.size()-1 if necessary
        if (izrangelist.size() == 0 ||
            izrangelist[izrangelist.size() - 1].GetEnd() <
                result->Redshifts.size() - 1)
        {
            if (izrangelist.size() > 0)
            {
                izmin = zindexesFullLstSquare[zindexesFullLstSquare.size() - 1] + 1;
            }
            UInt32 izmax = result->Redshifts.size() - 1;
            izrangelist.push_back(TInt32Range(izmin, izmax));
        }

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

    const TAxisSampleList & spectrumRebinedLambda = spectrumRebinedSpectralAxis.GetSamplesVector();
    const TAxisSampleList & spectrumRebinedFluxRaw = spectrumRebinedFluxAxis.GetSamplesVector();
    UInt32 nSpc = spectrumRebinedSpectralAxis.GetSamplesCount();
    if (verboseLogFitAllz)
    {
        Log.LogDebug("  Operator-TemplateFittingLog: FitAllz: spc lbda min=%f, max=%f",
                     spectrumRebinedLambda[0], spectrumRebinedLambda[nSpc - 1]);
    }

    const TAxisSampleList & tplRebinedLambdaGlobal = tplRebinedSpectralAxis.GetSamplesVector();
    const TAxisSampleList & tplRebinedFluxRawGlobal = tplRebinedFluxAxis.GetSamplesVector();
    TAxisSampleList tplRebinedLambda(tplRebinedLambdaGlobal.size());
    TAxisSampleList tplRebinedFluxRaw(tplRebinedLambdaGlobal.size());

    for (Int32 k = 0; k < nzranges; k++)
    {
        // prepare the zrange-result container
        std::shared_ptr<CTemplateFittingResult> subresult = std::shared_ptr<CTemplateFittingResult>(new CTemplateFittingResult());
        TFloat64Range zrange = TFloat64Range(result->Redshifts[izrangelist[k].GetBegin()],
                                             result->Redshifts[izrangelist[k].GetEnd()]);
        TInt32Range ilbda;
        if (enableIGM && result->Redshifts.size() > 1)
        {
            TFloat64List subRedshifts;
            for (Int32 kzsub = izrangelist[k].GetBegin();
                 kzsub <= izrangelist[k].GetEnd(); kzsub++)
            {
                subRedshifts.push_back(result->Redshifts[kzsub]);
            }
            subresult->Init(subRedshifts.size(),
                            std::max((Int32)ismEbmvCoeffs.size(), 1),
                            std::max((Int32)igmMeiksinCoeffs.size(), 1));
            subresult->Redshifts = subRedshifts;

            // slice the template
            Float64 redshiftStep_toBeDoneDifferently =
                result->Redshifts[izrangelist[k].GetBegin() + 1] -
                result->Redshifts[izrangelist[k].GetBegin()];
            ilbda = FindTplSpectralIndex(spectrumRebinedLambda, tplRebinedLambdaGlobal,
                                         zrange, redshiftStep_toBeDoneDifferently);
        } else
        {
            ilbda =
                TInt32Range(0, tplRebinedSpectralAxis.GetSamplesCount() - 1);
            subresult->Init(result->Redshifts.size(), std::max((Int32)ismEbmvCoeffs.size(), 1),
                            std::max((Int32)igmMeiksinCoeffs.size(), 1));
            subresult->Redshifts = result->Redshifts;
        }
        if (ilbda.GetBegin() == -1 || ilbda.GetEnd() == -1)
        {
            Log.LogError("  Operator-TemplateFittingLog: FitAllz: Can't find tpl indexes");
            return -1;
        }

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
                         0, tplRebinedSpectralAxis.GetSamplesCount()-1);
            Log.LogDebug("  Operator-TemplateFittingLog: FitAllz: tpl lbda "
                         "min*zmax=%f, max*zmin=%f",
                         tplRebinedLambdaGlobal[ilbda.GetBegin()] * (1.0 + zrange.GetEnd()),
                         tplRebinedLambdaGlobal[ilbda.GetEnd()] * (1.0 + zrange.GetBegin()));
        }

        UInt32 nTpl = ilbda.GetEnd() - ilbda.GetBegin() + 1;
        tplRebinedLambda.resize(nTpl);
        tplRebinedFluxRaw.resize(nTpl);
        for (UInt32 j = 0; j < nTpl; j++)
        {
            tplRebinedLambda[j] = tplRebinedLambdaGlobal[j + ilbda.GetBegin()];
            tplRebinedFluxRaw[j] = tplRebinedFluxRawGlobal[j + ilbda.GetBegin()];
        }

        //*
        Float64 dtd = 0.0;
        Float64 inv_err2 = 1.0;
        for (Int32 j = 0; j < nSpc; j++)
        {
            inv_err2 = 1.0 / (error[j] * error[j]);
            dtd += spectrumRebinedFluxRaw[j] * spectrumRebinedFluxRaw[j] * inv_err2;
        }

        FitRangez(spectrumRebinedLambda, spectrumRebinedFluxRaw, error,
                  tplRebinedLambda, tplRebinedFluxRaw, subresult,
                  igmMeiksinCoeffs, ismEbmvCoeffs);
        //*/

        // copy subresults into global results
        for (UInt32 isubz = 0; isubz < subresult->Redshifts.size(); isubz++)
        {
            UInt32 fullResultIdx = isubz + izrangelist[k].GetBegin();
            result->ChiSquare[fullResultIdx] = subresult->ChiSquare[isubz];
            result->FitAmplitude[fullResultIdx] = subresult->FitAmplitude[isubz];
            result->FitAmplitudeError[fullResultIdx] = subresult->FitAmplitudeError[isubz];
            result->FitAmplitudeNegative[fullResultIdx] = subresult->FitAmplitudeNegative[isubz];
            result->FitDtM[fullResultIdx] = subresult->FitDtM[isubz];
            result->FitMtM[fullResultIdx] = subresult->FitMtM[isubz];


            Float64 logprior = 0.;
            if(logpriorze.size()>0)
            {
                bool verbose_priorA = false;
                Int32 kism_best = -1;
                if(subresult->FitDustCoeff[isubz]==-1)
                {
                    kism_best=0;
                }else{
                    Int32 nISM = ismEbmvCoeffs.size();
                    for (Int32 kISM = 0; kISM < nISM; kISM++)
                    {
                        Float64 ebmv = m_ismCorrectionCalzetti->GetEbmvValue(kISM);
                        if(ebmv==subresult->FitDustCoeff[isubz])
                        {
                            kism_best = kISM;
                            break;
                        }
                    }
                }

                CPriorHelper::SPriorTZE &pTZE = logpriorze[fullResultIdx][kism_best];
                logprior += -2.0*pTZE.betaTE*pTZE.logprior_precompTE;
                logprior += -2.0*pTZE.betaA*pTZE.logprior_precompA;
                logprior += -2.0*pTZE.betaZ*pTZE.logprior_precompZ;


                if(pTZE.A_sigma>0.0 && pTZE.A_mean>0.0)
                {
                    //now update the amplitude if there is any constraints from the priors
                    Float64 ampl = result->FitAmplitude[fullResultIdx];
                    Float64 ampl_err = result->FitAmplitudeError[fullResultIdx];
                    Float64 ampl_neg = result->FitAmplitudeNegative[fullResultIdx];
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
                        Log.LogDebug("    TemplateFittingLog: update the amplitude (a_mean=%e, a_sigma=%e)",
                                     pTZE.A_mean,
                                     pTZE.A_sigma);
                        Log.LogDebug("    TemplateFittingLog: update the amplitude (ampl was = %e, updated to %e)",
                                     result->FitAmplitude[fullResultIdx],
                                     ampl);
                    }

                    // check negative amplitude
                    ampl_neg = ampl < -3*ampl_err ? 1 : 0;
                    
                    // force positivity
                    ampl = max(0., ampl);

                    result->FitAmplitude[fullResultIdx] = ampl;
                    result->FitAmplitudeError[fullResultIdx] = ampl_err;
                    result->FitAmplitudeNegative[fullResultIdx] = ampl_neg;
                    result->ChiSquare[fullResultIdx] = dtd + result->FitMtM[fullResultIdx]*ampl*ampl - 2.*ampl*result->FitDtM[fullResultIdx];

                    Float64 logPa = pTZE.betaA*(ampl-pTZE.A_mean)*(ampl-pTZE.A_mean)/(pTZE.A_sigma*pTZE.A_sigma);
                    if(std::isnan(logPa) || logPa!=logPa || std::isinf(logPa))
                    {
                        Log.LogError("    TemplateFittingLog: logPa is NAN (a_mean=%e, a_sigma=%e)",
                                     pTZE.A_mean,
                                     pTZE.A_sigma);
                        throw std::runtime_error("    TemplateFittingLog: logPa is NAN or inf, or invalid");
                    }
                    logprior += logPa;
                }else{
                    if(verbose_priorA)
                    {
                        Log.LogDebug("    TemplateFittingLog: NOT updating the amplitude (a_mean=%e, a_sigma=%e)",
                                     pTZE.A_mean,
                                     pTZE.A_sigma);
                    }
                }
                if(std::isnan(logprior) || logprior!=logprior || std::isinf(logprior))
                {
                    Log.LogError("    TemplateFittingLog: logPa is NAN (a_mean=%e, a_sigma=%e, precompA=%e)",
                                 pTZE.A_mean,
                                 pTZE.A_sigma,
                                 pTZE.logprior_precompA);
                    throw std::runtime_error("    TemplateFittingLog: logPrior is NAN or inf, or invalid");
                }
                result->ChiSquare[fullResultIdx] += logprior;
            }
            result->Overlap[fullResultIdx] = subresult->Overlap[isubz];
            result->LogPrior[fullResultIdx] = logprior;
            result->FitDustCoeff[fullResultIdx] = subresult->FitDustCoeff[isubz];
            result->FitMeiksinIdx[fullResultIdx] = subresult->FitMeiksinIdx[isubz];
            result->Status[fullResultIdx] = subresult->Status[isubz];

            for (Int32 kism = 0;
                 kism < result->ChiSquareIntermediate[fullResultIdx].size();
                 kism++)
            {
                for (Int32 kigm = 0;
                     kigm < result->ChiSquareIntermediate[fullResultIdx][kism].size();
                     kigm++)
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
 * @param igmMeiksinCoeffs
 * @param ismEbmvCoeffs
 * @return
 */
Int32 COperatorTemplateFittingLog::FitRangez(const TAxisSampleList & spectrumRebinedLambda,
                                             const TAxisSampleList & spectrumRebinedFluxRaw,
                                             const TAxisSampleList & error,
                                             const TAxisSampleList & tplRebinedLambda,
                                             const TAxisSampleList & tplRebinedFluxRaw,
                                             std::shared_ptr<CTemplateFittingResult> result,
                                             std::vector<Int32> igmMeiksinCoeffs,
                                             std::vector<Int32> ismEbmvCoeffs)
{
    Int32 nTpl = tplRebinedLambda.size(),
          nSpc = spectrumRebinedLambda.size();

   Float64 redshiftValueMeiksin = result->Redshifts[0];

   if (verboseLogFitFitRangez)
    {
        Log.LogDebug("  Operator-TemplateFittingLog: FitRangez: redshiftValueMeiksin = %f", redshiftValueMeiksin);
        Log.LogDebug("  Operator-TemplateFittingLog: FitRangez: spc[0] = %f", spectrumRebinedLambda[0]);
        Log.LogDebug("  Operator-TemplateFittingLog: FitRangez: spc[max] = %f", spectrumRebinedLambda[nSpc - 1]);
        Log.LogDebug("  Operator-TemplateFittingLog: FitRangez: tpl[0]*zmax = %f", tplRebinedLambda[0] * (1.0 + result->Redshifts[result->Redshifts.size() - 1]));
        Log.LogDebug("  Operator-TemplateFittingLog: FitRangez: tpl[max]*zmin = %f", tplRebinedLambda[nTpl - 1] * (1 + result->Redshifts[0]));
    }

    Int32 nshifts = nTpl - nSpc + 1;
    m_nPaddedSamples = nTpl * 2.0;
    /*
    //next power of two
    Int32 maxnsamples = (Int32)(nTpl*2);
    //Log.LogInfo("  Operator-TemplateFittingLog: FitRangez: maxnsamples = %d",
    maxnsamples); Int32 power = 1; while(power < maxnsamples)
    {
        power*=2;
    }
    m_nPaddedSamples = power;
    //*/

    Log.LogDetail("  Operator-TemplateFittingLog: Now fitting using the FFT on "
                  "nshifts=%d values, for Meiksin redshift=%f",
                  nshifts, redshiftValueMeiksin);

    if (verboseLogFitFitRangez)
    {
        Log.LogDebug("  Operator-TemplateFittingLog: FitRangez: initializing FFT with n = %d points",
                     m_nPaddedSamples);
    }
    InitFFT(m_nPaddedSamples);

    // prepare z array
    std::vector<Float64> z_vect(nshifts, 0.0);
    Int32 z_vect_size = z_vect.size();
    for (Int32 t = 0; t < z_vect_size; t++)
    {
        z_vect[t] = (spectrumRebinedLambda[0] - tplRebinedLambda[t]) /
                    tplRebinedLambda[t];
    }

    // Estimate DtD
    //    std::vector<Float64> dtd_vec;
    //    CSpectrumFluxAxis spcSquareFluxAxis(spectrumRebinedFluxRaw);
    //    for(Int32 k=0; k<spcSquareFluxAxis.GetSamplesCount(); k++)
    //    {
    //        spcSquareFluxAxis[k] =
    //        spectrumRebinedFluxRaw[k]*spectrumRebinedFluxRaw[k];
    //    }
    //    EstimateXtY(OneFluxAxis, spcSquareFluxAxis, dtd_vec);
    Float64 dtd = 0.0;
    Float64 inv_err2 = 1.0;
    for (Int32 j = 0; j < nSpc; j++)
    {
        inv_err2 = 1.0 / (error[j] * error[j]);
        dtd += spectrumRebinedFluxRaw[j] * spectrumRebinedFluxRaw[j] * inv_err2;
    }

    // prepare arrays
    TAxisSampleList spcRebinedFluxOverErr2(nSpc);
    TAxisSampleList oneSpcRebinedFluxOverErr2(nSpc);
    for (Int32 j = 0; j < nSpc; j++)
    {
        inv_err2 = 1.0 / (error[j] * error[j]);
        spcRebinedFluxOverErr2[j] = spectrumRebinedFluxRaw[j] * inv_err2;
        oneSpcRebinedFluxOverErr2[j] = inv_err2;
    }

    if (verboseExportFitRangez)
    {
        // save spcRebinedFluxOverErr2 input xty data
        FILE *f = fopen("loglbda_spcRebinedFluxOverErr2l_dbg.txt", "w+");
        for (Int32 j = 0; j < nSpc; j++)
        {
            fprintf(f, "%f\t%e\n", spectrumRebinedLambda[j],
                    spcRebinedFluxOverErr2[j]);
        }
        fclose(f);
    }

    TAxisSampleList tplRebinedFluxIgm = tplRebinedFluxRaw;
    TAxisSampleList tplRebinedFlux = tplRebinedFluxRaw;
    TAxisSampleList tpl2RebinedFlux(nTpl);

    for (Int32 j = 0; j < nTpl; j++)
    {
        tpl2RebinedFlux[j] = tplRebinedFlux[j] * tplRebinedFlux[j];
    }

    bool enableISM = true;
    Int32 nISM = ismEbmvCoeffs.size();
    if (nISM == 0)
    {
        nISM = 1;
        enableISM = false;
    }
    bool enableIGM = true;
    Int32 nIGM = igmMeiksinCoeffs.size();
    if (nIGM == 0)
    {
        nIGM = 1;
        enableIGM = false;
    }
    // disable IGM if the redshift range and lambda range do not make the IGM
    // wavelength appear
    Float64 lambdaMinTpl = tplRebinedLambda[0];
    Int32 overrideNIGMTobesaved = -1;
    if (lambdaMinTpl > 1216. && nIGM > 1)
    {
        overrideNIGMTobesaved = nIGM;
        nIGM = 1;
        enableIGM = false;
        if (verboseLogFitFitRangez)
        {
            Log.LogDebug("  Operator-TemplateFittingLog: FitRangez: IGM disabled, min-tpl-lbda=%f", lambdaMinTpl);
        }
    }

    // prepare best fit data buffer
    TFloat64List bestChi2(nshifts, DBL_MAX);
    TFloat64List bestFitAmp(nshifts, -1.0);
    TFloat64List bestFitAmpErr(nshifts, -1.0);
    TBoolList    bestFitAmpNeg(nshifts);
    TFloat64List bestFitDtm(nshifts, -1.0);
    TFloat64List bestFitMtm(nshifts, -1.0);
    TFloat64List bestISMCoeff(nshifts, -1.0);
    TFloat64List bestIGMIdx(nshifts, -1.0);

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

    Int32 errorWhileFitting = 0;
    //#pragma omp parallel for
    for (Int32 kIGM = 0; kIGM < nIGM; kIGM++)
    {

        if (verboseLogFitFitRangez && enableIGM)
        {
            Log.LogDebug("  Operator-TemplateFittingLog: FitRangez: IGM index=%d", kIGM);
        }

        if (enableIGM)
        {

            Int32 meiksinIdx = igmMeiksinCoeffs[kIGM];

            for (Int32 j = 0; j < nTpl; j++)
            {
                Float64 restLambda = tplRebinedLambda[j];
                Float64 coeffIGM = m_igmCorrectionMeiksin->getCoeff(meiksinIdx, redshiftValueMeiksin, restLambda);

                tplRebinedFluxIgm[j] = tplRebinedFluxRaw[j] * coeffIGM;
                tplRebinedFlux[j] = tplRebinedFluxIgm[j];
                tpl2RebinedFlux[j] = tplRebinedFlux[j] * tplRebinedFlux[j];
            }

        } else
        {
            for (Int32 j = 0; j < nTpl; j++)
            {
                tplRebinedFluxIgm[j] = tplRebinedFluxRaw[j];
                tplRebinedFlux[j] = tplRebinedFluxIgm[j];
                tpl2RebinedFlux[j] = tplRebinedFlux[j] * tplRebinedFlux[j];
            }
        }

        for (Int32 kISM = 0; kISM < nISM; kISM++)
        {
            if (verboseLogFitFitRangez && enableISM)
            {
                Log.LogDebug("  Operator-TemplateFittingLog: FitRangez: ISM index =%d", kISM);
            }

            if (enableISM)
            {
                Int32 kDustCalzetti = ismEbmvCoeffs[kISM];

                // correct tplRebinedFlux
                // correct tpl2RebinedFlux
                Float64 restLambda;
                Float64 ebmvDustCoeff = 1.0;
                for (Int32 j = 0; j < nTpl; j++)
                {
                    // apply ism dust correction from Calzetti
                    restLambda = tplRebinedLambda[j];
                    ebmvDustCoeff = m_ismCorrectionCalzetti->getDustCoeff(kDustCalzetti, restLambda);

                    tplRebinedFlux[j] = tplRebinedFluxIgm[j] * ebmvDustCoeff;

                    tpl2RebinedFlux[j] = tplRebinedFlux[j] * tplRebinedFlux[j];
                }
            }

            if (verboseExportFitRangez_model)
            {
                if ((enableISM && exportISMIdx == ismEbmvCoeffs[kISM]) ||
                    (!enableISM))
                {
                    if ((enableIGM && exportIGMIdx == igmMeiksinCoeffs[kIGM]) ||
                        (!enableIGM))
                    {

                        // save chi2 data
                        FILE *f = fopen("loglbda_model_dbg.txt", "w+");
                        for (Int32 j = 0; j < nTpl; j++)
                        {
                            fprintf(f, "%f\t%e\n", tplRebinedLambda[j], tplRebinedFlux[j]);
                        }
                        fclose(f);
                    }
                }
            }

            // Estimate DtM
            std::vector<Float64> dtm_vec;
            //*
            EstimateXtY(spcRebinedFluxOverErr2, tplRebinedFlux, nshifts, dtm_vec,
                        //-1);
                        0);
            //*/
            // EstimateXtYSlow(spcRebinedFluxOverErr2, tplRebinedFlux, nSpc,
            // nshifts, dtm_vec);
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

            // Estimate MtM
            std::vector<Float64> mtm_vec;
            //*
            EstimateXtY(oneSpcRebinedFluxOverErr2, tpl2RebinedFlux, nshifts, mtm_vec,
                        //-1);
                        1);
            //*/
            // EstimateXtYSlow(oneSpcRebinedFluxOverErr2, tpl2RebinedFlux, nSpc,
            // nshifts, mtm_vec); EstimateMtMFast(oneSpcRebinedFluxOverErr2,
            // tpl2RebinedFlux, nSpc, nshifts, mtm_vec);

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

            // return -1;

            // Estimate Chi2
            if (mtm_vec.size() != dtm_vec_size)
            {
                errorWhileFitting = 1;
                Log.LogError("  Operator-TemplateFittingLog: FitRangez: xty vectors sizes don't match: dtm size = %d, mtm size = %d", dtm_vec_size, mtm_vec.size());
                // return 2;
                continue;
            }
            //Log.LogDetail("  Operator-TemplateFittingLog: FitRangez: kISM = %d, kIGM = %d", kISM, kIGM);
            std::vector<Float64> chi2(dtm_vec_size, DBL_MAX);
            std::vector<Float64> amp(dtm_vec_size, DBL_MAX);
            TBoolList amp_neg(dtm_vec_size);
            std::vector<Float64> amp_err(dtm_vec_size, DBL_MAX);
            for (Int32 k = 0; k < dtm_vec_size; k++)
            {
                if (mtm_vec[k] == 0.0)
                {
                    amp[k] = 0.0;
                    amp_err[k] = 0.0;
                    amp_neg[k] = 0;
                    chi2[k] = dtd;
                } else
                {
                    amp[k] = dtm_vec[k] / mtm_vec[k];
                    amp_err[k] = sqrt(1./mtm_vec[k]);
                    amp_neg[k] = amp[k] < -3*amp_err[k] ? 1 : 0;
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
                    for (Int32 koigm = 1; koigm < overrideNIGMTobesaved;
                         koigm++)
                    {
                        intermediateChi2[k][kISM][koigm] = chi2[k];
                    }
                }

                if (bestChi2[k] > chi2[k])
                {
                    bestChi2[k] = chi2[k];
                    bestFitAmp[k] = amp[k];
                    bestFitAmpErr[k] = amp_err[k];
                    bestFitAmpNeg[k] = amp_neg[k];
                    bestFitDtm[k] = dtm_vec[k];
                    bestFitMtm[k] = mtm_vec[k];
                    bestISMCoeff[k] = enableISM ? m_ismCorrectionCalzetti->GetEbmvValue(ismEbmvCoeffs[kISM]) : -1;
                    bestIGMIdx[k] = enableIGM ? igmMeiksinCoeffs[kIGM] : -1;
                }
            }

            if (verboseExportFitRangez)
            {
                Log.LogDebug("  Operator-TemplateFittingLog: FitRangez: spc lbda 0 =%f", spectrumRebinedLambda[0]);
                Log.LogDebug("  Operator-TemplateFittingLog: FitRangez: tpl lbda 0 =%f", tplRebinedLambda[0]);
                Float64 z_O = (spectrumRebinedLambda[0] - tplRebinedLambda[0]) / tplRebinedLambda[0];
                Log.LogDebug("  Operator-TemplateFittingLog: FitRangez: z 0 =%f", z_O);

                //        Float64 logstep =
                //        log(spectrumRebinedLambda[1])-log(spectrumRebinedLambda[0]);
                //        Float64 logInitialStep =
                //        log(spectrumRebinedLambda[0])-log(tplRebinedLambda[0]);
                //        Log.LogInfo("  Operator-TemplateFittingLog: FitRangez:
                //        logstep=%f", logstep); Log.LogInfo("
                //        Operator-TemplateFittingLog: FitRangez: step=%f A",
                //        exp(logstep)); Log.LogInfo("  Operator-TemplateFittingLog:
                //        FitRangez: logInitialStep=%f", logInitialStep);
                //        Log.LogInfo("  Operator-TemplateFittingLog: FitRangez:
                //        InitialStep=%f A", exp(logInitialStep));
                //        std::vector<Float64> log1pz(dtm_vec.size(), 0.0);
                //        for( Int32 t=0;t<log1pz.size();t++)
                //        {
                //            log1pz[t] = t*logstep + logInitialStep;
                //        }

                // save chi2 data
                FILE *f_chi2 = fopen("loglbda_chi2output_dbg.txt", "w+");
                for (Int32 t = 0; t < chi2.size(); t++)
                {
                    fprintf(f_chi2, "%f\t%e\n",
                            z_vect[t] /*(Float64)(exp(log1pz[t])-1)*/, chi2[t]);
                }
                fclose(f_chi2);
            }
        }
    }

    if (errorWhileFitting != 0)
    {
      freeFFTPlans();
      throw runtime_error("Error while fitting");
        // return -2;
    }

    // interpolating on the initial z grid
    std::reverse(z_vect.begin(), z_vect.end());
    std::reverse(bestChi2.begin(), bestChi2.end());
    std::reverse(bestFitAmp.begin(), bestFitAmp.end());
    std::reverse(bestFitAmpErr.begin(), bestFitAmpErr.end());
    std::reverse(bestFitAmpNeg.begin(), bestFitAmpNeg.end());
    std::reverse(bestFitDtm.begin(), bestFitDtm.end());
    std::reverse(bestFitMtm.begin(), bestFitMtm.end());
    std::reverse(bestISMCoeff.begin(), bestISMCoeff.end());
    std::reverse(bestIGMIdx.begin(), bestIGMIdx.end());


    Int32 k = 0;
    Int32 klow = 0;
    for (Int32 iz = 0; iz < result->Redshifts.size(); iz++)
    {
        //* //NGP
        // Log.LogInfo("  Operator-TemplateFittingLog: FitRangez: interpolating
        // gsl-bsearch z result, , ztgt=%f, zcalc_min=%f, zcalc_max=%f",
        // result->Redshifts[iz], zreversed_array[0],
        // zreversed_array[z_vect_size-1]);
        k = gsl_interp_bsearch(&z_vect.at(0), result->Redshifts[iz], klow, z_vect_size - 1);
        klow = k;

        /*
        if(result->Redshifts[iz]==0.0)
        {
            Log.LogInfo("  Operator-TemplateFittingLog: FitRangez: interpolating z
        result, kshift=%f", k); Log.LogInfo("  Operator-TemplateFittingLog: FitRangez:
        interpolating z result, zcalc=%f, kfound=%d, zfound=%f",
        result->Redshifts[iz], k, z_vect[z_vect.size()-1-k]); Log.LogInfo("
        Operator-TemplateFittingLog: FitRangez: interpolating z result, bestChi2=%f",
        bestChi2[z_vect.size()-1-k]);
        }
        // closest value
        result->ChiSquare[iz] = bestChi2[z_vect.size()-1-k];
        //*/

        result->Overlap[iz] = 1.0;
        result->FitAmplitude[iz] = bestFitAmp[k];
        result->FitAmplitudeError[iz] = bestFitAmpErr[k];
        result->FitAmplitudeNegative[iz] = bestFitAmpNeg[k];
        result->FitDtM[iz] = bestFitDtm[k];
        result->FitMtM[iz] = bestFitMtm[k];
        result->FitDustCoeff[iz] = bestISMCoeff[k];
        result->FitMeiksinIdx[iz] = bestIGMIdx[k];
        result->Status[iz] = nStatus_OK;

        //*/
    }

    //*
    Log.LogDetail("  Operator-TemplateFittingLog: FitRangez: interpolating (lin) z result from n=%d (min=%f, max=%f) to n=%d (min=%f, max=%f)",
            z_vect_size,
            z_vect.front(),
            z_vect.back(),
            result->Redshifts.size(),
            result->Redshifts.front(),
            result->Redshifts.back());
    InterpolateResult(bestChi2,
                      z_vect,
                      result->Redshifts,
                      result->ChiSquare, DBL_MAX);
    //*/

    //*
    // Interpolating intermediate chisquare results
    TFloat64List intermChi2BufferReversed_array(intermediateChi2.size());
    TFloat64List intermChi2BufferRebinned_array(
        result->Redshifts.size(), boost::numeric::bounds<float>::highest());
    for (Int32 kism = 0; kism < nISM; kism++)
    {
        for (Int32 kigm = 0; kigm < nIGMFinal; kigm++)
        {
            for (Int32 t = 0; t < z_vect_size; t++)
            {
                intermChi2BufferReversed_array[t] =
                    intermediateChi2[z_vect_size - 1 - t][kism][kigm];
            }
            InterpolateResult(intermChi2BufferReversed_array, z_vect,
                              result->Redshifts,
                              intermChi2BufferRebinned_array, DBL_MAX);
            for (Int32 t = 0; t < result->Redshifts.size(); t++)
            {
                result->ChiSquareIntermediate[t][kism][kigm] =
                    intermChi2BufferRebinned_array[t];
            }
        }
    }

    freeFFTPlans();
    return 0;
}

Int32 COperatorTemplateFittingLog::InterpolateResult(const std::vector<Float64>& in,
                                                     std::vector<Float64>& inGrid,
                                                     const std::vector<Float64>& tgtGrid,
                                                     std::vector<Float64> &out,
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
            Log.LogError("Error while interpolating loglambda chi2 result : xrebin(%f) < x[0](%f)", tgtGrid[0], inGrid[0]);
            throw runtime_error("Error while interpolating loglambda chi2 result.");
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
            Log.LogError("Error while interpolating loglambda chi2 result : xrebin(%f) > x[n-1](%f)", tgtGrid[tgtn-1], inGrid[n-1]);
            throw runtime_error("Error while interpolating loglambda chi2 result.");
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
          Log.LogError("   Operator-TemplateFittingLog: InterpolateError: GSL code = %d, %s ; see file: %s at line: %d", status, gsl_strerror(status), __FILENAME__, __LINE__);
          throw std::runtime_error("GSL Error during the interpolation evaluation");
        }
        // out[j] = gsl_spline_eval (spline, Xrebin, accelerator); //spline
        // Log.LogInfo("  Operator-TemplateFittingLog: FitAllz: interpolating
        // gsl-spline z result, , ztgt=%f, rebinY=%f", tgtGrid[j], out[j]);
    }

    gsl_set_error_handler(gsl_error_handler_old);

    gsl_interp_free(interpolation);
    // gsl_spline_free (spline);
    gsl_interp_accel_free(accelerator);

    return 0;
}

TInt32Range COperatorTemplateFittingLog::FindTplSpectralIndex(
    const TAxisSampleList & spcLambda, const TAxisSampleList & tplLambda,
        TFloat64Range redshiftrange, Float64 redshiftStep)
{
    Float64 redshiftMargin =
        10.0 *
        redshiftStep; // this redshiftstep margin might be useless now that
                      // there are spclambdamargin and tpllambdamargin
    Float64 spcLambdaMargin = abs(spcLambda[1] - spcLambda[0]);
    Float64 tplLambdaMargin = abs(tplLambda[1] - tplLambda[0]);

    const UInt32 nTpl = tplLambda.size(),
                 nSpc = spcLambda.size();

    //Log.LogInfo("  Operator-TemplateFittingLog: FindTplSpectralIndex: redshiftMargin = %f", redshiftMargin);
    //Log.LogInfo("  Operator-TemplateFittingLog: FindTplSpectralIndex: spcLambdaMargin = %f", spcLambdaMargin);
    //Log.LogInfo("  Operator-TemplateFittingLog: FindTplSpectralIndex: tplLambdaMargin = %f", tplLambdaMargin);

    UInt32 ilbdamin = 0;
    //Log.LogInfo("  Operator-TemplateFittingLog: FindTplSpectralIndex: ilbdamin computation starts");
    for (UInt32 k = 0; k < nTpl; k++)
    {
        Float64 _spcLambda = spcLambda[0] - spcLambdaMargin;
        Float64 _tplLambda = tplLambda[k] + tplLambdaMargin;
        Float64 _z = (_spcLambda - _tplLambda) / _tplLambda;

	//Log.LogInfo("  Operator-TemplateFittingLog: FindTplSpectralIndex: _spcLambda = %f, spcLambda[%d]=%f", _spcLambda, 0, spcLambda[0]);
	//Log.LogInfo("  Operator-TemplateFittingLog: FindTplSpectralIndex: _tplLambda = %f, tplLambda[%d]=%f", _tplLambda, k, tplLambda[k]);
	//Log.LogInfo("  Operator-TemplateFittingLog: FindTplSpectralIndex: _z = %f", _z);

        if (_z > redshiftrange.GetEnd() + redshiftMargin)
        {
            ilbdamin = k;
        } else
        {
	    //Log.LogInfo("  Operator-TemplateFittingLog: FindTplSpectralIndex: ilbdamin computation breaks");
            break;
        }
    }

    UInt32 ilbdamax = nTpl - 1;
    //Log.LogInfo("  Operator-TemplateFittingLog: FindTplSpectralIndex: ilbdamax computation starts");
    for (UInt32 k = nTpl - 1; k > 0; k--)
    {
        Float64 _spcLambda = spcLambda[nSpc - 1] + spcLambdaMargin;
        Float64 _tplLambda = tplLambda[k] - tplLambdaMargin;
        Float64 _z = (_spcLambda - _tplLambda) / _tplLambda;

        //Log.LogInfo("  Operator-TemplateFittingLog: FindTplSpectralIndex: _spcLambda = %f, spcLambda[%d]=%f", _spcLambda, nSpc - 1, spcLambda[nSpc - 1]);
        //Log.LogInfo("  Operator-TemplateFittingLog: FindTplSpectralIndex: _tplLambda = %f, tplLambda[%d]=%f", _tplLambda, k, tplLambda[k]);
        //Log.LogInfo("  Operator-TemplateFittingLog: FindTplSpectralIndex: _z = %f", _z);

        if (_z < redshiftrange.GetBegin() - redshiftMargin)
        {
            ilbdamax = k;
        } else
        {
	    //Log.LogInfo("  Operator-TemplateFittingLog: FindTplSpectralIndex: ilbdamax computation breaks");
            break;
        }
    }

    if (ilbdamin < 0 || ilbdamin > nTpl - 1)
    {
        Log.LogError("  Operator-TemplateFittingLog: Problem with tpl indexes for "
                     "zranges, found lbdamin=%d",
                     ilbdamin);
        return TInt32Range(-1, -1);
    }
    if (ilbdamax < 0 || ilbdamax > nTpl - 1)
    {
        Log.LogError("  Operator-TemplateFittingLog: Problem with tpl indexes for "
                     "zranges, found ilbdamax=%d",
                     ilbdamax);
        return TInt32Range(-1, -1);
    }
    if (ilbdamin > ilbdamax)
    {
        Log.LogError("  Operator-TemplateFittingLog: Problem with tpl indexes for "
                     "zranges, found ilbdamin=%d > ilbdamax=%d",
                     ilbdamin, ilbdamax);
        return TInt32Range(-1, -1);
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
 **/
std::shared_ptr<COperatorResult> COperatorTemplateFittingLog::Compute(const CSpectrum &spectrum,
        const CTemplate &tpl,
        const TFloat64Range &lambdaRange,
        const TFloat64List &redshifts,
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
    Log.LogDetail("  Operator-TemplateFittingLog: starting computation for template: %s", tpl.GetName().c_str());

    if(redshifts.size()<2){
        Log.LogError("       Operator-TemplateFittingLog::Compute: Cannot compute on a redshift array %d <2", redshifts.size());
        throw runtime_error("Operator-TemplateFittingLog::Compute: Cannot compute on a redshift array <2");
    }
    if ((opt_dustFitting==-10 || opt_dustFitting>-1) && m_ismCorrectionCalzetti->calzettiInitFailed)
    {
        Log.LogError("  Operator-TemplateFittingLog: no calzetti calib. file loaded... aborting!");
        return NULL;
    }
    if( opt_dustFitting>-1 && opt_dustFitting>m_ismCorrectionCalzetti->GetNPrecomputedDustCoeffs()-1)
    {
        Log.LogError("  Operator-TemplateFittingLog: calzetti index overflow (opt=%d, while NPrecomputedDustCoeffs=%d)... aborting!",
                     opt_dustFitting,
                     m_ismCorrectionCalzetti->GetNPrecomputedDustCoeffs());
        return NULL;
    }

    if (opt_extinction && m_igmCorrectionMeiksin->meiksinInitFailed)
    {
        Log.LogError("  Operator-TemplateFittingLog: no meiksin calib. file loaded... aborting!");
        return NULL;
    }

    if (spectrum.GetSpectralAxis().IsInLinearScale() == false ||
        tpl.GetSpectralAxis().IsInLinearScale() == false)
    {
        Log.LogError("  Operator-TemplateFittingLog: input spectrum or template are not in log scale (ignored)");
        // return NULL;
    }

    // assuming that lambdarange is strictly included in the spectrum spectral
    // range
    if (std::max(spectrum.GetSpectralAxis()[0], lambdaRange.GetBegin()) !=
        lambdaRange.GetBegin())
    {
        Log.LogError("  Operator-TemplateFittingLog: lambdarange is not strictly included in the spectrum spectral range (min lambdarange)");
    }
    if (std::min(spectrum.GetSpectralAxis()[spectrum.GetSpectralAxis().GetSamplesCount() - 1],
                 lambdaRange.GetEnd()) != lambdaRange.GetEnd())
    {
        Log.LogError("  Operator-TemplateFittingLog: lambdarange is not strictly included in the spectrum spectral range (max lambdarange)");
    }


    // sort the redshifts and keep track of the indexes
    TFloat64List sortedRedshifts;
    TFloat64List sortedIndexes; // used for the correspondence between input redshifts
                                // (list) and input additionalMasks (list)
    // This is a vector of {value,index} pairs
    vector<pair<Float64, Int32>> vp;
    vp.reserve(redshifts.size());
    for (Int32 i = 0; i < redshifts.size(); i++)
    {
        vp.push_back(make_pair(redshifts[i], i));
    }
    std::sort(vp.begin(), vp.end());
    for (Int32 i = 0; i < vp.size(); i++)
    {
        sortedRedshifts.push_back(vp[i].first);
        sortedIndexes.push_back(vp[i].second);
    }

    /****************************************************************
     *                                                   
     *  Determine the loglambda step
     *   it will be used for: Spectrum, template and redshift grid
     *
     ****************************************************************/
    // Create/Retrieve the spectrum log-lambda spectral axis
    Int32 lbdaMinIdx = spectrum.GetSpectralAxis().GetIndexAtWaveLength(lambdaRange.GetBegin());
    Int32 lbdaMaxIdx = spectrum.GetSpectralAxis().GetIndexAtWaveLength(lambdaRange.GetEnd());
    if (lbdaMinIdx < 0 || lbdaMaxIdx < 0 || lbdaMaxIdx <= lbdaMinIdx)
    {
        Log.LogError("  Operator-TemplateFittingLog: problem while searching for lambdarange observed spectral indexes");
    }

    Float64 loglbdaStep_fromOriSpc = log(spectrum.GetSpectralAxis()[lbdaMaxIdx]) -
                          log(spectrum.GetSpectralAxis()[lbdaMaxIdx - 1]);

    //find loglbdaStep from tgt redshift grid
    Float64 tgtDzOnepzMin=DBL_MAX;
    for (Int32 k = 1; k < sortedRedshifts.size(); k++)
    {
        //Log.LogDetail("  Operator-TemplateFittingLog: Log-Rebin: sortedRedshifts[k] = %f", sortedRedshifts[k]);
        Float64 tgtDzOnepz=(sortedRedshifts[k]-sortedRedshifts[k-1])/(1+sortedRedshifts[k-1]);
        if(tgtDzOnepzMin>tgtDzOnepz)
        {
            tgtDzOnepzMin = tgtDzOnepz;
        }
    }
    Log.LogDetail("  Operator-TemplateFittingLog: Log-Rebin: tgtDzOnepzMin = %f", tgtDzOnepzMin);
    Float64 loglbdaStep_fromTgtZgrid = log(1.0+tgtDzOnepzMin);
    Log.LogDetail("  Operator-TemplateFittingLog: Log-Rebin: loglbdaStep_fromOriSpc = %f", loglbdaStep_fromOriSpc);
    Log.LogDetail("  Operator-TemplateFittingLog: Log-Rebin: loglbdaStep_fromTgtZgrid = %f", loglbdaStep_fromTgtZgrid);

    
    //Float64 loglbdaStep = std::max(loglbdaStep_fromOriSpc, loglbdaStep_fromTgtZgrid);
    //Force loglbdaStep to loglambdaStep_fromTgtZgrid
    Float64 loglbdaStep = loglbdaStep_fromTgtZgrid;
    Log.LogDetail("  Operator-TemplateFittingLog: Log-Rebin: loglbdaStep = %f", loglbdaStep);
    // check that the overlap is >1. for all sortedRedshifts
    if (overlapThreshold < 1.0)
    {
        Log.LogError("  Operator-TemplateFittingLog: overlap threshold can't be "
                     "lower than 1.0");
        return NULL;
    }

    // compute the effective zrange of the new redshift grid
    // set the min to the initial min
    // set the max to an interger number of log(z+1) steps 
    Float64 zmin_new = sortedRedshifts[0];
    Float64 zmax_new = sortedRedshifts[sortedRedshifts.size() - 1];
    {
        Float64 log_zmin_new_p1 = log(zmin_new + 1.);
        Float64 log_zmax_new_p1 = log(zmax_new + 1.);
        Int32 nb_z = Int32( ceil((log_zmax_new_p1 - log_zmin_new_p1)/loglbdaStep) );
        zmax_new = exp(log_zmin_new_p1 + nb_z*loglbdaStep) - 1.;
    }

    // Display the template coverage
    Log.LogDetail("  Operator-TemplateFittingLog: Log-Rebin: tpl[min] = %f, tpl[min]*(1 + zmax_new) = %f", tpl.GetSpectralAxis()[0], tpl.GetSpectralAxis()[0] * (1. + zmax_new));
    Log.LogDetail("  Operator-TemplateFittingLog: Log-Rebin: tpl[max] = %f, tpl[max]*(1 + zmin_new) = %f", tpl.GetSpectralAxis()[tpl.GetSpectralAxis().GetSamplesCount() - 1], tpl.GetSpectralAxis()[tpl.GetSpectralAxis().GetSamplesCount() - 1] * (1. + zmin_new));

    // Check the template coverage wrt the new effective redshift range
    Bool overlapFull = true;
    if (lambdaRange.GetBegin() < tpl.GetSpectralAxis()[0] * (1. + zmax_new))
    {
        overlapFull = false;
    }
    if (lambdaRange.GetEnd() >
        tpl.GetSpectralAxis()[tpl.GetSpectralAxis().GetSamplesCount() - 1] * (1. + zmin_new))
    {
        overlapFull = false;
    }
    if (!overlapFull)
    {
        Log.LogError("  Operator-TemplateFittingLog: overlap found to be lower than 1.0 for this redshift range");
        Log.LogError("  Operator-TemplateFittingLog: for zmin=%f, tpl.lbdamax is %f (should be >%f)",
                     zmin_new,
                     tpl.GetSpectralAxis()[tpl.GetSpectralAxis().GetSamplesCount() - 1],
                     lambdaRange.GetEnd() / (1 + zmin_new));
        Log.LogError("  Operator-TemplateFittingLog: for zmax=%f, tpl.lbdamin is %f (should be <%f)",
                     zmax_new,
                     tpl.GetSpectralAxis()[0],
                     lambdaRange.GetBegin() /(1 + zmax_new));
        return NULL;
    }

    /****************************************************************
     *                                                   
     *  Rebin the spectrum
     *  or check that the raw spectrum grid is in log-regular grid 
     *  with expected step
     *
     ****************************************************************/
    Float64 loglbdamin = log(lambdaRange.GetBegin());
    Float64 loglbdamax = log(lambdaRange.GetEnd());
    Int32 loglbdaCount = (Int32) floor((loglbdamax - loglbdamin) / loglbdaStep + 1);

    if(loglbdaCount<2){
        Log.LogError("       Operator-TemplateFittingLog::Compute: loglbdaCount = %d <2", loglbdaCount);
        throw runtime_error("Operator-TemplateFittingLog::Compute: loglbdaCount <2. Abort!");
    }
    // Allocate the Log-rebined spectrum and mask
    CSpectrumFluxAxis &spectrumRebinedFluxAxis = m_spectrumRebinedLog.GetFluxAxis();
    CSpectrumSpectralAxis &spectrumRebinedSpectralAxis =  m_spectrumRebinedLog.GetSpectralAxis();
    spectrumRebinedFluxAxis.SetSize(loglbdaCount);
    spectrumRebinedSpectralAxis.SetSize(loglbdaCount);
    m_mskRebinedLog.SetSize(loglbdaCount);

    if (m_opt_spcrebin)
    {

        Log.LogDetail("  Operator-TemplateFittingLog: Log-regular lambda resampling START");

        if (verboseLogRebin)
        {
            // Log.LogInfo("  Operator-TemplateFittingLog: Log-Rebin: zref = %f",
            // zRef); Log.LogInfo("  Operator-TemplateFittingLog: Log-Rebin: zStep =
            // %f", zStep);
            Log.LogDebug("  Operator-TemplateFittingLog: Log-Rebin: loglbdaStep = %f", loglbdaStep);
            Log.LogDebug("  Operator-TemplateFittingLog: Log-Rebin: loglbdamin=%f : loglbdamax=%f", loglbdamin, loglbdamax);
            Log.LogDebug("  Operator-TemplateFittingLog: Log-Rebin: lbdamin=%f : lbdamax=%f", exp(loglbdamin), exp(loglbdamax));
        }
        Log.LogDetail("  Operator-TemplateFittingLog: ORIGINAL grid spectrum count = %d", lbdaMaxIdx - lbdaMinIdx + 1);
        Log.LogDetail("  Operator-TemplateFittingLog: Log-Rebin: loglbdaCount = %d", loglbdaCount);

        CSpectrumSpectralAxis targetSpectralAxis; //define it earlier to use it in the different contexts
        targetSpectralAxis.SetSize(loglbdaCount);

        for (Int32 k = 0; k < loglbdaCount; k++)
        {
            targetSpectralAxis[k] = exp(loglbdamin + k * loglbdaStep);
        }
        // precision problem due to exp/log - beginning
        Float64 delta = targetSpectralAxis[0] - spectrum.GetSpectralAxis()[0];
        Float64 step = targetSpectralAxis[1] - targetSpectralAxis[0];
        if (delta < 0.0 && abs(delta) < step * 1e-5)
        {
            targetSpectralAxis[0] = spectrum.GetSpectralAxis()[0];
        }
        // precision problem due to exp/log - end
        delta = targetSpectralAxis[loglbdaCount - 1] -
                spectrum.GetSpectralAxis() [spectrum.GetSpectralAxis().GetSamplesCount() - 1];
        step = targetSpectralAxis[loglbdaCount - 1] -
               targetSpectralAxis[loglbdaCount - 2];
        if (delta > 0.0 && abs(delta) < step * 1e-5)
        {
            targetSpectralAxis[loglbdaCount - 1] =
                spectrum.GetSpectralAxis() [spectrum.GetSpectralAxis().GetSamplesCount() - 1];
        }

        if (verboseExportLogRebin)
        {
            // save rebinned data
            FILE *f_targetSpcAxis =
                fopen("loglbda_rebinlog_targetSpcAxis_dbg.txt", "w+");
            for (Int32 t = 0; t < targetSpectralAxis.GetSamplesCount(); t++)
            {
                fprintf(f_targetSpcAxis, "%f\n", targetSpectralAxis[t]);
            }
            fclose(f_targetSpcAxis);
        }

        bool enableVarianceWeightedRebin = true;
        std:string errorRebinMethod = "rebin";
        if (enableVarianceWeightedRebin)
            errorRebinMethod = "rebinVariance";

        // rebin the spectrum
        TFloat64Range spcLbdaRange(exp(loglbdamin - 0.5 * loglbdaStep),
                                   exp(loglbdamax + 0.5 * loglbdaStep));
        //linear
        spectrum.Rebin(spcLbdaRange, targetSpectralAxis,
                       m_spectrumRebinedLog,
                       m_mskRebinedLog, rebinMethod, errorRebinMethod);

        if (verboseExportLogRebin)
        {
            // save rebinned data
            FILE *f = fopen("loglbda_rebinlog_spclogrebin_dbg.txt", "w+");
            const TAxisSampleList & w = spectrumRebinedSpectralAxis.GetSamplesVector();
            const TAxisSampleList & F = spectrumRebinedFluxAxis.GetSamplesVector();
            for (Int32 t = 0; t < w.size(); t++)
            {
                fprintf(f, "%f\t%e\n", w[t], F[t]);
            }
            fclose(f);
        }
        if (verboseExportLogRebin)
        {
            // save error rebinned data
            FILE *f = fopen("loglbda_rebinlog_errorlogrebin_dbg.txt", "w+");
            const TAxisSampleList & w = spectrumRebinedSpectralAxis.GetSamplesVector();
            const TAxisSampleList & err = spectrumRebinedFluxAxis.GetError().GetSamplesVector();
            for (Int32 t = 0; t < w.size(); t++)
            {
                fprintf(f, "%f\t%e\n", w[t], err[t]);
            }
            fclose(f);
        }

        Log.LogDetail("  Operator-TemplateFittingLog: Log-regular lambda resampling FINISHED");
    } else // check that the raw spectrum grid is in log-regular grid
    {
        Log.LogDetail("  Operator-TemplateFittingLog: Log-regular lambda resampling OFF");

        if (verboseLogRebin)
        {
            Log.LogDebug("  Operator-TemplateFittingLog: lbda raw min_idx=%d, max_idx=%d", lbdaMinIdx, lbdaMaxIdx);
        }

        bool logRegLbdaCheck = true;
        Float64 relativeLogLbdaStepTol = 1e-1;
        Float64 maxAbsRelativeError = 0.0; //
        Float64 lbda1 = spectrum.GetSpectralAxis()[lbdaMinIdx];
        for (Int32 t = 1; t < loglbdaCount; t++)
        {
            Float64 lbda2 = spectrum.GetSpectralAxis()[lbdaMinIdx + t];
            Float64 _loglbdaStep = log(lbda2) - log(lbda1);

            Float64 relativeErrAbs = std::abs((_loglbdaStep - loglbdaStep) / loglbdaStep);
            // if(verboseLogRebin)
            // {
            //     Log.LogInfo("  Operator-TemplateFittingLog: _loglbdastep = %f, relativeErrAbs = %f", _loglbdaStep, relativeErrAbs);
            // }

            if (relativeErrAbs > maxAbsRelativeError)
            {
                maxAbsRelativeError = relativeErrAbs;
            }
            if (relativeErrAbs > relativeLogLbdaStepTol)
            {
                logRegLbdaCheck = false;
                break;
            }
            lbda1 = lbda2;
        }
        Log.LogDetail("  Operator-TemplateFittingLog: Log-regular lambda check: max Abs Relative Error (log lbda step)= %f", maxAbsRelativeError);
        Log.LogDetail("  Operator-TemplateFittingLog: Log-regular lambda check (for rel-tol=%f): (1=success, 0=fail) CHECK = %d", relativeLogLbdaStepTol, logRegLbdaCheck);
        if (!logRegLbdaCheck)
        {
            Log.LogError("  Operator-TemplateFittingLog: Log-regular lambda check FAILED");
        }
        // prepare data for log-regular computation
        for (Int32 t = 0; t < loglbdaCount; t++)
        {
            spectrumRebinedFluxAxis[t] = spectrum.GetFluxAxis()[lbdaMinIdx + t];
            spectrumRebinedSpectralAxis[t] =
                spectrum.GetSpectralAxis()[lbdaMinIdx + t];
            spectrumRebinedFluxAxis.GetError()[t] =
                spectrum.GetFluxAxis().GetError()[lbdaMinIdx + t];
            m_mskRebinedLog[t] = 1;
        }
        if (verboseExportLogRebin)
        {
            // save rebinned data
            FILE *f = fopen("loglbda_rebinlog_spclogrebin_dbg.txt", "w+");
            for (Int32 t = 0; t < m_spectrumRebinedLog.GetSampleCount(); t++)
            {
                fprintf(f, "%f\t%e\n",
                        m_spectrumRebinedLog.GetSpectralAxis()[t],
                        m_spectrumRebinedLog.GetFluxAxis()[t]);
            }
            fclose(f);
        }
    }

     /****************************************************************
     *                                                   
     *  Rebin the template
     *
     ****************************************************************/
    // Create the Template Log-Rebined spectral axis
    {

        // The template grid has to be aligned with the spectrum log-grid redshifted at all the new redshift grid
        // We align the rebined spectrum max lambda at min redshift.

        Float64 tpl_raw_loglbdamin = log( spectrumRebinedSpectralAxis[0]/(1.0 + zmax_new));
        Float64 tpl_raw_loglbdamax = log( spectrumRebinedSpectralAxis[loglbdaCount - 1]/ (1.0 + zmin_new));

        Float64 tpl_tgt_loglbdamax = tpl_raw_loglbdamax;
        Int32 tpl_loglbdaCount = Int32(ceil((tpl_raw_loglbdamax - tpl_raw_loglbdamin)/loglbdaStep)) + 1;
        Float64 tpl_tgt_loglbdamin = tpl_raw_loglbdamax - (tpl_loglbdaCount-1)*loglbdaStep;

	// Display lambda min and max for raw templates
	Log.LogDetail("  Operator-TemplateFittingLog: Log-Rebin: tpl raw loglbdamin=%f : raw loglbdamax=%f", tpl_raw_loglbdamin, tpl_raw_loglbdamax);

	// Display zmin and zmax used to rebin templates
	Log.LogDetail("  Operator-TemplateFittingLog: zmin_new = %f, tpl.lbdamax = %f", zmin_new, exp(tpl_raw_loglbdamax));
	Log.LogDetail("  Operator-TemplateFittingLog: zmax_new = %f, tpl.lbdamin = %f", zmax_new, exp(tpl_raw_loglbdamin));

	// Recheck the template coverage is larger,
	//  by construction only the min has to be re-checked,
	//  since the max is aligned to the max of the rebined spectrum at zmin which is
	// smaller to the max input lamdba range at min already checked 
	if (exp(tpl_tgt_loglbdamin) < tpl.GetSpectralAxis()[0] )
	{
	    Log.LogError("  Operator-TemplateFittingLog: overlap found to be lower than "
			 "1.0 for this redshift range");
            Log.LogError("  Operator-TemplateFittingLog: for zmax=%f, tpl.lbdamin is %f "
			 "(should be <%f)",
			 zmax_new,
			 tpl.GetSpectralAxis()[0],
			 exp(tpl_tgt_loglbdamin));
	    return NULL;
	}
	
        if (verboseLogRebin)
        {
            Log.LogDebug("  Operator-TemplateFittingLog: Log-Rebin: tpl (tgt) loglbdamin=%f "
                         ": loglbdamax=%f",
                         tpl_tgt_loglbdamin, tpl_tgt_loglbdamax);
            Log.LogDebug("  Operator-TemplateFittingLog: Log-Rebin: tpl (tgt) lbdamin=%f : "
                         "lbdamax=%f",
                         exp(tpl_tgt_loglbdamin), exp(tpl_tgt_loglbdamax));
            Log.LogDebug("  Operator-TemplateFittingLog: Log-Rebin: tpl loglbdaCount = %d",
                         tpl_loglbdaCount);
        }
        // todo: check that the coverage is ok with the current tgtTplAxis ?

        // rebin the template
        CSpectrumSpectralAxis tpl_targetSpectralAxis;
        tpl_targetSpectralAxis.SetSize(tpl_loglbdaCount);
        for (Int32 k = 0; k < tpl_loglbdaCount; k++)
        {
            tpl_targetSpectralAxis[k] = exp(tpl_tgt_loglbdamin + k * loglbdaStep);
                exp(tpl_tgt_loglbdamin + k * loglbdaStep);
        }
        // precision problem due to exp/log - beginning
        Float64 delta = tpl_targetSpectralAxis[0] - tpl.GetSpectralAxis()[0];
        Float64 step = tpl_targetSpectralAxis[1] - tpl_targetSpectralAxis[0];
        if (delta < 0.0 && abs(delta) < step * 1e-5)
        {
            tpl_targetSpectralAxis[0] = tpl.GetSpectralAxis()[0];
        }
        // precision problem due to exp/log - end
        delta = tpl_targetSpectralAxis[loglbdaCount - 1] -
                tpl.GetSpectralAxis()[tpl.GetSpectralAxis().GetSamplesCount() - 1];
        step = tpl_targetSpectralAxis[loglbdaCount - 1] -
               tpl_targetSpectralAxis[loglbdaCount - 2];
        if (delta > 0.0 && abs(delta) < step * 1e-5)
        {
            tpl_targetSpectralAxis[loglbdaCount - 1] = tpl.GetSpectralAxis()[tpl.GetSpectralAxis().GetSamplesCount() -1];
        }

        if (verboseExportLogRebin)
        {
            // save rebinned data
            FILE *f_targetTplAxis =
                fopen("loglbda_rebinlog_targetTplAxis_dbg.txt", "w+");
            for (Int32 t = 0; t < tpl_targetSpectralAxis.GetSamplesCount(); t++)
            {
                fprintf(f_targetTplAxis, "%f\n", tpl_targetSpectralAxis[t]);
            }
            fclose(f_targetTplAxis);
        }

        TFloat64Range tplLbdaRange(exp(tpl_tgt_loglbdamin - 0.5 * loglbdaStep),
                                   exp(tpl_tgt_loglbdamax + 0.5 * loglbdaStep));
        // precision problem due to exp/log - beginning
        delta = tpl_targetSpectralAxis[0] - tpl.GetSpectralAxis()[0];
        if (delta < 0.0)
        {
            Log.LogError("  Operator-TemplateFittingLog: Log-Rebin: tpl rebin error. Target "
                         "MIN lbda value=%f, input tpl min lbda value=%f",
                         tpl_targetSpectralAxis[0], tpl.GetSpectralAxis()[0]);
            Log.LogError("  Operator-TemplateFittingLog: Log-Rebin: tpl rebin error. "
                         "Extend your input template wavelength range or "
                         "modify the processing parameter <lambdarange>");
        }
        // precision problem due to exp/log - end
        delta = tpl_targetSpectralAxis[tpl_targetSpectralAxis.GetSamplesCount() -1] -
                tpl.GetSpectralAxis()[tpl.GetSpectralAxis().GetSamplesCount() - 1];
        if (delta > 0.0)
        {
            Log.LogError("  Operator-TemplateFittingLog: Log-Rebin: tpl rebin error. Target "
                         "MAX lbda value=%f, input tpl max lbda value=%f",
                         tpl_targetSpectralAxis[tpl_targetSpectralAxis.GetSamplesCount() - 1],
                         tpl.GetSpectralAxis()[tpl.GetSpectralAxis().GetSamplesCount() -1]);
            Log.LogError("  Operator-TemplateFittingLog: Log-Rebin: tpl rebin error. "
                         "Extend your input template wavelength range or "
                         "modify the processing parameter <lambdarange>");
        }

        tpl.Rebin(tplLbdaRange,
                  tpl_targetSpectralAxis,
                  m_templateRebinedLog,
                  m_mskRebinedLog);

        if (verboseExportLogRebin)
        {
            // save rebinned data
            FILE *f_tpl_tgtlbda =
                fopen("loglbda_rebinlog_tpllogrebin_dbg.txt", "w+");
            for (Int32 t = 0; t < m_templateRebinedLog.GetSampleCount(); t++)
            {
                fprintf(f_tpl_tgtlbda, "%f\t%e\n",
                        m_templateRebinedLog.GetSpectralAxis()[t],
                        m_templateRebinedLog.GetFluxAxis()[t]);
            }
            fclose(f_tpl_tgtlbda);
        }
    }

    //**************** Fitting at all redshifts ****************//
    // Optionally apply some IGM absorption
    std::vector<Int32> igmMeiksinCoeffs;
    if (opt_extinction)
    {
        // Log.LogError("  Operator-TemplateFittingLog: FitAllz opt_extinction not
        // validated yet...");
        Int32 nIGMCoeffs = m_igmCorrectionMeiksin->GetIdxCount();
        igmMeiksinCoeffs.resize(nIGMCoeffs);
        for (Int32 kigm = 0; kigm < nIGMCoeffs; kigm++)
        {
            igmMeiksinCoeffs[kigm] = kigm;
        }
    } else
    {
        igmMeiksinCoeffs.clear();
    }

    // Optionally apply some ISM attenuation
    std::vector<Int32> ismEbmvCoeffs;
    if (opt_dustFitting==-10)
    {
        // Log.LogError("  Operator-TemplateFittingLog: FitAllz opt_dustFitting not
        // validated yet...");
        Int32 nISMCoeffs = m_ismCorrectionCalzetti->GetNPrecomputedDustCoeffs();
        ismEbmvCoeffs.resize(nISMCoeffs);
        for (Int32 kism = 0; kism < nISMCoeffs; kism++)
        {
            ismEbmvCoeffs[kism] = kism;
        }
    }else if (opt_dustFitting>-1)
    {
        Int32 nISMCoeffs = 1;
        ismEbmvCoeffs.resize(nISMCoeffs);
        ismEbmvCoeffs[0] = opt_dustFitting;

    }else
    {
        ismEbmvCoeffs.clear();
    }

    std::shared_ptr<CTemplateFittingResult> result =
        std::shared_ptr<CTemplateFittingResult>(new CTemplateFittingResult());
    result->Init(sortedRedshifts.size(),
                 std::max((Int32)ismEbmvCoeffs.size(), 1),
                 std::max((Int32)igmMeiksinCoeffs.size(), 1));
    result->Redshifts = sortedRedshifts;

    // WARNING: no additional masks coded for use as of 2017-06-13
    if (additional_spcMasks.size() != 0)
    {
        Log.LogError("  Operator-TemplateFittingLog: No additional masks used. "
                     "Feature not coded for this log-lambda operator!)");
    }

    if(logpriorze.size()>0 && logpriorze.size()!=sortedRedshifts.size())
    {
        Log.LogError("  Operator-TemplateFittingLog: prior list size (%d) didn't match the input redshift-list (%d) !)", logpriorze.size(), sortedRedshifts.size());
        throw std::runtime_error("  Operator-TemplateFittingLog: prior list size didn't match the input redshift-list size");
    }

    //*
    Int32 retFit = FitAllz(lambdaRange, result, igmMeiksinCoeffs, ismEbmvCoeffs, CMask(), logpriorze);
    if (retFit != 0)
    {
        Log.LogError("  Operator-TemplateFittingLog: FitAllz failed with error %d", retFit);
    }
    //*/
    //**************** End Fitting at all redshifts ****************//

    // overlap warning
    Float64 overlapValidInfZ = -1;
    for (Int32 i = 0; i < sortedRedshifts.size(); i++)
    {
        if (result->Overlap[i] >= overlapThreshold && overlapValidInfZ == -1)
        {
            overlapValidInfZ = sortedRedshifts[i];
            break;
        }
    }
    Float64 overlapValidSupZ = -1;
    for (Int32 i = sortedRedshifts.size() - 1; i >= 0; i--)
    {
        if (result->Overlap[i] >= overlapThreshold && overlapValidSupZ == -1)
        {
            overlapValidSupZ = sortedRedshifts[i];
            break;
        }
    }
    if (overlapValidInfZ != sortedRedshifts[0] ||
        overlapValidSupZ != sortedRedshifts[sortedRedshifts.size() - 1])
    {
        Log.LogInfo("  Operator-TemplateFittingLog: overlap warning for %s: minz=%.3f, maxz=%.3f",
                    tpl.GetName().c_str(), overlapValidInfZ, overlapValidSupZ);
    }

    // estimate CstLog for PDF estimation
    result->CstLog = EstimateLikelihoodCstLog(
        m_spectrumRebinedLog, lambdaRange); // 0.0;//Todo: check how to estimate
                                            // that value for loglambda//

    return result;
}

/**
 * \brief this function estimates the likelihood_cstLog term withing the
 *wavelength range
 **/
Float64 COperatorTemplateFittingLog::EstimateLikelihoodCstLog(
    const CSpectrum &spectrum, const TFloat64Range &lambdaRange)
{
    const CSpectrumSpectralAxis &spcSpectralAxis = spectrum.GetSpectralAxis();
    const TFloat64List &error = spectrum.GetFluxAxis().GetError().GetSamplesVector();

    Int32 numDevs = 0;
    Float64 cstLog = 0.0;
    Float64 sumLogNoise = 0.0;

    Float64 imin = spcSpectralAxis.GetIndexAtWaveLength(lambdaRange.GetBegin());
    Float64 imax = spcSpectralAxis.GetIndexAtWaveLength(lambdaRange.GetEnd());
    for (UInt32 j = imin; j < imax; j++)
    {
        numDevs++;
        sumLogNoise += log(error[j]);
    }
    // Log.LogDebug( "CLineModelElementList::EstimateMTransposeM val = %f", mtm
    // );

    cstLog = -numDevs * 0.5 * log(2 * M_PI) - sumLogNoise;

    return cstLog;
}

void COperatorTemplateFittingLog::enableSpcLogRebin(Bool enable)
{
    m_opt_spcrebin = enable;
    Log.LogDetail("  Operator-TemplateFittingLog: Spectrum REBIN-LOG enabled=%d", m_opt_spcrebin);
}
