#include <RedshiftLibrary/operator/chisquareloglambda.h>

#include <RedshiftLibrary/common/mask.h>
#include <RedshiftLibrary/common/quicksort.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/operator/chisquareresult.h>
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

COperatorChiSquareLogLambda::COperatorChiSquareLogLambda(
    std::string calibrationPath)
{
    m_opt_spcrebin = true;

    // ISM
    m_ismCorrectionCalzetti = new CSpectrumFluxCorrectionCalzetti();
    m_ismCorrectionCalzetti->Init(calibrationPath, 0.0, 0.1, 10);

    // IGM
    m_igmCorrectionMeiksin = new CSpectrumFluxCorrectionMeiksin();
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

COperatorChiSquareLogLambda::~COperatorChiSquareLogLambda()
{
    delete m_ismCorrectionCalzetti;
    delete m_igmCorrectionMeiksin;
    freeFFTPlans();
}

Int32 COperatorChiSquareLogLambda::EstimateXtYSlow(const Float64 *X,
                                                   const Float64 *Y, UInt32 nX,
                                                   UInt32 nShifts,
                                                   std::vector<Float64> &XtY)
{
    XtY.resize(nShifts);

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
Int32 COperatorChiSquareLogLambda::EstimateMtMFast(const Float64 *X,
                                                   const Float64 *Y, UInt32 nX,
                                                   UInt32 nShifts,
                                                   std::vector<Float64> &XtY)
{
    XtY.resize(nShifts);

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

Int32 COperatorChiSquareLogLambda::EstimateXtY(const Float64 *X,
                                               const Float64 *Y, UInt32 nx,
                                               UInt32 ny, UInt32 nshifts,
                                               std::vector<Float64> &XtY,
                                               Int32 precomputedFFT)
{
    // Processing the FFT
    Int32 nSpc = nx;
    Int32 nTpl = ny;
    Int32 nPadded = m_nPaddedSamples;

    Int32 nPadBeforeSpc = nPadded - nSpc; //(Int32)nPadded/2.0;
    Int32 nPadBeforeTpl = 0;

    if (verboseLogXtYFFT)
    {
        Log.LogInfo("  Operator-ChisquareLog: FitAllz: Processing spc-fft with "
                    "n=%d, padded to n=%d",
                    nSpc, nPadded);
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
            Log.LogInfo("  Operator-ChisquareLog: FitAllz: Processing spc-fft "
                        "with nPadBeforeSpc=%d",
                        nPadBeforeSpc);
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
            Log.LogInfo("  Operator-ChisquareLog: FitAllz: spc-fft done");
        }
        if (verboseExportXtYFFT)
        {
            // save spc-fft data
            FILE *f_fftoutput =
                fopen("loglbda_fitallz_xfftOutput_dbg.txt", "w+");
            for (Int32 t = 0; t < nPadded; t++)
            {
                fprintf(f_fftoutput, "%f\t%e\n", (Float64)t,
                        outSpc[t][0] * outSpc[t][0] +
                            outSpc[t][1] * outSpc[t][1]);
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
                Log.LogError("  Operator-ChisquareLog: InitFFT: Unable to "
                             "allocate precomputedFFT_spcFluxOverErr2");
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
                Log.LogError("  Operator-ChisquareLog: InitFFT: Unable to "
                             "allocate precomputedFFT_spcOneOverErr2");
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
        Log.LogInfo("  Operator-ChisquareLog: FitAllz: Processing tpl-fft with "
                    "n=%d, padded to n=%d",
                    nTpl, nPadded);
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
        Log.LogInfo("  Operator-ChisquareLog: FitAllz: tpl-fft done");
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
        outCombined[k][0] =
            (outTpl[k][0] * outSpc[k][0] - outTpl[k][1] * outSpc[k][1]);
        outCombined[k][1] =
            (outTpl[k][0] * outSpc[k][1] + outTpl[k][1] * outSpc[k][0]);
        // Y conjugate
        // outCombined[k][0] =
        // (outTpl[k][0]*outSpc[k][0]+outTpl[k][1]*outSpc[k][1]);
        // outCombined[k][1] =
        // (outTpl[k][0]*outSpc[k][1]-outTpl[k][1]*outSpc[k][0]);
    }

    fftw_execute(pBackward);
    if (verboseLogXtYFFT)
    {
        Log.LogInfo("  Operator-ChisquareLog: FitAllz: backward-fft done");
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
        FILE *f_fftoutput =
            fopen("loglbda_fitallz_xtyfft-output_dbg.txt", "w+");
        for (Int32 t = 0; t < nshifts; t++)
        {
            fprintf(f_fftoutput, "%f\t%e\n", (Float64)t, XtY[t]);
        }
        fclose(f_fftoutput);
    }

    return 0;
}

Int32 COperatorChiSquareLogLambda::InitFFT(Int32 nPadded)
{
    freeFFTPlans();

    inSpc = (Float64 *)fftw_malloc(sizeof(Float64) * nPadded);
    outSpc = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nPadded);
    pSpc = fftw_plan_dft_r2c_1d(nPadded, inSpc, outSpc, FFTW_ESTIMATE);
    if (inSpc == 0)
    {
        Log.LogError(
            "  Operator-ChisquareLog: InitFFT: Unable to allocate inSpc");
        return -1;
    }
    if (outSpc == 0)
    {
        Log.LogError(
            "  Operator-ChisquareLog: InitFFT: Unable to allocate inSpc");
        return -1;
    }

    inTpl_padded = (Float64 *)fftw_malloc(sizeof(Float64) * nPadded);
    inTpl = (Float64 *)fftw_malloc(sizeof(Float64) * nPadded);
    outTpl = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nPadded);
    pTpl = fftw_plan_dft_r2c_1d(nPadded, inTpl, outTpl, FFTW_ESTIMATE);
    if (inTpl_padded == 0)
    {
        Log.LogError("  Operator-ChisquareLog: InitFFT: Unable to allocate "
                     "inTpl_padded");
        return -1;
    }
    if (inTpl == 0)
    {
        Log.LogError(
            "  Operator-ChisquareLog: InitFFT: Unable to allocate inTpl");
        return -1;
    }
    if (outTpl == 0)
    {
        Log.LogError(
            "  Operator-ChisquareLog: InitFFT: Unable to allocate outTpl");
        return -1;
    }

    outCombined = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nPadded);
    inCombined = (Float64 *)fftw_malloc(sizeof(Float64) * nPadded);
    pBackward =
        fftw_plan_dft_c2r_1d(nPadded, outCombined, inCombined, FFTW_ESTIMATE);
    if (outCombined == 0)
    {
        Log.LogError(
            "  Operator-ChisquareLog: InitFFT: Unable to allocate outCombined");
        return -1;
    }
    if (inCombined == 0)
    {
        Log.LogError(
            "  Operator-ChisquareLog: InitFFT: Unable to allocate inCombined");
        return -1;
    }

    // reinit fft precomputed buffers
    freeFFTPrecomputedBuffers();

    return 0;
}

void COperatorChiSquareLogLambda::freeFFTPrecomputedBuffers()
{
    if (!precomputedFFT_spcFluxOverErr2 == 0)
    {
        fftw_free(precomputedFFT_spcFluxOverErr2);
        precomputedFFT_spcFluxOverErr2 = 0;
    }
    if (!precomputedFFT_spcOneOverErr2 == 0)
    {
        fftw_free(precomputedFFT_spcOneOverErr2);
        precomputedFFT_spcOneOverErr2 = 0;
    }
}

void COperatorChiSquareLogLambda::freeFFTPlans()
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
 * @brief COperatorChiSquareLogLambda::FitAllz
 *
 * @param lambdaRange
 * @param result
 * @param opt_extinction
 * @param opt_dustFitting
 * @param spcMaskAdditional
 * @return
 */
Int32 COperatorChiSquareLogLambda::FitAllz(
    const TFloat64Range &lambdaRange, std::shared_ptr<CChisquareResult> result,
    std::vector<Int32> igmMeiksinCoeffs, std::vector<Int32> ismEbmvCoeffs,
    CMask spcMaskAdditional)
{
    bool verboseLogFitAllz = false;

    CSpectrumFluxAxis &spectrumRebinedFluxAxis =
        m_spectrumRebinedLog.GetFluxAxis();
    Float64 *error = m_errorRebinedLog.GetSamples();
    CSpectrumSpectralAxis &spectrumRebinedSpectralAxis =
        m_spectrumRebinedLog.GetSpectralAxis();
    CSpectrumFluxAxis &tplRebinedFluxAxis = m_templateRebinedLog.GetFluxAxis();
    CSpectrumSpectralAxis &tplRebinedSpectralAxis =
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
                    // Log.LogInfo("  Operator-ChisquareLog:
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
            Log.LogInfo("  Operator-ChisquareLog: FitAllz: indexes for full "
                        "LstSquare calculation, count = %d",
                        zindexesFullLstSquare.size());
            for (Int32 k = 0; k < zindexesFullLstSquare.size(); k++)
            {
                Log.LogInfo("  Operator-ChisquareLog: FitAllz: indexes ranges: "
                            "for i=%d, zindexesFullLstSquare=%d",
                            k, zindexesFullLstSquare[k]);
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
                izmin =
                    zindexesFullLstSquare[zindexesFullLstSquare.size() - 1] + 1;
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
        Log.LogInfo("  Operator-ChisquareLog: FitAllz: indexes - izrangelist "
                    "calculation, count = %d",
                    izrangelist.size());
        for (Int32 k = 0; k < nzranges; k++)
        {
            Log.LogInfo("  Operator-ChisquareLog: FitAllz: indexes ranges: for "
                        "i=%d, zmin=%f, zmax=%f",
                        k, result->Redshifts[izrangelist[k].GetBegin()],
                        result->Redshifts[izrangelist[k].GetEnd()]);
        }
    }

    Float64 *spectrumRebinedLambda = spectrumRebinedSpectralAxis.GetSamples();
    Float64 *spectrumRebinedFluxRaw = spectrumRebinedFluxAxis.GetSamples();
    UInt32 nSpc = spectrumRebinedSpectralAxis.GetSamplesCount();
    if (verboseLogFitAllz)
    {
        Log.LogInfo("  Operator-ChisquareLog: FitAllz: spc lbda min=%f, max=%f",
                    spectrumRebinedLambda[0], spectrumRebinedLambda[nSpc - 1]);
    }

    Float64 *tplRebinedLambdaGlobal = tplRebinedSpectralAxis.GetSamples();
    Float64 *tplRebinedFluxRawGlobal = tplRebinedFluxAxis.GetSamples();
    //    UInt32 nTpl = tplRebinedSpectralAxis.GetSamplesCount();
    Float64 *tplRebinedLambda =
        new Float64[(int)tplRebinedSpectralAxis.GetSamplesCount()]();
    Float64 *tplRebinedFluxRaw =
        new Float64[(int)tplRebinedSpectralAxis.GetSamplesCount()]();

    for (Int32 k = 0; k < nzranges; k++)
    {
        // prepare the zrange-result container
        std::shared_ptr<CChisquareResult> subresult =
            std::shared_ptr<CChisquareResult>(new CChisquareResult());
        TFloat64Range zrange =
            TFloat64Range(result->Redshifts[izrangelist[k].GetBegin()],
                          result->Redshifts[izrangelist[k].GetEnd()]);
        TInt32Range ilbda;
        if (result->Redshifts.size() > 1)
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
            ilbda = FindTplSpectralIndex(
                spectrumRebinedLambda, tplRebinedLambdaGlobal,
                spectrumRebinedSpectralAxis.GetSamplesCount(),
                tplRebinedSpectralAxis.GetSamplesCount(), zrange,
                redshiftStep_toBeDoneDifferently);
        } else
        {
            ilbda =
                TInt32Range(0, tplRebinedSpectralAxis.GetSamplesCount() - 1);
            subresult->Init(1, std::max((Int32)ismEbmvCoeffs.size(), 1),
                            std::max((Int32)igmMeiksinCoeffs.size(), 1));
            subresult->Redshifts = result->Redshifts;
        }
        if (ilbda.GetBegin() == -1 || ilbda.GetEnd() == -1)
        {
            Log.LogError(
                "  Operator-ChisquareLog: FitAllz: Can't find tpl indexes");
            return -1;
        }

        if (verboseLogFitAllz)
        {
            Log.LogInfo(
                "  Operator-ChisquareLog: FitAllz: zrange min=%f, max=%f",
                zrange.GetBegin(), zrange.GetEnd());
            Log.LogInfo(
                "  Operator-ChisquareLog: FitAllz: full zmin=%f, full zmax=%f",
                result->Redshifts[0],
                result->Redshifts[result->Redshifts.size() - 1]);
            Log.LogInfo("  Operator-ChisquareLog: FitAllz: indexes tpl crop: "
                        "lbda min=%d, max=%d",
                        ilbda.GetBegin(), ilbda.GetEnd());
            Log.LogInfo("  Operator-ChisquareLog: FitAllz: indexes tpl full: "
                        "lbda min=%d, max=%d",
                        0, tplRebinedSpectralAxis.GetSamplesCount());
            Log.LogInfo("  Operator-ChisquareLog: FitAllz: tpl lbda "
                        "min*zmax=%f, max*zmin=%f",
                        tplRebinedLambdaGlobal[ilbda.GetBegin()] *
                            (1.0 + zrange.GetEnd()),
                        tplRebinedLambdaGlobal[ilbda.GetEnd()] *
                            (1.0 + zrange.GetBegin()));
        }

        UInt32 nTpl = ilbda.GetEnd() - ilbda.GetBegin() + 1;
        for (UInt32 j = 0; j < nTpl; j++)
        {
            tplRebinedLambda[j] = tplRebinedLambdaGlobal[j + ilbda.GetBegin()];
            tplRebinedFluxRaw[j] =
                tplRebinedFluxRawGlobal[j + ilbda.GetBegin()];
        }

        //*
        FitRangez(spectrumRebinedLambda, spectrumRebinedFluxRaw, error,
                  tplRebinedLambda, tplRebinedFluxRaw, nSpc, nTpl, subresult,
                  igmMeiksinCoeffs, ismEbmvCoeffs);
        //*/

        // copy subresults into global results
        for (UInt32 isubz = 0; isubz < subresult->Redshifts.size(); isubz++)
        {
            UInt32 fullResultIdx = isubz + izrangelist[k].GetBegin();
            result->ChiSquare[fullResultIdx] = subresult->ChiSquare[isubz];
            result->Overlap[fullResultIdx] = subresult->Overlap[isubz];
            result->FitAmplitude[fullResultIdx] =
                subresult->FitAmplitude[isubz];
            result->FitDtM[fullResultIdx] = subresult->FitDtM[isubz];
            result->FitMtM[fullResultIdx] = subresult->FitMtM[isubz];
            result->FitDustCoeff[fullResultIdx] =
                subresult->FitDustCoeff[isubz];
            result->FitMeiksinIdx[fullResultIdx] =
                subresult->FitMeiksinIdx[isubz];
            result->Status[fullResultIdx] = subresult->Status[isubz];

            for (Int32 kism = 0;
                 kism < result->ChiSquareIntermediate[fullResultIdx].size();
                 kism++)
            {
                for (Int32 kigm = 0;
                     kigm <
                     result->ChiSquareIntermediate[fullResultIdx][kism].size();
                     kigm++)
                {
                    result->ChiSquareIntermediate[fullResultIdx][kism][kigm] =
                        subresult->ChiSquareIntermediate[isubz][kism][kigm];
                }
            }
        }
    }

    delete[] tplRebinedLambda;
    delete[] tplRebinedFluxRaw;

    return 0;
}

/**
  // TODO : many vectors allocated in this function. Check if the allocation
 time is significant, and eventually use preallocated member buffers...
 * @brief COperatorChiSquareLogLambda::FitRangez
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
Int32 COperatorChiSquareLogLambda::FitRangez(
    Float64 *spectrumRebinedLambda, Float64 *spectrumRebinedFluxRaw,
    Float64 *error, Float64 *tplRebinedLambda, Float64 *tplRebinedFluxRaw,
    UInt32 nSpc, UInt32 nTpl, std::shared_ptr<CChisquareResult> result,
    std::vector<Int32> igmMeiksinCoeffs, std::vector<Int32> ismEbmvCoeffs)
{
    Float64 redshiftValueMeiksin = result->Redshifts[0];
    if (verboseLogFitFitRangez)
    {
        Log.LogInfo(
            "  Operator-ChisquareLog: FitRangez: redshiftValueMeiksin = %f",
            redshiftValueMeiksin);
        Log.LogInfo("  Operator-ChisquareLog: FitRangez: spc[0] = %f",
                    spectrumRebinedLambda[0]);
        Log.LogInfo("  Operator-ChisquareLog: FitRangez: spc[max] = %f",
                    spectrumRebinedLambda[nSpc - 1]);
        Log.LogInfo(
            "  Operator-ChisquareLog: FitRangez: tpl[0]*zmax = %f",
            tplRebinedLambda[0] *
                (1.0 + result->Redshifts[result->Redshifts.size() - 1]));
        Log.LogInfo("  Operator-ChisquareLog: FitRangez: tpl[max]*zmin = %f",
                    tplRebinedLambda[nTpl - 1] * (1 + result->Redshifts[0]));
    }

    //
    Int32 nshifts = nTpl - nSpc;
    // Int32 nshifts = nTpl*2.0;
    m_nPaddedSamples = (Int32)nTpl * 2.0;
    /*
    //next power of two
    Int32 maxnsamples = (Int32)(nTpl*2);
    //Log.LogInfo("  Operator-ChisquareLog: FitRangez: maxnsamples = %d",
    maxnsamples); Int32 power = 1; while(power < maxnsamples)
    {
        power*=2;
    }
    m_nPaddedSamples = power;
    //*/

    Log.LogDetail("  Operator-ChisquareLog: Now fitting using the FFT on "
                  "nshifts=%d values, for Meiksin redshift=%f",
                  nshifts, redshiftValueMeiksin);

    if (verboseLogFitFitRangez)
    {
        Log.LogInfo("  Operator-ChisquareLog: FitRangez: initializing FFT with "
                    "n = %d points",
                    m_nPaddedSamples);
    }
    InitFFT(m_nPaddedSamples);

    // prepare z array
    Int32 zoff = 1;
    std::vector<Float64> z_vect(nshifts, 0.0);
    for (Int32 t = 0; t < z_vect.size(); t++)
    {
        z_vect[t] = (spectrumRebinedLambda[0] - tplRebinedLambda[t + zoff]) /
                    tplRebinedLambda[t + zoff];
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
    Float64 *spcRebinedFluxOverErr2 = new Float64[nSpc]();
    Float64 *oneSpcRebinedFluxOverErr2 = new Float64[nSpc]();
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

    Float64 *tplRebinedFluxIgm = new Float64[nTpl]();
    Float64 *tplRebinedFlux = new Float64[nTpl]();
    Float64 *tpl2RebinedFlux = new Float64[nTpl]();
    for (Int32 j = 0; j < nTpl; j++)
    {
        tplRebinedFluxIgm[j] = tplRebinedFluxRaw[j];
        tplRebinedFlux[j] = tplRebinedFluxRaw[j];
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
            Log.LogInfo("  Operator-ChisquareLog: FitRangez: IGM disabled, "
                        "min-tpl-lbda=%f",
                        lambdaMinTpl);
        }
    }

    // prepare best fit data buffer
    std::vector<Float64> bestChi2(nshifts, DBL_MAX);
    std::vector<Float64> bestFitAmp(nshifts, -1.0);
    std::vector<Float64> bestFitDtm(nshifts, -1.0);
    std::vector<Float64> bestFitMtm(nshifts, -1.0);
    std::vector<Float64> bestISMCoeff(nshifts, -1.0);
    std::vector<Float64> bestIGMIdx(nshifts, -1.0);

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
            Log.LogInfo("  Operator-ChisquareLog: FitRangez: IGM index=%d",
                        kIGM);
        }

        if (enableIGM)
        {

            Int32 meiksinIdx = igmMeiksinCoeffs[kIGM];

            for (Int32 j = 0; j < nTpl; j++)
            {
                Float64 restLambda = tplRebinedLambda[j];
                Float64 coeffIGM = m_igmCorrectionMeiksin->getCoeff(
                    meiksinIdx, redshiftValueMeiksin, restLambda);

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
                Log.LogInfo("  Operator-ChisquareLog: FitRangez: ISM index =%d",
                            kISM);
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
                    ebmvDustCoeff = m_ismCorrectionCalzetti->getDustCoeff(
                        kDustCalzetti, restLambda);

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
                            fprintf(f, "%f\t%e\n", tplRebinedLambda[j],
                                    tplRebinedFlux[j]);
                        }
                        fclose(f);
                    }
                }
            }

            // Estimate DtM
            std::vector<Float64> dtm_vec;
            //*
            EstimateXtY(spcRebinedFluxOverErr2, tplRebinedFlux, nSpc, nTpl,
                        nshifts, dtm_vec,
                        //-1);
                        0);
            //*/
            // EstimateXtYSlow(spcRebinedFluxOverErr2, tplRebinedFlux, nSpc,
            // nshifts, dtm_vec);

            if (verboseExportFitRangez)
            {
                // save chi2 data
                FILE *f = fopen("loglbda_dtm_dbg.txt", "w+");
                for (Int32 t = 0; t < dtm_vec.size(); t++)
                {
                    fprintf(f, "%f\t%e\n", (Float64)t, dtm_vec[t]);
                }
                fclose(f);
            }

            // Estimate MtM
            std::vector<Float64> mtm_vec;
            //*
            EstimateXtY(oneSpcRebinedFluxOverErr2, tpl2RebinedFlux, nSpc, nTpl,
                        nshifts, mtm_vec,
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
                Log.LogInfo("  Operator-ChisquareLog: FitRangez: dtd = %e",
                            dtd);
            }

            // return -1;

            // Estimate Chi2
            if (mtm_vec.size() != dtm_vec.size())
            {
                errorWhileFitting = 1;
                Log.LogError("  Operator-ChisquareLog: FitRangez: xty vectors "
                             "sizes don't match: dtm size = %d, mtm size = %d",
                             dtm_vec.size(), mtm_vec.size());
                // return 2;
                continue;
            }
            std::vector<Float64> chi2(dtm_vec.size(), DBL_MAX);
            std::vector<Float64> amp(dtm_vec.size(), DBL_MAX);
            for (Int32 k = 0; k < dtm_vec.size(); k++)
            {
                if (mtm_vec[k] == 0.0)
                {
                    amp[k] = 0.0;
                    chi2[k] = dtd;
                } else
                {
                    amp[k] = max(0.0, dtm_vec[k] / mtm_vec[k]);
                    chi2[k] = dtd - dtm_vec[k] * amp[k];
                }
                // chi2[k] = dtm_vec[k];
            }

            for (Int32 k = 0; k < dtm_vec.size(); k++)
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
                    bestFitDtm[k] = dtm_vec[k];
                    bestFitMtm[k] = mtm_vec[k];
                    if (enableISM)
                    {
                        Int32 kDustCalzetti = ismEbmvCoeffs[kISM];
                        bestISMCoeff[k] = m_ismCorrectionCalzetti->GetEbmvValue(
                            kDustCalzetti);
                    } else
                    {
                        bestISMCoeff[k] = -1;
                    }
                    if (enableIGM)
                    {
                        bestIGMIdx[k] = igmMeiksinCoeffs[kIGM];
                    } else
                    {
                        bestIGMIdx[k] = -1;
                    }
                }
            }

            if (verboseExportFitRangez)
            {
                Log.LogInfo(
                    "  Operator-ChisquareLog: FitRangez: spc lbda 0 =%f",
                    spectrumRebinedLambda[0]);
                Log.LogInfo(
                    "  Operator-ChisquareLog: FitRangez: tpl lbda 0 =%f",
                    tplRebinedLambda[0]);
                Float64 z_O = (spectrumRebinedLambda[0] - tplRebinedLambda[0]) /
                              tplRebinedLambda[0];
                Log.LogInfo("  Operator-ChisquareLog: FitRangez: z 0 =%f", z_O);

                //        Float64 logstep =
                //        log(spectrumRebinedLambda[1])-log(spectrumRebinedLambda[0]);
                //        Float64 logInitialStep =
                //        log(spectrumRebinedLambda[0])-log(tplRebinedLambda[0]);
                //        Log.LogInfo("  Operator-ChisquareLog: FitRangez:
                //        logstep=%f", logstep); Log.LogInfo("
                //        Operator-ChisquareLog: FitRangez: step=%f A",
                //        exp(logstep)); Log.LogInfo("  Operator-ChisquareLog:
                //        FitRangez: logInitialStep=%f", logInitialStep);
                //        Log.LogInfo("  Operator-ChisquareLog: FitRangez:
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
      delete[] spcRebinedFluxOverErr2;
      delete[] oneSpcRebinedFluxOverErr2;
      delete[] tplRebinedFlux;
      delete[] tpl2RebinedFlux;
      delete[] tplRebinedFluxIgm;

      freeFFTPlans();
      throw runtime_error("Error while fitting");
        // return -2;
    }

    // interpolating on the regular z grid
    Float64 *zreversed_array = new Float64[(int)z_vect.size()]();
    Float64 *chi2reversed_array = new Float64[(int)bestChi2.size()]();
    Float64 *ampreversed_array = new Float64[(int)bestFitAmp.size()]();
    Float64 *dtmreversed_array = new Float64[(int)bestFitDtm.size()]();
    Float64 *mtmreversed_array = new Float64[(int)bestFitMtm.size()]();
    Float64 *ismCoeffreversed_array = new Float64[(int)bestISMCoeff.size()]();
    Float64 *igmIdxreversed_array = new Float64[(int)bestIGMIdx.size()]();
    for (Int32 t = 0; t < z_vect.size(); t++)
    {
        zreversed_array[t] = z_vect[z_vect.size() - 1 - t];
        chi2reversed_array[t] = bestChi2[z_vect.size() - 1 - t];
        ampreversed_array[t] = bestFitAmp[z_vect.size() - 1 - t];
        dtmreversed_array[t] = bestFitDtm[z_vect.size() - 1 - t];
        mtmreversed_array[t] = bestFitMtm[z_vect.size() - 1 - t];
        ismCoeffreversed_array[t] = bestISMCoeff[z_vect.size() - 1 - t];
        igmIdxreversed_array[t] = bestIGMIdx[z_vect.size() - 1 - t];
        // Log.LogInfo("  Operator-ChisquareLog: FitRangez: interpolating z
        // result, for z=%f, chi2reversed_array=%f", zreversed_array[t],
        // chi2reversed_array[t]);
    }

    Int32 k = 0;
    Int32 klow = 0;
    for (Int32 iz = 0; iz < result->Redshifts.size(); iz++)
    {
        //* //NGP
        // Log.LogInfo("  Operator-ChisquareLog: FitRangez: interpolating
        // gsl-bsearch z result, , ztgt=%f, zcalc_min=%f, zcalc_max=%f",
        // result->Redshifts[iz], zreversed_array[0],
        // zreversed_array[z_vect.size()-1]);
        k = gsl_interp_bsearch(zreversed_array, result->Redshifts[iz], klow,
                               z_vect.size() - 1);
        klow = k;

        /*
        if(result->Redshifts[iz]==0.0)
        {
            Log.LogInfo("  Operator-ChisquareLog: FitRangez: interpolating z
        result, kshift=%f", k); Log.LogInfo("  Operator-ChisquareLog: FitRangez:
        interpolating z result, zcalc=%f, kfound=%d, zfound=%f",
        result->Redshifts[iz], k, z_vect[z_vect.size()-1-k]); Log.LogInfo("
        Operator-ChisquareLog: FitRangez: interpolating z result, bestChi2=%f",
        bestChi2[z_vect.size()-1-k]);
        }
        // closest value
        result->ChiSquare[iz] = bestChi2[z_vect.size()-1-k];
        //*/

        result->Overlap[iz] = 1.0;
        result->FitAmplitude[iz] = ampreversed_array[k];
        result->FitDtM[iz] = dtmreversed_array[k];
        result->FitMtM[iz] = mtmreversed_array[k];
        result->FitDustCoeff[iz] = ismCoeffreversed_array[k];
        result->FitMeiksinIdx[iz] = igmIdxreversed_array[k];
        result->Status[iz] = nStatus_OK;

        //*/
    }

    //*
    Log.LogDetail("  Operator-ChisquareLog: FitRangez: interpolating (lin) z "
                  "result from n=%d (min=%f, max=%f) to n=%d (min=%f, max=%f)",
                  z_vect.size(), zreversed_array[0],
                  zreversed_array[z_vect.size() - 1], result->Redshifts.size(),
                  result->Redshifts[0],
                  result->Redshifts[result->Redshifts.size() - 1]);
    InterpolateResult(chi2reversed_array, zreversed_array,
                      &result->Redshifts.front(), z_vect.size(),
                      result->Redshifts.size(), result->ChiSquare, DBL_MAX);
    //*/

    //*
    // Interpolating intermediate chisquare results
    Float64 *intermChi2BufferReversed_array =
        new Float64[(int)intermediateChi2.size()]();
    std::vector<Float64> intermChi2BufferRebinned_array(
        result->Redshifts.size(), boost::numeric::bounds<float>::highest());
    for (Int32 kism = 0; kism < nISM; kism++)
    {
        for (Int32 kigm = 0; kigm < nIGMFinal; kigm++)
        {
            for (Int32 t = 0; t < z_vect.size(); t++)
            {
                intermChi2BufferReversed_array[t] =
                    intermediateChi2[z_vect.size() - 1 - t][kism][kigm];
            }
            InterpolateResult(intermChi2BufferReversed_array, zreversed_array,
                              &result->Redshifts.front(), z_vect.size(),
                              result->Redshifts.size(),
                              intermChi2BufferRebinned_array, DBL_MAX);
            for (Int32 t = 0; t < result->Redshifts.size(); t++)
            {
                result->ChiSquareIntermediate[t][kism][kigm] =
                    intermChi2BufferRebinned_array[t];
            }
        }
    }

    delete[] intermChi2BufferReversed_array;

    delete[] zreversed_array;
    delete[] chi2reversed_array;
    delete[] ampreversed_array;
    delete[] dtmreversed_array;
    delete[] mtmreversed_array;
    delete[] ismCoeffreversed_array;
    delete[] igmIdxreversed_array;

    delete[] spcRebinedFluxOverErr2;
    delete[] oneSpcRebinedFluxOverErr2;
    delete[] tplRebinedFlux;
    delete[] tpl2RebinedFlux;
    delete[] tplRebinedFluxIgm;

    freeFFTPlans();
    return 0;
}

Int32 COperatorChiSquareLogLambda::InterpolateResult(
    const Float64 *in, const Float64 *inGrid, const Float64 *tgtGrid, Int32 n,
    Int32 tgtn, std::vector<Float64> &out, Float64 defaultValue)
{
    out.resize(tgtn);
    for (Int32 j = 0; j < tgtn; j++)
    {
        out[j] = defaultValue;
    }

    //* // GSL method spline
    // initialise and allocate the gsl objects
    // lin
    gsl_interp *interpolation = gsl_interp_alloc(gsl_interp_linear, n);
    gsl_interp_init(interpolation, inGrid, in, n);
    gsl_interp_accel *accelerator = gsl_interp_accel_alloc();
    // spline
    // gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n);
    // gsl_spline_init (spline, inGrid, in, n);
    // gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();

    for (Int32 j = 0; j < tgtn; j++)
    {
        Float64 Xrebin = tgtGrid[j];
        out[j] = gsl_interp_eval(interpolation, inGrid, in, Xrebin,
                                 accelerator); // lin
        // out[j] = gsl_spline_eval (spline, Xrebin, accelerator); //spline
        // Log.LogInfo("  Operator-ChisquareLog: FitAllz: interpolating
        // gsl-spline z result, , ztgt=%f, rebinY=%f", tgtGrid[j], out[j]);
    }

    gsl_interp_free(interpolation);
    // gsl_spline_free (spline);
    gsl_interp_accel_free(accelerator);

    return 0;
}

TInt32Range COperatorChiSquareLogLambda::FindTplSpectralIndex(
    const Float64 *spcLambda, const Float64 *tplLambda, UInt32 nSpc,
    UInt32 nTpl, TFloat64Range redshiftrange, Float64 redshiftStep)
{
    Float64 redshiftMargin =
        10.0 *
        redshiftStep; // this redshiftstep margin might be useless now that
                      // there are spclambdamargin and tpllambdamargin
    Float64 spcLambdaMargin = abs(spcLambda[1] - spcLambda[0]);
    Float64 tplLambdaMargin = abs(tplLambda[1] - tplLambda[0]);

    UInt32 ilbdamin = 0;
    for (UInt32 k = 0; k < nTpl; k++)
    {
        Float64 _spcLambda = spcLambda[0] - spcLambdaMargin;
        Float64 _tplLambda = tplLambda[k] + tplLambdaMargin;
        Float64 _z = (_spcLambda - _tplLambda) / _tplLambda;
        if (_z > redshiftrange.GetEnd() + redshiftMargin)
        {
            ilbdamin = k;
        } else
        {
            break;
        }
    }
    UInt32 ilbdamax = nTpl - 1;

    for (UInt32 k = nTpl - 1; k > 0; k--)
    {
        Float64 _spcLambda = spcLambda[nSpc - 1] + spcLambdaMargin;
        Float64 _tplLambda = tplLambda[k] - tplLambdaMargin;
        Float64 _z = (_spcLambda - _tplLambda) / _tplLambda;
        if (_z < redshiftrange.GetBegin() - redshiftMargin)
        {
            ilbdamax = k;
        } else
        {
            break;
        }
    }

    if (ilbdamin < 0 || ilbdamin > nTpl - 1)
    {
        Log.LogError("  Operator-ChisquareLog: Problem with tpl indexes for "
                     "zranges, found lbdamin=%d",
                     ilbdamin);
        return TInt32Range(-1, -1);
    }
    if (ilbdamax < 0 || ilbdamax > nTpl - 1)
    {
        Log.LogError("  Operator-ChisquareLog: Problem with tpl indexes for "
                     "zranges, found ilbdamax=%d",
                     ilbdamax);
        return TInt32Range(-1, -1);
    }
    if (ilbdamin > ilbdamax)
    {
        Log.LogError("  Operator-ChisquareLog: Problem with tpl indexes for "
                     "zranges, found ilbdamin=%d > ilbdamax=%d",
                     ilbdamin, ilbdamax);
        return TInt32Range(-1, -1);
    }

    return TInt32Range(ilbdamin, ilbdamax);
}

/**
 * \brief COperatorChiSquareLogLambda::Compute
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
 **/
std::shared_ptr<COperatorResult> COperatorChiSquareLogLambda::Compute(
    const CSpectrum &spectrum, const CTemplate &tpl,
    const TFloat64Range &lambdaRange, const TFloat64List &redshifts,
    Float64 overlapThreshold, std::vector<CMask> additional_spcMasks,
    std::string opt_interp, Int32 opt_extinction, Int32 opt_dustFitting)
{
    Log.LogInfo(
        "  Operator-ChisquareLog: starting computation for template: %s",
        tpl.GetName().c_str());

    if (opt_dustFitting && m_ismCorrectionCalzetti->calzettiInitFailed)
    {
        Log.LogError("  Operator-ChisquareLog: no calzetti calib. file "
                     "loaded... aborting!");
        return NULL;
    }

    if (opt_extinction && m_igmCorrectionMeiksin->meiksinInitFailed)
    {
        Log.LogError("  Operator-ChisquareLog: no meiksin calib. file "
                     "loaded... aborting!");
        return NULL;
    }

    if (spectrum.GetSpectralAxis().IsInLinearScale() == false ||
        tpl.GetSpectralAxis().IsInLinearScale() == false)
    {
        Log.LogError("  Operator-ChisquareLog: input spectrum or template are "
                     "not in log scale (ignored)");
        // return NULL;
    }

    // assuming that lambdarange is strictly included in the spectrum spectral
    // range
    if (std::max(spectrum.GetSpectralAxis()[0], lambdaRange.GetBegin()) !=
        lambdaRange.GetBegin())
    {
        Log.LogError(
            "  Operator-ChisquareLog: lambdarange is not strictly included in "
            "the spectrum spectral range (min lambdarange)");
    }
    if (std::min(spectrum.GetSpectralAxis()
                     [spectrum.GetSpectralAxis().GetSamplesCount() - 1],
                 lambdaRange.GetEnd()) != lambdaRange.GetEnd())
    {
        Log.LogError(
            "  Operator-ChisquareLog: lambdarange is not strictly included in "
            "the spectrum spectral range (max lambdarange)");
    }

    // sort the redshift and keep track of the indexes
    TFloat64List sortedRedshifts;
    TFloat64List
        sortedIndexes; // used for the correspondence between input redshifts
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

    if (overlapThreshold < 1.0)
    {
        Log.LogError("  Operator-ChisquareLog: overlap threshold can't be "
                     "lower than 1.0");
        return NULL;
    }
    // check that the overlap is >1. for all sortedRedshifts
    Bool overlapFull = true;
    if (lambdaRange.GetBegin() <
        tpl.GetSpectralAxis()[0] *
            (1 + sortedRedshifts[sortedRedshifts.size() - 1]))
    {
        overlapFull = false;
    }
    if (lambdaRange.GetEnd() >
        tpl.GetSpectralAxis()[tpl.GetSpectralAxis().GetSamplesCount() - 1] *
            (1 + sortedRedshifts[0]))
    {
        overlapFull = false;
    }
    if (!overlapFull)
    {
        Log.LogError("  Operator-ChisquareLog: overlap found to be lower than "
                     "1.0 for this redshift range");
        Log.LogError(
            "  Operator-ChisquareLog: for zmin=%f, tpl.lbdamax is %f (should "
            "be >%f)",
            sortedRedshifts[0],
            tpl.GetSpectralAxis()[tpl.GetSpectralAxis().GetSamplesCount() - 1],
            lambdaRange.GetEnd() / (1 + sortedRedshifts[0]));
        Log.LogError("  Operator-ChisquareLog: for zmax=%f, tpl.lbdamin is %f "
                     "(should be <%f)",
                     sortedRedshifts[sortedRedshifts.size() - 1],
                     tpl.GetSpectralAxis()[0],
                     lambdaRange.GetBegin() /
                         (1 + sortedRedshifts[sortedRedshifts.size() - 1]));
        return NULL;
    }

    // Create/Retrieve the spectrum log-lambda spectral axis
    Int32 lbdaMinIdx =
        spectrum.GetSpectralAxis().GetIndexAtWaveLength(lambdaRange.GetBegin());
    Int32 lbdaMaxIdx =
        spectrum.GetSpectralAxis().GetIndexAtWaveLength(lambdaRange.GetEnd());
    if (lbdaMinIdx < 0 || lbdaMaxIdx < 0 || lbdaMaxIdx <= lbdaMinIdx)
    {
        Log.LogError("  Operator-ChisquareLog: problem while searching for "
                     "lambdarange observed spectral indexes");
    }

    Float64 loglbdaStep = log(spectrum.GetSpectralAxis()[lbdaMaxIdx]) -
                          log(spectrum.GetSpectralAxis()[lbdaMaxIdx - 1]);
    // Float64 loglbdaStep =
    // log(spectrum.GetSpectralAxis()[1])-log(spectrum.GetSpectralAxis()[0]);
    Float64 loglbdamin = log(lambdaRange.GetBegin());
    Float64 loglbdamax = log(lambdaRange.GetEnd());
    Int32 loglbdaCount = int((loglbdamax - loglbdamin) / loglbdaStep + 1);

    // Allocate the Log-rebined spectrum and mask
    CSpectrumFluxAxis &spectrumRebinedFluxAxis =
        m_spectrumRebinedLog.GetFluxAxis();
    CSpectrumSpectralAxis &spectrumRebinedSpectralAxis =
        m_spectrumRebinedLog.GetSpectralAxis();
    spectrumRebinedFluxAxis.SetSize(loglbdaCount);
    spectrumRebinedSpectralAxis.SetSize(loglbdaCount);
    m_mskRebinedLog.SetSize(loglbdaCount);
    m_errorRebinedLog.SetSize(m_spectrumRebinedLog.GetSampleCount());

    if (m_opt_spcrebin)
    {

        Log.LogInfo(
            "  Operator-ChisquareLog: Log-regular lambda resampling ON");

        if (verboseLogRebin)
        {
            Log.LogInfo(
                "  Operator-ChisquareLog: ORIGINAL grid spectrum count = %d",
                lbdaMaxIdx - lbdaMinIdx + 1);
            // Log.LogInfo("  Operator-ChisquareLog: Log-Rebin: zref = %f",
            // zRef); Log.LogInfo("  Operator-ChisquareLog: Log-Rebin: zStep =
            // %f", zStep);
            Log.LogInfo("  Operator-ChisquareLog: Log-Rebin: loglbdaStep = %f",
                        loglbdaStep);
            Log.LogInfo("  Operator-ChisquareLog: Log-Rebin: loglbdamin=%f : "
                        "loglbdamax=%f",
                        loglbdamin, loglbdamax);
            Log.LogInfo(
                "  Operator-ChisquareLog: Log-Rebin: lbdamin=%f : lbdamax=%f",
                exp(loglbdamin), exp(loglbdamax));
            Log.LogInfo("  Operator-ChisquareLog: Log-Rebin: loglbdaCount = %d",
                        loglbdaCount);
        }
        CSpectrumSpectralAxis targetSpectralAxis;
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
                spectrum.GetSpectralAxis()
                    [spectrum.GetSpectralAxis().GetSamplesCount() - 1];
        step = targetSpectralAxis[loglbdaCount - 1] -
               targetSpectralAxis[loglbdaCount - 2];
        if (delta > 0.0 && abs(delta) < step * 1e-5)
        {
            targetSpectralAxis[loglbdaCount - 1] =
                spectrum.GetSpectralAxis()
                    [spectrum.GetSpectralAxis().GetSamplesCount() - 1];
        }
        if (verboseLogRebin)
        {
            // zstep lower lbda range for this sampling
            Float64 dlambda_begin =
                (targetSpectralAxis[1] - targetSpectralAxis[0]);
            Float64 dz_zrangemin_begin =
                dlambda_begin /
                (targetSpectralAxis[0] / (1. + sortedRedshifts[0]));
            Float64 dz_zrangemax_begin =
                dlambda_begin /
                (targetSpectralAxis[0] /
                 (1. + sortedRedshifts[sortedRedshifts.size() - 1]));
            Log.LogInfo("  Operator-ChisquareLog: Log-Rebin: "
                        "dz_zrangemin_begin=%f : dz_zrangemax_begin=%f",
                        dz_zrangemin_begin, dz_zrangemax_begin);

            Float64 dlambda_end = (targetSpectralAxis[loglbdaCount - 1] -
                                   targetSpectralAxis[loglbdaCount - 2]);
            Float64 dz_zrangemin_end =
                dlambda_end / (targetSpectralAxis[loglbdaCount - 2] /
                               (1. + sortedRedshifts[0]));
            Float64 dz_zrangemax_end =
                dlambda_end /
                (targetSpectralAxis[loglbdaCount - 2] /
                 (1. + sortedRedshifts[sortedRedshifts.size() - 1]));
            Log.LogInfo("  Operator-ChisquareLog: Log-Rebin: "
                        "dz_zrangemin_end=%f : dz_zrangemax_end=%f",
                        dz_zrangemin_end, dz_zrangemax_end);
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
        if (enableVarianceWeightedRebin)
        {
            const TFloat64List &error = spectrum.GetFluxAxis().GetError();
            CSpectrumFluxAxis errorFluxAxis(spectrum.GetSampleCount());
            for (Int32 t = 0; t < spectrum.GetSampleCount(); t++)
            {
                errorFluxAxis[t] = error[t];
            }
            CSpectrumFluxAxis::RebinVarianceWeighted(
                spectrum.GetFluxAxis(), spectrum.GetSpectralAxis(),
                errorFluxAxis, targetSpectralAxis, spectrumRebinedFluxAxis,
                spectrumRebinedSpectralAxis, m_errorRebinedLog, "lin");
        } else
        {
            // rebin the spectrum
            TFloat64Range spcLbdaRange(exp(loglbdamin - 0.5 * loglbdaStep),
                                       exp(loglbdamax + 0.5 * loglbdaStep));
            Float64 *pfgBuffer_unused = 0;
            Float64 redshift_unused = 0.0;
            CSpectrumFluxAxis::Rebin2(
                spcLbdaRange, spectrum.GetFluxAxis(), pfgBuffer_unused,
                redshift_unused, spectrum.GetSpectralAxis(), targetSpectralAxis,
                spectrumRebinedFluxAxis, spectrumRebinedSpectralAxis,
                m_mskRebinedLog, rebinMethod);

            // rebin the variance
            const TFloat64List &error = spectrum.GetFluxAxis().GetError();
            CSpectrumFluxAxis errorFluxAxis(spectrum.GetSampleCount());
            for (Int32 t = 0; t < spectrum.GetSampleCount(); t++)
            {
                errorFluxAxis[t] = error[t];
            }
            CSpectrumSpectralAxis errorRebinedSpectralAxis;
            errorRebinedSpectralAxis.SetSize(loglbdaCount);
            CMask error_mskRebinedLog;
            error_mskRebinedLog.SetSize(loglbdaCount);

            CSpectrumFluxAxis::Rebin2(
                spcLbdaRange, errorFluxAxis, pfgBuffer_unused, redshift_unused,
                spectrum.GetSpectralAxis(), targetSpectralAxis,
                m_errorRebinedLog, errorRebinedSpectralAxis,
                error_mskRebinedLog, rebinMethod);
        }

        // put the error in the spectrum object: needed for cstLog estimation
        // (at least)
        for (Int32 t = 0; t < m_spectrumRebinedLog.GetSampleCount(); t++)
        {
            m_spectrumRebinedLog.GetFluxAxis().GetError()[t] =
                m_errorRebinedLog[t];
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
        if (verboseExportLogRebin)
        {
            // save rebinned data
            FILE *f = fopen("loglbda_rebinlog_errorlogrebin_dbg.txt", "w+");
            for (Int32 t = 0; t < m_errorRebinedLog.GetSamplesCount(); t++)
            {
                fprintf(f, "%f\t%e\n",
                        m_spectrumRebinedLog.GetSpectralAxis()[t],
                        m_errorRebinedLog[t]);
            }
            fclose(f);
        }

    } else // check that the raw spectrum grid is in log-regular grid
    {
        Log.LogInfo(
            "  Operator-ChisquareLog: Log-regular lambda resampling OFF");

        if (verboseLogRebin)
        {
            Log.LogInfo(
                "  Operator-ChisquareLog: lbda raw min_idx=%d, max_idx=%d",
                lbdaMinIdx, lbdaMaxIdx);
        }

        bool logRegLbdaCheck = true;
        Float64 relativeLogLbdaStepTol = 1e-1;
        Float64 maxAbsRelativeError = 0.0; //
        Float64 lbda1 = spectrum.GetSpectralAxis()[lbdaMinIdx];
        for (Int32 t = 1; t < loglbdaCount; t++)
        {
            Float64 lbda2 = spectrum.GetSpectralAxis()[lbdaMinIdx + t];
            Float64 _loglbdaStep = log(lbda2) - log(lbda1);

            Float64 relativeErrAbs =
                std::abs((_loglbdaStep - loglbdaStep) / loglbdaStep);
            // if(verboseLogRebin)
            //{
            //    Log.LogInfo("  Operator-ChisquareLog: _loglbdastep = %f,
            //    relativeErrAbs = %f", _loglbdaStep, relativeErrAbs);
            //}

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
        Log.LogInfo("  Operator-ChisquareLog: Log-regular lambda check: max "
                    "Abs Relative Error (log lbda step)= %f",
                    maxAbsRelativeError);
        Log.LogInfo("  Operator-ChisquareLog: Log-regular lambda check (for "
                    "rel-tol=%f): (1=success, 0=fail) CHECK = %d",
                    relativeLogLbdaStepTol, logRegLbdaCheck);
        if (!logRegLbdaCheck)
        {
            Log.LogError(
                "  Operator-ChisquareLog: Log-regular lambda check FAILED");
        }
        // prepare data for log-regular computation
        for (Int32 t = 0; t < loglbdaCount; t++)
        {
            spectrumRebinedFluxAxis[t] = spectrum.GetFluxAxis()[lbdaMinIdx + t];
            spectrumRebinedSpectralAxis[t] =
                spectrum.GetSpectralAxis()[lbdaMinIdx + t];
            m_errorRebinedLog[t] =
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

    // Create the Template Log-Rebined spectral axis
    {

        // The template grid has to be aligned with the spectrum log-grid (will
        // use the spc first element)
        // Float64 tpl_raw_loglbdamin = log(tpl.GetSpectralAxis()[0]); //full
        // tpl lambda range Float64 tpl_raw_loglbdamax =
        // log(tpl.GetSpectralAxis()[tpl.GetSpectralAxis().GetSamplesCount()-1]);
        // //full tpl lambda range
        Float64 tplloglbdaStep =
            log(tpl.GetSpectralAxis()[1]) -
            log(tpl.GetSpectralAxis()[0]); // warning : considering constant
                                           // dlambda for the input template
        Float64 roundingErrorsMargin = tplloglbdaStep / 2.0; //
        Float64 tpl_raw_loglbdamin =
            log(lambdaRange.GetBegin() /
                (1.0 + sortedRedshifts[sortedRedshifts.size() - 1])) -
            roundingErrorsMargin; // lambdarange cropped to useful range given
                                  // zmin
        Float64 tpl_raw_loglbdamax =
            log(lambdaRange.GetEnd() / (1.0 + sortedRedshifts[0])) +
            roundingErrorsMargin; // lambdarange cropped to useful range given
                                  // zmin

        Float64 tpl_tgt_loglbdamin = loglbdamin;
        if (tpl_raw_loglbdamin <= loglbdamin)
        {
            while (tpl_tgt_loglbdamin - loglbdaStep >= tpl_raw_loglbdamin)
            {
                tpl_tgt_loglbdamin -= loglbdaStep;
            }
        } else
        {
            while (tpl_tgt_loglbdamin + loglbdaStep <= tpl_raw_loglbdamin)
            {
                tpl_tgt_loglbdamin += loglbdaStep;
            }
        }
        Float64 tpl_tgt_loglbdamax = loglbdamax;
        if (tpl_raw_loglbdamax <= loglbdamax)
        {
            while (tpl_tgt_loglbdamax - loglbdaStep >= tpl_raw_loglbdamax)
            {
                tpl_tgt_loglbdamax -= loglbdaStep;
            }
        } else
        {
            while (tpl_tgt_loglbdamax + loglbdaStep <= tpl_raw_loglbdamax)
            {
                tpl_tgt_loglbdamax += loglbdaStep;
            }
        }
        Int32 tpl_loglbdaCount =
            int((tpl_tgt_loglbdamax - tpl_tgt_loglbdamin) / loglbdaStep + 1);
        if (verboseLogRebin)
        {
            Log.LogInfo("  Operator-ChisquareLog: Log-Rebin: tpl loglbdamin=%f "
                        ": loglbdamax=%f",
                        tpl_tgt_loglbdamin, tpl_tgt_loglbdamax);
            Log.LogInfo("  Operator-ChisquareLog: Log-Rebin: tpl lbdamin=%f : "
                        "lbdamax=%f",
                        exp(tpl_tgt_loglbdamin), exp(tpl_tgt_loglbdamax));
            Log.LogInfo(
                "  Operator-ChisquareLog: Log-Rebin: tpl loglbdaCount = %d",
                tpl_loglbdaCount);
        }
        // todo: check that the coverage is ok with teh current tgtTplAxis ?

        // rebin the template
        CSpectrumSpectralAxis tpl_targetSpectralAxis;
        tpl_targetSpectralAxis.SetSize(tpl_loglbdaCount);
        for (Int32 k = 0; k < tpl_loglbdaCount; k++)
        {
            tpl_targetSpectralAxis[k] =
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
        delta =
            tpl_targetSpectralAxis[loglbdaCount - 1] -
            tpl.GetSpectralAxis()[tpl.GetSpectralAxis().GetSamplesCount() - 1];
        step = tpl_targetSpectralAxis[loglbdaCount - 1] -
               tpl_targetSpectralAxis[loglbdaCount - 2];
        if (delta > 0.0 && abs(delta) < step * 1e-5)
        {
            tpl_targetSpectralAxis[loglbdaCount - 1] =
                tpl.GetSpectralAxis()[tpl.GetSpectralAxis().GetSamplesCount() -
                                      1];
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
            Log.LogError(
                "  Operator-ChisquareLog: Log-Rebin: tpl rebin error. Target "
                "MIN lbda value=%f, input tpl min lbda value=%f",
                tpl_targetSpectralAxis[0], tpl.GetSpectralAxis()[0]);
            Log.LogError("  Operator-ChisquareLog: Log-Rebin: tpl rebin error. "
                         "Extend your input template wavelength range or "
                         "modify the processing parameter <lambdarange>");
        }
        // precision problem due to exp/log - end
        delta =
            tpl_targetSpectralAxis[tpl_targetSpectralAxis.GetSamplesCount() -
                                   1] -
            tpl.GetSpectralAxis()[tpl.GetSpectralAxis().GetSamplesCount() - 1];
        if (delta > 0.0)
        {
            Log.LogError(
                "  Operator-ChisquareLog: Log-Rebin: tpl rebin error. Target "
                "MAX lbda value=%f, input tpl max lbda value=%f",
                tpl_targetSpectralAxis
                    [tpl_targetSpectralAxis.GetSamplesCount() - 1],
                tpl.GetSpectralAxis()[tpl.GetSpectralAxis().GetSamplesCount() -
                                      1]);
            Log.LogError("  Operator-ChisquareLog: Log-Rebin: tpl rebin error. "
                         "Extend your input template wavelength range or "
                         "modify the processing parameter <lambdarange>");
        }

        Float64 *pfgBuffer_unused = 0;
        Float64 redshift_unused = 0.0;
        CSpectrumFluxAxis &templateRebinedFluxAxis =
            m_templateRebinedLog.GetFluxAxis();
        CSpectrumSpectralAxis &templateRebinedSpectralAxis =
            m_templateRebinedLog.GetSpectralAxis();
        templateRebinedFluxAxis.SetSize(tpl_loglbdaCount);
        templateRebinedSpectralAxis.SetSize(tpl_loglbdaCount);
        CMask tpl_mskRebinedLog;
        tpl_mskRebinedLog.SetSize(tpl_loglbdaCount);
        CSpectrumFluxAxis::Rebin2(
            tplLbdaRange, tpl.GetFluxAxis(), pfgBuffer_unused, redshift_unused,
            tpl.GetSpectralAxis(), tpl_targetSpectralAxis,
            templateRebinedFluxAxis, templateRebinedSpectralAxis,
            tpl_mskRebinedLog, rebinMethod);
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
        // Log.LogError("  Operator-ChisquareLog: FitAllz opt_extinction not
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
    if (opt_dustFitting)
    {
        // Log.LogError("  Operator-ChisquareLog: FitAllz opt_dustFitting not
        // validated yet...");
        Int32 nISMCoeffs = m_ismCorrectionCalzetti->GetNPrecomputedDustCoeffs();
        ismEbmvCoeffs.resize(nISMCoeffs);
        for (Int32 kism = 0; kism < nISMCoeffs; kism++)
        {
            ismEbmvCoeffs[kism] = kism;
        }
    } else
    {
        ismEbmvCoeffs.clear();
    }

    std::shared_ptr<CChisquareResult> result =
        std::shared_ptr<CChisquareResult>(new CChisquareResult());
    result->Init(sortedRedshifts.size(),
                 std::max((Int32)ismEbmvCoeffs.size(), 1),
                 std::max((Int32)igmMeiksinCoeffs.size(), 1));
    result->Redshifts = sortedRedshifts;

    // WARNING: no additional masks coded for use as of 2017-06-13
    if (additional_spcMasks.size() != 0)
    {
        Log.LogError("  Operator-ChisquareLog: No additional masks used. "
                     "Feature not coded for this log-lambda operator!)");
    }

    //*
    Int32 retFit =
        FitAllz(lambdaRange, result, igmMeiksinCoeffs, ismEbmvCoeffs);
    if (retFit != 0)
    {
        Log.LogError("  Operator-ChisquareLog: FitAllz failed with error %d",
                     retFit);
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
        Log.LogInfo("  Operator-ChisquareLog: overlap warning for %s: "
                    "minz=%.3f, maxz=%.3f",
                    tpl.GetName().c_str(), overlapValidInfZ, overlapValidSupZ);
    }

    // estimate CstLog for PDF estimation
    result->CstLog = EstimateLikelihoodCstLog(
        m_spectrumRebinedLog, lambdaRange); // 0.0;//Todo: check how to estimate
                                            // that value for loglambda//

    // extrema
    Int32 extremumCount = 10;
    if (result->Redshifts.size() > extremumCount)
    {
        TPointList extremumList;
        TFloat64Range redshiftsRange(
            result->Redshifts[0],
            result->Redshifts[result->Redshifts.size() - 1]);
        CExtremum extremum(redshiftsRange, extremumCount, true);
        extremum.Find(result->Redshifts, result->ChiSquare, extremumList);

        //*
        // Refine Extremum with a second maximum search around the z candidates:
        // This corresponds to the finer xcorrelation in EZ Pandora (in
        // standard_DP fctn in SolveKernel.py)
        Float64 radius = 0.001;
        for (Int32 i = 0; i < extremumList.size(); i++)
        {
            Float64 x = extremumList[i].X;
            Float64 left_border = max(redshiftsRange.GetBegin(), x - radius);
            Float64 right_border = min(redshiftsRange.GetEnd(), x + radius);

            TPointList extremumListFine;
            TFloat64Range rangeFine = TFloat64Range(left_border, right_border);
            CExtremum extremumFine(rangeFine, 1, true);
            extremumFine.Find(result->Redshifts, result->ChiSquare,
                              extremumListFine);
            if (extremumListFine.size() > 0)
            {
                extremumList[i] = extremumListFine[0];
            }
        }
        //*/
        // store extrema results
        result->Extrema.resize(extremumCount);
        for (Int32 i = 0; i < extremumList.size(); i++)
        {

            result->Extrema[i] = extremumList[i].X;
        }

    } else
    {
        // store extrema results
        result->Extrema.resize(result->Redshifts.size());
        TFloat64List tmpX;
        TFloat64List tmpY;
        for (Int32 i = 0; i < result->Redshifts.size(); i++)
        {
            tmpX.push_back(result->Redshifts[i]);
            tmpY.push_back(result->ChiSquare[i]);
        }
        // sort the results by merit
        CQuickSort<Float64> sort;
        vector<Int32> sortedIndexes(result->Redshifts.size());
        sort.SortIndexes(tmpY.data(), sortedIndexes.data(),
                         sortedIndexes.size());
        for (Int32 i = 0; i < result->Redshifts.size(); i++)
        {
            result->Extrema[i] = tmpX[sortedIndexes[i]];
        }
    }

    return result;
}

/**
 * \brief this function estimates the likelihood_cstLog term withing the
 *wavelength range
 **/
Float64 COperatorChiSquareLogLambda::EstimateLikelihoodCstLog(
    const CSpectrum &spectrum, const TFloat64Range &lambdaRange)
{
    const CSpectrumSpectralAxis &spcSpectralAxis = spectrum.GetSpectralAxis();
    const TFloat64List &error = spectrum.GetFluxAxis().GetError();

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

const Float64 *COperatorChiSquareLogLambda::getDustCoeff(Float64 dustCoeff,
                                                         Float64 maxLambda)
{
    return m_ismCorrectionCalzetti->getDustCoeff(dustCoeff, maxLambda);
}

const Float64 *COperatorChiSquareLogLambda::getMeiksinCoeff(Int32 meiksinIdx,
                                                            Float64 redshift,
                                                            Float64 maxLambda)
{
    return m_igmCorrectionMeiksin->getMeiksinCoeff(meiksinIdx, redshift,
                                                   maxLambda);
}

void COperatorChiSquareLogLambda::enableSpcLogRebin(Bool enable)
{
    m_opt_spcrebin = enable;
    Log.LogInfo("  Operator-ChisquareLog: Spectrum REBIN-LOG enabled=%d",
                m_opt_spcrebin);
}
