#include <RedshiftLibrary/operator/chisquareloglambda.h>

#include <RedshiftLibrary/spectrum/axis.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/tools.h>
#include <RedshiftLibrary/common/mask.h>
#include <RedshiftLibrary/operator/chisquareresult.h>
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
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <fftw3.h>

#include <assert.h>

#define NOT_OVERLAP_VALUE NAN
#include <stdio.h>

namespace bfs = boost::filesystem;
using namespace NSEpic;
using namespace std;

COperatorChiSquareLogLambda::COperatorChiSquareLogLambda( std::string calibrationPath )
{   
    //ISM
    m_ismCorrectionCalzetti = new CSpectrumFluxCorrectionCalzetti();
    m_ismCorrectionCalzetti->Init(calibrationPath);

    //IGM
    m_igmCorrectionMeiksin = new CSpectrumFluxCorrectionMeiksin();
    m_igmCorrectionMeiksin->Init(calibrationPath);

}

COperatorChiSquareLogLambda::~COperatorChiSquareLogLambda()
{
}


Int32 COperatorChiSquareLogLambda::EstimateXtYSlow(const Float64* X, const Float64* Y, UInt32 nX, UInt32 nShifts, std::vector<Float64>& XtY)
{
    XtY.resize(nShifts);

    Float64 xty = 0.0;
    for(Int32 k=0; k<nShifts; k++)
    {
        xty = 0.0;
        for(Int32 j=0; j<nX; j++)
        {
            xty += X[j]*Y[j+k];
        }
        XtY[k]=xty;
    }
    return 0;
}

//only works for mtm, Y=model^2, X=1.
Int32 COperatorChiSquareLogLambda::EstimateMtMFast(const Float64 *X, const Float64 *Y, UInt32 nX, UInt32 nShifts, std::vector<Float64>& XtY)
{
    XtY.resize(nShifts);

    Float64 xty = 0.0;
    for(Int32 j=0; j<nX; j++)
    {
        xty += X[j]*Y[j];
    }
    XtY[0
            ]=xty;

    for(Int32 k=1; k<nShifts; k++)
    {
        xty = XtY[k-1];
        xty -= X[0]*Y[0+k-1];
        xty += X[nX-1]*Y[nX-1+k-1];

        XtY[k]=xty;
    }
    return 0;
}

Int32 COperatorChiSquareLogLambda::EstimateXtY(const Float64* X, const Float64* Y, UInt32 nx, UInt32 ny, UInt32 nshifts, std::vector<Float64>& XtY)
{
    bool verbose = false;
    //Processing the FFT
    Int32 nSpc = nx;
    Int32 nTpl = ny;
    //Int32 nPadded = (Int32)std::max((Float64)(ny*2.0-1), (Float64)(nx*2.0-1));
    //Int32 nPadded = nshifts;
    Int32 nPadded = (Int32)(nTpl*2.0);

    /*
    //next power of two
    Float64 maxnsamples = (Int32)std::max((Float64)(ny), (Float64)(nx));;
    int power = 1;
    while(std::pow(2, power) < maxnsamples)
    {
        power*=2;
    }
    Int32 nPadded = std::pow(2, power);
    //*/
    Int32 nPadBeforeSpc = nPadded-nSpc;//(Int32)nPadded/2.0;
    Int32 nPadBeforeTpl = 0;

    Log.LogInfo("ChisquareLog, FitAllz: Processing spc-fft with n=%d, padded to n=%d", nSpc, nPadded);
    fftw_complex *inSpc;
    fftw_complex *outSpc;
    inSpc = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nPadded);
    outSpc = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nPadded);
    fftw_plan pSpc = fftw_plan_dft_1d(nPadded, inSpc, outSpc, FFTW_FORWARD, FFTW_ESTIMATE);
    if(inSpc==0)
    {
        Log.LogError("ChisquareLog, FitAllz: Unable to allocate inSpc");
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

    for(Int32 k=0; k<nPadBeforeSpc; k++)
    {
        inSpc[k][0] = 0.0;
        inSpc[k][1] = 0.0;
    }
    for(Int32 k=nPadBeforeSpc; k<nPadBeforeSpc+nSpc; k++)
    {
        inSpc[k][0] = X[nSpc-(k-nPadBeforeSpc)];//X[k-nPadBeforeSpc];
        inSpc[k][1] = 0.0;
    }
    for(Int32 k=nPadBeforeSpc+nSpc; k<nPadded; k++)
    {
        inSpc[k][0] = 0.0;
        inSpc[k][1] = 0.0;
    }
    if(verbose)
    {
        // save spc-input data
        FILE* f_fftinput = fopen( "loglbda_fitallz_xfftInput_dbg.txt", "w+" );
        for( Int32 t=0;t<nPadded;t++)
        {
            fprintf( f_fftinput, "%f\t%e\n", (Float64)t, inSpc[t][0]);
        }
        fclose( f_fftinput );
    }
    if(outSpc==0)
    {
        Log.LogError("ChisquareLog, FitAllz: Unable to allocate outSpc");
    }
    fftw_execute(pSpc);
    Log.LogInfo("ChisquareLog, FitAllz: spc-fft done");
    if(verbose)
    {
        // save spc-fft data
        FILE* f_fftoutput = fopen( "loglbda_fitallz_xfftOutput_dbg.txt", "w+" );
        for( Int32 t=0;t<nPadded;t++)
        {
            fprintf( f_fftoutput, "%f\t%e\n", (Float64)t, outSpc[t][0]*outSpc[t][0]+outSpc[t][1]*outSpc[t][1]);
        }
        fclose( f_fftoutput );
    }


    Log.LogInfo("ChisquareLog, FitAllz: Processing tpl-fft with n=%d, padded to n=%d", nTpl, nPadded);
    fftw_complex *inTpl;
    fftw_complex *inTpl_padded;
    fftw_complex *outTpl;
    inTpl_padded = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nPadded);
    inTpl = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nPadded);
    outTpl = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nPadded);
    fftw_plan pTpl = fftw_plan_dft_1d(nPadded, inTpl, outTpl, FFTW_FORWARD, FFTW_ESTIMATE);
    if(inTpl_padded==0)
    {
        Log.LogError("ChisquareLog, FitAllz: Unable to allocate inTpl_padded");
    }
    if(inTpl==0)
    {
        Log.LogError("ChisquareLog, FitAllz: Unable to allocate inTpl");
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
    for(Int32 k=0; k<nPadBeforeTpl; k++)
    {
        inTpl_padded[k][0] = 0.0;
        inTpl_padded[k][1] = 0.0;
    }
    for(Int32 k=nPadBeforeTpl; k<nPadBeforeTpl+nTpl; k++)
    {
        inTpl_padded[k][0] = Y[k-nPadBeforeTpl];
        inTpl_padded[k][1] = 0.0;
    }
    for(Int32 k=nPadBeforeTpl+nTpl; k<nPadded; k++)
    {
        inTpl_padded[k][0] = 0.0;
        inTpl_padded[k][1] = 0.0;
    }

    for(Int32 k=0; k<nPadded; k++)
    {
        inTpl[k][0] = inTpl_padded[k][0];//inTpl_padded[nPadded-1-k][0];//
        inTpl[k][1] = 0.0;
    }
    if(outTpl==0)
    {
        Log.LogError("ChisquareLog, FitAllz: Unable to allocate outTpl");
    }
    if(verbose)
    {
        // save tpl input data
        FILE* f_fftinput = fopen( "loglbda_fitallz_yfftInput_dbg.txt", "w+" );
        for( Int32 t=0;t<nPadded;t++)
        {
            fprintf( f_fftinput, "%f\t%e\n", (Float64)t, inTpl[t][0]);
        }
        fclose( f_fftinput );
    }
    fftw_execute(pTpl);
    Log.LogInfo("ChisquareLog, FitAllz: tpl-fft done");
    if(verbose)
    {
        // save tpl-fft data
        FILE* f_fftoutput = fopen( "loglbda_fitallz_yfftOutput_dbg.txt", "w+" );
        for( Int32 t=0;t<nPadded;t++)
        {
            fprintf( f_fftoutput, "%f\t%e\n", (Float64)t, outTpl[t][0]*outTpl[t][0]+outTpl[t][1]*outTpl[t][1]);
        }
        fclose( f_fftoutput );
    }


    //Multiplying the FFT outputs
    fftw_complex* outCombined = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nPadded);
    fftw_complex* inCombined = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nPadded);
    fftw_plan pBackward = fftw_plan_dft_1d(nPadded, outCombined, inCombined, FFTW_BACKWARD, FFTW_ESTIMATE);
    if(outCombined==0)
    {
        Log.LogError("ChisquareLog, FitAllz: Unable to allocate outCombined");
    }
    for(Int32 k=0; k<nPadded; k++)
    {
        outCombined[k][0] = (outTpl[k][0]*outSpc[k][0]-outTpl[k][1]*outSpc[k][1]);
        outCombined[k][1] = (outTpl[k][0]*outSpc[k][1]+outTpl[k][1]*outSpc[k][0]);
        //Y conjugate
        //outCombined[k][0] = (outTpl[k][0]*outSpc[k][0]+outTpl[k][1]*outSpc[k][1]);
        //outCombined[k][1] = (outTpl[k][0]*outSpc[k][1]-outTpl[k][1]*outSpc[k][0]);
    }
    if(inCombined==0)
    {
        Log.LogError("ChisquareLog, FitAllz: Unable to allocate inCombined");
    }
    fftw_execute(pBackward);
    Log.LogInfo("ChisquareLog, FitAllz: backward-fft done");
    XtY.resize(nshifts);
    Int32 offsetSamples = 0.0;//nPadded/2.0;//nSpc;//(nTpl+nSpc);
    for (Int32 k=0; k<nshifts; k++)
    {
        //XtY[k] = inCombined[k+(Int32)(nPadded/2.0)][0]/nPadded;
        Int32 IndexCyclic = k+offsetSamples;
        if(IndexCyclic<0.0)
        {
            IndexCyclic+=nPadded;
        }else if(IndexCyclic>=nPadded)
        {
            IndexCyclic-=nPadded;
        }

        XtY[k] = inCombined[IndexCyclic][0]/nPadded;
    }

    if(verbose)
    {
        // save rebinned data
        FILE* f_fftoutput = fopen( "loglbda_fitallz_xtyfft-output_dbg.txt", "w+" );
        for( Int32 t=0;t<nshifts;t++)
        {
            fprintf( f_fftoutput, "%f\t%e\n", (Float64)t, XtY[t]);
        }
        fclose( f_fftoutput );
    }



    fftw_destroy_plan(pSpc);
    fftw_free(inSpc); fftw_free(outSpc);
    fftw_destroy_plan(pTpl);
    fftw_free(inTpl_padded);
    fftw_free(inTpl); fftw_free(outTpl);
    fftw_destroy_plan(pBackward);

    return 0;
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
Int32 COperatorChiSquareLogLambda::FitAllz(const TFloat64Range& lambdaRange,
                                           std::shared_ptr<CChisquareResult> result,
                                           std::vector<Int32> igmMeiksinCoeffs,
                                           std::vector<Int32> ismEbmvCoeffs,
                                           CMask spcMaskAdditional)
{
    bool verboseExportFitAllz = false;
    CSpectrumFluxAxis& spectrumRebinedFluxAxis = m_spectrumRebinedLog.GetFluxAxis();
    const Float64* error = spectrumRebinedFluxAxis.GetError();
    CSpectrumSpectralAxis& spectrumRebinedSpectralAxis = m_spectrumRebinedLog.GetSpectralAxis();
    CSpectrumFluxAxis& tplRebinedFluxAxis = m_templateRebinedLog.GetFluxAxis();
    CSpectrumSpectralAxis& tplRebinedSpectralAxis = m_templateRebinedLog.GetSpectralAxis();



    //prepare list of redshifts that need a full lbdarange lst-square calculation
    std::vector<Int32> zindexesFullLstSquare;
    zindexesFullLstSquare.push_back(0); //first index is always a mandatory full Lstsq Calculation case
    std::vector<Float64> zlistsegments = m_igmCorrectionMeiksin->GetSegmentsStartRedshiftList();
    for(Int32 k=0; k<zlistsegments.size(); k++)
    {
        for (Int32 i=0;i<result->Redshifts.size()-1;i++)
        {
            if(zlistsegments[k]>=result->Redshifts[i] && zlistsegments[k]<result->Redshifts[i+1])
            {
                zindexesFullLstSquare.push_back(i);
                //Log.LogInfo("ChisquareLog, zindexesFullLstSquare: index found for zlistsegments_k = %f, redshift_i = %f, i = %d", zlistsegments[k], i);
                break; //index found, go to next zlistsegments item
            }
        }
    }
    zindexesFullLstSquare.erase( std::unique( zindexesFullLstSquare.begin(), zindexesFullLstSquare.end() ), zindexesFullLstSquare.end() );
    if(verboseExportFitAllz)
    {
        Log.LogInfo("ChisquareLog, FitAllz: indexes for full LstSquare calculation, count = %d", zindexesFullLstSquare.size());

        Log.LogInfo("ChisquareLog, FitAllz: spc[0] = %f", spectrumRebinedSpectralAxis[0]);
        Log.LogInfo("ChisquareLog, FitAllz: spc[max] = %f", spectrumRebinedSpectralAxis[spectrumRebinedSpectralAxis.GetSamplesCount()-1]);
        Log.LogInfo("ChisquareLog, FitAllz: tpl[0]*zmax = %f", tplRebinedSpectralAxis[0]*(1.0+result->Redshifts[result->Redshifts.size()-1]));
        Log.LogInfo("ChisquareLog, FitAllz: tpl[max]*zmin = %f", tplRebinedSpectralAxis[tplRebinedSpectralAxis.GetSamplesCount()-1]*(1+result->Redshifts[0]));
    }

    //
    Int32 nshifts = tplRebinedSpectralAxis.GetSamplesCount()-spectrumRebinedSpectralAxis.GetSamplesCount();
    //Int32 nshifts = tplRebinedSpectralAxis.GetSamplesCount()*2.0;

    //prepare z array
    std::vector<Float64> z_vect(nshifts, 0.0);
    for( Int32 t=0;t<z_vect.size();t++)
    {
        z_vect[t] = (spectrumRebinedSpectralAxis[0]-tplRebinedSpectralAxis[t])/tplRebinedSpectralAxis[t];
    }

    //best fit data
    std::vector<Float64> bestChi2(nshifts, 0.0);


    // Estimate DtD
    //    std::vector<Float64> dtd_vec;
    //    CSpectrumFluxAxis spcSquareFluxAxis(spectrumRebinedFluxAxis);
    //    for(Int32 k=0; k<spcSquareFluxAxis.GetSamplesCount(); k++)
    //    {
    //        spcSquareFluxAxis[k] = spectrumRebinedFluxAxis[k]*spectrumRebinedFluxAxis[k];
    //    }
    //    EstimateXtY(OneFluxAxis, spcSquareFluxAxis, dtd_vec);
    Float64 dtd = 0.0;
    Float64 inv_err2 = 1.0;
    for(Int32 j=0; j<spectrumRebinedFluxAxis.GetSamplesCount(); j++)
    {
        inv_err2 = 1.0/(error[j]*error[j]);
        dtd += spectrumRebinedFluxAxis[j]*spectrumRebinedFluxAxis[j]*inv_err2;
    }

    //prepare arrays
    Float64* spcRebinedFlux = spectrumRebinedFluxAxis.GetSamples();
    Float64* spcRebinedFluxOverErr2 = new Float64 [(int)spectrumRebinedFluxAxis.GetSamplesCount()]();
    Float64* oneSpcRebinedFluxOverErr2 = new Float64 [(int)spectrumRebinedFluxAxis.GetSamplesCount()]();
    for(Int32 j=0; j<spectrumRebinedFluxAxis.GetSamplesCount(); j++)
    {
        inv_err2 = 1.0/(error[j]*error[j]);
        spcRebinedFluxOverErr2[j] = spectrumRebinedFluxAxis[j]*inv_err2;
        oneSpcRebinedFluxOverErr2[j] = inv_err2;
    }

    Float64* tplRebinedFluxRaw = tplRebinedFluxAxis.GetSamples();
    Float64* tplRebinedFlux = new Float64 [(int)tplRebinedFluxAxis.GetSamplesCount()]();
    Float64* tpl2RebinedFlux = new Float64 [(int)tplRebinedFluxAxis.GetSamplesCount()]();
    for(Int32 j=0; j<tplRebinedFluxAxis.GetSamplesCount(); j++)
    {
        tplRebinedFlux[j] = tplRebinedFluxRaw[j];
        tpl2RebinedFlux[j] = tplRebinedFlux[j]*tplRebinedFlux[j];
    }


    bool enableIGM = true;
    UInt32 nIGM = igmMeiksinCoeffs.size();
    if(nIGM==0)
    {
        nIGM=1;
        enableIGM = false;
    }
    for(Int32 kIGM=0; kIGM<nIGM; kIGM++)
    {
//        if(enableIGM)
//        {
//            Int32 kmeiksin = igmMeiksinCoeffs[kIGM];
//        }

        bool enableISM = true;
        UInt32 nISM = ismEbmvCoeffs.size();
        if(nISM==0)
        {
            nISM=1;
            enableISM = false;
        }
        for(Int32 kISM=0; kISM<nISM; kISM++)
        {
            if(verboseExportFitAllz && enableIGM)
            {
                Log.LogInfo("ChisquareLog, FitAllz: IGM index=%d", kIGM);
            }
            if(verboseExportFitAllz && enableISM)
            {
                Log.LogInfo("ChisquareLog, FitAllz: ISM index =%d", kISM);
            }


            if(enableISM)
            {
                Int32 kDustCalzetti = ismEbmvCoeffs[kISM];

                //correct tplRebinedFlux
                //correct tpl2RebinedFlux
                Float64 restLambda;
                Float64 ebmvDustCoeff=1.0;
                for(Int32 j=0; j<tplRebinedFluxAxis.GetSamplesCount(); j++)
                {
                    if(enableISM)
                    {
                        //apply ism dust correction from Calzetti
                        restLambda = tplRebinedSpectralAxis[j];
                        ebmvDustCoeff = m_ismCorrectionCalzetti->getDustCoeff( kDustCalzetti, restLambda);
                    }

                    tplRebinedFlux[j] = tplRebinedFluxRaw[j]*ebmvDustCoeff;

                    tpl2RebinedFlux[j] = tplRebinedFlux[j]*tplRebinedFlux[j];
                }
            }

            // Estimate DtM
            std::vector<Float64> dtm_vec;
            EstimateXtY(spcRebinedFluxOverErr2, tplRebinedFlux, spectrumRebinedFluxAxis.GetSamplesCount(), tplRebinedFluxAxis.GetSamplesCount(), nshifts, dtm_vec);
            //EstimateXtYSlow(spcRebinedFluxOverErr2, tplRebinedFlux, spectrumRebinedFluxAxis.GetSamplesCount(), nshifts, dtm_vec);

            if(verboseExportFitAllz)
            {
                // save chi2 data
                FILE* f = fopen( "loglbda_dtm_dbg.txt", "w+" );
                for( Int32 t=0;t<dtm_vec.size();t++)
                {
                    fprintf( f, "%f\t%e\n", (Float64)t, dtm_vec[t]);
                }
                fclose( f );
            }

            // Estimate MtM
            std::vector<Float64> mtm_vec;
            EstimateXtY(oneSpcRebinedFluxOverErr2, tpl2RebinedFlux, spectrumRebinedFluxAxis.GetSamplesCount(),tplRebinedFluxAxis.GetSamplesCount(), nshifts, mtm_vec);
            //EstimateXtYSlow(oneSpcRebinedFluxOverErr2, tpl2RebinedFlux, spectrumRebinedFluxAxis.GetSamplesCount(), nshifts, mtm_vec);
            //EstimateMtMFast(oneSpcRebinedFluxOverErr2, tpl2RebinedFlux, spectrumRebinedFluxAxis.GetSamplesCount(), nshifts, mtm_vec);

            if(verboseExportFitAllz)
            {
                // save chi2 data
                FILE* f = fopen( "loglbda_mtm_dbg.txt", "w+" );
                for( Int32 t=0;t<mtm_vec.size();t++)
                {
                    fprintf( f, "%f\t%e\n", (Float64)t, mtm_vec[t]);
                }
                fclose( f );
            }



            if(verboseExportFitAllz)
            {
                Log.LogInfo("ChisquareLog, FitAllz: dtd = %e", dtd);
            }

            //return -1;

            // Estimate Chi2
            if( mtm_vec.size()!=dtm_vec.size())
            {
                Log.LogError("ChisquareLog, FitAllz: xty vectors sizes don't match: dtm size = %d, mtm size = %d", dtm_vec.size(), mtm_vec.size());
                return 2;
            }
            std::vector<Float64> chi2(dtm_vec.size(), DBL_MAX);
            for(Int32 k=0; k<dtm_vec.size(); k++)
            {
                if(mtm_vec[k]==0.0)
                {
                    chi2[k] = dtd;
                }else{
                    Float64 ampl = max(0.0, dtm_vec[k] / mtm_vec[k]);
                    chi2[k] = dtd - dtm_vec[k]*ampl;
                }
                //chi2[k] = dtm_vec[k];
            }

            for(Int32 k=0; k<dtm_vec.size(); k++)
            {
                bestChi2[k] = chi2[k];
            }

            if(verboseExportFitAllz)
            {
                Log.LogInfo("ChisquareLog, FitAllz: spc lbda 0 =%f", spectrumRebinedSpectralAxis[0]);
                Log.LogInfo("ChisquareLog, FitAllz: tpl lbda 0 =%f", tplRebinedSpectralAxis[0]);
                Float64 z_O = (spectrumRebinedSpectralAxis[0]-tplRebinedSpectralAxis[0])/tplRebinedSpectralAxis[0];
                Log.LogInfo("ChisquareLog, FitAllz: z 0 =%f", z_O);

                //        Float64 logstep = log(spectrumRebinedSpectralAxis[1])-log(spectrumRebinedSpectralAxis[0]);
                //        Float64 logInitialStep = log(spectrumRebinedSpectralAxis[0])-log(tplRebinedSpectralAxis[0]);
                //        Log.LogInfo("ChisquareLog, FitAllz: logstep=%f", logstep);
                //        Log.LogInfo("ChisquareLog, FitAllz: step=%f A", exp(logstep));
                //        Log.LogInfo("ChisquareLog, FitAllz: logInitialStep=%f", logInitialStep);
                //        Log.LogInfo("ChisquareLog, FitAllz: InitialStep=%f A", exp(logInitialStep));
                //        std::vector<Float64> log1pz(dtm_vec.size(), 0.0);
                //        for( Int32 t=0;t<log1pz.size();t++)
                //        {
                //            log1pz[t] = t*logstep + logInitialStep;
                //        }

                // save chi2 data
                FILE* f_chi2 = fopen( "loglbda_chi2output_dbg.txt", "w+" );
                for( Int32 t=0;t<chi2.size();t++)
                {
                    fprintf( f_chi2, "%f\t%e\n", z_vect[t]/*(Float64)(exp(log1pz[t])-1)*/, chi2[t]);
                }
                fclose( f_chi2 );
            }
        }
    }


    //interpolating on the regular z grid
    Float64* zreversed_array = new Float64 [(int)z_vect.size()]();
    Float64* chi2reversed_array = new Float64 [(int)bestChi2.size()]();
    for( Int32 t=0;t<z_vect.size();t++)
    {
        zreversed_array[t] = z_vect[z_vect.size()-1-t];
        chi2reversed_array[t] = bestChi2[z_vect.size()-1-t];
    }

    Int32 k = 0;
    Int32 klow = 0;
    for (Int32 iz=0;iz<result->Redshifts.size();iz++)
    {
        /* //NGP
        k = gsl_interp_bsearch (zreversed_array, result->Redshifts[iz], klow, z_vect.size()-1);
        klow = k;

//        if(result->Redshifts[iz]<1.1)
//        {
//            Log.LogInfo("ChisquareLog, FitAllz: interpolating z result, zcalc=%f, kfound=%d, zfound=%f", result->Redshifts[iz], k, z_vect[z_vect.size()-1-k]);
//        }
        // closest value
        result->ChiSquare[iz] = bestChi2[z_vect.size()-1-k];
        //*/

        result->Overlap[iz] = 1.0;
        result->FitAmplitude[iz] = 1.0;
        result->FitDtM[iz] = 1.0;
        result->FitMtM[iz] = 1.0;
        result->FitDustCoeff[iz] = 1.0;
        result->FitMeiksinIdx[iz] = 1.0;
        result->Status[iz] = nStatus_OK;

        //*/
    }

    //*
    Log.LogInfo("ChisquareLog, FitAllz: interpolating z result (linear)");
    Int32 interpRet = InterpolateResult(
                chi2reversed_array,
                zreversed_array,
                &result->Redshifts.front(),
                z_vect.size(),
                result->Redshifts.size(),
                result->ChiSquare,
                DBL_MAX);
    //*/


    delete[] zreversed_array;
    delete[] chi2reversed_array;
    delete[] spcRebinedFluxOverErr2;
    delete[] oneSpcRebinedFluxOverErr2;
    delete[] tplRebinedFlux;
    delete[] tpl2RebinedFlux;

    return 0;
}

Int32 COperatorChiSquareLogLambda::InterpolateResult(const Float64* in, const Float64* inGrid, const Float64* tgtGrid, Int32 n, Int32 tgtn, std::vector<Float64> &out, Float64 defaultValue)
{
    out.resize(tgtn);
    for(Int32 j=0; j<tgtn; j++)
    {
        out[j] = defaultValue;
    }

    //* // GSL method spline
    //initialise and allocate the gsl objects
    //lin
    gsl_interp *interpolation = gsl_interp_alloc (gsl_interp_linear,n);
    gsl_interp_init(interpolation, inGrid, in, n);
    gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();
    //spline
    //gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n);
    //gsl_spline_init (spline, inGrid, in, n);
    //gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();

    for( Int32 j=0; j<tgtn; j++ )
    {
        Float64 Xrebin = tgtGrid[j];
        out[j] = gsl_interp_eval(interpolation, inGrid, in, Xrebin, accelerator); //lin
        //out[j] = gsl_spline_eval (spline, Xrebin, accelerator); //spline
        //Log.LogInfo("ChisquareLog, FitAllz: interpolating gsl-spline z result, , ztgt=%f, rebinY=%f", tgtGrid[j], out[j]);
    }

    gsl_interp_free(interpolation);
    //gsl_spline_free (spline);
    gsl_interp_accel_free (accelerator);

    return 0;
}


/**
 * \brief COperatorChiSquareLogLambda::Compute
 *
 * This method computes the log_likelihood for the input spc and the tpl on a given redshift range (should be a regular grid):
 * 0. checks :
 *      - is the redshift list a regular grid ? (assumed in the rest of the method)
 *      - is the overlap always >100% in the given redshift range ?
 * 1. resample the input spectrum/tpl on a loglambda regular grid (always resampling as of 2017-06-13, option todo: use the already log-sampled input spectrum grid)
 * 2.
 * input: if additional_spcMasks size is 0, no additional mask will be used, otherwise its size should match the redshifts list size
 *
 **/
std::shared_ptr<COperatorResult> COperatorChiSquareLogLambda::Compute(const CSpectrum& spectrum, const CTemplate& tpl,
                          const TFloat64Range& lambdaRange, const TFloat64List& redshifts,
                          Float64 overlapThreshold , std::vector<CMask> additional_spcMasks, std::string opt_interp, Int32 opt_extinction, Int32 opt_dustFitting)
{

    if( opt_dustFitting && m_ismCorrectionCalzetti->calzettiInitFailed)
    {
        Log.LogError("ChisquareLog, no calzetti calib. file loaded... aborting!");
        return NULL;
    }

    if( spectrum.GetSpectralAxis().IsInLinearScale() == false || tpl.GetSpectralAxis().IsInLinearScale() == false )
    {
        Log.LogError("ChisquareLog, input spectrum or template are not in log scale (ignored)");
        //return NULL;
    }

    //assuming that lambdarange is strictly included in the spectrum spectral range
    if(std::max(spectrum.GetSpectralAxis()[0], lambdaRange.GetBegin()) != lambdaRange.GetBegin())
    {
        Log.LogError("ChisquareLog, lambdarange is not strictly included in the spectrum spectral range (min lambdarange)");
    }
    if(std::min(spectrum.GetSpectralAxis()[spectrum.GetSpectralAxis().GetSamplesCount()-1], lambdaRange.GetEnd())!=lambdaRange.GetEnd())
    {
        Log.LogError("ChisquareLog, lambdarange is not strictly included in the spectrum spectral range (max lambdarange)");
    }

    //sort the redshift and keep track of the indexes
    TFloat64List sortedRedshifts;
    TFloat64List sortedIndexes; //used for the correspondence between input redshifts (list) and input additionalMasks (list)
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

    if(overlapThreshold<1.0)
    {
        Log.LogError("ChisquareLog, overlap threshold can't be lower than 1.0");
        return NULL;
    }
    //check that the overlap is >1. for all sortedRedshifts
    Bool overlapFull=true;
    if(lambdaRange.GetBegin()<tpl.GetSpectralAxis()[0]*(1+sortedRedshifts[sortedRedshifts.size()-1]))
    {
        overlapFull=false;
    }
    if(lambdaRange.GetEnd()>tpl.GetSpectralAxis()[tpl.GetSpectralAxis().GetSamplesCount()-1]*(1+sortedRedshifts[0]))
    {
        overlapFull=false;
    }
    if(!overlapFull)
    {
        Log.LogError("ChisquareLog, overlap found to be lower than 1.0 for this redshift range");
        Log.LogError("ChisquareLog, for zmin=%f, tpl.lbdamax is %f (should be >%f)",
                     sortedRedshifts[0],
                     tpl.GetSpectralAxis()[tpl.GetSpectralAxis().GetSamplesCount()-1],
                     lambdaRange.GetEnd()/(1+sortedRedshifts[0]));
        Log.LogError("ChisquareLog, for zmax=%f, tpl.lbdamin is %f (should be <%f)",
                     sortedRedshifts[sortedRedshifts.size()-1],
                     tpl.GetSpectralAxis()[0],
                     lambdaRange.GetBegin()/(1+sortedRedshifts[sortedRedshifts.size()-1]));
        return NULL;
    }


    std::string rebinMethod = "spline";
    Int32 enableLogRebin = 1;
    Bool verboseLogRebin = 1;
    Bool verboseExportLogRebin = 1;
    if(enableLogRebin==1){

        // Create the Spectrum Log-Rebined spectral axis
        //Float64 zRef = sortedRedshifts[sortedRedshifts.size()-1];           //max redshift
        //Float64 zRef = sortedRedshifts[int(sortedRedshifts.size()/2.0)];  //middle redshift value
        //Float64 zRef = sortedRedshifts[0];                          //min redshift value
        Float64 zStep = sortedRedshifts[1]-sortedRedshifts[0];
        //Float64 loglbdaStep = log(1.0+zStep/(1.0+zRef));
        Int32 lbdaMinIdx = spectrum.GetSpectralAxis().GetIndexAtWaveLength(lambdaRange.GetBegin());
        Int32 lbdaMaxIdx = spectrum.GetSpectralAxis().GetIndexAtWaveLength(lambdaRange.GetEnd());
        if(lbdaMinIdx<0 || lbdaMaxIdx<0 || lbdaMaxIdx<=lbdaMinIdx)
        {
            Log.LogError("ChisquareLog, problem while searching for lambdarange observed spectral indexes");
        }
        Float64 loglbdaStep = log(spectrum.GetSpectralAxis()[lbdaMaxIdx])-log(spectrum.GetSpectralAxis()[lbdaMaxIdx-1]);
        //Float64 loglbdaStep = log(spectrum.GetSpectralAxis()[1])-log(spectrum.GetSpectralAxis()[0]);
        Float64 loglbdamin = log(lambdaRange.GetBegin());
        Float64 loglbdamax = log(lambdaRange.GetEnd());
        Int32 loglbdaCount = int((loglbdamax-loglbdamin)/loglbdaStep+1);
        if(verboseLogRebin){
            Log.LogInfo("ChisquareLog, ORIGINAL grid spectrum count = %d", lbdaMaxIdx-lbdaMinIdx+1);
            //Log.LogInfo("ChisquareLog, Log-Rebin: zref = %f", zRef);
            //Log.LogInfo("ChisquareLog, Log-Rebin: zStep = %f", zStep);
            Log.LogInfo("ChisquareLog, Log-Rebin: loglbdaStep = %f", loglbdaStep);
            Log.LogInfo("ChisquareLog, Log-Rebin: loglbdamin=%f : loglbdamax=%f", loglbdamin, loglbdamax);
            Log.LogInfo("ChisquareLog, Log-Rebin: lbdamin=%f : lbdamax=%f", exp(loglbdamin), exp(loglbdamax));
            Log.LogInfo("ChisquareLog, Log-Rebin: loglbdaCount = %d", loglbdaCount);

        }
        CSpectrumSpectralAxis targetSpectralAxis;
        targetSpectralAxis.SetSize(loglbdaCount);

        for(Int32 k=0; k<loglbdaCount; k++)
        {
            targetSpectralAxis[k] = exp(loglbdamin+k*loglbdaStep);
        }
        //precision problem due to exp/log - beginning
        Float64 delta = targetSpectralAxis[0]-spectrum.GetSpectralAxis()[0];
        Float64 step = targetSpectralAxis[1]-targetSpectralAxis[0];
        if( delta < 0.0 && abs(delta)<step*1e-5)
        {
            targetSpectralAxis[0] = spectrum.GetSpectralAxis()[0];
        }
        //precision problem due to exp/log - end
        delta = targetSpectralAxis[loglbdaCount-1]-spectrum.GetSpectralAxis()[spectrum.GetSpectralAxis().GetSamplesCount()-1];
        step = targetSpectralAxis[loglbdaCount-1]-targetSpectralAxis[loglbdaCount-2];
        if( delta > 0.0 && abs(delta)<step*1e-5)
        {
            targetSpectralAxis[loglbdaCount-1] = spectrum.GetSpectralAxis()[spectrum.GetSpectralAxis().GetSamplesCount()-1];
        }
        if(verboseLogRebin){
            //zstep lower lbda range for this sampling
            Float64 dlambda_begin = (targetSpectralAxis[1]-targetSpectralAxis[0]);
            Float64 dz_zrangemin_begin = dlambda_begin/(targetSpectralAxis[0]/(1.+sortedRedshifts[0]));
            Float64 dz_zrangemax_begin = dlambda_begin/(targetSpectralAxis[0]/(1.+sortedRedshifts[sortedRedshifts.size()-1]));
            Log.LogInfo("ChisquareLog, Log-Rebin: dz_zrangemin_begin=%f : dz_zrangemax_begin=%f", dz_zrangemin_begin, dz_zrangemax_begin);

            Float64 dlambda_end = (targetSpectralAxis[loglbdaCount-1]-targetSpectralAxis[loglbdaCount-2]);
            Float64 dz_zrangemin_end = dlambda_end/(targetSpectralAxis[loglbdaCount-2]/(1.+sortedRedshifts[0]));
            Float64 dz_zrangemax_end = dlambda_end/(targetSpectralAxis[loglbdaCount-2]/(1.+sortedRedshifts[sortedRedshifts.size()-1]));
            Log.LogInfo("ChisquareLog, Log-Rebin: dz_zrangemin_end=%f : dz_zrangemax_end=%f", dz_zrangemin_end, dz_zrangemax_end);
        }

        if(verboseExportLogRebin)
        {
            // save rebinned data
            FILE* f_targetSpcAxis = fopen( "loglbda_rebinlog_targetSpcAxis_dbg.txt", "w+" );
            for( Int32 t=0;t<targetSpectralAxis.GetSamplesCount();t++)
            {
                fprintf( f_targetSpcAxis, "%f\n", targetSpectralAxis[t]);
            }
            fclose( f_targetSpcAxis );
        }

        // Allocate the Log-rebined spectrum and mask
        CSpectrumFluxAxis& spectrumRebinedFluxAxis = m_spectrumRebinedLog.GetFluxAxis();
        CSpectrumSpectralAxis& spectrumRebinedSpectralAxis = m_spectrumRebinedLog.GetSpectralAxis();
        spectrumRebinedFluxAxis.SetSize(loglbdaCount);
        spectrumRebinedSpectralAxis.SetSize(loglbdaCount);
        m_mskRebinedLog.SetSize(loglbdaCount);

        //rebin the spectrum
        TFloat64Range spcLbdaRange( exp(loglbdamin-0.5*loglbdaStep), exp(loglbdamax+0.5*loglbdaStep) );
        Float64* pfgBuffer_unused = 0;
        Float64 redshift_unused = 0.0;
        CSpectrumFluxAxis::Rebin2( spcLbdaRange,
                                   spectrum.GetFluxAxis(),
                                   pfgBuffer_unused,
                                   redshift_unused,
                                   spectrum.GetSpectralAxis(),
                                   targetSpectralAxis,
                                   spectrumRebinedFluxAxis,
                                   spectrumRebinedSpectralAxis,
                                   m_mskRebinedLog,
                                   rebinMethod );

        if(verboseExportLogRebin)
        {
            // save rebinned data
            FILE* f = fopen( "loglbda_rebinlog_spclogrebin_dbg.txt", "w+" );
            for( Int32 t=0;t<m_spectrumRebinedLog.GetSampleCount();t++)
            {
                fprintf( f, "%f\t%e\n", m_spectrumRebinedLog.GetSpectralAxis()[t], m_spectrumRebinedLog.GetFluxAxis()[t]);
            }
            fclose( f );
        }

        // Create the Template Log-Rebined spectral axis
        // The template grid has to be aligned with the spectrum log-grid (will use the spc first element)
        //Float64 tpl_raw_loglbdamin = log(tpl.GetSpectralAxis()[0]); //full tpl lambda range
        //Float64 tpl_raw_loglbdamax = log(tpl.GetSpectralAxis()[tpl.GetSpectralAxis().GetSamplesCount()-1]); //full tpl lambda range
        Float64 tplloglbdaStep = log(tpl.GetSpectralAxis()[1])-log(tpl.GetSpectralAxis()[0]); //warning : considering constant dlambda for the input template
        Float64 roundingErrorsMargin = tplloglbdaStep/10.0; //
        Float64 tpl_raw_loglbdamin = log(lambdaRange.GetBegin()/(1.0+sortedRedshifts[sortedRedshifts.size()-1]))-roundingErrorsMargin; // lambdarange cropped to useful range given zmin
        Float64 tpl_raw_loglbdamax = log(lambdaRange.GetEnd()/(1.0+sortedRedshifts[0]))+roundingErrorsMargin; // lambdarange cropped to useful range given zmin

        Float64 tpl_tgt_loglbdamin = loglbdamin;
        if(tpl_raw_loglbdamin<=loglbdamin)
        {
            while( tpl_tgt_loglbdamin-loglbdaStep >= tpl_raw_loglbdamin )
            {
                tpl_tgt_loglbdamin-=loglbdaStep;
            }
        }else
        {
            while( tpl_tgt_loglbdamin+loglbdaStep <= tpl_raw_loglbdamin )
            {
                tpl_tgt_loglbdamin+=loglbdaStep;
            }
        }
        Float64 tpl_tgt_loglbdamax = loglbdamax;
        if(tpl_raw_loglbdamax<=loglbdamax)
        {
            while( tpl_tgt_loglbdamax-loglbdaStep >= tpl_raw_loglbdamax )
            {
                tpl_tgt_loglbdamax-=loglbdaStep;
            }
        }else
        {
            while( tpl_tgt_loglbdamax+loglbdaStep <= tpl_raw_loglbdamax )
            {
                tpl_tgt_loglbdamax+=loglbdaStep;
            }
        }
        Int32 tpl_loglbdaCount = int((tpl_tgt_loglbdamax-tpl_tgt_loglbdamin)/loglbdaStep+1);
        if(verboseLogRebin){
            Log.LogInfo("ChisquareLog, Log-Rebin: tpl loglbdamin=%f : loglbdamax=%f", tpl_tgt_loglbdamin, tpl_tgt_loglbdamax);
            Log.LogInfo("ChisquareLog, Log-Rebin: tpl lbdamin=%f : lbdamax=%f", exp(tpl_tgt_loglbdamin), exp(tpl_tgt_loglbdamax));
            Log.LogInfo("ChisquareLog, Log-Rebin: tpl loglbdaCount = %d", tpl_loglbdaCount);

        }

        //rebin the template
        CSpectrumSpectralAxis tpl_targetSpectralAxis;
        tpl_targetSpectralAxis.SetSize(tpl_loglbdaCount);
        for(Int32 k=0; k<tpl_loglbdaCount; k++)
        {
            tpl_targetSpectralAxis[k] = exp(tpl_tgt_loglbdamin+k*loglbdaStep);
        }
        //precision problem due to exp/log - beginning
        delta = tpl_targetSpectralAxis[0]-tpl.GetSpectralAxis()[0];
        step = tpl_targetSpectralAxis[1]-tpl_targetSpectralAxis[0];
        if( delta < 0.0 && abs(delta)<step*1e-5)
        {
            tpl_targetSpectralAxis[0] = tpl.GetSpectralAxis()[0];
        }
        //precision problem due to exp/log - end
        delta = tpl_targetSpectralAxis[loglbdaCount-1]-tpl.GetSpectralAxis()[tpl.GetSpectralAxis().GetSamplesCount()-1];
        step = tpl_targetSpectralAxis[loglbdaCount-1]-tpl_targetSpectralAxis[loglbdaCount-2];
        if( delta > 0.0 && abs(delta)<step*1e-5)
        {
            tpl_targetSpectralAxis[loglbdaCount-1] = tpl.GetSpectralAxis()[tpl.GetSpectralAxis().GetSamplesCount()-1];
        }

        if(verboseExportLogRebin)
        {
            // save rebinned data
            FILE* f_targetTplAxis = fopen( "loglbda_rebinlog_targetTplAxis_dbg.txt", "w+" );
            for( Int32 t=0;t<tpl_targetSpectralAxis.GetSamplesCount();t++)
            {
                fprintf( f_targetTplAxis, "%f\n", tpl_targetSpectralAxis[t]);
            }
            fclose( f_targetTplAxis );
        }

        TFloat64Range tplLbdaRange( exp(tpl_tgt_loglbdamin-0.5*loglbdaStep), exp(tpl_tgt_loglbdamax+0.5*loglbdaStep) );
        //precision problem due to exp/log - beginning
        delta = tpl_targetSpectralAxis[0]-tpl.GetSpectralAxis()[0];
        if( delta < 0.0 )
        {
            Log.LogError("ChisquareLog, Log-Rebin: tpl rebin error. Target MIN lbda value=%f, input tpl min lbda value=%f",
                         tpl_targetSpectralAxis[0], tpl.GetSpectralAxis()[0]);
            Log.LogError("ChisquareLog, Log-Rebin: tpl rebin error. Extend your input template wavelength range or modify the processing parameter <lambdarange>");
        }
        //precision problem due to exp/log - end
        delta = tpl_targetSpectralAxis[tpl_targetSpectralAxis.GetSamplesCount()-1]-tpl.GetSpectralAxis()[tpl.GetSpectralAxis().GetSamplesCount()-1];
        if( delta > 0.0 )
        {
            Log.LogError("ChisquareLog, Log-Rebin: tpl rebin error. Target MAX lbda value=%f, input tpl max lbda value=%f",
                         tpl_targetSpectralAxis[tpl_targetSpectralAxis.GetSamplesCount()-1], tpl.GetSpectralAxis()[tpl.GetSpectralAxis().GetSamplesCount()-1]);
            Log.LogError("ChisquareLog, Log-Rebin: tpl rebin error. Extend your input template wavelength range or modify the processing parameter <lambdarange>");
        }

        //Float64* pfgBuffer_unused = 0;
        //Float64 redshift_unused = 0.0;
        CSpectrumFluxAxis& templateRebinedFluxAxis = m_templateRebinedLog.GetFluxAxis();
        CSpectrumSpectralAxis& templateRebinedSpectralAxis = m_templateRebinedLog.GetSpectralAxis();
        templateRebinedFluxAxis.SetSize(tpl_loglbdaCount);
        templateRebinedSpectralAxis.SetSize(tpl_loglbdaCount);
        CMask   tpl_mskRebinedLog;
        tpl_mskRebinedLog.SetSize(tpl_loglbdaCount);
        CSpectrumFluxAxis::Rebin2( tplLbdaRange,
                                   tpl.GetFluxAxis(),
                                   pfgBuffer_unused,
                                   redshift_unused,
                                   tpl.GetSpectralAxis(),
                                   tpl_targetSpectralAxis,
                                   templateRebinedFluxAxis,
                                   templateRebinedSpectralAxis,
                                   tpl_mskRebinedLog,
                                   rebinMethod );
        if(verboseExportLogRebin)
        {
            // save rebinned data
            FILE* f_tpl_tgtlbda = fopen( "loglbda_rebinlog_tpllogrebin_dbg.txt", "w+" );
            for( Int32 t=0;t<m_templateRebinedLog.GetSampleCount();t++)
            {
                fprintf( f_tpl_tgtlbda, "%f\t%e\n", m_templateRebinedLog.GetSpectralAxis()[t], m_templateRebinedLog.GetFluxAxis()[t]);
            }
            fclose( f_tpl_tgtlbda );
        }
    }


    std::shared_ptr<CChisquareResult> result = std::shared_ptr<CChisquareResult>( new CChisquareResult() );
    result->Init(sortedRedshifts.size());
    result->Redshifts = sortedRedshifts;

    //WARNING: no additional masks coded for use as of 2017-06-13
    if(additional_spcMasks.size()!=0)
    {
        Log.LogError("ChisquareLog, No additional masks used. Feature not coded for this log-lambda operator!)");
    }


    //**************** Fitting at all redshifts ****************//
    //Optionally apply some IGM absorption
    std::vector<Int32> igmMeiksinCoeffs;
    if(opt_extinction)
    {
        Log.LogError("ChisquareLog, FitAllz opt_extinction not implemented yet...");
        Int32 nIGMCoeffs = m_igmCorrectionMeiksin->GetIdxCount();
        igmMeiksinCoeffs.resize(nIGMCoeffs);
        for(Int32 kigm=0; kigm<nIGMCoeffs; kigm++)
        {
            igmMeiksinCoeffs[kigm] = kigm;
        }
    }else{
        igmMeiksinCoeffs.clear();
    }

    //Optionally apply some ISM attenuation
    std::vector<Int32> ismEbmvCoeffs;
    if(opt_dustFitting)
    {
        Log.LogError("ChisquareLog, FitAllz opt_dustFitting not validated yet...");
        Int32 nISMCoeffs = m_ismCorrectionCalzetti->GetNPrecomputedDustCoeffs();
        ismEbmvCoeffs.resize(nISMCoeffs);
        for(Int32 kism=0; kism<nISMCoeffs; kism++)
        {
            ismEbmvCoeffs[kism] = kism;
        }
    }else{
        ismEbmvCoeffs.clear();
    }


    Int32 retFit = FitAllz(lambdaRange, result, igmMeiksinCoeffs, ismEbmvCoeffs);
    if(retFit!=0)
    {
        Log.LogError("ChisquareLog, FitAllz failed with error %d", retFit);
    }
    //**************** End Fitting at all redshifts ****************//


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
        Log.LogInfo("ChisquareLog, overlap warning for %s: minz=%.3f, maxz=%.3f", tpl.GetName().c_str(), overlapValidInfZ, overlapValidSupZ);
    }

    //estimate CstLog for PDF estimation
    result->CstLog = 0.0;//Todo: check how to estimate that value for loglambda// EstimateLikelihoodCstLog(spectrum, lambdaRange);

    // extrema
    Int32 extremumCount = 10;
    if(result->Redshifts.size()>extremumCount)
    {
        TPointList extremumList;
        TFloat64Range redshiftsRange(result->Redshifts[0], result->Redshifts[result->Redshifts.size()-1]);
        CExtremum extremum( redshiftsRange, extremumCount, true);
        extremum.Find( result->Redshifts, result->ChiSquare, extremumList );

        //*
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
        //*/
        // store extrema results
        result->Extrema.resize( extremumCount );
        for( Int32 i=0; i<extremumList.size(); i++ )
        {

            result->Extrema[i] = extremumList[i].X;
        }

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
    }

    return result;

}

const Float64*  COperatorChiSquareLogLambda::getDustCoeff(Float64 dustCoeff, Float64 maxLambda)
{
    return m_ismCorrectionCalzetti->getDustCoeff(dustCoeff, maxLambda);
}

const Float64*  COperatorChiSquareLogLambda::getMeiksinCoeff(Int32 meiksinIdx, Float64 redshift, Float64 maxLambda)
{
    return m_igmCorrectionMeiksin->getMeiksinCoeff(meiksinIdx, redshift, maxLambda);
}
