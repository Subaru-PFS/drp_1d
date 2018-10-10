#include <RedshiftLibrary/operator/chisquare2.h>

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

#include <assert.h>

#define NOT_OVERLAP_VALUE NAN
#include <stdio.h>

namespace bfs = boost::filesystem;
using namespace NSEpic;
using namespace std;

COperatorChiSquare2::COperatorChiSquare2( std::string calibrationPath )
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

COperatorChiSquare2::~COperatorChiSquare2()
{
    delete[] m_YtplRawBuffer;
    delete m_ismCorrectionCalzetti;
    delete m_igmCorrectionMeiksin;
}


void COperatorChiSquare2::BasicFit(const CSpectrum& spectrum, const CTemplate& tpl, Float64* pfgTplBuffer,
                                const TFloat64Range& lambdaRange, Float64 redshift, Float64 overlapThreshold,
                                   Float64& overlapRate,
                                   Float64& chiSquare,
                                   Float64& fittingAmplitude,
                                   Float64& fittingDtM,
                                   Float64& fittingMtM,
                                   Float64 &fittingDustCoeff,
                                   Float64 &fittingMeiksinIdx,
                                   EStatus& status,
                                   std::vector<TFloat64List>& ChiSquareInterm,
                                   std::string opt_interp, Float64 forcedAmplitude, Int32 opt_extinction, Int32 opt_dustFitting, CMask spcMaskAdditional)
{
    Bool amplForcePositive=true;
    chiSquare = boost::numeric::bounds<float>::highest();
    bool status_chisquareSetAtLeastOnce = false;

    for(Int32 kism=0; kism<ChiSquareInterm.size(); kism++)
    {
        for(Int32 kigm=0; kigm<ChiSquareInterm[kism].size(); kigm++)
        {
            ChiSquareInterm[kism][kigm] = boost::numeric::bounds<float>::highest();
        }
    }
    fittingAmplitude = -1.0;
    overlapRate = 0.0;
    status = nStatus_DataError;

    //CSpectrumSpectralAxis shiftedTplSpectralAxis( tpl.GetSampleCount(), false );

    const CSpectrumSpectralAxis& spcSpectralAxis = spectrum.GetSpectralAxis();
    const CSpectrumFluxAxis& spcFluxAxis = spectrum.GetFluxAxis();

    const CSpectrumSpectralAxis& tplSpectralAxis = tpl.GetSpectralAxis();
    const CSpectrumFluxAxis& tplFluxAxis = tpl.GetFluxAxis();


    if(spcMaskAdditional.GetMasksCount()!=spcFluxAxis.GetSamplesCount())
    {
        Log.LogInfo("  Operator-Chisquare2: spcMaskAdditional does not have the same size as the spectrum flux vector... (%d vs %d), aborting!", spcMaskAdditional.GetMasksCount(), spcFluxAxis.GetSamplesCount());
        status = nStatus_DataError;
        return ;
    }

    // Compute clamped lambda range over spectrum
    TFloat64Range spcLambdaRange;
    spcSpectralAxis.ClampLambdaRange( lambdaRange, spcLambdaRange );

    // Compute shifted template
    Float64 onePlusRedshift = 1.0 + redshift;
    m_shiftedTplSpectralAxis_bf.ShiftByWaveLength( tplSpectralAxis, onePlusRedshift, CSpectrumSpectralAxis::nShiftForward );
    TFloat64Range intersectedLambdaRange( 0.0, 0.0 );

    // Compute clamped lambda range over template
    TFloat64Range tplLambdaRange;
    m_shiftedTplSpectralAxis_bf.ClampLambdaRange( lambdaRange, tplLambdaRange );

    // if there is any intersection between the lambda range of the spectrum and the lambda range of the template
    // Compute the intersected range
    TFloat64Range::Intersect( tplLambdaRange, spcLambdaRange, intersectedLambdaRange );

    //UInt32 tgtn = spcSpectralAxis.GetSamplesCount() ;
    CSpectrumFluxAxis& itplTplFluxAxis = m_templateRebined_bf.GetFluxAxis();
    CSpectrumSpectralAxis& itplTplSpectralAxis = m_templateRebined_bf.GetSpectralAxis();
    CMask& itplMask = m_mskRebined_bf;

    //CSpectrumFluxAxis::Rebin( intersectedLambdaRange, tplFluxAxis, shiftedTplSpectralAxis, spcSpectralAxis, itplTplFluxAxis, itplTplSpectralAxis, itplMask );
    CSpectrumFluxAxis::Rebin2( intersectedLambdaRange, tplFluxAxis, pfgTplBuffer, redshift, m_shiftedTplSpectralAxis_bf, spcSpectralAxis, itplTplFluxAxis, itplTplSpectralAxis, itplMask, opt_interp );


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
    overlapRate = spcSpectralAxis.IntersectMaskAndComputeOverlapRate( lambdaRange, itplMask );

    // Check for overlap rate
    if( overlapRate < overlapThreshold || overlapRate<=0.0 )
    {
        status = nStatus_NoOverlap;
        return ;
    }


    const Float64* Xtpl = itplTplSpectralAxis.GetSamples();
    Float64* Ytpl = itplTplFluxAxis.GetSamples();
    const Float64* Xspc = spcSpectralAxis.GetSamples();
    const Float64* Yspc = spcFluxAxis.GetSamples();
    TFloat64Range logIntersectedLambdaRange( log( intersectedLambdaRange.GetBegin() ), log( intersectedLambdaRange.GetEnd() ) );
    //the spectral axis should be in the same scale
    TFloat64Range currentRange = logIntersectedLambdaRange;
    if( spcSpectralAxis.IsInLinearScale() != tplSpectralAxis.IsInLinearScale() )
    {
        Log.LogError( "    chisquare operator: data and model not in the same scale (lin/log) ! Aborting.");
        status = nStatus_DataError;
        return;
    }
    if(spcSpectralAxis.IsInLinearScale()){
        currentRange = intersectedLambdaRange;
    }



    /*//debug:
    // save rebinedtpl
    if(redshift==0.0){
        FILE* f = fopen( "chisquare2_template_rebined.txt", "w+" );
        for(Int32 m=0; m<itplTplSpectralAxis.GetSamplesCount(); m++){
            fprintf( f, "%e %e\n", Xtpl[m], Ytpl[m]);
        }
        fclose( f );
    }
    //*/


    //save Tpl Flux without dust or any other weighting
    bool ytpl_modified = false;
    for(Int32 k=0; k<itplTplSpectralAxis.GetSamplesCount(); k++)
    {
        m_YtplRawBuffer[k] = Ytpl[k];
    }

    // Optionally Apply some Calzetti Extinction for DUST
    bool opt_dust_calzetti = opt_dustFitting;

    Int32 nDustCoeffs=1;
    if(!opt_dust_calzetti)
    {
        nDustCoeffs = 1;
    }else{
        nDustCoeffs = m_ismCorrectionCalzetti->GetNPrecomputedDustCoeffs();
        if(m_YtplRawBufferMaxBufferSize<itplTplSpectralAxis.GetSamplesCount())
        {
            Log.LogError( "chisquare operator: rebinned tpl size > buffer size for dust-fit ! Aborting.");
            status = nStatus_DataError;
            return ;
        }
    }

    //Optionally apply some IGM absorption
    Int32 nIGMCoeffs=1;
    if(opt_extinction)
    {
        nIGMCoeffs = m_igmCorrectionMeiksin->GetIdxCount();
    }

    Bool option_igmFastProcessing = true; //todo: find a way to unit-test this acceleration
    //Prepare the wavelengthRange Limits
    Float64 lbda_min = currentRange.GetBegin();
    Float64 lbda_max = currentRange.GetEnd();
    std::vector<Float64> sumCross_outsideIGM(nDustCoeffs, 0.0);
    std::vector<Float64>  sumT_outsideIGM(nDustCoeffs, 0.0);
    std::vector<Float64>  sumS_outsideIGM(nDustCoeffs, 0.0);


    //Loop on the meiksin Idx
    Bool igmLoopUseless_WavelengthRange = false;
    for(Int32 kMeiksin=0; kMeiksin<nIGMCoeffs; kMeiksin++)
    {
        if(igmLoopUseless_WavelengthRange)
        {
            //Now copy from the already calculated k>0 igm values
            for(Int32 kism=0; kism<ChiSquareInterm.size(); kism++)
            {
                for(Int32 kigm=1; kigm<ChiSquareInterm[kism].size(); kigm++)
                {
                    ChiSquareInterm[kism][kigm] = ChiSquareInterm[kism][0];
                }
            }
            break;
        }
        Int32 meiksinIdx = kMeiksin; //index for the Meiksin curve (0-6; 3 being the median extinction value)

        if(option_igmFastProcessing && kMeiksin>0)
        {
            lbda_max = m_igmCorrectionMeiksin->GetLambdaMax()*(1+redshift);
            if(lbda_max>currentRange.GetEnd())
            {
                lbda_max = currentRange.GetEnd();
            }

        }else{
            lbda_max = currentRange.GetEnd();
        }

        //find samples limits
        Int32 kStart = -1;
        Int32 kEnd = -1;
        for(Int32 k=0; k<spcSpectralAxis.GetSamplesCount(); k++)
        {
            if(Xspc[k] >= lbda_min && kStart==-1){
                kStart=k;
            }
            if(Xspc[k] <= lbda_max){
                kEnd=k;
            }

        }
        if(kStart==-1 || kEnd==-1)
        {
            Log.LogDebug( "chisquare operator: kStart=%d, kEnd=%d ! Aborting.", kStart, kEnd);
            break;
        }


        //Loop on the EBMV dust coeff
        for(Int32 kDust=0; kDust<nDustCoeffs; kDust++)
        {
            if(ytpl_modified)
            {
                //re-init flux tpl without dust and other weightin
                for(Int32 k=kStart; k<=kEnd; k++)
                {
                    Ytpl[k] = m_YtplRawBuffer[k];
                }
            }

            Float64 coeffEBMV = m_ismCorrectionCalzetti->GetEbmvValue(kDust);
            //Log.LogInfo("  Operator-Chisquare2: fitting with dust coeff value: %f", coeffEBMV);

            Float64 z = redshift;
            for(Int32 k=kStart; k<=kEnd; k++)
            {
                Float64 restLambda = Xtpl[k]/(1.0+z);
                Float64 coeffDust = m_ismCorrectionCalzetti->getDustCoeff( kDust, restLambda);

                Ytpl[k] *= coeffDust;
            }
            ytpl_modified = true;

            //Meiksin IGM extinction
            if(opt_extinction)
            {
                Int32 redshiftIdx = m_igmCorrectionMeiksin->GetRedshiftIndex(redshift); //index for IGM Meiksin redshift range

                Float64 coeffIGM = 1.0;
                Float64 z = redshift;
                Bool igmCorrectionAppliedOnce = false;
                for(Int32 k=kStart; k<=kEnd; k++)
                {

                    Float64 restLambda = Xtpl[k]/(1.0+z);
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
                        Ytpl[k] *= coeffIGM;
                        igmCorrectionAppliedOnce = true;
                    }
                }
                if(!igmCorrectionAppliedOnce)
                {
                    igmLoopUseless_WavelengthRange = true;
                }
            }
            //*/


            /*//debug:
            // save final ISM/IGM template model
            if(redshift>=2.4 && redshift<2.4001 && meiksinIdx==6){
                FILE* f = fopen( "chisquare2_template_finalModel.txt", "w+" );
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
                        Log.LogDebug("  Operator-Chisquare2: noise invalid=%e for i=%d", error[j], j);
                        Log.LogDebug("  Operator-Chisquare2: noise invalid, w=%e for i=%d", spcSpectralAxis[j], j);
                        Log.LogDebug("  Operator-Chisquare2: noise invalid,samples count=%d", spcSpectralAxis.GetSamplesCount(), j);
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
		       	Log.LogError("  Operator-Chisquare2: found invalid inverse variance : err2=%e, for index=%d at wl=%f", err2, j, spcSpectralAxis[j]);
			status = nStatus_InvalidProductsError;
            		return ; 
		    }


		    if( std::isinf(sumS) || std::isnan(sumS) || sumS!=sumS ){
		       	Log.LogError("  Operator-Chisquare2: found invalid dtd : dtd=%e, for index=%d at wl=%f", sumS, j, spcSpectralAxis[j]);
		       	Log.LogError("  Operator-Chisquare2: found invalid dtd : Yspc=%e, for index=%d at wl=%f", Yspc[j], j, spcSpectralAxis[j]);
		       	Log.LogError("  Operator-Chisquare2: found invalid dtd : err2=%e, for index=%d at wl=%f", err2, j, spcSpectralAxis[j]);
		       	Log.LogError("  Operator-Chisquare2: found invalid dtd : error=%e, for index=%d at wl=%f", error[j], j, spcSpectralAxis[j]);
			status = nStatus_InvalidProductsError;
            		return ; 
		    }
                    
		    if( std::isinf(sumT) || std::isnan(sumT) ){
		       	Log.LogError("  Operator-Chisquare2: found invalid mtm : mtm=%e, for index=%d at wl=%f", sumT, j, spcSpectralAxis[j]);
			status = nStatus_InvalidProductsError;
            		return ; 
		    }
                    

		    if(option_igmFastProcessing && kMeiksin==0)
                    {
                        //store intermediate sums for IGM range
                        if(sumsIgmSaved==0)
                        {
                            if(Xspc[j]>=m_igmCorrectionMeiksin->GetLambdaMax()*(1+redshift))
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
            if(option_igmFastProcessing && kMeiksin==0 && sumsIgmSaved==1)
            {
                sumCross_outsideIGM[kDust] = sumCross-sumCross_IGM;
                sumT_outsideIGM[kDust] = sumT-sumT_IGM;
                sumS_outsideIGM[kDust] = sumS-sumS_IGM;
            }
            if(option_igmFastProcessing && kMeiksin>0)
            {
                sumCross += sumCross_outsideIGM[kDust];
                sumT += sumT_outsideIGM[kDust];
                sumS += sumS_outsideIGM[kDust];
            }

            if ( numDevs==0 )
            {
                status = nStatus_DataError;
                return;
            }

            Float64 ampl = 0.0;
            if( sumT==0 )
            {
                ampl = 0.0;
                fit = sumS;
                //status = nStatus_DataError;
                //return;
            }else{
                if(amplForcePositive)
                {
                    ampl = max(0.0, sumCross / sumT);
                }else{
                    ampl = sumCross / sumT;
                }
                if(forcedAmplitude !=-1){
                    ampl = forcedAmplitude;
                }
                //* //1. fast method: D. Vibert, Amazed methods improvements, 10/06/2015
                fit = sumS - sumCross*ampl;
            }

            /*
            Log.LogDebug("  Operator-Chisquare2: fit=%e", fit);
            Log.LogDebug("  Operator-Chisquare2: sumT=%e", sumT);
            Log.LogDebug("  Operator-Chisquare2: sumS=%e", sumS);
            Log.LogDebug("  Operator-Chisquare2: sumCross=%e", sumCross);
            //*/

            /*
            if(redshift==0.0)
            {
                Log.LogInfo( "chisquare2 operator: for z=%f dtd = %e", redshift, sumS);
                Log.LogInfo( "chisquare2 operator: for z=%f dtm = %e", redshift, sumCross);
                Log.LogInfo( "chisquare2 operator: for z=%f mtm = %e", redshift, sumT);
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

            ChiSquareInterm[kDust][kMeiksin] = fit;

            if(fit<chiSquare)
            {
                chiSquare = fit;
                fittingDustCoeff = coeffEBMV;
                fittingMeiksinIdx = meiksinIdx;
                fittingAmplitude = ampl;
                fittingDtM = sumCross;
                fittingMtM = sumT;

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
 **/
std::shared_ptr<COperatorResult> COperatorChiSquare2::Compute(const CSpectrum& spectrum, const CTemplate& tpl,
                          const TFloat64Range& lambdaRange, const TFloat64List& redshifts,
                          Float64 overlapThreshold , std::vector<CMask> additional_spcMasks, std::string opt_interp, Int32 opt_extinction, Int32 opt_dustFitting)
{
    Log.LogInfo("  Operator-Chisquare2: starting computation for template: %s", tpl.GetName().c_str());


    if( opt_dustFitting && m_ismCorrectionCalzetti->calzettiInitFailed)
    {
        Log.LogError("  Operator-Chisquare2: no calzetti calib. file loaded... aborting!");
        return NULL;
    }

    if( opt_extinction && m_igmCorrectionMeiksin->meiksinInitFailed)
    {
        Log.LogError("  Operator-Chisquare2: no meiksin calib. file loaded... aborting!");
        return NULL;
    }

    if( spectrum.GetSpectralAxis().IsInLinearScale() == false || tpl.GetSpectralAxis().IsInLinearScale() == false )
    {
        Log.LogError("  Operator-Chisquare2: input spectrum or template are not in log scale (ignored)");
        //return NULL;
    }

    //**************************************************
//    Bool enableInputSpcCheck = false;
//    if(enableInputSpcCheck)
//    {
//        //Check if the Spectrum is valid on the lambdarange
//        const Float64 lmin = lambdaRange.GetBegin();
//        const Float64 lmax = lambdaRange.GetEnd();
//        if( !spectrum.IsFluxValid( lmin, lmax ) ){
//            Log.LogError( "    Operator-Chisquare2: debug - Failed to validate spectrum flux: %s, on wavelength range (%.1f ; %.1f)", spectrum.GetName().c_str(), lmin, lmax );
//            return NULL;
//        }else{
//            Log.LogDetail( "    Operator-Chisquare2: debug - Successfully validated spectrum flux: %s, on wavelength range (%.1f ; %.1f)", spectrum.GetName().c_str(), lmin, lmax );
//        }
//        if( !spectrum.IsNoiseValid( lmin, lmax ) ){
//            Log.LogError( "    Operator-Chisquare2: debug - Failed to validate noise from spectrum: %s, on wavelength range (%.1f ; %.1f)", spectrum.GetName().c_str(), lmin, lmax );
//            return NULL;
//        }else{
//            Log.LogDetail( "    Operator-Chisquare2: debug - Successfully validated noise from spectrum: %s, on wavelength range (%.1f ; %.1f)", spectrum.GetName().c_str(), lmin, lmax );
//        }
//    }


    // Pre-Allocate the rebined template and mask with regard to the spectrum size
    m_templateRebined_bf.GetSpectralAxis().SetSize(spectrum.GetSampleCount());
    m_templateRebined_bf.GetFluxAxis().SetSize(spectrum.GetSampleCount());
    m_mskRebined_bf.SetSize(spectrum.GetSampleCount());
    m_shiftedTplSpectralAxis_bf.SetSize( tpl.GetSampleCount());

    Float64* precomputedFineGridTplFlux;
    if(opt_interp=="precomputedfinegrid"){
        //*/
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

    std::shared_ptr<CChisquareResult> result = std::shared_ptr<CChisquareResult>( new CChisquareResult() );
    Int32 nDustCoeffs=1;
    if(opt_dustFitting)
    {
        nDustCoeffs = m_ismCorrectionCalzetti->GetNPrecomputedDustCoeffs();
    }
    Int32 nIGMCoeffs=1;
    if(opt_extinction)
    {
        nIGMCoeffs = m_igmCorrectionMeiksin->GetIdxCount();
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
        Log.LogError("  Operator-Chisquare2: using default mask, masks-list size (%d) didn't match the input redshift-list (%d) !)", additional_spcMasks.size(), sortedRedshifts.size());
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

        Float64 redshift = result->Redshifts[i];

        BasicFit( spectrum,
                  tpl,
                  precomputedFineGridTplFlux,
                  lambdaRange,
                  redshift,
                  overlapThreshold,
                  result->Overlap[i],
                  result->ChiSquare[i],
                  result->FitAmplitude[i],
                  result->FitDtM[i],
                  result->FitMtM[i],
                  result->FitDustCoeff[i],
                  result->FitMeiksinIdx[i],
                  result->Status[i],
                  result->ChiSquareIntermediate[i],
                  opt_interp,
                  -1,
                  opt_extinction,
                  opt_dustFitting,
                  additional_spcMask);

	if(result->Status[i]==nStatus_InvalidProductsError)
        {
            Log.LogError("  Operator-Chisquare2: found invalid chisquare products for z=%f. Now breaking z loop.", redshift);
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
        Log.LogInfo("  Operator-Chisquare2: overlap warning for %s: minz=%.3f, maxz=%.3f", tpl.GetName().c_str(), overlapValidInfZ, overlapValidSupZ);
    }

    //only bad status warning
    Int32 oneValidStatusFoundIndex = -1;
    for (Int32 i=0;i<sortedRedshifts.size();i++)
    {
        if(result->Status[i]==nStatus_OK)
        {
            oneValidStatusFoundIndex=i;
            Log.LogDebug("  Operator-Chisquare2: STATUS VALID found for %s: at least at index=%d", tpl.GetName().c_str(), i);
            break;
        }
    }if(oneValidStatusFoundIndex==-1)
    {
        Log.LogWarning("  Operator-Chisquare2: STATUS WARNING for %s: Not even one single valid fit/merit value found", tpl.GetName().c_str());
    }


    //loop error status warning
    Int32 loopErrorStatusFoundIndex = -1;
    for (Int32 i=0;i<sortedRedshifts.size();i++)
    {
        if(result->Status[i]==nStatus_LoopError)
        {
            loopErrorStatusFoundIndex=i;
            Log.LogDebug("  Operator-Chisquare2: STATUS Loop Error found for %s: at least at index=%d", tpl.GetName().c_str(), i);
            break;
        }
    }if(loopErrorStatusFoundIndex!=-1)
    {
        Log.LogWarning("    chisquare2-operator: Loop Error - chisquare values not set even once");
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

        Log.LogDebug("  Operator-Chisquare2: EXTREMA found n=%d", extremumList.size());
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
        Log.LogDebug("  Operator-Chisquare2: EXTREMA forced n=%d", result->Extrema.size());
    }

    if(opt_interp=="precomputedfinegrid"){
        free(precomputedFineGridTplFlux);
    }
    return result;

}


/* @brief COperatorChiSquare2::getDustCoeff: get the dust coeff at a fixed resolution of 1A
* @param dustCoeff
* @param maxLambda
* @return
*/
const Float64*  COperatorChiSquare2::getDustCoeff(Float64 dustCoeff, Float64 maxLambda)
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
 * @brief COperatorChiSquare2::getMeiksinCoeff: get the IGM Meiksin coeff at a fixed resolution of 1A
 * @param dustCoeff
 * @param maxLambda
 * @return
 */
const Float64*  COperatorChiSquare2::getMeiksinCoeff(Int32 meiksinIdx, Float64 redshift, Float64 maxLambda)
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


const COperatorResult* COperatorChiSquare2::ExportChi2versusAZ(const CSpectrum& spectrum, const CTemplate& tpl,
                          const TFloat64Range& lambdaRange, const TFloat64List& redshifts,
                          Float64 overlapThreshold )
{
    if( spectrum.GetSpectralAxis().IsInLinearScale() == false || tpl.GetSpectralAxis().IsInLinearScale() == false )
    {
        Log.LogError("  Operator-Chisquare2: input spectrum or template are not in log scale (ignored)");
        //return NULL;
    }

    // Pre-Allocate the rebined template and mask with regard to the spectrum size
    m_templateRebined_bf.GetSpectralAxis().SetSize(spectrum.GetSampleCount());
    m_templateRebined_bf.GetFluxAxis().SetSize(spectrum.GetSampleCount());
    m_mskRebined_bf.SetSize(spectrum.GetSampleCount());
    m_shiftedTplSpectralAxis_bf.SetSize( tpl.GetSampleCount());


    //*/
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
    Float64* precomputedFineGridTplFlux = (Float64*)malloc(nTgt*sizeof(Float64));
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
    //*/

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

    CChisquareResult* result = new CChisquareResult();
    result->ChiSquare.resize( sortedRedshifts.size() );
    result->FitAmplitude.resize( sortedRedshifts.size() );
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
            BasicFit( spectrum, tpl, precomputedFineGridTplFlux, lambdaRange, result->Redshifts[i], overlapThreshold,
                      result->Overlap[i],
                      result->ChiSquare[i],
                      result->FitAmplitude[i],
                      result->FitDtM[i],
                      result->FitMtM[i],
                      result->FitDustCoeff[i],
                      result->FitMeiksinIdx[i],
                      result->Status[i],
                      result->ChiSquareIntermediate[i], "lin", ampl );

            fprintf( f, "%.15e", result->ChiSquare[i]);
            if(j<sortedAmplitudes.size()-1){
                fprintf( f, "\t");
            }
        }
        fprintf( f, "\n");
    }
    fclose( f );



    free(precomputedFineGridTplFlux);
    return result;

}

/**
 * \brief this function estimates the likelihood_cstLog term withing the wavelength range
 **/
Float64 COperatorChiSquare2::EstimateLikelihoodCstLog(const CSpectrum& spectrum, const TFloat64Range& lambdaRange)
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
