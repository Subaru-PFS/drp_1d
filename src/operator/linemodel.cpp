#include <epic/redshift/operator/linemodel.h>

#include <epic/redshift/spectrum/axis.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/spectrum/tools.h>
#include <epic/redshift/common/mask.h>
#include <epic/redshift/operator/chisquareresult.h>
#include <epic/redshift/extremum/extremum.h>

#include <epic/redshift/spectrum/io/fitswriter.h>
#include <epic/core/log/log.h>

#include <boost/numeric/conversion/bounds.hpp>

#include <epic/redshift/gaussianfit/gaussianfitsimple.h>

#include <math.h>
#include <float.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <algorithm>    // std::sort

#include <assert.h>

#define NOT_OVERLAP_VALUE NAN
#include <stdio.h>

using namespace NSEpic;
using namespace std;


IMPLEMENT_MANAGED_OBJECT(COperatorLineModel)

COperatorLineModel::COperatorLineModel()
{

}

COperatorLineModel::~COperatorLineModel()
{

}


const COperatorResult* COperatorLineModel::Compute(const CSpectrum& spectrum, const CRayCatalog& restraycatalog,
                          const TFloat64Range& lambdaRange, const TFloat64List& redshifts)
{

    if( spectrum.GetSpectralAxis().IsInLinearScale() == false)
    {
        Log.LogError("Line Model, input spectrum is not in log scale (ignored)");
        //return NULL;
    }

    TFloat64List sortedRedshifts = redshifts;
    std::sort(sortedRedshifts.begin(), sortedRedshifts.end());

    Int32 typeFilter = CRay::nType_Emission;
    Int32 forceFilter = -1;//CRay::nForce_Strong;
    CRayCatalog::TRayVector restRayList = restraycatalog.GetFilteredList(typeFilter, forceFilter);

    CLineModelResult* result = new CLineModelResult();
    result->ChiSquare.resize( sortedRedshifts.size() );
    result->Redshifts.resize( sortedRedshifts.size() );
    result->Redshifts = sortedRedshifts;
    result->Status.resize( sortedRedshifts.size() );
    result->restRayList = restRayList;
    result->LineModelSolutions.resize( sortedRedshifts.size() );

    CLineModelElementList model(spectrum, restRayList);

    for (Int32 i=0;i<sortedRedshifts.size();i++)
    {
        ModelFit( spectrum, model, result->restRayList, lambdaRange, result->Redshifts[i], result->ChiSquare[i], result->LineModelSolutions[i]);
    }

    /* //saving the model for viewing
    CSpectrumFluxAxis& sfluxAxisPtr = model.GetFluxAxis();
    sfluxAxisPtr = modelFluxAxis;
    CSpectrumIOFitsWriter writer;
    Bool retVal = writer.Write( "model.fits", model );
    //*/

    // extrema
    Int32 extremumCount = 15;
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
    result->Area.resize( extremumCount );
    for( Int32 i=0; i<extremumList.size(); i++ )
    {
        result->Extrema[i] = extremumList[i].X;
        result->Area[i] = -1.0;
    }
    ComputeGaussAreaForExtrema(result);


    /*
    //  //saving the best model for viewing
    if(result->Extrema.size()>0){
        Float64 _chi=0.0;
        CLineModelResult::SLineModelSolution _lineModelSolution;
        ModelFit( spectrum, model, result->restRayList, lambdaRange, result->Extrema[0], _chi, _lineModelSolution);
        //*
        //CSpectrumFluxAxis& sfluxAxisPtr = model.GetFluxAxis();
        //CSpectrumFluxAxis& modelFluxAxis = model.GetFluxAxis();
        //sfluxAxisPtr = modelFluxAxis;
        CSpectrum spcmodel = model.GetModelSpectrum();

        CSpectrumIOFitsWriter writer;
        Bool retVal = writer.Write( "model.fits",  spcmodel);
    }
    //*/
    return result;

}

void COperatorLineModel::ComputeGaussAreaForExtrema(CLineModelResult* results)
{
    //prepare p
    Float64 maxp = DBL_MIN;
    CSpectrum pspc;
    CSpectrumFluxAxis& spcFluxAxis = pspc.GetFluxAxis();
    spcFluxAxis.SetSize( results->Redshifts.size() );
    CSpectrumSpectralAxis& spcSpectralAxis = pspc.GetSpectralAxis();
    spcSpectralAxis.SetSize( results->Redshifts.size()  );
    for( Int32 i2=0; i2<results->Redshifts.size(); i2++ )
    {
        if( maxp < results->ChiSquare[i2]){
            maxp = results->ChiSquare[i2];
        }
    }
    for( Int32 i2=0; i2<results->Redshifts.size(); i2++ )
    {
        spcFluxAxis[i2] = exp(-(results->ChiSquare[i2]-maxp)/2.0)-1.0;
        spcSpectralAxis[i2] = results->Redshifts[i2];
    }

    /*//debug:
    FILE* f = fopen( "getbestredshiftbayes_dbg.txt", "w+" );
    for( Int32 i=0; i<spcFluxAxis.GetSamplesCount(); i++ )
    {
        if( spcFluxAxis[i] < 0.0001 ){
            fprintf( f, "%e %e\n", spcSpectralAxis[i], spcFluxAxis[i]);
        }else{
            fprintf( f, "%f %f\n", spcSpectralAxis[i], spcFluxAxis[i]);
        }
    }
    fclose( f );
    //*/


    Float64 winsize = 0.0025;
    Float64 inclusionThresRatio = 0.25;
    Int32 iz0=0;
    for( Int32 i=0; i<results->Extrema.size(); i++ )
    {
        //find iz, izmin and izmax
        Int32 izmin= -1;
        Int32 iz= -1;
        Int32 izmax= -1;
        for( Int32 i2=0; i2<results->Redshifts.size(); i2++ )
        {
            if(iz == -1 && (results->Extrema[i]) <= results->Redshifts[i2]){
                iz = i2;
                if(i==0){
                    iz0=iz;
                }
            }
            if(izmin == -1 && (results->Extrema[i] - winsize/2.0) <= results->Redshifts[i2]){
                izmin = i2;
            }
            if(izmax == -1 && (results->Extrema[i] + winsize/2.0) <= results->Redshifts[i2]){
                izmax = i2;
                break;
            }
        }
        Float64 di = abs(results->ChiSquare[iz]-maxp) ;
        Float64 d0 = abs(results->ChiSquare[iz0]-maxp) ;
        if( di < inclusionThresRatio*d0){
            continue;
        }

        /*
        CGaussianFitSimple fitter;
        CGaussianFitSimple::EStatus status = fitter.Compute( pspc, TInt32Range( izmin, izmax ) );
        if(status!=NSEpic::CGaussianFitSimple::nStatus_Success){
            continue;
        }

        Float64 gaussAmp;
        Float64 gaussPos;
        Float64 gaussWidth;
        fitter.GetResults( gaussAmp, gaussPos, gaussWidth );
        Float64 gaussAmpErr;
        Float64 gaussPosErr;
        Float64 gaussWidthErr;
        fitter.GetResultsError( gaussAmpErr, gaussPosErr, gaussWidthErr );
        */
        Float64 gaussWidth = FitBayesWidth( spcSpectralAxis, spcFluxAxis, results->Extrema[i], izmin, izmax);
        Float64 gaussAmp = spcFluxAxis[iz];

        Float64 intsize = 0.001;
        Float64 area=0.0;
        for( Int32 i2=izmin; i2<izmax; i2++ )
        {
            Float64 x = spcSpectralAxis[i2];
            Float64 Yi = gaussAmp * exp (-1.*(x-results->Extrema[i])*(x-results->Extrema[i])/(2*gaussWidth*gaussWidth));
            area += Yi;
        }

        //Float64 area = gaussAmp*gaussWidth*sqrt(2.0*3.141592654);
        results->Area[i] = area;
    }
}

Void COperatorLineModel::ModelFit(const CSpectrum& spectrum, CLineModelElementList& model, const CRayCatalog::TRayVector& restRayList,
                                const TFloat64Range& lambdaRange, Float64 redshift,
                                   Float64& chiSquare, CLineModelResult::SLineModelSolution& modelSolution)
{
    chiSquare = boost::numeric::bounds<float>::highest();

    model.fit(redshift, modelSolution);

    const CSpectrumSpectralAxis& spcSpectralAxis = spectrum.GetSpectralAxis();
    const CSpectrumFluxAxis& spcFluxAxis = spectrum.GetFluxAxis();
    const CSpectrumFluxAxis& modelFluxAxis = model.GetModelSpectrum().GetFluxAxis();

    Int32 numDevs = 0;
    Float64 fit = 0;
    const Float64* error = spcFluxAxis.GetError();
    const Float64* Ymodel = modelFluxAxis.GetSamples();
    const Float64* Yspc = spcFluxAxis.GetSamples();
    for( UInt32 j=0; j<spcSpectralAxis.GetSamplesCount(); j++ )
    {
        numDevs++;
        // fit
        fit += pow( Yspc[j] - Ymodel[j] , 2.0 ) / pow( error[j], 2.0 );
        //fit += pow( Yspc[j] - Ymodel[j] , 2.0 );
    }
    fit /= numDevs;

    chiSquare = fit;
    return;
}


Float64 COperatorLineModel::FitAmplitude( const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis, Float64 lambda, Float64 width, Int32 start, Int32 end)
{
    const Float64* flux = fluxAxis.GetSamples();
    const Float64* spectral = spectralAxis.GetSamples();
    const Float64* error = fluxAxis.GetError();
    Float64 mu = lambda;
    Float64 c = width;

    Float64 y = 0.0;
    Float64 x = 0.0;
    Float64 yg = 0.0;

    Float64 sumCross = 0.0;
    Float64 sumGauss = 0.0;
    Float64 err2 = 0.0;
    Int32 num = 0;

    //A estimation
    for ( Int32 i = start; i < end; i++)
    {
        y = flux[i];
        x = spectral[i];
        yg = exp (-1.*(x-mu)*(x-mu)/(2*c*c));

        num++;
        err2 = 1.0 / (error[i] * error[i]);
        sumCross += yg*y*err2;
        sumGauss += yg*yg*err2;
    }

    if ( num==0 || sumCross==0 || sumGauss==0 )
    {
        return 0.0;
    }

    Float64 A = max(0.0, sumCross / sumGauss);

    /*
    //SNR estimation
    Float64 sumErr  = 0.0;
    Float64 sumGaussA = 0.0;
    for ( Int32 i = start; i < end; i++)
    {
        x = spectral[i];
        sumGaussA += A*exp (-1.*(x-mu)*(x-mu)/(2*c*c));
        sumErr += (error[i] * error[i]);
    }
    Float64 SNRThres = 0.0001;
    if(sumGaussA/sumErr < SNRThres){
        A = 0.0;
    }
    */

    return A;
}


Float64 COperatorLineModel::FitBayesWidth( CSpectrumSpectralAxis& spectralAxis, CSpectrumFluxAxis& fluxAxis, Float64 z, Int32 start, Int32 end)
{
    Float64 A = boost::numeric::bounds<float>::lowest();
    const Float64* flux = fluxAxis.GetSamples();
    const Float64* spectral = spectralAxis.GetSamples();
    //const Float64* error = fluxAxis.GetError();

    //A = max, good guess ?
    for ( Int32 i = start; i < end; i++)
    {
        Float64 y = flux[i];
        if(y>A){
            A = y;
        }
    }

    if(A<=0){
        return 0.0;
    }
    //c fitting iteration loop
    Float64 mu = z;
    Float64 c = 0.0001;
    Float64 cmax = 0.05;
    Int32 maxIteration = 500;
    Float64 cstepup = (cmax-c)/((Float64)(maxIteration+1));
    Float64 sum2 = boost::numeric::bounds<float>::highest();
    Float64 minsum2 = boost::numeric::bounds<float>::highest();
    Float64 minc = c;

    Int32 icmpt = 0;
    while( icmpt<maxIteration){
        sum2 = 0.0;
        for ( Int32 i = start; i < end; i++)
        {
            Float64 x = spectral[i];
            Float64 Yi = A * exp (-1.*(x-mu)*(x-mu)/(2*c*c));
            sum2 += pow( Yi - flux[i] , 2.0 );
            //sum2 += pow( Yi - flux[i] , 2.0 ) / pow( error[i], 2.0 );
        }
        if(sum2<minsum2){
            minc = c;
            minsum2 = sum2;
        }
        icmpt++;
        c = c+cstepup;
    }

    if(minc<0){
        minc=0;
    }
    return minc;
}

Float64 COperatorLineModel::FitAmplitudeIterative( const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis, Float64 lambda, Float64 width, Int32 start, Int32 end)
{
    Float64 A = boost::numeric::bounds<float>::lowest();
    const Float64* flux = fluxAxis.GetSamples();
    const Float64* spectral = spectralAxis.GetSamples();
    const Float64* error = fluxAxis.GetError();

    //A first guess
    for ( Int32 i = start; i < end; i++)
    {
        Float64 y = flux[i];
        if(y>A){
            A = y;
        }
    }

    if(A<=0){
        return 0.0;
    }
    //A fitting iteration loop
    A = A*1.5;
    Float64 mu = lambda;
    Float64 c = width;
    Float64 thres = 1e-5;
    Int32 maxIteration = 100;
    Float64 AstepDown = A/((Float64)(maxIteration+1));
    Float64 sum2 = boost::numeric::bounds<float>::highest();
    Float64 sum2prev = boost::numeric::bounds<float>::highest();
    Int32 icmpt = 0;
    while( sum2prev>=sum2 && sum2>thres && icmpt<maxIteration){
        sum2prev = sum2;
        sum2 = 0.0;
        for ( Int32 i = start; i < end; i++)
        {
            Float64 x = spectral[i];
            Float64 Yi = A * exp (-1.*(x-mu)*(x-mu)/(2*c*c));
            //sum2 += Yi-flux[i];
            sum2 += pow( Yi - flux[i] , 2.0 ) / pow( error[i], 2.0 );
        }
        //sum2 /= (Float64)(end-start+1);
        icmpt++;
        A = A-AstepDown;
    }

    if(A<0){
        A=0;
    }
    return A;
}


Void COperatorLineModel::Apply2LinesAmplitudeRule(const CRayCatalog::TRayVector& restRayList, std::vector<Float64>& Amplitudes, std::vector<Bool> outsidelambdarange,
                                                  std::string lineA, std::string lineB, Float64 coeff )
{
    Int32 iA = -1;
    for( UInt32 iRestRay=0; iRestRay<restRayList.size(); iRestRay++ )
    {
        std::string name = restRayList[iRestRay].GetName();
        std::size_t foundstra = name.find(lineA.c_str());
        if (foundstra!=std::string::npos){
            if(!outsidelambdarange[iRestRay]){
                iA = iRestRay;
            }
            break;
        }
    }

    if(iA == -1){
        return;
    }

    for( UInt32 iRestRay=0; iRestRay<restRayList.size(); iRestRay++ )
    {
        if(iA == iRestRay){
            continue;
        }
        std::string name = restRayList[iRestRay].GetName();
        std::size_t foundstra = name.find(lineB.c_str());
        if (foundstra!=std::string::npos){
            if(Amplitudes[iRestRay] > coeff*Amplitudes[iA]){
                Amplitudes[iRestRay] = coeff*Amplitudes[iA];
            }
        }
    }

}
