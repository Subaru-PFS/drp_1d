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

#include <math.h>
#include <float.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <algorithm>    // std::sort
#include <gsl/gsl_multifit.h>
#include <assert.h>

#define NOT_OVERLAP_VALUE NAN
#include <stdio.h>

using namespace NSEpic;
using namespace std;


IMPLEMENT_MANAGED_OBJECT(COperatorLineModel)

COperatorLineModel::COperatorLineModel()
{
    mSumLogErr=0.0;
}

COperatorLineModel::~COperatorLineModel()
{

}


const COperatorResult* COperatorLineModel::Compute(const CSpectrum& spectrum, const CSpectrum& spectrumContinuum, const CRayCatalog& restraycatalog,
                          const TFloat64Range& lambdaRange, const TFloat64List& redshifts, Int32 lineWidthType)
{

    if( spectrum.GetSpectralAxis().IsInLinearScale() == false)
    {
        Log.LogError("Line Model, input spectrum is not in log scale (ignored)");
        //return NULL;
    }

    TFloat64List sortedRedshifts = redshifts;
    std::sort(sortedRedshifts.begin(), sortedRedshifts.end());

    Int32 typeFilter = -1;//CRay::nType_Absorption;//CRay::nType_Emission;
    Int32 forceFilter = -1;//CRay::nForce_Strong;
    CRayCatalog::TRayVector restRayList = restraycatalog.GetFilteredList(typeFilter, forceFilter);

    CLineModelResult* result = new CLineModelResult();
    result->ChiSquare.resize( sortedRedshifts.size() );
    result->Redshifts.resize( sortedRedshifts.size() );
    result->Redshifts = sortedRedshifts;
    result->Status.resize( sortedRedshifts.size() );
    result->restRayList = restRayList;
    result->LineModelSolutions.resize( sortedRedshifts.size() );


    CLineModelElementList model(spectrum, spectrumContinuum, restRayList, lineWidthType);
    //model.LoadContinuum(); //in order to use a fit with continuum

    PrecomputeLogErr(spectrum);
    for (Int32 i=0;i<sortedRedshifts.size();i++)
    {
        ModelFit( spectrum, model, result->restRayList, lambdaRange, result->Redshifts[i], result->ChiSquare[i], result->LineModelSolutions[i]);
    }

    // extrema
    Int32 extremumCount = 15;
    TPointList extremumList;
    TFloat64Range redshiftsRange(result->Redshifts[0], result->Redshifts[result->Redshifts.size()-1]);
    CExtremum extremum( redshiftsRange, extremumCount, true, 1);
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
    extremumCount = extremumList.size();
    result->Extrema.resize( extremumCount );
    result->LogArea.resize( extremumCount );
    result->LogAreaCorrectedExtrema.resize( extremumCount );
    result->SigmaZ.resize( extremumCount );
    result->bic.resize( extremumCount );
    Int32 start = spectrum.GetSpectralAxis().GetIndexAtWaveLength(lambdaRange.GetBegin());
    Int32 end = spectrum.GetSpectralAxis().GetIndexAtWaveLength(lambdaRange.GetEnd());
    Int32 nsamples = end - start + 1;
    for( Int32 i=0; i<extremumList.size(); i++ )
    {
        Float64 z = extremumList[i].X;
        Float64 m = extremumList[i].Y;

        //find the index in the zaxis results
        Int32 idx=0;
        for ( UInt32 i2=0; i2<result->Redshifts.size(); i2++)
        {
            if(result->Redshifts[i2] == z){
                idx = i2;
                break;
            }
        }

        result->Extrema[i] = z;
        result->LogArea[i] = -DBL_MAX;
        result->LogAreaCorrectedExtrema[i] = -1.0;


        Int32 nddl = model.GetNElements(); //get the total number of elements in the model
        nddl = result->LineModelSolutions[idx].nDDL; //override nddl by the actual number of elements in the fitted model

        //result->bic[i] = m + nddl*log(nsamples); //BIC
        Float64 aic = m + 2*nddl; //AIC
        result->bic[i] = aic;
        //result->bic[i] = aic + (2*nddl*(nddl+1) )/(nsamples-nddl-1);  //AICc, better when nsamples small
    }
    ComputeArea2(result);


    if(result->Extrema.size()>0){
        Log.LogInfo( "LineModel Solution: best z found = %.5f", result->Extrema[0]);
    }else{
        Log.LogInfo( "LineModel Solution: no extrema found...");
    }

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
        Bool retVal1 = writer.Write( "model.fits",  spcmodel);

        if(retVal1){
            CSpectrum s(spectrum);
            Bool retVal2 = writer.Write( "spectrum.fits",  s);
        }
    }
    //*/
    return result;

}

void COperatorLineModel::ComputeArea1(CLineModelResult* results)
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
        results->LogArea[i] = area;
    }
}

///
/// \brief COperatorLineModel::ComputeArea2
/// computes the Laplace approx for a given Chi2 result around the N best extrema
///
void COperatorLineModel::ComputeArea2(CLineModelResult* results)
{
    Float64 maxp = DBL_MIN;
    for( Int32 i2=0; i2<results->Redshifts.size(); i2++ )
    {
        if( maxp < results->ChiSquare[i2]){
            maxp = results->ChiSquare[i2];
        }
    }
    Float64 winsize = 0.001;
    Float64 inclusionThresRatio = 0.01;
    Int32 iz0=0;
    for( Int32 indz=0; indz<results->Extrema.size(); indz++ )
    {
        //find iz, izmin and izmax
        Int32 izmin= -1;
        Int32 iz= -1;
        Int32 izmax= -1;
        for( Int32 i2=0; i2<results->Redshifts.size(); i2++ )
        {
            if(iz == -1 && (results->Extrema[indz]) <= results->Redshifts[i2]){
                iz = i2;
                if(indz==0){
                    iz0=iz;
                }
            }
            if(izmin == -1 && (results->Extrema[indz] - winsize/2.0) <= results->Redshifts[i2]){
                izmin = i2;
            }
            if(izmax == -1 && (results->Extrema[indz] + winsize/2.0) <= results->Redshifts[i2]){
                izmax = i2;
                break;
            }
        }
        Float64 di = abs(results->ChiSquare[iz]-maxp) ;
        Float64 d0 = abs(results->ChiSquare[iz0]-maxp) ;
        if( di < inclusionThresRatio*d0){
            continue;
        }
        if(izmin == -1 || izmax == -1){
            continue;
        }

        //quadratic fit
        int i, n;
        double xi, yi, ei, chisq;
        gsl_matrix *X, *cov;
        gsl_vector *y, *w, *c;

        n = izmax - izmin +1;

        X = gsl_matrix_alloc (n, 3);
        y = gsl_vector_alloc (n);
        w = gsl_vector_alloc (n);

        c = gsl_vector_alloc (3);
        cov = gsl_matrix_alloc (3, 3);

        double x0 = results->Extrema[indz];
        for (i = 0; i < n; i++)
        {
            xi = results->Redshifts[i+izmin];
            yi = results->ChiSquare[i+izmin];
            ei = 1.0; //todo, estimate weighting ?
            gsl_matrix_set (X, i, 0, 1.0);
            gsl_matrix_set (X, i, 1, xi-x0);
            gsl_matrix_set (X, i, 2, (xi-x0)*(xi-x0));

            gsl_vector_set (y, i, yi);
            gsl_vector_set (w, i, 1.0/(ei*ei));
        }

        {
          gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, 3);
          gsl_multifit_wlinear (X, w, y, c, cov, &chisq, work);
          gsl_multifit_linear_free (work);
        }

#define C(i) (gsl_vector_get(c,(i)))
#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))

        double zcorr = x0-C(1)/(2.0*C(2));
        double sigma = sqrt(1.0/C(2));
        Float64 a = (Float64)(C(0));
        Float64 b2sur4c = (Float64)(C(1)*C(1)/((Float64)(4.0*C(2))));
        Float64 logK = ( -(a - b2sur4c)/2.0 );
        Float64 logarea = log(sigma) + logK + log(2.0*M_PI);
        if(0){
            Log.LogInfo("Extrema: %g", results->Extrema[indz]);
            Log.LogInfo("# best fit: Y = %g + %g X + %g X^2", C(0), C(1), C(2));
            //Log.LogInfo("# covariance matrix:\n");
            //Log.LogInfo("[ %+.5e, %+.5e, %+.5e  \n", COV(0,0), COV(0,1), COV(0,2));
            //Log.LogInfo("  %+.5e, %+.5e, %+.5e  \n", COV(1,0), COV(1,1), COV(1,2));
            //Log.LogInfo("  %+.5e, %+.5e, %+.5e ]\n", COV(2,0), COV(2,1), COV(2,2));
            Log.LogInfo("# chisq/n = %g", chisq/n);
            Log.LogInfo("# zcorr = %g", zcorr);
            Log.LogInfo("# sigma = %g", sigma);
            Log.LogInfo("# logarea = %g", logarea);
            Log.LogInfo("\n");
        }

        gsl_matrix_free (X);
        gsl_vector_free (y);
        gsl_vector_free (w);
        gsl_vector_free (c);
        gsl_matrix_free (cov);

        results->LogArea[indz] = logarea;
        results->SigmaZ[indz] = sigma;
        results->LogAreaCorrectedExtrema[indz] = zcorr;
    }
}

Void COperatorLineModel::ModelFit(const CSpectrum& spectrum, CLineModelElementList& model, const CRayCatalog::TRayVector& restRayList,
                                const TFloat64Range& lambdaRange, Float64 redshift,
                                   Float64& chiSquare, CLineModelResult::SLineModelSolution& modelSolution)
{
    chiSquare = boost::numeric::bounds<float>::highest();

    model.fit(redshift, lambdaRange, modelSolution);

    Float64 fit = model.getLeastSquareMerit(lambdaRange);

    chiSquare = fit;// + mSumLogErr;
    return;
}

Float64 COperatorLineModel::FitBayesWidth( CSpectrumSpectralAxis& spectralAxis, CSpectrumFluxAxis& fluxAxis, Float64 z, Int32 start, Int32 end)
{
    Float64 A = boost::numeric::bounds<float>::lowest();
    const Float64* flux = fluxAxis.GetSamples();
    const Float64* spectral = spectralAxis.GetSamples();
    //const Float64* error = fluxAxis.GetError();

    //A = max, good value ?
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

Float64 COperatorLineModel::PrecomputeLogErr(const CSpectrum& spectrum)
{
    const CSpectrumSpectralAxis& spcSpectralAxis = spectrum.GetSpectralAxis();
    const CSpectrumFluxAxis& spcFluxAxis = spectrum.GetFluxAxis();

    Int32 numDevs = 0;
    Float64 logerrsum = 0.0;
    const Float64* error = spcFluxAxis.GetError();
    for( UInt32 j=0; j<spcSpectralAxis.GetSamplesCount(); j++ )
    {
        numDevs++;
        logerrsum += log(error[j]);
    }
    logerrsum *= 2.0;
    logerrsum += (Float64)numDevs*log(2*M_PI);

    mSumLogErr = logerrsum;
}
