#include <epic/redshift/operator/linemodel.h>
#include <epic/redshift/spectrum/axis.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/spectrum/tools.h>
#include <epic/redshift/common/mask.h>
#include <epic/redshift/operator/chisquareresult.h>
#include <epic/redshift/extremum/extremum.h>
#include <epic/redshift/linemodel/modelspectrumresult.h>
#include <epic/redshift/linemodel/modelfittingresult.h>
#include <epic/redshift/spectrum/io/fitswriter.h>
#include <epic/core/log/log.h>

#include <boost/numeric/conversion/bounds.hpp>
#include "boost/format.hpp"

#include <epic/redshift/processflow/datastore.h>

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

/**
 * \brief Attributes 0 to mSumLogErr.
 **/
COperatorLineModel::COperatorLineModel()
{
    mSumLogErr=0.0;
}

/**
 * \brief Empty destructor.
 **/
COperatorLineModel::~COperatorLineModel()
{

}

/**
 * \brief Call the Linemodel to fit models and select the parts of the results that are relevant.
 * Print an error if the spectral axis is not in linear scale.
 * Sort the list of redshifts.
 * Get the list of rest lines filtered by type and force.
 * Create a CLineModelResult object, resize its members and populate it with the filtered list.
 * Create an ElementList object.
 * Call PrecomputeLogErr on the argument spectrum.
 * For each sorted redshift value, call ModelFit.
 * Find the extrema in the set of per-redshift results.
 * Refine Extremum with a second maximum search around the z candidates.
 * Filter out peaks outside intervals centered on the refined extrema.
 * Save the results.
 **/
std::shared_ptr<COperatorResult> COperatorLineModel::Compute( CDataStore &dataStore,
							      const CSpectrum& spectrum,
							      const CSpectrum& spectrumContinuum,
							      const CRayCatalog& restraycatalog,
							      const std::string& opt_lineTypeFilter,
							      const std::string& opt_lineForceFilter,
							      const TFloat64Range& lambdaRange,
							      const TFloat64List& redshifts,
							      const Int32 opt_extremacount,
							      const std::string& opt_fittingmethod,
							      const std::string& opt_continuumcomponent,
							      const std::string& opt_lineWidthType,
							      const Float64 opt_resolution,
							      const Float64 opt_velocityEmission,
							      const Float64 opt_velocityAbsorption,
							      const std::string& opt_continuumreest,
							      const std::string& opt_rules )
{
  if( spectrum.GetSpectralAxis().IsInLinearScale()==false )
    {
      Log.LogError( "Line Model, input spectrum is not in linear scale (ignored)." );
    }
  TFloat64List sortedRedshifts = redshifts;
  std::sort( sortedRedshifts.begin(), sortedRedshifts.end() );
  Int32 typeFilter = -1;
  if( opt_lineTypeFilter == "A" )
    {
      typeFilter = CRay::nType_Absorption;
    }
  if( opt_lineTypeFilter == "E" )
    {
      typeFilter = CRay::nType_Emission;
    }
  Int32 forceFilter = -1;
  if( opt_lineForceFilter == "S" )
    {
        forceFilter = CRay::nForce_Strong;
    }
  CRayCatalog::TRayVector restRayList = restraycatalog.GetFilteredList(typeFilter, forceFilter);

  auto result = std::shared_ptr<CLineModelResult>( new CLineModelResult() );
  result->ChiSquare.resize( sortedRedshifts.size() );
  result->Redshifts.resize( sortedRedshifts.size() );
  result->Redshifts = sortedRedshifts;
  result->Status.resize( sortedRedshifts.size() );
  result->restRayList = restRayList;
  result->LineModelSolutions.resize( sortedRedshifts.size() );

  CLineModelElementList model( spectrum, spectrumContinuum, restRayList, opt_fittingmethod, opt_continuumcomponent, opt_lineWidthType, opt_resolution, opt_velocityEmission, opt_velocityAbsorption, opt_rules );
  //model.LoadContinuum(); //in order to use a fit with continuum
  result->nSpcSamples = model.getSpcNSamples(lambdaRange);

  PrecomputeLogErr( spectrum );

  Int32 contreest_iterations = 0;
  if( opt_continuumreest == "always" )
    {
      contreest_iterations = 1;
    }
  else
    {
      contreest_iterations = 0;
    }

  for( Int32 i=0; i<sortedRedshifts.size(); i++ )
    {
      ModelFit( model, lambdaRange, result->Redshifts[i], result->ChiSquare[i], result->LineModelSolutions[i], contreest_iterations);
    }

  // extrema
  Int32 extremumCount = opt_extremacount;
  TPointList extremumList;
  TFloat64Range redshiftsRange( result->Redshifts[0], result->Redshifts[result->Redshifts.size()-1] );
  CExtremum extremum( redshiftsRange, extremumCount, true, 2 );
  extremum.Find( result->Redshifts, result->ChiSquare, extremumList );

  // Refine Extremum with a second maximum search around the z candidates:
  // This corresponds to the finer xcorrelation in EZ Pandora (in standard_DP fctn in SolveKernel.py)
  Float64 radius = 0.001;
  for( Int32 i=0; i<extremumList.size(); i++ )
    {
      Float64 x = extremumList[i].X;
      Float64 left_border = max( redshiftsRange.GetBegin(), x-radius );
      Float64 right_border = min( redshiftsRange.GetEnd(), x+radius );
      TPointList extremumListFine;
      TFloat64Range rangeFine = TFloat64Range( left_border, right_border );
      CExtremum extremumFine( rangeFine , 1, true );
      extremumFine.Find( result->Redshifts, result->ChiSquare, extremumListFine );
      if( extremumListFine.size()>0 )
	{
	  extremumList[i] = extremumListFine[0];
        }
    }

  // extend z around the extrema
  Float64 extensionradius = 0.01;
  for( Int32 i=0; i<extremumList.size(); i++ )
    {
      Float64 x = extremumList[i].X;
      Float64 left_border = max( redshiftsRange.GetBegin(), x-extensionradius);
      Float64 right_border = min( redshiftsRange.GetEnd(), x+extensionradius);
      for( Int32 i=0; i<result->Redshifts.size(); i++ )
        {
	  if( result->Redshifts[i]>=left_border && result->Redshifts[i]<=right_border )
	    {
	      result->ExtremaExtendedRedshifts.push_back(result->Redshifts[i]);
	    }
        }
    }

  //todo: remove duplicate redshifts from the extended extrema list

  // store extrema results
  extremumCount = extremumList.size();
  result->Extrema.resize( extremumCount );
  result->Posterior.resize( extremumCount );
  result->LogArea.resize( extremumCount );
  result->LogAreaCorrectedExtrema.resize( extremumCount );
  result->SigmaZ.resize( extremumCount );
  result->bic.resize( extremumCount );
  Int32 start = spectrum.GetSpectralAxis().GetIndexAtWaveLength( lambdaRange.GetBegin() );
  Int32 end = spectrum.GetSpectralAxis().GetIndexAtWaveLength( lambdaRange.GetEnd() );
  Int32 nsamples = end - start + 1;
  Int32 savedModels = 0;
  for( Int32 i=0; i<extremumList.size(); i++ )
    {
      Float64 z = extremumList[i].X;
      Float64 m = extremumList[i].Y;

      //find the index in the zaxis results
      Int32 idx=-1;
      for ( UInt32 i2=0; i2<result->Redshifts.size(); i2++ )
        {
	  if( result->Redshifts[i2] == z )
	    {
	      idx = i2;
	      break;
            }
        }
      if( idx==-1 )
	{
	  Log.LogInfo( "Problem. could not find extrema solution index..." );
	  continue;
	}

      // reestimate the model (eventually with continuum reestimation) on the extrema selected
      if( opt_continuumreest == "always" || opt_continuumreest == "onlyextrema" )
	{
	  contreest_iterations = 1;
        }
      else
	{
	  contreest_iterations  = 0;
        }
      ModelFit( model, lambdaRange, result->Redshifts[idx], result->ChiSquare[idx], result->LineModelSolutions[idx], contreest_iterations);
      m = result->ChiSquare[idx];


      //save the model result
      static Int32 maxModelSave = 5;
      if( savedModels<maxModelSave )
	{
	  // CModelSpectrumResult
	  std::shared_ptr<CModelSpectrumResult> resultspcmodel = std::shared_ptr<CModelSpectrumResult>( new CModelSpectrumResult(model.GetModelSpectrum()) );
	  std::string fname_spc = ( boost::format("linemodel_spc_extrema_%1%" ) % savedModels).str();
	  dataStore.StoreScopedGlobalResult( fname_spc.c_str(), resultspcmodel );

	  // CModelFittingResult
	  std::shared_ptr<CModelFittingResult> resultfitmodel = std::shared_ptr<CModelFittingResult>( new CModelFittingResult( result->LineModelSolutions[idx], result->Redshifts[idx], result->ChiSquare[idx], result->restRayList ) );
	  std::string fname_fit = ( boost::format( "linemodel_fit_extrema_%1%" ) % savedModels ).str();
	  dataStore.StoreScopedGlobalResult( fname_fit.c_str(), resultfitmodel );
	  savedModels++;
        }

      result->Extrema[i] = z;
      static Float64 cutThres = 5.0;
      Int32 nValidLines = result->GetNLinesOverCutThreshold(i, cutThres, cutThres);
      result->Posterior[i] = m/Float64(1+nValidLines);
      result->LogArea[i] = -DBL_MAX;
      result->LogAreaCorrectedExtrema[i] = -1.0;

      Int32 nddl = model.GetNElements(); //get the total number of elements in the model
      nddl = result->LineModelSolutions[idx].nDDL; //override nddl by the actual number of elements in the fitted model

      Float64 aic = m + 2*nddl; //AIC
      result->bic[i] = aic;
    }

  if( result->Extrema.size()>0 )
    {
      Log.LogInfo( "LineModel Solution: best z found = %.5f", result->Extrema[0] );
    }
  else
    {
      Log.LogInfo( "LineModel Solution: no extrema found...");
    }

    return result;
}

/**
 * \brief Populates LogArea with the area of the spectrum that obeys the arguments' criteria.
 * Prepare p.
 * For each extremum,
 *   Find the first redshift where the extremum is not greater than the redshift.   
 *   Ignore results that have abdolute difference between the chisquare and the current maximum below a threshold.
 *   Call FitsBayesWidth with found criteria.
 *   Accumulate the width inside the interval defined by the index limits.
 *   Set the LogArea[i] to this sum.
 **/
void COperatorLineModel::ComputeArea1( CLineModelResult& results )
{
  //prepare p
  Float64 maxp = DBL_MIN;
  CSpectrum pspc;
  CSpectrumFluxAxis& spcFluxAxis = pspc.GetFluxAxis();
  spcFluxAxis.SetSize( results.Redshifts.size() );
  CSpectrumSpectralAxis& spcSpectralAxis = pspc.GetSpectralAxis();
  spcSpectralAxis.SetSize( results.Redshifts.size() );
  for( Int32 i2=0; i2<results.Redshifts.size(); i2++ )
    {
      if( maxp < results.ChiSquare[i2] )
	{
	  maxp = results.ChiSquare[i2];
        }
    }
  for( Int32 i2=0; i2<results.Redshifts.size(); i2++ )
    {
      spcFluxAxis[i2] = exp(-(results.ChiSquare[i2]-maxp)/2.0)-1.0;
      spcSpectralAxis[i2] = results.Redshifts[i2];
    }

  //debug:
  if ( false )
    {
      FILE* f = fopen( "getbestredshiftbayes_dbg.txt", "w+" );
      for( Int32 i=0; i<spcFluxAxis.GetSamplesCount(); i++ )
	{
	  if( spcFluxAxis[i] < 0.0001 )
	    {
	      fprintf( f, "%e %e\n", spcSpectralAxis[i], spcFluxAxis[i]);
	    }
	  else
	    {
	      fprintf( f, "%f %f\n", spcSpectralAxis[i], spcFluxAxis[i]);
	    }
	}
      fclose( f );
    }

  Float64 winsize = 0.0025;
  Float64 inclusionThresRatio = 0.25;
  Int32 iz0=0;
  for( Int32 i=0; i<results.Extrema.size(); i++ )
    {
      //find iz, izmin and izmax
      Int32 izmin= -1;
      Int32 iz= -1;
      Int32 izmax= -1;
      for( Int32 i2=0; i2<results.Redshifts.size(); i2++ )
        {
	  if( iz==-1 && (results.Extrema[i])<=results.Redshifts[i2] )
	    {
	      iz = i2;
	      if(i==0)
		{
		  iz0=iz;
                }
            }
	  if( izmin==-1 && (results.Extrema[i]-winsize/2.0)<=results.Redshifts[i2] )
	    {
	      izmin = i2;
            }
	  if( izmax==-1 && (results.Extrema[i]+winsize/2.0)<=results.Redshifts[i2] )
	    {
	      izmax = i2;
	      break;
            }
        }
      Float64 di = abs( results.ChiSquare[iz]-maxp );
      Float64 d0 = abs( results.ChiSquare[iz0]-maxp );
      if( di < inclusionThresRatio*d0 )
	{
	  continue;
        }
      Float64 gaussWidth = FitBayesWidth( spcSpectralAxis, spcFluxAxis, results.Extrema[i], izmin, izmax );
      Float64 gaussAmp = spcFluxAxis[iz];
      Float64 intsize = 0.001;
      Float64 area=0.0;
      for( Int32 i2=izmin; i2<izmax; i2++ )
        {
	  Float64 x = spcSpectralAxis[i2];
	  Float64 Yi = gaussAmp * exp (-1.*(x-results.Extrema[i])*(x-results.Extrema[i])/(2*gaussWidth*gaussWidth));
	  area += Yi;
        }
      results.LogArea[i] = area;
    }
}

/**
 * \brief Computes the Laplace approx for a given Chi2 result around the N best extrema
 **/
void COperatorLineModel::ComputeArea2(CLineModelResult& results)
{
    Float64 maxp = DBL_MIN;
    for( Int32 i2=0; i2<results.Redshifts.size(); i2++ )
    {
        if( maxp < results.ChiSquare[i2]){
            maxp = results.ChiSquare[i2];
        }
    }
    Float64 winsize = 0.001;
    Float64 inclusionThresRatio = 0.01;
    Int32 iz0=0;
    for( Int32 indz=0; indz<results.Extrema.size(); indz++ )
    {
        //find iz, izmin and izmax
        Int32 izmin= -1;
        Int32 iz= -1;
        Int32 izmax= -1;
        for( Int32 i2=0; i2<results.Redshifts.size(); i2++ )
        {
            if(iz == -1 && (results.Extrema[indz]) <= results.Redshifts[i2]){
                iz = i2;
                if(indz==0){
                    iz0=iz;
                }
            }
            if(izmin == -1 && (results.Extrema[indz] - winsize/2.0) <= results.Redshifts[i2]){
                izmin = i2;
            }
            if(izmax == -1 && (results.Extrema[indz] + winsize/2.0) <= results.Redshifts[i2]){
                izmax = i2;
                break;
            }
        }
        Float64 di = abs(results.ChiSquare[iz]-maxp) ;
        Float64 d0 = abs(results.ChiSquare[iz0]-maxp) ;
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

        double x0 = results.Extrema[indz];
        for (i = 0; i < n; i++)
        {
            xi = results.Redshifts[i+izmin];
            yi = results.ChiSquare[i+izmin];
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
            Log.LogInfo("Extrema: %g", results.Extrema[indz]);
            Log.LogInfo("# best fit: Y = %g + %g X + %g X^2", C(0), C(1), C(2));
	    if( false ) //debug
	      {
		Log.LogInfo("# covariance matrix:\n");
		Log.LogInfo("[ %+.5e, %+.5e, %+.5e  \n", COV(0,0), COV(0,1), COV(0,2));
		Log.LogInfo("  %+.5e, %+.5e, %+.5e  \n", COV(1,0), COV(1,1), COV(1,2));
		Log.LogInfo("  %+.5e, %+.5e, %+.5e ]\n", COV(2,0), COV(2,1), COV(2,2));
	      }
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

        results.LogArea[indz] = logarea;
        results.SigmaZ[indz] = sigma;
        results.LogAreaCorrectedExtrema[indz] = zcorr;
    }
}

/**
 * \brief Calls model.fit() and sets the chisquare to the return value of that method.
 **/
Void COperatorLineModel::ModelFit(CLineModelElementList& model, const TFloat64Range& lambdaRange, Float64 redshift,
				  Float64& chiSquare, CLineModelResult::SLineModelSolution& modelSolution, Int32 contreest_iterations)
{
    chiSquare = boost::numeric::bounds<float>::highest();
    Float64 fit = model.fit( redshift, lambdaRange, modelSolution, contreest_iterations );
    chiSquare = fit;
}

/**
 * \brief Returns a non-negative value for the width that yields the least squared difference between the flux and a exponentially decayed maximum amplitude.
 * Find the maximum flux amplitude. If this not greater than zero, return zero.
 * For each value of c within the range:
 *   Sum the squared difference between the flux and the maximum amplitude with a exponential decay parameterized by c.
 *   Save the minimal result.
 * If the result is not greater than zero, return zero.
 * Return the result.
 **/
Float64 COperatorLineModel::FitBayesWidth( CSpectrumSpectralAxis& spectralAxis, CSpectrumFluxAxis& fluxAxis, Float64 z, Int32 start, Int32 end)
{
    Float64 A = boost::numeric::bounds<float>::lowest();
    const Float64* flux = fluxAxis.GetSamples();
    const Float64* spectral = spectralAxis.GetSamples();
    for ( Int32 i = start; i < end; i++)
    {
        Float64 y = flux[i];
        if( y>A )
	  {
            A = y;
	  }
    }

    if( A<=0 )
      {
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
    while( icmpt<maxIteration )
      {
        sum2 = 0.0;
        for ( Int32 i = start; i < end; i++)
        {
            Float64 x = spectral[i];
            Float64 Yi = A * exp (-1.*(x-mu)*(x-mu)/(2*c*c));
            sum2 += pow( Yi - flux[i] , 2.0 );
        }
        if( sum2<minsum2 )
	  {
            minc = c;
            minsum2 = sum2;
	  }
        icmpt++;
        c = c+cstepup;
      }

    if( minc<0 )
      {
        minc=0;
      }
    return minc;
}

/**
 * \brief Returns the result of an error sum calulation.
 **/
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

    return mSumLogErr;
}
