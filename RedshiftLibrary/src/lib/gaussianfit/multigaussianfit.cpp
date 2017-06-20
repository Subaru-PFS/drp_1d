#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/gaussianfit/multigaussianfit.h>
#include <RedshiftLibrary/spectrum/spectrum.h>

#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_multimin.h>

using namespace NSEpic;

/**
 * Stack attribution constructor.
 */
CMultiGaussianFit::CMultiGaussianFit( ) :
    m_PolyOrder( 2 ),
    m_AbsTol( 0.0 ),
    m_RelTol( 1e-2 ),
    m_Amplitude( 0.0 ),
    m_AmplitudeErr( 0.0 ),
    m_Mu( 0.0 ),
    m_MuErr( 0.0 ),
    m_C( 0.0 ),
    m_CErr( 0.0 ),
    m_coeff0( 0.0)
{

}

/**
 * Empty destructor.
 */
CMultiGaussianFit::~CMultiGaussianFit()
{

}


/**
 * 
 */
Int32 CMultiGaussianFit::Compute( CLineModelElementList model )
{
    SUserData userData;
    userData.model = &model;

    Int32 n = model.GetModelValidElementsNDdl();
    if( n<1 ){
        return -1;
    }
    std::vector<Int32> modelIdx = model.GetModelValidElementsIndexes();
    userData.modelIdx = modelIdx;
    userData.nddl = n;

     const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
     gsl_multimin_fminimizer *s = NULL;
     gsl_vector *ss, *x;
     gsl_multimin_function minex_func;

     size_t iter = 0;
     Int32 maxIter = 200;
     Float64 minSizeSimplex = 0.0001;
     int status;
     double size;

     /* Starting point */
     x = gsl_vector_alloc (n);
     for(Int32 j=0; j<n; j++)
       {
         gsl_vector_set( x, j, model.GetElementAmplitude( modelIdx[j] ) );
       }

     // Set initial step sizes
     ss = gsl_vector_alloc( n );
     for( Int32 j=0; j<n; j++ )
       {
         gsl_vector_set( ss, j, model.GetElementAmplitude( modelIdx[j] )/10.0 );
       }

     /* Initialize method and iterate */
     minex_func.n = n;
     minex_func.f = &my_f;
     minex_func.params = &userData;

     s = gsl_multimin_fminimizer_alloc( T, n );
     gsl_multimin_fminimizer_set( s, &minex_func, x, ss );

     do
       {
         iter++;
         status = gsl_multimin_fminimizer_iterate( s );

         if( status )
           break;

         size = gsl_multimin_fminimizer_size( s );
         status = gsl_multimin_test_size( size, minSizeSimplex );

         if( status == GSL_SUCCESS )
	   {
	     printf( "converged to minimum at\n" );
	   }

         if( 1 )
	   {
             if( n>3 )
	       {
                 printf( "%5d a0=%10.3e a1=%10.3e a2=%10.3e a3=%10.3e f() = %7.3f size = %.12f\n",
                         iter,
                         gsl_vector_get( s->x, 0 ),
                         gsl_vector_get( s->x, 1 ),
                         gsl_vector_get( s->x, 2 ),
                         gsl_vector_get( s->x, 3 ),
                         s->fval, size );
	       }
	     else if( n==2 )
	       {
                 printf( "%5d a0=%10.3e a1=%10.3e f() = %7.3f size = %.12f\n",
                         iter,
                         gsl_vector_get( s->x, 0 ),
                         gsl_vector_get( s->x, 1 ),
                         s->fval, size );
	       } 
	     else
	       {
                 printf( "%5d a0=%10.3e f() = %7.3f size = %.12f\n",
                         iter,
                         gsl_vector_get( s->x, 0 ),
                         s->fval, size );
	       }
	   }
       }
     while( status==GSL_CONTINUE && iter<maxIter );

     if( status==GSL_SUCCESS || iter>=maxIter )
       {
         for( Int32 j=0; j<n; j++ )
	   {
             Float64 a = gsl_vector_get( s->x, j );
             model.SetElementAmplitude( modelIdx[j], a, 0.0 ); //todo: add support for sigma error, now=0.0
	   }
         model.refreshModel();
       }

     printf( "ITER, END: %5d size = %.12f\n",
             iter, size );
     gsl_vector_free( x );
     gsl_vector_free( ss );
     gsl_multimin_fminimizer_free( s );

     return status;
}

/**
 * function used to model a fit of multi Gaussians
 *
 * mu = c1
 * Y[i] = c0*exp(-1*(x-mu)**2/c2**2)
 *
 */
Float64 CMultiGaussianFit::my_f( const gsl_vector *v, void *data )
{
  SUserData* userData = (SUserData*)data;
  Int32 n = userData->nddl;

  Float64 pond = 1.0;
  for( Int32 j=0; j<n; j++ )
    {
      Float64 a = (gsl_vector_get( v, j ));
      if( a<0 )
	{
          pond = pond*(1+std::abs( a ));
	}
      userData->model->SetElementAmplitude( userData->modelIdx[j], a, 0.0 );
    }
  userData->model->refreshModel();

  Float64 fit = userData->model->getLeastSquareMeritUnderElements();
  fit = fit*pond;
  return (fit);
}
