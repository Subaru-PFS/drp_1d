#include <epic/redshift/common/gaussianfit.h>

#include <epic/core/debug/assert.h>

#include <math.h>
/*
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
*/
using namespace __NS__;

CGaussianFit::CGaussianFit()
{
}

CGaussianFit::~CGaussianFit()
{

}



Bool CGaussianFit::Compute( const Float64* x, const Float64* y, Int32 n, Int32 polyOrder )
{
    /*
    Int32 size = 3 + polyOrder +1;

    std::vector<Float64> coeffs;
    std::vector<Float64> coeffErrs;
    std::vector<Float64> covar;
    std::vector<Float64> chisq;

    std::vector<Float64> err;
    err.resize( n );
*/
    //pndMathGaussFit( x, y, err.data(), n, polyOrder, 0, NULL, coeffs.data(), coeffErrs.data(), 3 + polyOrder + 1, -1, -1);

    return true;
}
/*
struct datafit {
  Int32 n;
  Float64 * x;
  Float64 * y;
  Float64 * err;
};



// Function definition
int _gauss_f (const gsl_vector *param, void *data, gsl_vector * f)
{
    int n = (int) ((struct datafit *)data)->n;
    double *x = ((struct datafit *)data)->x;
    double *y = ((struct datafit *)data)->y;
    double *err = ((struct datafit *) data)->err;

    double A = gsl_vector_get (param, 0);
    double mu = gsl_vector_get (param, 1);
    double c = gsl_vector_get (param, 2);

    // Polynomial order
    int order = (int) param->size - 3 - 1;

    int i, k;

    for (i = 0; i < n; i++) {
        // Polynomial term
        double Pi = 0;
        for(k = 0; k <= order; k++)
            Pi += gsl_vector_get(param, 3 + k) * pow(x[i] - mu, k);

        // Add gaussina term to polynomial term
        double Yi = A * exp (-1.*(x[i]-mu)*(x[i]-mu)/(c*c)) + Pi;
        gsl_vector_set (f, i, (Yi - y[i])/err[i]);
    }

  return GSL_SUCCESS;
}


// Jacobian function definition
int _gauss_df (const gsl_vector * param, void *data, gsl_matrix * J) {

    int n = (int) ((struct datafit *)data)->n;
    double *x = ((struct datafit *)data)->x;
    double *err = ((struct datafit *) data)->err;

    double A = gsl_vector_get (param, 0);
    double mu = gsl_vector_get (param, 1);
    double c = gsl_vector_get (param, 2);

    double A_d = 2*A/(c*c);

    // Polynomial order
    int order = (int) param->size - 3 - 1;

    int i, k;

    // Number of parameters
    for (i = 0; i < n; i++) {
        // Jacobian matrix J(i,j) = dfi / dxj,
        // where fi = (Yi - yi)/err[i],
        //       Yi = A * exp(-1*(xi-mu)**2/c**2) + P(n, x-mu)

        // Exponential term
        double e = exp(-1.*(x[i]-mu)*(x[i]-mu)/(c*c))/err[i];

        // Gaussian term
        gsl_matrix_set (J, i, 0, e);
        double P_mu = 0;
        for(k=1; k <= order; k++)
            P_mu+=(-1.*(k*pow(x[i]-mu, k-1)*gsl_vector_get(param, 3+k)));
        gsl_matrix_set (J, i, 1, A_d*(x[i]-mu)*e + P_mu/err[i]);
        gsl_matrix_set (J, i, 2, A_d*(x[i]-mu)*(x[i]-mu)*e/c);

        // Polynomial term
        for(k=0; k <= order; k++)
            gsl_matrix_set (J, i, 3+k, pow(x[i] - mu, k)/err[i]);
    }
    return GSL_SUCCESS;

}


int _gauss_fdf (const gsl_vector * param, void *data, gsl_vector * f, gsl_matrix * J)
{
    _gauss_f (param, data, f);
    _gauss_df (param, data, J);

    return GSL_SUCCESS;

}


// Function definition
int _gauss_origin_f (const gsl_vector *param, void *data, gsl_vector * f)
{

    int n = (int) ((struct datafit *)data)->n;
    double *x = ((struct datafit *)data)->x;
    double *y = ((struct datafit *)data)->y;
    double *err = ((struct datafit *) data)->err;

    double A = gsl_vector_get (param, 0);
    double mu = gsl_vector_get (param, 1);
    double c = gsl_vector_get (param, 2);

    // Polynomial order
    int order = (int) param->size - 3 - 1;

    int i, k;

    for (i = 0; i < n; i++) {
        // Polynomial term
        double Pi = 0;
        for(k = 0; k <= order; k++)
            Pi += gsl_vector_get(param, 3 + k) * pow(x[i], k);

        // Add gaussina term to polynomial term
        double Yi = A * exp (-1.*(x[i]-mu)*(x[i]-mu)/(c*c)) + Pi;
        gsl_vector_set (f, i, (Yi - y[i])/err[i]);
    }

  return GSL_SUCCESS;

}


// Jacobian function definition
int _gauss_origin_df (const gsl_vector * param, void *data, gsl_matrix * J) {

    int n = (int) ((struct datafit *)data)->n;
    double *x = ((struct datafit *)data)->x;
    double *err = ((struct datafit *) data)->err;

    double A = gsl_vector_get (param, 0);
    double mu = gsl_vector_get (param, 1);
    double c = gsl_vector_get (param, 2);

    double A_d = 2*A/(c*c);

    // Polynomial order
    int order = (int) param->size - 3 - 1;

    int i, k;

    // Number of parameters
    for (i = 0; i < n; i++) {
        // Jacobian matrix J(i,j) = dfi / dxj,
        // where fi = (Yi - yi)/err[i],
        //       Yi = A * exp(-1*(xi-mu)**2/c**2) + P(n, x-mu)

        // Exponential term
        double e = exp(-1.*(x[i]-mu)*(x[i]-mu)/(c*c))/err[i];

        // Gaussian term
        gsl_matrix_set (J, i, 0, e);
        double P_mu = 0;
        for(k=1; k <= order; k++)
            P_mu+=(-1.*(k*pow(x[i], k-1)*gsl_vector_get(param, 3+k)));
        gsl_matrix_set (J, i, 1, A_d*(x[i]-mu)*e + P_mu);
        gsl_matrix_set (J, i, 2, A_d*(x[i]-mu)*(x[i]-mu)*e/c);

        // Polynomial term
        for(k=0; k <= order; k++)
            gsl_matrix_set (J, i, 3+k, pow(x[i], k)/err[i]);
    }
    return GSL_SUCCESS;

}


int _gauss_origin_fdf (const gsl_vector * param, void *data, gsl_vector * f, gsl_matrix * J)
{
    _gauss_origin_f (param, data, f);
    _gauss_origin_df (param, data, J);

    return GSL_SUCCESS;

}


Bool pndMathGaussFit( const Float64* x, const Float64* y, Float64* err, Int32 n, Int32 poly_order, Int32 origin,
        Float64* p_init, Float64* p_out, Float64* p_err_out, Int32 np, Float64 abs_toll, Float64 rel_toll)
{

    const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
    gsl_multifit_fdfsolver *s;
    int status;
    int i, iter = 0;

    gsl_vector_view p;

    double *p_scratch = NULL;
    double *tmp_err = NULL;
    struct datafit d;

    // Fill data struct
    d.n = n;
    d.x = (Float64*)x;
    d.y = (Float64*)y;

    DebugAssert( err != NULL );

    for(i=0;i<n;i++)
        err[i] = 1.;

    d.err = err;

    gsl_matrix *covar = gsl_matrix_alloc (np, np);

    gsl_multifit_function_fdf f;

    if (origin == 0)
    {
        // Use polinomial part in the form P(x) = a0 + a1*(x-mu) + ... + an*(x-mu)^n
        // mu is the center of the gaussian
        f.f = &_gauss_f;
        f.df = &_gauss_df;
        f.fdf = &_gauss_fdf;
    }
    else
    {
        // Use polinomial part in the form P(x) = a0 + a1*x + ... + an*x^n
        f.f = &_gauss_origin_f;
        f.df = &_gauss_origin_df;
        f.fdf = &_gauss_origin_fdf;
    }
    f.n = n;
    f.p = np;
    f.params = &d;

    // Retrieve first guesses
    if (p_init)
    {
        p = gsl_vector_view_array (p_init, np);
    }
    else
    {
        // Create suitable first guess
        // To BE improved, trivial version now
        p_scratch = (double*) calloc (np, sizeof(double));
        p = gsl_vector_view_array (p_scratch, np);

        double *v = (double*) calloc(n ,sizeof(double));
        for(i=0;i<n;i++)
            v[i] = y[i];
        gsl_sort (v, 1, n);
        double y_median = gsl_stats_median_from_sorted_data (v, 1, n);

        // Gaussian peak +/- and peak position
        for(i=0;i<n;i++)
            v[i] =y[i] - y_median;
        double max = gsl_stats_max(v, 1, n);
        double min = gsl_stats_min(v, 1, n);

        if (fabs(max) > fabs(min))
        {
            // Peak value
            p_scratch[0] = max;
            // Peak position: mu
            p_scratch[1] = x[(int) gsl_stats_max_index(v, 1, n)];
        }
        else
        {
            // Peak value
            p_scratch[0] = min;
            // Peak position: mu
            p_scratch[1] = x[(int) gsl_stats_min_index(v, 1, n)];
        }

        // Gaussian amplitude
        double std = 0;
        for(i=0; i<n;i++)
            std+=(y[i] - std)*(y[i] - std);
        std/=(n-1);
        int count = 0;
        for (i=0; i<n; i++)
        {
            if (y[i]> y_median + 3.*std) count++;
        }
        if (count > 2)
        {
            p_scratch[2] = (x[count] -x[0])/6.;
        }
        else
        {
            p_scratch[2] = (x[n-1] - x[1]) / 6.;
        }
    }

    if (abs_toll < 0)
        abs_toll = 0.;
    if (rel_toll < 0)
        rel_toll = 1e-5;

    // Setup solver
    if (n < np)
        return false;

     s = gsl_multifit_fdfsolver_alloc (T, n, np);
     gsl_multifit_fdfsolver_set (s, &f, &p.vector);

     // Iterate
     do
     {
         iter++;
         status = gsl_multifit_fdfsolver_iterate (s);

         if (status)
             break;

         status = gsl_multifit_test_delta (s->dx, s->x, abs_toll, rel_toll);
     }
     while (status == GSL_CONTINUE && iter < 500);

     // Set values and errors
     gsl_multifit_covar (s->J, 0.0, covar);

     {
         double chi = gsl_blas_dnrm2(s->f);
         double dof = n - np;
         double c = chi / sqrt(dof);

         for(i=0;i<np;i++)
         {
             p_out[i] = gsl_vector_get(s->x, i);
             p_err_out[i] = c*sqrt(gsl_matrix_get(covar,i,i));
         }

     // Set amplitude > 0 by default
         p_out[2] = fabs(p_out[2]);

     }

     Bool error;

     if (status)
         error = false;
     else
         error = true;

     gsl_multifit_fdfsolver_free (s);
     gsl_matrix_free (covar);

     if (tmp_err)
         free(tmp_err);
     if (p_scratch)
         free(p_scratch);


    return error;
}
*/

