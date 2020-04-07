#include <RedshiftLibrary/statistics/deltaz.h>

#include <RedshiftLibrary/log/log.h>

using namespace NSEpic;
using namespace std;
#include <fstream>

#include <gsl/gsl_multifit.h>

CDeltaz::CDeltaz()
{

}

CDeltaz::~CDeltaz()
{

}

/**
 * @brief CDeltaz::Compute
 * @param merits
 * @param redshifts
 * @param redshift
 * @param redshiftRange
 * @param, output: sigma
 * @return 0: success, 1:problem with indexes,
 */
Int32 CDeltaz::Compute(TFloat64List merits, TFloat64List redshifts, Float64 redshift, TFloat64Range redshiftRange, Float64& sigma)
{
    Bool verbose = false;

    sigma = -1.0;

    //find indexes: iz, izmin and izmax
    Int32 izmin= -1;
    Int32 iz= -1;
    Int32 izmax= -1;
    for( Int32 i2=0; i2<redshifts.size(); i2++ )
    {
        if(iz == -1 && redshift <= redshifts[i2]){
            iz = i2;
        }
        if(izmin == -1 && (redshiftRange.GetBegin()) <= redshifts[i2]){
            izmin = i2;
        }
        if(izmax == -1 && (redshiftRange.GetEnd()) <= redshifts[i2]){
            izmax = i2;
            break;
        }
    }
    if(izmin == -1 || izmax == -1 || iz == -1){
        return 1;
    }

    //quadratic fit
    Int32 i, n;
    Float64 xi, yi, ei, chisq;
    gsl_matrix *X, *cov;
    gsl_vector *y, *w, *c;

    n = izmax - izmin +1;

    X = gsl_matrix_alloc (n, 1);
    y = gsl_vector_alloc (n);
    w = gsl_vector_alloc (n);

    c = gsl_vector_alloc (1);
    cov = gsl_matrix_alloc (1, 1);

    Float64 x0 = redshifts[iz];
    Float64 y0 = merits[iz];
    for (i = 0; i < n; i++)
    {
        xi = redshifts[i+izmin];
        yi = merits[i+izmin]-y0;
        if(verbose){
            fprintf (stderr, "  x = %+.5e,  y = %+.5e\n",xi, yi);
        }
        ei = 1.0; //todo, estimate weighting ?
        gsl_matrix_set (X, i, 0, (xi-x0)*(xi-x0));

        gsl_vector_set (y, i, yi);
        gsl_vector_set (w, i, 1.0/(ei*ei));
    }

    {
        gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, 1);
        gsl_multifit_wlinear (X, w, y, c, cov, &chisq, work);
        gsl_multifit_linear_free (work);
    }

#define C(i) (gsl_vector_get(c,(i)))
#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))

    sigma = sqrt(1.0/C(0));

    if(verbose){
        Log.LogInfo("Center Redshift: %g", x0);
        Log.LogInfo("# best fit: Y = %g + %g X + %g X^2\n", y0, 0.0, C(0));
        fprintf (stderr, "# best fit: Y = %g + %g X + %g X^2\n", y0, 0.0, C(0));

        Log.LogInfo("# covariance matrix:\n");
        Log.LogInfo("[ %+.5e \n", COV(0,0));

        fprintf (stderr, "# covariance matrix:\n");
        fprintf (stderr, "[ %+.5e \n", COV(0,0));

        Log.LogInfo("# chisq/n = %g", chisq/n);
        Log.LogInfo("# sigma = %g", sigma);
        Log.LogInfo("\n");

        fprintf (stderr, "# sigma = %g\n", sigma);
        fprintf (stderr, "# chisq/n = %g\n", chisq/n);
    }

    gsl_matrix_free (X);
    gsl_vector_free (y);
    gsl_vector_free (w);
    gsl_vector_free (c);
    gsl_matrix_free (cov);

    return 0;
}

//todo : merge with previous function instead of duplicating code...
Int32 CDeltaz::Compute3ddl(TFloat64List merits, TFloat64List redshifts, Float64 redshift, TFloat64Range redshiftRange, Float64& sigma)
{
    sigma = -1.0; //default value
    Bool verbose = false;

    //find iz, izmin and izmax
    Int32 izmin= -1;
    Int32 iz= -1;
    Int32 izmax= -1;
    for( Int32 i2=0; i2<redshifts.size(); i2++ )
    {
        if(iz == -1 && redshift <= redshifts[i2]){
            iz = i2;
        }
        if(izmin == -1 && (redshiftRange.GetBegin()) <= redshifts[i2]){
            izmin = i2;
        }
        if(izmax == -1 && (redshiftRange.GetEnd()) <= redshifts[i2]){
            izmax = i2;
            break;
        }
    }

    if(izmin == -1 || izmax == -1){
        return 1;
    }

    //quadratic fit
    Int32 i, n;
    Float64 xi, yi, ei, chisq;
    gsl_matrix *X, *cov;
    gsl_vector *y, *w, *c;

    n = izmax - izmin +1;
    if(n<3)
    {
        return 1;
    }

    X = gsl_matrix_alloc (n, 3);
    y = gsl_vector_alloc (n);
    w = gsl_vector_alloc (n);

    c = gsl_vector_alloc (3);
    cov = gsl_matrix_alloc (3, 3);

    double x0 = redshift;
    for (i = 0; i < n; i++)
    {
        xi = redshifts[i+izmin];
        yi = merits[i+izmin];
        if(verbose){
            fprintf (stderr, "  y = %+.5e ]\n", yi);
            fprintf (stderr, "  x = %+.5e ]\n", xi);
        }
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
    Float64 cc_ = C(2);
    sigma = sqrt(1.0/cc_);
    //sigma = sqrt(1.0/C(2));

    //Float64 a = (Float64)(C(0));
    //Float64 b2sur4c = (Float64)(C(1)*C(1)/((Float64)(4.0*C(2))));
    //Float64 logK = ( -(a - b2sur4c)/2.0 );
    //Float64 logarea = log(sigma) + logK + log(2.0*M_PI);
    if(verbose){
        Log.LogInfo("Center Redshift: %g", x0);
        Log.LogInfo("# best fit: Y = %g + %g X + %g X^2", C(0), C(1), C(2));
        fprintf (stderr, "# best fit: Y = %g + %g X + %g X^2\n", C(0), C(1), C(2));
        if(1) //debug
        {
            Log.LogInfo("# covariance matrix:\n");
            Log.LogInfo("[ %+.5e, %+.5e, %+.5e  \n", COV(0,0), COV(0,1), COV(0,2));
            Log.LogInfo("  %+.5e, %+.5e, %+.5e  \n", COV(1,0), COV(1,1), COV(1,2));
            Log.LogInfo("  %+.5e, %+.5e, %+.5e ]\n", COV(2,0), COV(2,1), COV(2,2));

            fprintf (stderr, "# covariance matrix:\n");
            fprintf (stderr, "[ %+.5e, %+.5e, %+.5e  \n", COV(0,0), COV(0,1), COV(0,2));
            fprintf (stderr, "  %+.5e, %+.5e, %+.5e  \n", COV(1,0), COV(1,1), COV(1,2));
            fprintf (stderr, "  %+.5e, %+.5e, %+.5e ]\n", COV(2,0), COV(2,1), COV(2,2));
        }
        Log.LogInfo("# chisq/n = %g", chisq/n);
        Log.LogInfo("# zcorr = %g", zcorr);
        Log.LogInfo("# sigma = %g", sigma);
        //Log.LogInfo("# logarea = %g", logarea);
        Log.LogInfo("\n");

        fprintf (stderr, "# sigma = %g\n", sigma);
        fprintf (stderr, "# chisq/n = %g\n", chisq/n);
    }

    gsl_matrix_free (X);
    gsl_vector_free (y);
    gsl_vector_free (w);
    gsl_vector_free (c);
    gsl_matrix_free (cov);

    //results.LogArea[indz] = logarea;
    //results.SigmaZ[indz] = sigma;
    //results.LogAreaCorrectedExtrema[indz] = zcorr;

    return 0;
}
