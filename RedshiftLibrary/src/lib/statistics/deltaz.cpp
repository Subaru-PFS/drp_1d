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

Int32 CDeltaz::Compute(TFloat64List merits, TFloat64List redshifts, Float64 redshift, TFloat64Range redshiftRange, Float64& sigma)
{
    sigma = -1.0;
    TFloat64Range refRange = TFloat64Range(redshifts[0], redshifts[redshifts.size()-1]);
    bool ret = TFloat64Range::Intersect(redshiftRange, refRange, redshiftRange);
    if(!ret){
        Log.LogError("  Deltaz: Deltaz range for candidate %f is outside the redshift range", redshift);
        throw runtime_error("Deltaz: Candidate is outside range!");
    }

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
    if( iz==-1 || izmin == -1 || izmax == -1 ){
        return 1;
    }

    Int32 n = 0;
    n = izmax - izmin +1;
    Float64 x0 = redshifts[iz];
    Float64 y0 = merits[iz];
    Float64 xi2, yi, c0; 
    Float64 sum = 0, sum2 = 0;
    for (Int32 i = 0; i < n; i++)
    {
        xi2 = redshifts[i+izmin]-x0;
        xi2 *= xi2;
        yi = y0 - merits[i+izmin]; //pour re-caler les y pour que le centre soit à zero pour x0
        if(yi<0){
            Log.LogError("    CDeltaz::  Xi2 value of zi = %f should be greater than Xi2 of Zcand = %f \n", redshifts[i+izmin], redshift);
            //two possible reasons: 1) deltaz range is too large and a neighboring peak can fall within this range and
            //the second pass updates the candidate value but didnt check that the new Zcand is a true peak on the deltaz range
            return -1;
        }
        sum += xi2*yi;
        sum2 += xi2*xi2;
    }
    c0 = sum/sum2; 
    if(c0<=0){
        return -1;
    }
    sigma = sqrt(1.0/c0);
    return 0;
}

//todo : merge with previous function instead of duplicating code...
Int32 CDeltaz::Compute3ddl(TFloat64List merits, TFloat64List redshifts, Float64 redshift, TFloat64Range redshiftRange, Float64& sigma)
{
    sigma = -1.0; //default value
    Bool verbose = false;
    
    TFloat64Range refRange = TFloat64Range(redshifts[0], redshifts[redshifts.size()-1]);
    bool ret = TFloat64Range::Intersect(redshiftRange, refRange, redshiftRange);
    if(!ret){
        Log.LogError("  Deltaz: Deltaz range for candidate %f is completely outside the redshift range", redshift);
        throw runtime_error("Deltaz: Candidate is outside range!");
    }
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
    Float64 c2 = C(2);
    if(c2>0)
        sigma = sqrt(1.0/c2);
    
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
    if(c2<=0) 
        return -1;
    return 0;
}
