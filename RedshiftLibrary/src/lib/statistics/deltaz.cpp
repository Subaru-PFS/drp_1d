#include "RedshiftLibrary/statistics/deltaz.h"
#include "RedshiftLibrary/log/log.h"

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

Float64 CDeltaz::GetDeltaz(const TFloat64List & redshifts, const TFloat64List & pdf, const Float64 z, const Int32 gslfit)
{
    Float64 dz = -1;
    if (!redshifts.size() )
        return 0;
    Int32 ret = -1, deltaz_i = 0, maxIter = 2;
    while(ret == -1 && deltaz_i < maxIter){//iterate only twice
            Int32 izmin= -1;
            Int32 iz= -1;
            Int32 izmax= -1;

            // Float64 zRangeHalf = 0.002/(deltaz_i+1); 
            // Log.LogInfo("  Deltaz: Deltaz computation nb %i with zRangeHalf %f", deltaz_i, zRangeHalf);
            // TFloat64Range range = TFloat64Range(z - zRangeHalf*(1+z), z + zRangeHalf*(1+z));
            // ret = GetRangeIndices(redshifts, z, range, iz, izmin, izmax );
            
            Int32 half_samples_nb = 5/(deltaz_i+1);
            ret = GetIndices(redshifts, z, half_samples_nb, iz, izmin, izmax);

            if(gslfit)
                ret = Compute3ddl(pdf, redshifts, iz, izmin, izmax, dz);
            else{
                ret = Compute(pdf, redshifts, iz, izmin, izmax, dz);
            }          
            if (ret == -1){
                //Log.LogWarning("  Deltaz: Deltaz computation failed for half range %f", zRangeHalf);
                Log.LogWarning("  Deltaz: Deltaz computation failed for half range %i samples", half_samples_nb);
                deltaz_i++; 
            }
    }
    if(dz == -1){    
        Log.LogError("  Deltaz: Deltaz for candidate %f couldnt be calculated", z);
        dz = 0.001; //default value
    }
    return dz;
}

Int32 CDeltaz::GetIndices(const TFloat64List & redshifts, const Float64 redshift, const Int32 HalfNbSamples, 
                          Int32 & iz, Int32 & izmin, Int32 & izmax )
{
    //find indexes: iz, izmin and izmax
    izmin= -1;
    TFloat64List::const_iterator iiz = std::lower_bound(redshifts.begin(),redshifts.end(),redshift); 
    izmax= -1;

    iz = iiz - redshifts.begin();
    if( *iiz !=redshift ){
        Log.LogError("  Deltaz: Impossible to get redshift index %f (%d)",
                     redshift, iz);
        throw runtime_error("Deltaz: impossible to get redshift indices!");
    }

    izmin = max(iz - HalfNbSamples, 0);
    izmax = min(iz + HalfNbSamples, Int32(redshifts.size()-1));

    return 0;
}

Int32 CDeltaz::GetRangeIndices(const TFloat64List & redshifts, const Float64 redshift, const TFloat64Range & redshiftRange, 
                               Int32 & iz, Int32 & izmin, Int32 & izmax )
{
    TFloat64Range effectiveRange = TFloat64Range(redshifts.front(), redshifts.back());
    bool ret = effectiveRange.IntersectWith(redshiftRange);
    if(!ret){
        Log.LogError("  Deltaz: Deltaz range for candidate %f is outside the redshift range", redshift);
        throw runtime_error("Deltaz: Candidate is outside range!");
    }
    
    //find indexes: iz, izmin and izmax
    bool ok = effectiveRange.getEnclosingIntervalIndices(const_cast<TFloat64List&>(redshifts),izmin,izmax);
    iz= std::lower_bound(redshifts.begin(),redshifts.end(),redshift) - redshifts.begin();
    
    
    if( !ok || iz == -1 ){
        Log.LogError("  Deltaz: Impossible to get redshift index %f (%d) or redshift range indices %f,%f (%d,%d)",
                     redshift, iz, effectiveRange.GetBegin(), effectiveRange.GetEnd(), izmin, izmax);
        throw runtime_error("Deltaz: impossible to get redshift indices!");
    }
    return 0;
}


Int32 CDeltaz::Compute(const TFloat64List & merits, const TFloat64List & redshifts, const Int32 iz, const Int32 izmin, const Int32 izmax, Float64& sigma)
 {
    Float64 x0 = redshifts[iz];
    Float64 y0 = merits[iz];
    Float64 xi2, yi, c0; 
    Float64 sum = 0, sum2 = 0;

    sigma = -1.0;

    for (Int32 i = izmin; i < izmax+1; i++)
    {
        xi2 = redshifts[i]-x0;
        xi2 *= xi2;
        yi = y0 - merits[i]; //pour re-caler les y pour que le centre soit à zero pour x0
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
Int32 CDeltaz::Compute3ddl(const TFloat64List &merits, const TFloat64List &redshifts, const Int32 iz, const Int32 izmin, const Int32 izmax, Float64& sigma)
{
    sigma = -1.0; //default value
    Bool verbose = false;
    
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

    double x0 = redshifts[iz];
    double y0 = merits[iz];
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

        gsl_vector_set (y, i, y0-yi);
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
        Log.LogDebug("Center Redshift: %g", x0);
        Log.LogDebug("# best fit: Y = %g + %g X + %g X^2", C(0), C(1), C(2));
        fprintf (stderr, "# best fit: Y = %g + %g X + %g X^2\n", C(0), C(1), C(2));
        if(1) //debug
        {
            Log.LogDebug("# covariance matrix:\n");
            Log.LogDebug("[ %+.5e, %+.5e, %+.5e  \n", COV(0,0), COV(0,1), COV(0,2));
            Log.LogDebug("  %+.5e, %+.5e, %+.5e  \n", COV(1,0), COV(1,1), COV(1,2));
            Log.LogDebug("  %+.5e, %+.5e, %+.5e ]\n", COV(2,0), COV(2,1), COV(2,2));

            fprintf (stderr, "# covariance matrix:\n");
            fprintf (stderr, "[ %+.5e, %+.5e, %+.5e  \n", COV(0,0), COV(0,1), COV(0,2));
            fprintf (stderr, "  %+.5e, %+.5e, %+.5e  \n", COV(1,0), COV(1,1), COV(1,2));
            fprintf (stderr, "  %+.5e, %+.5e, %+.5e ]\n", COV(2,0), COV(2,1), COV(2,2));
        }
        Log.LogDebug("# chisq/n = %g", chisq/n);
        Log.LogDebug("# zcorr = %g", zcorr);
        Log.LogDebug("# sigma = %g", sigma);
        //Log.LogDebug("# logarea = %g", logarea);
        Log.LogDebug("\n");

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
