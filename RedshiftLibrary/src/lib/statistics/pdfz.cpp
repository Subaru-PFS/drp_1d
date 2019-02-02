#include <RedshiftLibrary/statistics/pdfz.h>

#include <RedshiftLibrary/log/log.h>

using namespace NSEpic;
using namespace std;
#include <fstream>
#include <RedshiftLibrary/extremum/extremum.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlin.h>

//** gaussian fit **//
struct pdfz_lmfitdata {
    size_t n;
    Float64 *y;
    Float64 *z;
    Float64 zcenter;
};

int pdfz_lmfit_f(const gsl_vector *x, void *data, gsl_vector *f)
{
    Int32 verbose = 0;
    size_t n = ((struct pdfz_lmfitdata *)data)->n;
    Float64 *y = ((struct pdfz_lmfitdata *)data)->y;
    Float64 *z = ((struct pdfz_lmfitdata *)data)->z;
    Float64 zcenter = ((struct pdfz_lmfitdata *)data)->zcenter;

    double a = gsl_vector_get(x, 0);
    double sigma = gsl_vector_get(x, 1);
    if (verbose)
    {
        Log.LogInfo("Pdfz: Pdfz computation: pdfz_lmfit_f : a=%e, sigma=%e, "
                    "zcenter=%.5f ",
                    a, sigma, zcenter);
    }

    for (Int32 i = 0; i < n; i++)
    {
        Float64 t = z[i] - zcenter;
        const Float64 xsurc = t / sigma;
        Float64 Yi = a * exp(-0.5 * xsurc * xsurc);
        if (0 && verbose)
        {
            Log.LogInfo(
                "Pdfz: Pdfz computation: pdfz_lmfit_f : for i=%d, Yi=%e", i,
                Yi);
        }
        gsl_vector_set(f, i, Yi - y[i]);
    }

    return GSL_SUCCESS;
}

int pdfz_lmfit_df(const gsl_vector *x, void *data, gsl_matrix *J)
{
    size_t n = ((struct pdfz_lmfitdata *)data)->n;
    Float64 *z = ((struct pdfz_lmfitdata *)data)->z;
    Float64 zcenter = ((struct pdfz_lmfitdata *)data)->zcenter;

    double a = gsl_vector_get(x, 0);
    double sigma = gsl_vector_get(x, 1);

    size_t i;

    for (i = 0; i < n; i++)
    {
        Float64 t = z[i] - zcenter;
        const Float64 xsurc = t / sigma;

        gsl_matrix_set(J, i, 0, exp(-0.5 * xsurc * xsurc));
        gsl_matrix_set(J, i, 1,
                       t * t / (sigma * sigma * sigma) * a *
                           exp(-0.5 * xsurc * xsurc));
    }
    return GSL_SUCCESS;
}
//** gaussian fit end**//

CPdfz::CPdfz() {}

CPdfz::~CPdfz() {}

/**
 * @brief CPdfz::Compute
 * @param merits
 * @param redshifts
 * ...
 * NB-2018-02-19 : this method works with IRREGULAR z grid. No need to have a
 * regular grid z-pdf anymore
 * @return 0: success, 1:problem, 3 not enough z values, 4: zPrior not valid
 */
Int32 CPdfz::Compute(TFloat64List merits, TFloat64List redshifts,
                     Float64 cstLog, TFloat64List logZPrior,
                     TFloat64List &logPdf, Float64 &logEvidence)
{
    Bool verbose = true;
    Int32 sumMethod = 1; // 0=rect, 1=trapez
    logPdf.clear();

    if (verbose)
    {
        Float64 meritmax = -DBL_MAX;
        Float64 meritmin = DBL_MAX;
        for (UInt32 k = 0; k < redshifts.size(); k++)
        {
            if (meritmax < merits[k])
            {
                meritmax = merits[k];
            }
            if (meritmin > merits[k])
            {
                meritmin = merits[k];
            }
        }
        Log.LogDetail("Pdfz: Pdfz computation: using merit min=%e", meritmin);
        Log.LogDetail("Pdfz: Pdfz computation: using merit max=%e", meritmax);
    }

    //
    //    //check if the z step is constant. If not, pdf cannot be estimated by
    //    the current method. Float64 reldzThreshold = 0.05; //relative
    //    difference accepted bool constantdz = true; Float64 mindz = DBL_MAX;
    //    Float64 maxdz = -DBL_MAX;
    //    for ( UInt32 k=1; k<redshifts.size(); k++)
    //    {
    //        Float64 diff = redshifts[k]-redshifts[k-1];
    //        if(mindz > diff)
    //        {
    //            mindz = diff;
    //        }
    //        if(maxdz < diff)
    //        {
    //            maxdz = diff;
    //        }
    //    }
    //    Float64 zstep_mean = (maxdz+mindz)/2.0;
    //    if(abs(maxdz-mindz)/zstep_mean>reldzThreshold)
    //    {
    //        constantdz = false;
    //        return 2;
    //    }
    //

    // deactivate logPrior just to see...
    //    for ( UInt32 k=0; k<redshifts.size(); k++)
    //    {
    //        logZPrior[k] = 0;
    //    }

    // check if there is at least 1 redshifts values
    if (redshifts.size() == 1) // consider this as a success
    {
        logPdf.resize(redshifts.size());
        logPdf[0] = 1.0;
        logEvidence = 1.0;
        return 0;
    } else if (redshifts.size() < 1)
    {
        return 2;
    }

    // check that the zPrior is size-compatible
    if (logZPrior.size() != redshifts.size())
    {
        return 4;
        // Float64 logPrior = log(1.0/redshifts.size()); //log(1.0)
    }
    if (verbose)
    {
        Float64 logZPriorMax = -DBL_MAX;
        Float64 logZPriorMin = DBL_MAX;
        for (UInt32 k = 0; k < redshifts.size(); k++)
        {
            if (logZPriorMax < logZPrior[k])
            {
                logZPriorMax = logZPrior[k];
            }
            if (logZPriorMin > logZPrior[k])
            {
                logZPriorMin = logZPrior[k];
            }
        }
        Log.LogDetail("Pdfz: Pdfz computation: using logZPrior min=%e",
                      logZPriorMax);
        Log.LogDetail("Pdfz: Pdfz computation: using logZPrior max=%e",
                      logZPriorMin);
    }
    //    std::vector<Float64> logZPrior(zPrior.size(), 1.0);
    //    for(UInt32 kz=0; kz<zPrior.size(); kz++)
    //    {
    //        logZPrior[kz] = log(zPrior[kz]);
    //    }

    logPdf.resize(redshifts.size());

    /* ------------------------------------------------------------------
     * NOTE: this uses the LOG-SUM-EXP trick originally suggested by S. Jamal
     * ------------------------------------------------------------------  */

    // prepare logLikelihood and LogEvidence
    Float64 maxi = -DBL_MAX;
    std::vector<Float64> mchi2Sur2(redshifts.size(), 0.0);
    std::vector<Float64> smallVALUES(redshifts.size(), 0.0);
    for (UInt32 k = 0; k < redshifts.size(); k++)
    {
        mchi2Sur2[k] = -0.5 * merits[k];
        Float64 zstep;
        // here, using rect. approx. (enough precision for maxi estimation)
        if (k == 0)
        {
            zstep = (redshifts[k + 1] - redshifts[k]) * 0.5;
        } else if (k == redshifts.size() - 1)
        {
            zstep = (redshifts[k] - redshifts[k - 1]) * 0.5;
        } else
        {
            zstep = (redshifts[k + 1] + redshifts[k]) * 0.5 -
                    (redshifts[k] + redshifts[k - 1]) * 0.5;
        }
        smallVALUES[k] = mchi2Sur2[k] + logZPrior[k];
        if (maxi < smallVALUES[k] + log(zstep))
        {
            maxi = smallVALUES[k] +
                   log(zstep); // maxi will be used to avoid underflows when
                               // summing exponential of small values
        }
    }
    Log.LogDebug(
        "Pdfz: Pdfz computation: using maxi value for log-sum-exp trick=%e",
        maxi);

    Float64 sumModifiedExp = 0.0;
    if (sumMethod == 0)
    {
        Log.LogDebug(
            "Pdfz: Pdfz computation: summation method option = RECTANGLES");
        for (UInt32 k = 0; k < redshifts.size(); k++)
        {
            Float64 modifiedEXPO = exp(smallVALUES[k] - maxi);
            Float64 area = (modifiedEXPO);
            Float64 zstep;
            if (k == 0)
            {
                zstep = (redshifts[k + 1] - redshifts[k]) * 0.5;
            } else if (k == redshifts.size() - 1)
            {
                zstep = (redshifts[k] - redshifts[k - 1]) * 0.5;
            } else
            {
                zstep = (redshifts[k + 1] + redshifts[k]) * 0.5 -
                        (redshifts[k] + redshifts[k - 1]) * 0.5;
            }
            area *= zstep;
            sumModifiedExp += area;
        }
    } else if (sumMethod == 1)
    {
        Log.LogDebug(
            "Pdfz: Pdfz computation: summation method option = TRAPEZOID");
        Float64 modifiedEXPO_previous = exp(smallVALUES[0] - maxi);
        for (UInt32 k = 1; k < redshifts.size(); k++)
        {
            Float64 modifiedEXPO = exp(smallVALUES[k] - maxi);
            Float64 trapezArea = (modifiedEXPO + modifiedEXPO_previous) / 2.0;
            trapezArea *= (redshifts[k] - redshifts[k - 1]);
            sumModifiedExp += trapezArea;
            modifiedEXPO_previous = modifiedEXPO;
        }
    } else
    {
        Log.LogError(
            "Pdfz: Pdfz computation: unable to parse summation method option");
    }
    logEvidence = cstLog + maxi + log(sumModifiedExp);

    if (verbose)
    {
        Log.LogDetail("Pdfz: Pdfz computation: using cstLog=%e", cstLog);
        Log.LogDetail("Pdfz: Pdfz computation: using logEvidence=%e",
                      logEvidence);
        // Log.LogInfo("Pdfz: Pdfz computation: using logPrior=%f",  logPrior);
        // //logPrior can be variable with z
    }

    for (UInt32 k = 0; k < redshifts.size(); k++)
    {
        logPdf[k] = logZPrior[k] + (mchi2Sur2[k] + cstLog) - logEvidence;
    }

    if (verbose)
    {
        Float64 pdfmax = -DBL_MAX;
        Float64 pdfmin = DBL_MAX;
        for (UInt32 k = 0; k < redshifts.size(); k++)
        {
            if (pdfmax < logPdf[k])
            {
                pdfmax = logPdf[k];
            }
            if (pdfmin > logPdf[k])
            {
                pdfmin = logPdf[k];
            }
        }
        Log.LogDetail("Pdfz: Pdfz computation: found pdf min=%e", pdfmin);
        Log.LogDetail("Pdfz: Pdfz computation: found pdf max=%e", pdfmax);
    }

    return 0;
}

Float64 CPdfz::getSumTrapez(std::vector<Float64> redshifts,
                            std::vector<Float64> valprobalog)
{
    Float64 sum = 0.0;
    if(redshifts.size()==0)
    {
        return sum;
    }
    if(redshifts.size()!=valprobalog.size()) //this should raise an exception ? or return some error values ?
    {
        return sum;
    }

    // prepare LogEvidence
    Float64 maxi = -DBL_MAX;
    std::vector<Float64> smallVALUES(redshifts.size(), 0.0);
    for (UInt32 k = 0; k < redshifts.size(); k++)
    {
        // find the smallest zstep in order to use most penalizing case fot the
        // log-sum-exp trick
        Float64 zstepPrevious = -1.;
        if (k > 0)
        {
            zstepPrevious = (redshifts[k] - redshifts[k - 1]);
        } else if (k < redshifts.size() - 1)
        {
            zstepPrevious = (redshifts[k + 1] - redshifts[k]);
        }
        Float64 zstepNext = -1.;
        if (k < redshifts.size() - 1)
        {
            zstepNext = (redshifts[k + 1] - redshifts[k]);
        } else if (k > 0)
        {
            zstepNext = (redshifts[k] - redshifts[k - 1]);
        }
        Float64 zstepCurrent = min(zstepPrevious, zstepNext);
        smallVALUES[k] = valprobalog[k];
        if (maxi < smallVALUES[k] + log(zstepCurrent))
        {
            maxi =
                smallVALUES[k] +
                log(zstepCurrent); // maxi will be used to avoid underflows when
                                   // summing exponential of small values
        }
    }

    Float64 sumModifiedExp = 0.0;
    Float64 modifiedEXPO_previous = exp(smallVALUES[0] - maxi);
    for (UInt32 k = 1; k < redshifts.size(); k++)
    {
        Float64 modifiedEXPO = exp(smallVALUES[k] - maxi);
        Float64 trapezArea = (modifiedEXPO + modifiedEXPO_previous) / 2.0;
        trapezArea *= (redshifts[k] - redshifts[k - 1]);
        sumModifiedExp += trapezArea;
        modifiedEXPO_previous = modifiedEXPO;
    }
    Float64 logSum = maxi + log(sumModifiedExp);

    sum = exp(logSum);

    return sum;
}

Float64 CPdfz::getSumRect(std::vector<Float64> redshifts,
                          std::vector<Float64> valprobalog)
{
    Float64 sum = 0.0;

    // prepare LogEvidence
    Float64 maxi = -DBL_MAX;
    std::vector<Float64> smallVALUES(redshifts.size(), 0.0);
    for (UInt32 k = 0; k < redshifts.size(); k++)
    {
        Float64 zstep;
        if (k == 0)
        {
            zstep = (redshifts[k + 1] - redshifts[k]) * 0.5;
        } else if (k == redshifts.size() - 1)
        {
            zstep = (redshifts[k] - redshifts[k - 1]) * 0.5;
        } else
        {
            zstep = (redshifts[k + 1] + redshifts[k]) * 0.5 -
                    (redshifts[k] + redshifts[k - 1]) * 0.5;
        }
        smallVALUES[k] = valprobalog[k];
        if (maxi < smallVALUES[k] + log(zstep))
        {
            maxi = smallVALUES[k] +
                   log(zstep); // maxi will be used to avoid underflows when
                               // summing exponential of small values
        }
    }

    Log.LogDebug("  pdfz: getSumRect - found maxi = %e", maxi);
    Float64 sumModifiedExp = 0.0;
    for (UInt32 k = 0; k < redshifts.size() - 1; k++)
    {
        Float64 modifiedEXPO = exp(smallVALUES[k] - maxi);
        Float64 area = modifiedEXPO;
        Float64 zstep;
        if (k == 0)
        {
            zstep = (redshifts[k + 1] - redshifts[k]) * 0.5;
        } else if (k == redshifts.size() - 1)
        {
            zstep = (redshifts[k] - redshifts[k - 1]) * 0.5;
        } else
        {
            zstep = (redshifts[k + 1] + redshifts[k]) * 0.5 -
                    (redshifts[k] + redshifts[k - 1]) * 0.5;
        }
        area *= zstep;
        sumModifiedExp += area;
    }
    Float64 logSum = maxi + log(sumModifiedExp);

    sum = exp(logSum);

    return sum;
}

/**
 * @brief CPdfz::getCandidateSumTrapez
 * @param redshifts
 * @param valprobalog
 * @param zcandidate
 * @param zwidth
 * @return -1 if error, else sum around the candidate
 */
Float64 CPdfz::getCandidateSumTrapez(std::vector<Float64> redshifts,
                                     std::vector<Float64> valprobalog,
                                     Float64 zcandidate, Float64 zwidth)
{
    // check that redshifts are sorted
    for (UInt32 k = 1; k < redshifts.size(); k++)
    {
        if (redshifts[k] < redshifts[k - 1])
        {
            Log.LogError("    CPdfz::getCandidateSumTrapez - redshifts are not "
                         "sorted for (at least) index={}",
                         k);
            return -1.0;
        }
    }

    // find indexes kmin, kmax so that zmin and zmax are inside [
    // redshifts[kmin]:redshifts[kmax] ]
    Int32 kmin = 0;
    Int32 kmax = redshifts.size() - 1;
    Float64 halfzwidth = zwidth / 2.0;
    for (UInt32 k = 0; k < redshifts.size(); k++)
    {
        if (redshifts[k] < (zcandidate - halfzwidth))
        {
            kmin = k;
        }
    }
    Log.LogDebug("    CPdfz::getCandidateSumTrapez - kmin index=%d", kmin);

    for (UInt32 k = redshifts.size() - 1; k > 0; k--)
    {
        if (redshifts[k] > (zcandidate + halfzwidth))
        {
            kmax = k;
        }
    }
    Log.LogDebug("    CPdfz::getCandidateSumTrapez - kmax index=%d", kmax);

    // initialize the LOG-SUM-EXP trick
    Float64 maxi = -DBL_MAX;
    std::vector<Float64> smallVALUES(redshifts.size(), 0.0);
    for (UInt32 k = kmin; k <= kmax; k++)
    {
        // estimate zstep for the log-sum-exp trick
        Float64 zstep;
        if (k == 0)
        {
            zstep = (redshifts[k + 1] - redshifts[k]);
        } else if (k == redshifts.size() - 1)
        {
            zstep = (redshifts[k] - redshifts[k - 1]);
        } else
        {
            zstep = (redshifts[k + 1] + redshifts[k]) * 0.5 -
                    (redshifts[k] + redshifts[k - 1]) * 0.5;
        }
        smallVALUES[k] = valprobalog[k];
        if (maxi < smallVALUES[k] + log(zstep))
        {
            maxi = smallVALUES[k] +
                   log(zstep); // maxi will be used to avoid underflows when
                               // summing exponential of small values
        }
    }

    // for now the sum is estimated between kmin and kmax.
    // todo: INTERPOLATE (linear) in order to start exactly at zmin and stop at
    // zmax
    Float64 sum = 0.0;
    Float64 sumModifiedExp = 0.0;
    Float64 modifiedEXPO_previous = exp(smallVALUES[kmin] - maxi);
    for (UInt32 k = kmin + 1; k <= kmax; k++)
    {
        Float64 modifiedEXPO = exp(smallVALUES[k] - maxi);
        Float64 trapezArea = (modifiedEXPO + modifiedEXPO_previous) / 2.0;
        trapezArea *= (redshifts[k] - redshifts[k - 1]);
        sumModifiedExp += trapezArea;
        modifiedEXPO_previous = modifiedEXPO;
    }
    Float64 logSum = maxi + log(sumModifiedExp);

    sum = exp(logSum);

    return sum;
}

Int32 CPdfz::getCandidateRobustGaussFit(std::vector<Float64> redshifts,
                                        std::vector<Float64> valprobalog,
                                        Float64 zcandidate, Float64 zwidth,
                                        Float64 &gaussAmp, Float64 &gaussAmpErr,
                                        Float64 &gaussSigma,
                                        Float64 &gaussSigmaErr)
{
    Int32 fitSuccessful = false;
    Int32 nTry = 5;
    Int32 iTry = 0;

    Float64 current_zwidth = zwidth;
    while (!fitSuccessful && iTry < nTry)
    {
        Int32 retFit = getCandidateGaussFit(
            redshifts, valprobalog, zcandidate, current_zwidth, gaussAmp,
            gaussAmpErr, gaussSigma, gaussSigmaErr);
        if (!retFit && gaussSigma < zwidth * 2.0 &&
            std::abs(gaussSigma / gaussSigmaErr) > 1e-2)
        {
            fitSuccessful = true;
        } else
        {
            Log.LogDebug("    CPdfz::getCandidateRobustGaussFit - iTry=%d",
                         iTry);
            Log.LogDebug("    CPdfz::getCandidateRobustGaussFit -    for "
                         "zcandidate=%.5f",
                         zcandidate);
            Log.LogDebug("    CPdfz::getCandidateRobustGaussFit -       found "
                         "gaussAmp=%e",
                         gaussAmp);
            Log.LogDebug("    CPdfz::getCandidateRobustGaussFit -       found "
                         "gaussSigma=%e",
                         gaussSigma);
            Log.LogDebug("    CPdfz::getCandidateRobustGaussFit -       now "
                         "going to retry w. different parameters",
                         gaussSigma);
        }
        current_zwidth /= 2.0;
        iTry++;
    }

    if (fitSuccessful)
    {
        return 0;
    } else
    {
        return -1;
    }
}

Int32 CPdfz::getCandidateGaussFit(std::vector<Float64> redshifts,
                                  std::vector<Float64> valprobalog,
                                  Float64 zcandidate, Float64 zwidth,
                                  Float64 &gaussAmp, Float64 &gaussAmpErr,
                                  Float64 &gaussSigma, Float64 &gaussSigmaErr)
{
    Int32 verbose = 0;
    Int32 ret = 0;
    Log.LogDebug("    CPdfz::getCandidateSumGaussFit - Starting pdf peaks "
                 "gaussian fitting");

    // check that redshifts are sorted
    for (UInt32 k = 1; k < redshifts.size(); k++)
    {
        if (redshifts[k] < redshifts[k - 1])
        {
            Log.LogError("    CPdfz::getCandidateSumGaussFit - redshifts are "
                         "not sorted for (at least) index={}",
                         k);
            return -1.0;
        }
    }

    // find indexes kmin, kmax so that zmin and zmax are inside [
    // redshifts[kmin]:redshifts[kmax] ]
    Int32 kmin = 0;
    Int32 kmax = redshifts.size() - 1;
    Float64 halfzwidth = zwidth / 2.0;
    for (UInt32 k = 0; k < redshifts.size(); k++)
    {
        if (redshifts[k] < (zcandidate - halfzwidth))
        {
            kmin = k;
        }
    }
    if (verbose)
    {
        Log.LogInfo("    CPdfz::getCandidateSumGaussFit - kmin index=%d", kmin);
    }

    for (UInt32 k = redshifts.size() - 1; k > 0; k--)
    {
        if (redshifts[k] > (zcandidate + halfzwidth))
        {
            kmax = k;
        }
    }
    if (verbose)
    {
        Log.LogInfo("    CPdfz::getCandidateSumGaussFit - kmax index=%d", kmax);
    }

    // initialize GSL
    const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
    gsl_multifit_fdfsolver *s;
    int status, info;
    size_t i;
    size_t n =
        kmax - kmin +
        2; // n samples on the support, /* number of data points to fit */
    size_t p = 2; // DOF = 1.amplitude + 2.width

    if (verbose)
    {
        Log.LogInfo("    CPdfz::getCandidateSumGaussFit - n=%d, p=%d", n, p);
    }
    if (n < p)
    {
        Log.LogError("    CPdfz::getCandidateSumGaussFit - LMfit not enough "
                     "samples on support");
        return -1;
    }

    gsl_matrix *J = gsl_matrix_alloc(n, p);
    gsl_matrix *covar = gsl_matrix_alloc(p, p);
    double y[n], weights[n], z[n];
    Float64 zc = zcandidate;
    gsl_multifit_function_fdf f;

    Float64 *x_init = (Float64 *)calloc(p, sizeof(Float64));
    if (x_init == 0)
    {
        Log.LogError(
            "    CPdfz::getCandidateSumGaussFit - Unable to allocate x_init");
        return -1;
    }

    // initialize lmfit with previously estimated values ?
    //
    Float64 maxP = 0.0;
    for (i = 0; i < n; i++)
    {
        Float64 idx = i + kmin;
        Float64 _p = exp(valprobalog[idx]);
        if (_p > maxP)
        {
            maxP = _p;
        }
    }
    Float64 normFactor = maxP;
    if (normFactor <= 0.)
    {
        normFactor = 1.0;
    }
    if (maxP > 0.0)
    {
        x_init[0] = maxP / normFactor;
    } else
    {
        x_init[0] = 1.0;
    }
    x_init[1] = zwidth / 2.0;
    if (verbose)
    {
        Log.LogError("    CPdfz::getCandidateSumGaussFit - init a=%e",
                     x_init[0]);
        Log.LogError("    CPdfz::getCandidateSumGaussFit - init sigma=%e",
                     x_init[1]);
    }

    gsl_vector_view x = gsl_vector_view_array(x_init, p);
    //    if(x.vector==0){
    //        Log.LogError( "    CPdfz::getCandidateSumGaussFit - Unable to
    //        allocate x"); return -1;
    //    }
    gsl_vector_view w = gsl_vector_view_array(weights, n);
    //    if(w.vector==0){
    //        Log.LogError( "    CPdfz::getCandidateSumGaussFit - Unable to
    //        allocate w"); return -1;
    //    }
    gsl_vector *res_f;
    double chi, chi0;

    const double xtol = 1e-8;
    const double gtol = 1e-8;
    const double ftol = 1 - 8;
    Int32 maxIterations = 500;

    // This is the data to be fitted;
    for (i = 0; i < n; i++)
    {
        Float64 idx = i + kmin;
        weights[i] = 1.0; // no weights
        y[i] = exp(valprobalog[idx]) / normFactor;
        z[i] = redshifts[idx];
    }

    struct pdfz_lmfitdata d = {n, y, z, zc};
    f.f = &pdfz_lmfit_f;
    f.df = &pdfz_lmfit_df;
    f.n = n;
    f.p = p;
    f.params = &d;

    Log.LogDebug("    CPdfz::getCandidateSumGaussFit - LMfit data ready");

    s = gsl_multifit_fdfsolver_alloc(T, n, p);
    if (s == 0)
    {
        Log.LogError("    CPdfz::getCandidateSumGaussFit - Unable to allocate "
                     "the multifit solver s");
        return -1;
    }

    /* initialize solver with starting point and weights */
    gsl_multifit_fdfsolver_wset(s, &f, &x.vector, &w.vector);

    /* compute initial residual norm */
    res_f = gsl_multifit_fdfsolver_residual(s);
    chi0 = gsl_blas_dnrm2(res_f);

    /* solve the system with a maximum of maxIterations iterations */
    status = gsl_multifit_fdfsolver_driver(s, maxIterations, xtol, gtol, ftol,
                                           &info);

    gsl_multifit_fdfsolver_jac(s, J);
    gsl_multifit_covar(J, 0.0, covar);

    /* compute final residual norm */
    chi = gsl_blas_dnrm2(res_f);

    double dof = n - p;
    double c = GSL_MAX_DBL(1, chi / sqrt(dof));
    if (verbose)
    {
#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar, i, i))

        Log.LogInfo("summary from method '%s'", gsl_multifit_fdfsolver_name(s));
        Log.LogInfo("number of iterations: %zu",
                    gsl_multifit_fdfsolver_niter(s));
        Log.LogInfo("function evaluations: %zu", f.nevalf);
        Log.LogInfo("Jacobian evaluations: %zu", f.nevaldf);
        Log.LogInfo("reason for stopping: %s",
                    (info == 1)
                        ? "small step size"
                        : (info == 2) ? "small gradient" : "small change in f");
        Log.LogInfo("initial |f(x)| = %g", chi0);
        Log.LogInfo("final   |f(x)| = %g", chi);

        {

            Log.LogInfo("chisq/dof = %g", pow(chi, 2.0) / dof);

            for (Int32 k = 0; k < p; k++)
            {
                if (FIT(k) < 1e-3)
                {
                    Log.LogInfo("A %d     = %.3e +/- %.8f", k, FIT(k),
                                c * ERR(k));
                } else
                {
                    Log.LogInfo("A %d     = %.5f +/- %.8f", k, FIT(k),
                                c * ERR(k));
                }
            }
        }
        Log.LogInfo("status = %s (%d)", gsl_strerror(status), status);
    }

    gaussAmp = gsl_vector_get(s->x, 0) * normFactor;
    gaussAmpErr = c * ERR(0) * normFactor;
    gaussSigma = abs(gsl_vector_get(s->x, 1));
    gaussSigmaErr = c * ERR(1);

    return ret;
}

Int32   CPdfz::getPmis(std::vector<Float64> redshifts,
                       std::vector<Float64> valprobalog,
                       Float64 zbest,
                       std::vector<Float64> zcandidates,
                       Float64 zwidth,
                       Float64 &pmis)
{
    Float64 halfzwidth = zwidth/2.;
    pmis = -1;

    //definition of the lines for mismatch
    std::vector<Float64> lambda_mis;
    lambda_mis.push_back(9533.2); //SIII9530
    lambda_mis.push_back(9071.1); //SIII9068

    lambda_mis.push_back(7753.2); //ArIII7751

    lambda_mis.push_back(7332.2); //OII7330
    lambda_mis.push_back(7321); //OII7319

    lambda_mis.push_back(6564.61); //halpha
    lambda_mis.push_back(3729.88); //[OII]3729
    lambda_mis.push_back(3727.09); //[OII]3726
    lambda_mis.push_back(5008.24); //[OIII](doublet-1)
    lambda_mis.push_back(4960.29); //[OIII](doublet-1/3)
    lambda_mis.push_back(4862.72); //Hbeta


    //definition of relzerr_mis
    std::vector<Float64> relzerr_mis;
    relzerr_mis.clear();
    for (Int32 kl1 = 0; kl1<lambda_mis.size(); kl1++)
    {
        for (Int32 kl2 = 0; kl2<lambda_mis.size(); kl2++)
        {
            if(kl1 != kl2)
            {
                Float64 relzerr = (lambda_mis[kl1]/lambda_mis[kl2])-1.;
                relzerr_mis.push_back(relzerr);
                Log.LogDetail( "pdfz: Found l1=%f, l2=%f, relzerr=%f", lambda_mis[kl1], lambda_mis[kl2], relzerr);
            }
        }
    }

    Int32 n_relzerr_mis = relzerr_mis.size();
    Log.LogDetail( "pdfz: Found n relzerr for mismatch =%d", n_relzerr_mis);

    //override amazed-zcandidates by all pdf peaks found
    Int32 maxpeakscount = 1000;
    TFloat64Range redshiftsRange(
        redshifts[0],
        redshifts[redshifts.size() - 1]);
    CExtremum extremum(redshiftsRange, maxpeakscount, false, 2);
    TPointList extremumList;
    extremum.Find(redshifts, valprobalog, extremumList);
    zcandidates.clear();
    for (Int32 i = 0; i < extremumList.size(); i++)
    {
        Float64 x = extremumList[i].X;
        zcandidates.push_back(x);
    }
    Log.LogDetail( "pdfz: Found n candidates for mismatch calc.=%d", zcandidates.size());


    std::vector<Int32> indexes_candidates_selected;
    //for(Int32 k_zmap=0; k_zmap<zcandidates.size(); k_zmap++)
    {
        //Int32 zmap = zcandidates[k_zmap];
        Float64 zmap = zbest;

        //find zcandidates corresponding to a relzerr value wrt to zmap
        for(Int32 k_cand=0; k_cand<zcandidates.size(); k_cand++)
        {
            for(Int32 k_relzerr_mis=0; k_relzerr_mis<n_relzerr_mis; k_relzerr_mis++)
            {
                if(true)
                    //if(relzerr_mis[k_relzerr_mis] > halfzwidth*2.0) //avoid zmap candidates selection
                {
                    //convert relzerr in z value for zmap
                    Float64 z_mis = relzerr_mis[k_relzerr_mis]*(1+zmap)+zmap;
                    //Log.LogDetail( "pdfz: Found zmis=%f for relzerr=%f", z_mis, relzerr_mis[k_relzerr_mis]);

                    if(z_mis>(zcandidates[k_cand]-halfzwidth) && z_mis<(zcandidates[k_cand]+halfzwidth))
                    {
                        indexes_candidates_selected.push_back(k_cand);
                        Log.LogDetail( "pdfz: Found mismatch for relzerr=%f, for zmap=%f: candidate at z=%f", relzerr_mis[k_relzerr_mis], zmap, zcandidates[k_cand]);
                    }
                }
            }
        }
    }
    Log.LogDetail( "pdfz: Found n indexes_candidates_selected=%d", indexes_candidates_selected.size());

    //find redshift indexes corresponding to selected candidates surrounding
    std::vector<Int32> indexes_zsamples_selected;
    Int32 n_indexes_candidates_selected = indexes_candidates_selected.size();//min((Int32)1, (Int32)indexes_candidates_selected.size());
    for(Int32 k_cand_sel=0; k_cand_sel<n_indexes_candidates_selected; k_cand_sel++)
    {
        Int32 k_cand = indexes_candidates_selected[k_cand_sel];

        //find zcandidate index
        for(Int32 kz=0; kz<redshifts.size(); kz++)
        {
            if(redshifts[kz]>(zcandidates[k_cand]-halfzwidth) && redshifts[kz]<(zcandidates[k_cand]+halfzwidth))
            {
                indexes_zsamples_selected.push_back(kz);
            }
        }
    }

    Log.LogDetail( "pdfz: Found n NON-UNIQUE indexes_zsamples_selected=%d", indexes_zsamples_selected.size());
    //remove duplicate indexes
    std::sort(indexes_zsamples_selected.begin(), indexes_zsamples_selected.end());
    indexes_zsamples_selected.erase( std::unique( indexes_zsamples_selected.begin(), indexes_zsamples_selected.end() ), indexes_zsamples_selected.end() );

    Log.LogDetail( "pdfz: Found n UNIQUE indexes_zsamples_selected=%d", indexes_zsamples_selected.size());


    //now integrate over the selected unique z indexes
    //warning, integrating exp. values directly. supposing that no overflow will occur if pdf is normalized
    Float64 pmis_raw = 0;
    std::vector<Float64> valprobalog_selected;
    std::vector<Float64> redshifts_selected;
    for(Int32 kz=0; kz<indexes_zsamples_selected.size(); kz++)
    {
        redshifts_selected.push_back(redshifts[indexes_zsamples_selected[kz]]);
        valprobalog_selected.push_back(valprobalog[indexes_zsamples_selected[kz]]);
    }
    pmis_raw = getSumTrapez(redshifts_selected, valprobalog_selected);

    //estimate zcalc intg proba
    Float64 pzcalc = getCandidateSumTrapez( redshifts, valprobalog, zbest, zwidth);

    Log.LogInfo("pdfz: <pmisraw><%.6e>", pmis_raw);
    Log.LogInfo("pdfz: <pmap><%.6e>", pzcalc);

    pmis = pmis_raw/(1.-pzcalc);
    Log.LogInfo("pdfz: <pmis><%.6e>", pmis);

    return 0;
}

std::vector<Float64> CPdfz::GetConstantLogZPrior(UInt32 nredshifts)
{
    std::vector<Float64> zPrior(nredshifts, 1.0);
    for (UInt32 kz = 0; kz < nredshifts; kz++)
    {
        zPrior[kz] = 1.0 / nredshifts;
    }

    // switch to log
    std::vector<Float64> logzPrior(nredshifts, 0.0);
    for (UInt32 kz = 0; kz < nredshifts; kz++)
    {
        logzPrior[kz] = log(zPrior[kz]);
    }

    return logzPrior;
}

std::vector<Float64>
CPdfz::GetStrongLinePresenceLogZPrior(std::vector<bool> linePresence,
                                      Float64 penalization_factor)
{
    Float64 probaPresent = 1.0;
    Float64 probaAbsent = penalization_factor;
    std::vector<Float64> zPrior(linePresence.size(), probaAbsent);
    Float64 sum = 0.0;
    for (UInt32 kz = 0; kz < linePresence.size(); kz++)
    {
        if (linePresence[kz])
        {
            zPrior[kz] = probaPresent;
        } else
        {
            zPrior[kz] = probaAbsent;
        }
        sum += zPrior[kz];
    }

    if (sum > 0)
    {
        for (UInt32 kz = 0; kz < linePresence.size(); kz++)
        {
            zPrior[kz] /= sum;
        }
    }

    // switch to log
    std::vector<Float64> logzPrior(linePresence.size(), probaAbsent);
    for (UInt32 kz = 0; kz < linePresence.size(); kz++)
    {
        logzPrior[kz] = log(zPrior[kz]);
    }

    return logzPrior;
}

std::vector<Float64>
CPdfz::GetNLinesSNRAboveCutLogZPrior(std::vector<Int32> nlinesAboveSNR,
                                      Float64 penalization_factor)
{
    Int32 nz = nlinesAboveSNR.size();
    Int32 nlinesThres = 2;
    Float64 probaPresent = 1.0;
    Float64 probaAbsent = penalization_factor;
    std::vector<Float64> zPrior(nz, probaAbsent);
    Float64 sum = 0.0;
    for (UInt32 kz = 0; kz < nz; kz++)
    {
        if (nlinesAboveSNR[kz] >= nlinesThres)
        {
            zPrior[kz] = probaPresent;
            //Log.LogDetail("Pdfz: Prior: nlinesAboveSNR[kz] >= nlinesThres for kz=%d", kz);
        } else
        {
            zPrior[kz] = probaAbsent;
        }
        sum += zPrior[kz];
    }

    if (sum > 0)
    {
        for (UInt32 kz = 0; kz < nz; kz++)
        {
            zPrior[kz] /= sum;
        }
    }

    // switch to log
    std::vector<Float64> logzPrior(nz, probaAbsent);
    for (UInt32 kz = 0; kz < nz; kz++)
    {
        logzPrior[kz] = log(zPrior[kz]);
    }

    return logzPrior;
}

std::vector<Float64> CPdfz::GetEuclidNhaLogZPrior(std::vector<Float64> redshifts, Float64 aCoeff)
{
    if(aCoeff<=0)
    {
        Log.LogError( "    CPdfz::GetEuclidNhaLogZPrior: problem found aCoeff<=0: aCoeff=%f", aCoeff);
        throw std::runtime_error("    CPdfz::GetEuclidNhaLogZPrior: problem found aCoeff<=0");
    }
    std::vector<Float64> zPrior(redshifts.size(), 0.0);

    Float64 maxP = -DBL_MAX;
    Float64 minP = DBL_MAX;
    for (UInt32 kz = 0; kz < redshifts.size(); kz++)
    {
        Float64 z = redshifts[kz];
        Float64 z2 = z*z;
        Float64 z3 = z2*z;
        Float64 z4 = z3*z;
        Float64 z5 = z4*z;
        Float64 z6 = z5*z;

        //poly reg pozzetti model 1 at FHa=1e-16
        zPrior[kz] = (- 54.7422088727874*z6
                      + 1203.94994364807*z5
                      - 10409.6716744981*z4
                      + 44240.3837462642*z3
                      - 92914.84430357*z2
                      + 79004.76406*z
                      - 2288.98457865 );

        //shape prior at low z, left of the bell
        Bool enable_low_z_flat = true;
        if(enable_low_z_flat && z<0.7204452872044528){
            zPrior[kz]=20367.877916402278;
        }else{
            if(zPrior[kz]<0){
                zPrior[kz]=DBL_MIN;
            }
        }
        //apply strength
        zPrior[kz] = pow(zPrior[kz], aCoeff);

        if (zPrior[kz] > maxP)
        {
            maxP = zPrior[kz];
        }
        if (zPrior[kz] < minP)
        {
            minP = zPrior[kz];
        }
    }

    Log.LogDetail("Pdfz: zPrior: using HalphaZPrior min=%e", minP);
    Log.LogDetail("Pdfz: zPrior: using HalphaZPrior max=%e", maxP);
    Float64 dynamicCut = 1e12;
    if (maxP > 0)
    {
        for (UInt32 kz = 0; kz < redshifts.size(); kz++)
        {
            zPrior[kz] /= maxP;
            if (zPrior[kz] < 1. / dynamicCut)
            {
                zPrior[kz] = 1. / dynamicCut;
            }
        }
    }

    Float64 sum = 0.0;
    for (UInt32 kz = 0; kz < redshifts.size(); kz++)
    {
        sum += zPrior[kz];
    }
    if (sum > 0)
    {
        for (UInt32 kz = 0; kz < redshifts.size(); kz++)
        {
            zPrior[kz] /= sum;
        }
    }

    // switch to log
    std::vector<Float64> logzPrior(redshifts.size(), 0.0);
    for (UInt32 kz = 0; kz < redshifts.size(); kz++)
    {
        logzPrior[kz] = log(zPrior[kz]);
    }

    return logzPrior;
}

/**
 * @brief CombineLogZPrior
 * returns a vector with log(prior1*prior2) for each z, with normalization so
 * that sumPrior=1
 * @param logprior1
 * @param logprior2
 * @return
 */
std::vector<Float64> CPdfz::CombineLogZPrior(std::vector<Float64> logprior1,
                                             std::vector<Float64> logprior2)
{
    bool normalizePrior=true;
    std::vector<Float64> logzPriorCombined;
    if (logprior1.size() != logprior2.size())
    {
        return logzPriorCombined;
    }
    Int32 n = logprior1.size();

    Float64 maxi = DBL_MIN;
    for (UInt32 k = 0; k < n; k++)
    {
        Float64 val = logprior1[k] + logprior2[k];
        if (maxi < val)
        {
            maxi = val;
        }
    }
    std::vector<Float64> zPrior_modif(n, 0.0);
    for (UInt32 k = 0; k < n; k++)
    {
        zPrior_modif[k] = logprior1[k] + logprior2[k] - maxi;
    }

    Float64 sumExpModif = 0.0;
    for (UInt32 k = 0; k < n; k++)
    {
        sumExpModif += exp(zPrior_modif[k]);
    }

    Float64 logSum = log(exp(maxi)) * sumExpModif;
    Float64 lognormterm = 0.;
    if(normalizePrior)
    {
        lognormterm=logSum;
    }

    logzPriorCombined.resize(n);
    for (Int32 k = 0; k < n; k++)
    {
        logzPriorCombined[k] = logprior1[k] + logprior2[k] - lognormterm;
    }

    return logzPriorCombined;
}

Int32 CPdfz::Marginalize(TFloat64List redshifts,
                         std::vector<TFloat64List> meritResults,
                         std::vector<TFloat64List> zPriors, Float64 cstLog,
                         std::shared_ptr<CPdfMargZLogResult> postmargZResult,
                         std::vector<Float64> modelPriors)
{
    bool verbose = false;

    if (meritResults.size() != zPriors.size())
    {
        Log.LogError("Pdfz: Pdfz marginalize problem. merit.size (%d) != "
                     "prior.size (%d)",
                     meritResults.size(), zPriors.size());
        return -9;
    }
    if (meritResults.size() < 1 || zPriors.size() < 1 || redshifts.size() < 1)
    {
        Log.LogError("Pdfz: Pdfz marginalize problem. merit.size (%d), "
                     "prior.size (%d), or redshifts.size (%d) is zero !",
                     meritResults.size(), zPriors.size(), redshifts.size());
        return -99;
    }

    // check merit curves. Maybe this should be assert stuff ?
    for (Int32 km = 0; km < meritResults.size(); km++)
    {
        TFloat64List _merit = meritResults[km];
        Bool invalidFound = false;
        for (Int32 kz = 0; kz < _merit.size(); kz++)
        {
            if (_merit[kz] != _merit[kz])
            {
                Log.LogError("    CPdfz::Marginalize - merit result #%d has at "
                             "least one nan or invalid value at index=%d",
                             km, kz);
                invalidFound = true;
            }
            if (invalidFound)
            {
                break;
            }
        }
    }

    Bool initPostMarg = false;
    std::vector<UInt32> nSum;

    Float64 MaxiLogEvidence = -DBL_MAX;
    TFloat64List LogEvidencesWPriorM;
    Float64 sumModifiedEvidences = 0;

    std::vector<Float64> logPriorModel;
    if (/*false &&*/ modelPriors.size() != meritResults.size())
    {
        Float64 priorModelCst = 1.0 / ((Float64)meritResults.size());
        Log.LogInfo(
            "Pdfz: Marginalize: no priors loaded, using constant priors (=%f)",
            priorModelCst);
        for (Int32 km = 0; km < meritResults.size(); km++)
        {
            logPriorModel.push_back(log(priorModelCst));
        }
    } else
    {
        /*
        //override modelPriors with pypelid 10 knn templates priors
        logPriorModel.push_back(log(0.1490));
        logPriorModel.push_back(log(0.0794));
        logPriorModel.push_back(log(0.0744));
        logPriorModel.push_back(log(0.0836));
        logPriorModel.push_back(log(0.1089));
        logPriorModel.push_back(log(0.1124));
        logPriorModel.push_back(log(0.0786));
        logPriorModel.push_back(log(0.1563));
        logPriorModel.push_back(log(0.0509));
        logPriorModel.push_back(log(0.1060));
        //*/
        for (Int32 km = 0; km < meritResults.size(); km++)
        {
            logPriorModel.push_back(log(modelPriors[km]));
        }

        // logging priors used
        Float64 sumPriors = 0.;
        for (Int32 km = 0; km < meritResults.size(); km++)
        {
            sumPriors += exp(logPriorModel[km]);
            Log.LogDetail("Pdfz: Marginalize: for model k=%d, using prior=%f", km,
                        exp(logPriorModel[km]));
        }
        Log.LogInfo("Pdfz: Marginalize: sumPriors=%f", sumPriors);
        if (sumPriors > 1.1 || sumPriors < 0.9)
        {
            Log.LogError("Pdfz: sumPriors should be close to 1... !!!");
        }
    }

    for (Int32 km = 0; km < meritResults.size(); km++)
    {
        // Todo: Check if the status is OK ?
        // meritResult->Status[i] == COperator::nStatus_OK

        CPdfz pdfz;
        TFloat64List logProba;
        Float64 logEvidence;
        Int32 retPdfz = pdfz.Compute(meritResults[km], redshifts, cstLog,
                                     zPriors[km], logProba, logEvidence);
        if (retPdfz != 0)
        {
            Log.LogError("Pdfz: Pdfz computation - compute logEvidence: failed "
                         "for result km=%d",
                         km);
            return -1;
        } else
        {
            //            if(verbose)
            //            {
            //                Log.LogInfo("Pdfz: Marginalize: for km=%d,
            //                logEvidence=%e", km, MaxiLogEvidence);
            //            }
            Float64 logEvidenceWPriorM = logEvidence + logPriorModel[km];

            LogEvidencesWPriorM.push_back(logEvidenceWPriorM);
            if (MaxiLogEvidence < logEvidenceWPriorM)
            {
                MaxiLogEvidence = logEvidenceWPriorM;
            }
        }
    }
    if (verbose)
    {
        Log.LogInfo("Pdfz: Marginalize: MaxiLogEvidence=%e", MaxiLogEvidence);
    }

    // Using computational trick to sum the evidences
    for (Int32 k = 0; k < LogEvidencesWPriorM.size(); k++)
    {
        sumModifiedEvidences += exp(LogEvidencesWPriorM[k] - MaxiLogEvidence);
    }
    Float64 logSumEvidence = MaxiLogEvidence + log(sumModifiedEvidences);
    if (verbose)
    {
        Log.LogInfo("Pdfz: Marginalize: logSumEvidence=%e", logSumEvidence);
    }

    for (Int32 km = 0; km < meritResults.size(); km++)
    {
        if (verbose)
        {
            Log.LogInfo("Pdfz: Marginalize: processing chi2-result km=%d", km);
        }

        // Todo: Check if the status is OK ?
        // meritResult->Status[i] == COperator::nStatus_OK

        CPdfz pdfz;
        TFloat64List logProba;
        Float64 logEvidence;
        Int32 retPdfz = pdfz.Compute(meritResults[km], redshifts, cstLog,
                                     zPriors[km], logProba, logEvidence);
        if (retPdfz != 0)
        {
            Log.LogError("Pdfz: Pdfz computation failed for result km=%d", km);
            return -1;
        } else
        {
            if (!initPostMarg)
            {
                nSum.resize(redshifts.size());
                postmargZResult->countTPL =
                    redshifts.size(); // assumed 1 model per z
                postmargZResult->Redshifts.resize(redshifts.size());
                postmargZResult->valProbaLog.resize(redshifts.size());
                for (UInt32 k = 0; k < redshifts.size(); k++)
                {
                    postmargZResult->Redshifts[k] = redshifts[k];
                    postmargZResult->valProbaLog[k] = -INFINITY;
                    nSum[k] = 0;
                }
                initPostMarg = true;
            } else
            {
                // check if the redshift bins are the same
                for (UInt32 k = 0; k < redshifts.size(); k++)
                {
                    if (postmargZResult->Redshifts[k] != redshifts[k])
                    {
                        Log.LogError("pdfz: Pdfz computation (z-bins "
                                     "comparison) failed for result km=%d",
                                     km);
                        break;
                    }
                }
            }

            postmargZResult->valEvidenceLog = logSumEvidence;
            for (UInt32 k = 0; k < redshifts.size(); k++)
            {
                if (true /*meritResult->Status[k]== COperator::nStatus_OK*/) // todo: check (temporarily considers status is always OK for linemodel tplshape)
                {
                    Float64 logValProba = postmargZResult->valProbaLog[k];
                    Float64 logValProbaAdd = logProba[k] + logPriorModel[km] +
                                             logEvidence - logSumEvidence;
                    Float64 maxP = logValProba;
                    if (maxP < logValProbaAdd)
                    {
                        maxP = logValProbaAdd;
                    }
                    Float64 valExp =
                        exp(logValProba - maxP) + exp(logValProbaAdd - maxP);
                    postmargZResult->valProbaLog[k] = maxP + log(valExp);
                    nSum[k]++;
                }
            }
        }
    }

    // THIS DOES NOT ALLOW Marginalization with coverage<100% for ALL templates
    for (UInt32 k = 0; k < postmargZResult->Redshifts.size(); k++)
    {
        if (nSum[k] != meritResults.size())
        {
            postmargZResult->valProbaLog[k] = NAN;
            Log.LogError("Pdfz: Pdfz computation failed. For z=%f, nSum=%d",
                         postmargZResult->Redshifts[k], nSum[k]);
            Log.LogError("Pdfz: Pdfz computation failed. For z=%f, "
                         "meritResults.size()=%d",
                         postmargZResult->Redshifts[k], meritResults.size());
            Log.LogError("Pdfz: Pdfz computation failed. Not all templates "
                         "have 100 percent coverage for all redshifts!");
        }
    }

    return 0;
}

// This mathematically does not correspond to any valid method for combining
// PDFs.
// TODO: problem while estimating best proba. is it best proba for each z ? In
// that case: what about sum_z P = 1 ?
// TODO: this methid should be replaced/modified to correspond to the MaxPDF
// technique.
Int32 CPdfz::BestProba(TFloat64List redshifts,
                       std::vector<TFloat64List> meritResults,
                       std::vector<TFloat64List> zPriors, Float64 cstLog,
                       std::shared_ptr<CPdfMargZLogResult> postmargZResult)
{
    Log.LogError("Pdfz: Pdfz-bestproba computation ! This method is currently "
                 "not working !! It will produce bad results as is....");
    bool verbose = false;

    if (meritResults.size() != zPriors.size())
    {
        Log.LogError(
            "Pdfz: Pdfz-bestproba problem. merit.size (%d) != prior.size (%d)",
            meritResults.size(), zPriors.size());
        return -9;
    }
    if (meritResults.size() < 1 || zPriors.size() < 1 || redshifts.size() < 1)
    {
        Log.LogError("Pdfz: Pdfz-bestproba problem. merit.size (%d), "
                     "prior.size (%d), or redshifts.size (%d) is zero !",
                     meritResults.size(), zPriors.size(), redshifts.size());
        return -99;
    }

    Bool initPostMarg = false;

    for (Int32 km = 0; km < meritResults.size(); km++)
    {
        if (verbose)
        {
            Log.LogInfo("Pdfz:-bestproba: processing chi2-result km=%d", km);
        }

        // Todo: Check if the status is OK ?
        // meritResult->Status[i] == COperator::nStatus_OK

        CPdfz pdfz;
        TFloat64List logProba;
        Float64 logEvidence;
        Int32 retPdfz = pdfz.Compute(meritResults[km], redshifts, cstLog,
                                     zPriors[km], logProba, logEvidence);
        if (retPdfz != 0)
        {
            Log.LogError(
                "Pdfz: Pdfz-bestproba computation failed for result km=%d", km);
            return -1;
        } else
        {
            if (!initPostMarg)
            {
                postmargZResult->countTPL =
                    redshifts.size(); // assumed 1 model per z
                postmargZResult->Redshifts.resize(redshifts.size());
                postmargZResult->valProbaLog.resize(redshifts.size());
                for (UInt32 k = 0; k < redshifts.size(); k++)
                {
                    postmargZResult->Redshifts[k] = redshifts[k];
                    postmargZResult->valProbaLog[k] = -DBL_MAX;
                }
                initPostMarg = true;
            } else
            {
                // check if the redshift bins are the same
                for (UInt32 k = 0; k < redshifts.size(); k++)
                {
                    if (postmargZResult->Redshifts[k] != redshifts[k])
                    {
                        Log.LogError("pdfz: Pdfz-bestproba computation (z-bins "
                                     "comparison) failed for result km=%d",
                                     km);
                        break;
                    }
                }
            }
            for (UInt32 k = 0; k < redshifts.size(); k++)
            {
                if (true /*meritResult->Status[k]== COperator::nStatus_OK*/) // todo: check (temporarily considers status is always OK for linemodel tplshape)
                {
                    if (logProba[k] > postmargZResult->valProbaLog[k])
                    {
                        postmargZResult->valProbaLog[k] = logProba[k];
                    }
                }
            }
        }
    }

    // normalize: sum_z P = 1
    // 1. check if the z step is constant. If not, pdf cannot be estimated by
    // the current method.
    Float64 reldzThreshold = 0.05; // relative difference accepted
    Float64 mindz = DBL_MAX;
    Float64 maxdz = -DBL_MAX;
    for (UInt32 k = 1; k < redshifts.size(); k++)
    {
        Float64 diff = redshifts[k] - redshifts[k - 1];
        if (mindz > diff)
        {
            mindz = diff;
        }
        if (maxdz < diff)
        {
            maxdz = diff;
        }
    }
    Float64 zstep = (maxdz + mindz) / 2.0;
    if (abs(maxdz - mindz) / zstep > reldzThreshold)
    {
        return 2;
    }

    // 2. prepare LogEvidence
    Float64 maxi = -DBL_MAX;
    std::vector<Float64> smallVALUES(redshifts.size(), 0.0);
    for (UInt32 k = 0; k < redshifts.size(); k++)
    {
        smallVALUES[k] = postmargZResult->valProbaLog[k];
        if (maxi < smallVALUES[k])
        {
            maxi = smallVALUES[k]; // maxi will be used to avoid underflows when
                                   // summing exponential of small values
        }
    }

    Float64 sumModifiedExp = 0.0;
    for (UInt32 k = 0; k < redshifts.size(); k++)
    {
        Float64 modifiedEXPO = exp(smallVALUES[k] - maxi);
        sumModifiedExp += modifiedEXPO;
    }
    Float64 logEvidence = maxi + log(sumModifiedExp) + log(zstep);

    if (verbose)
    {
        Log.LogInfo("Pdfz: Pdfz-bestproba computation: using logEvidence=%e",
                    logEvidence);
        Log.LogInfo("Pdfz: Pdfz-bestproba computation: using log(zstep)=%e",
                    log(zstep));
    }

    for (UInt32 k = 0; k < redshifts.size(); k++)
    {
        postmargZResult->valProbaLog[k] =
            postmargZResult->valProbaLog[k] - logEvidence;
    }

    return 0;
}

/**
 * @brief CPdfz::BestChi2
 * Computes the combined pdf by taking the best chi2
 * WARNING: as long as the prior on the models are constant, it is equivalent to
 * compute the bestchi2 and the MaxPDF. If this prior is not constant any more,
 * the mas search has to be modified.
 * @param redshifts
 * @param meritResults
 * @param zPriors
 * @param cstLog
 * @param postmargZResult
 * @return
 */
Int32 CPdfz::BestChi2(TFloat64List redshifts,
                      std::vector<TFloat64List> meritResults,
                      std::vector<TFloat64List> zPriors, Float64 cstLog,
                      std::shared_ptr<CPdfMargZLogResult> postmargZResult)
{
    bool verbose = false;

    if (meritResults.size() != zPriors.size())
    {
        Log.LogError(
            "Pdfz: Pdfz-bestchi2 problem. merit.size (%d) != prior.size (%d)",
            meritResults.size(), zPriors.size());
        return -9;
    }
    if (meritResults.size() < 1 || zPriors.size() < 1 || redshifts.size() < 1)
    {
        Log.LogError("Pdfz: Pdfz-bestchi2 problem. merit.size (%d), prior.size "
                     "(%d), or redshifts.size (%d) is zero !",
                     meritResults.size(), zPriors.size(), redshifts.size());
        return -99;
    }

    // build best chi2 vector
    std::vector<Float64> chi2Min(redshifts.size(), DBL_MAX);
    for (Int32 km = 0; km < meritResults.size(); km++)
    {
        if (verbose)
        {
            Log.LogInfo("Pdfz:-bestchi2: processing chi2-result km=%d", km);
        }

        // Todo: Check if the status is OK ?
        // meritResult->Status[i] == COperator::nStatus_OK

        // Todo: use the priors for the min chi2 search ?
        for (UInt32 k = 0; k < redshifts.size(); k++)
        {
            if (true /*meritResult->Status[k]== COperator::nStatus_OK*/) // todo:
                                                                         // check
                                                                         // (temporarily
                                                                         // considers
                                                                         // status
                                                                         // is
                                                                         // always
                                                                         // OK)
            {
                if (meritResults[km][k] < chi2Min[k])
                {
                    chi2Min[k] = meritResults[km][k];
                }
            }
        }
    }

    // estimate Posterior on the best chi2
    CPdfz pdfz;
    TFloat64List logProba;
    Float64 logEvidence;
    TFloat64List zprior;
    zprior = pdfz.GetConstantLogZPrior(redshifts.size());
    Int32 retPdfz = pdfz.Compute(chi2Min, redshifts, cstLog, zprior, logProba, logEvidence);
    if (retPdfz != 0)
    {
        Log.LogError("Pdfz: Pdfz-bestchi2: Pdfz computation failed");
        return -1;
    } else
    {
        postmargZResult->valEvidenceLog = logEvidence;
        postmargZResult->countTPL = redshifts.size(); // assumed 1 model per z
        postmargZResult->Redshifts.resize(redshifts.size());
        postmargZResult->valProbaLog.resize(redshifts.size());
        for (UInt32 k = 0; k < redshifts.size(); k++)
        {
            postmargZResult->Redshifts[k] = redshifts[k];
            postmargZResult->valProbaLog[k] = logProba[k];
        }
    }

    return 0;
}
