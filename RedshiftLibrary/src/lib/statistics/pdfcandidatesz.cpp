#include "RedshiftLibrary/statistics/pdfcandidatesz.h"
#include "RedshiftLibrary/statistics/pdfcandidateszresult.h"
#include "RedshiftLibrary/operator/pdfz.h"
#include "RedshiftLibrary/statistics/deltaz.h"
#include "RedshiftLibrary/extremum/extremum.h"
#include "RedshiftLibrary/log/log.h"

using namespace NSEpic;
using namespace std;

#include <fstream>
#include <iostream>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlin.h>


CPdfCandidatesZ::CPdfCandidatesZ(const TCandidateZbyID & candidates):
    m_candidates(candidates)
{}

CPdfCandidatesZ::CPdfCandidatesZ(const TFloat64List & redshifts)
{
    for (UInt32 i =0; i<redshifts.size(); ++i)
    {
        const std::string Id = "EXT" + to_string(i);
        m_candidates[Id].Redshift = redshifts[i];
    }
}


/**
 * 1) Fix Deltaz problems: none is passed or none could be compute --> use default values
 * 2) Check if integration windows overlap, mostly for close candidates
 *      1) if no, keep deltaz value
 *      2) if overlapping, update the half-width of the left and right sides of the integration windows
 * Note: Output range includes the multiplication by (1+z).
 * Returns a list of identified very close candidates, at 2*1E-4
*/
TStringList CPdfCandidatesZ::SetIntegrationWindows(const TFloat64Range PdfZRange, TCandidateZRangebyID & ranges)
{
    Bool nodz = false;
    Int32 n = m_candidates.size();

    ranges.clear();

    for(auto & c: m_candidates){
        const std::string & Id = c.first;
        TCandidateZ & cand = c.second;

        //check cases where deltaz couldnt be computed or wasnt set--> use default value,
        if(cand.Deltaz == -1 || nodz)
            cand.Deltaz = m_dzDefault*(1+cand.Redshift);

        const Float64 halfWidth = 3*cand.Deltaz;

        //initialize range boundaries for each candidate]
        TFloat64Range range = {cand.Redshift - halfWidth, cand.Redshift + halfWidth};
        ranges.emplace(Id, std::move(range));
        ranges[Id].IntersectWith(PdfZRange);
    };

    // sort candidate keys (Ids) by decreasing redshifts
    std::vector<std::string> Ids;
    for (const auto & c : m_candidates){
        Ids.push_back(c.first); // keys = ids
    }
    TCandidateZbyID & c = m_candidates;    
    std::sort(Ids.rbegin(), Ids.rend(),
        [&c](std::string Id1, std::string Id2) {return c[Id1].Redshift < c[Id2].Redshift;});

    // trim overlapping ranges
    TStringList b = {}; 
    for(auto Id_it = Ids.begin(); Id_it != Ids.end()-1; ++Id_it){
        std::string & Id_h = *Id_it;
        std::string & Id_l = *(Id_it + 1);
        Float64 & z_h  = c[Id_h].Redshift;
        Float64 & z_l  = c[Id_l].Redshift;
        Float64 overlap = ranges[Id_h].GetBegin() - ranges[Id_l].GetEnd();
        if(overlap < 0){
            //in the case of duplicates, trim completely the range of the second cand
            if((z_h -  z_l)>2*1E-4){
                Log.LogDebug("    CPdfCandidatesZ::SetIntegrationWindows: integration supports overlap for %f and %f", z_h, z_l );
                ranges[Id_h].SetBegin(( max(z_l, ranges[Id_h].GetBegin()) + min(z_h, ranges[Id_l].GetEnd()) )/2);
            }else{
                Log.LogInfo(" CPdfCandidatesZ::SetIntegrationWindows: very close candidates are identified %f and %f", z_h,  z_l);
                b.push_back(Id_l);
            }
            ranges[Id_l].SetEnd(ranges[Id_h].GetBegin() - 1E-4);
         }
    }

    //iterate over computed ranges and check that corresponding zcandidates belong to that range, otherwise throw error
    for(const auto & Id: Ids){
        //if currend index belongs to the duplicates b vector, skip testing it and only test the others
        if(find(b.begin(), b.end(), Id)!= b.end())
            continue;
        if( c[Id].Redshift>= ranges[Id].GetBegin() && 
            c[Id].Redshift<= ranges[Id].GetEnd()){
            continue;
        }else{
            Log.LogError("CPdfCandidatesZ::SetIntegrationWindows: Failed to identify a range including the candidate %f", c[Id].Redshift);
            throw runtime_error("CPdfCandidatesZ::SetIntegrationWindows: Failed to identify a range including the candidate! Aborting");
        }
    } 
    return b; 
}

/**
 * @brief CPdfCandidatesZ::Compute
 */
std::shared_ptr<PdfCandidatesZResult> CPdfCandidatesZ::Compute(TRedshiftList const & PdfRedshifts, TFloat64List const & PdfProbaLog)
{
    if(m_optMethod==0)
    {
        Log.LogInfo("    CPdfCandidatesZ::Compute pdf peaks info (method=direct integration)" );
    }else{
        Log.LogInfo("    CPdfCandidatesZ::Compute pdf peaks info (method=gaussian fitting)" );
    }

    // compute deltaz
    CDeltaz deltaz_op;
    for (auto & c: m_candidates)
    {
        c.second.Deltaz = deltaz_op.GetDeltaz(PdfRedshifts, PdfProbaLog, c.second.Redshift);
    }    

    TCandidateZRangebyID zranges;
    TStringList duplicates = SetIntegrationWindows( TFloat64Range(PdfRedshifts), zranges);
    for(auto & c: m_candidates)
    {
        const std::string & Id = c.first;
        TCandidateZ & cand = c.second;
        if(m_optMethod==0)
        {
            //check if current candidate belongs to the identified duplicates list
            //if yes, force its pdf value to 0 and avoid callling getCandidateSumTrapez
            if(find(duplicates.begin(), duplicates.end(),Id)!=duplicates.end())
                cand.ValSumProba = 0;
            else
                getCandidateSumTrapez( PdfRedshifts, PdfProbaLog, zranges[Id], cand);
        }else
        {
            //TODO: this requires further check ?...
            if(find(duplicates.begin(), duplicates.end(),Id)!=duplicates.end()){
                cand.ValSumProba = 0;
                continue;
            }            
            Bool GaussFitok = getCandidateRobustGaussFit( PdfRedshifts, PdfProbaLog, zranges[Id], cand);
            if(GaussFitok)
            {
                cand.ValSumProba = cand.GaussAmp*cand.GaussSigma*sqrt(2*M_PI);
            }else{
                cand.ValSumProba = -1.;
            }
        }
    }

    std::shared_ptr<PdfCandidatesZResult> result = std::make_shared<PdfCandidatesZResult>(m_optMethod);

    SortByValSumProbaInt(result->m_ranked_candidates);

    return result;
}



void CPdfCandidatesZ::SortByValSumProbaInt(TCandidateZbyRank & ranked_candidates) const
{

    // sort m_candidates keys (Ids) by deacreasing integ proba
    std::vector<std::string> Ids;
    for (const auto & c : m_candidates){
        Ids.push_back(c.first); // keys = ids
    }
    const TCandidateZbyID & c = m_candidates; 
    std::stable_sort(Ids.rbegin(), Ids.rend(),
        [&c](std::string Id1, std::string Id2) {return c.at(Id1).ValSumProba < c.at(Id2).ValSumProba;});
    
    ranked_candidates.clear();
    for (const auto & Id: Ids)
        ranked_candidates.emplace_back(Id, m_candidates.at(Id));
}

/**
 * @brief CPdfCandidatesZ::getCandidateSumTrapez
 * @param redshifts
 * @param valprobalog
 * @param zcandidate
 * @param zrange corresponds to the concerned range
 * @return -1.0 if error, else sum around the candidate
 */
Bool CPdfCandidatesZ::getCandidateSumTrapez(const TRedshiftList & redshifts,
                                     const TFloat64List & valprobalog,
                                     const TFloat64Range & zrange,
                                     TCandidateZ & candidate
                                     ) const
{
  //TODO use a std function and throw exception
  //TODO check that this is really necessary and not just a debug functionnality
  //     -> use a DEBUG directive ? 

    // check that redshifts are sorted
    for (UInt32 k = 1; k < redshifts.size(); k++)
    {
        if (redshifts[k] < redshifts[k - 1])
        {
            Log.LogError("    CPdfCandidatesZ::getCandidateSumTrapez - redshifts are not sorted for (at least) index={}", k);
            throw runtime_error("CPdfCandidatesZ::getCandidateSumTrapez - redshifts are not sorted");
        }
    }

    // find indexes kmin, kmax so that zmin and zmax are inside [
    // redshifts[kmin]:redshifts[kmax] ]
    Int32 kmin = -1;
    Int32 kmax = -1;
    bool ok = zrange.getEnclosingIntervalIndices(redshifts, candidate.Redshift, kmin, kmax);

    if(!ok || kmin==-1 || kmax==-1){
        Log.LogError("CPdfCandidatesZ::getCandidateSumTrapez could not find enclosing interval"); 
        throw runtime_error("CPdfCandidatesZ::getCandidateSumTrapez could not find enclosing interval");
    }
        
    TFloat64List ZinRange = TFloat64List(redshifts.begin()+kmin, redshifts.begin()+kmax+1);
    candidate.ValSumProbaZmin = ZinRange.front();
    candidate.ValSumProbaZmax = ZinRange.back();
    TFloat64List valprobainRange = TFloat64List(valprobalog.begin()+kmin, valprobalog.begin()+kmax+1);

    Int32 sumMethod = 1;
    Float64 logSum = COperatorPdfz::logSumExpTrick( valprobainRange, ZinRange, sumMethod);
    candidate.ValSumProba = exp(logSum);
    
    return true;
}


//TODO: this requires a deeper check to include the updates window support range
Bool CPdfCandidatesZ::getCandidateRobustGaussFit(const TRedshiftList & redshifts,
                                        const TFloat64List & valprobalog,
                                        const TFloat64Range & zrange,
                                        TCandidateZ & candidate) const
{
    Int32 fitSuccessful = false;
    Int32 nTry = 5;
    Int32 iTry = 0;

    TFloat64Range current_zrange = zrange;
    Float64 zwidth_max = std::max(candidate.Redshift-zrange.GetBegin(), zrange.GetEnd() - candidate.Redshift);
    while (!fitSuccessful && iTry < nTry)
    {
        Int32 retFit = getCandidateGaussFit(
            redshifts, valprobalog, current_zrange, candidate);
        if (!retFit && candidate.GaussSigma < zwidth_max * 2.0 &&
            std::abs(candidate.GaussSigma / candidate.GaussSigmaErr) > 1e-2)
        {
            fitSuccessful = true;
        }
        else
        {
            Log.LogDebug("    CPdfCandidatesZ::getCandidateRobustSumGaussFit - iTry=%d", iTry);
            Log.LogDebug("    CPdfCandidatesZ::getCandidateRobustSumGaussFit -    for zcandidate=%.5f", candidate.Redshift);
            Log.LogDebug("    CPdfCandidatesZ::getCandidateRobustSumGaussFit -       found gaussAmp=%e", candidate.GaussAmp);
            Log.LogDebug("    CPdfCandidatesZ::getCandidateRobustSumGaussFit -       found gaussSigma=%e", candidate.GaussSigma);
            Log.LogDebug("    CPdfCandidatesZ::getCandidateRobustSumGaussFit -       now going to retry w. different parameters");
        }
        zwidth_max /= 2.0;
        current_zrange.IntersectWith(TFloat64Range(candidate.Redshift - zwidth_max, 
                                                    candidate.Redshift + zwidth_max));
        iTry++;
    }

    return fitSuccessful;
}


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
        Log.LogDebug("Pdfz: Pdfz computation: pdfz_lmfit_f : a=%e, sigma=%e, zcenter=%.5f ", a, sigma, zcenter);
    }

    for (Int32 i = 0; i < n; i++)
    {
        Float64 t = z[i] - zcenter;
        const Float64 xsurc = t / sigma;
        Float64 Yi = a * exp(-0.5 * xsurc * xsurc);
        if (0 && verbose)
        {
            Log.LogDebug("Pdfz: Pdfz computation: pdfz_lmfit_f : for i=%d, Yi=%e", i, Yi);
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
        gsl_matrix_set(J, i, 1, t * t / (sigma * sigma * sigma) * a * exp(-0.5 * xsurc * xsurc));
    }
    return GSL_SUCCESS;
}
//** gaussian fit end**//


Bool CPdfCandidatesZ::getCandidateGaussFit(const TRedshiftList & redshifts,
                                  const TFloat64List & valprobalog,
                                  const TFloat64Range & zrange,
                                  TCandidateZ & candidate) const
{
    Int32 verbose = 0;
    Log.LogDebug("    CPdfCandidatesZ::getCandidateSumGaussFit - Starting pdf peaks gaussian fitting");

    // check that redshifts are sorted
    for (UInt32 k = 1; k < redshifts.size(); k++)
    {
        if (redshifts[k] < redshifts[k - 1])
        {
            Log.LogError("    CPdfCandidatesZ::getCandidateSumGaussFit - redshifts are not sorted for (at least) index={}", k);
            throw runtime_error("CPdfCandidatesZ::getCandidateSumGaussFit - redshifts are not sorted");
            return false;
        }
    }

    // find indexes kmin, kmax so that zmin and zmax are inside [
    // redshifts[kmin]:redshifts[kmax] ]
    Int32 kmin = -1;
    Int32 kmax = - 1;
    bool ok = zrange.getEnclosingIntervalIndices(redshifts, candidate.Redshift, kmin, kmax);

    if(!ok || kmin==-1 || kmax==-1){
        Log.LogError("CPdfCandidatesZ::getCandidateSumGaussFit could not find enclosing interval"); 
        throw runtime_error("CPdfCandidatesZ::getCandidateSumGaussFit could not find enclosing interval");
    }

    if (verbose)
    {
        Log.LogDebug("    CPdfCandidatesZ::getCandidateSumGaussFit - kmax index=%d", kmax);
    }

    // initialize GSL
    const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
    gsl_multifit_fdfsolver *s;
    int status, info;
    size_t i;
    size_t n = kmax - kmin + 2; // n samples on the support, /* number of data points to fit */
    size_t p = 2; // DOF = 1.amplitude + 2.width

    if (verbose)
    {
        Log.LogDebug("    CPdfCandidatesZ::getCandidateSumGaussFit - n=%d, p=%d", n, p);
    }
    if (n < p)
    {
        Log.LogError("    CPdfCandidatesZ::getCandidateSumGaussFit - LMfit not enough samples on support");
        return false;
    }

    gsl_matrix *J = gsl_matrix_alloc(n, p);
    gsl_matrix *covar = gsl_matrix_alloc(p, p);
    double y[n], weights[n], z[n];
    const Float64 & zc = candidate.Redshift;
    gsl_multifit_function_fdf f;

    Float64 *x_init = (Float64 *)calloc(p, sizeof(Float64));
    if (x_init == 0)
    {
        Log.LogError("    CPdfCandidatesZ::getCandidateSumGaussFit - Unable to allocate x_init");
        return false;
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
    if (normFactor <= 0.0)
    {
        normFactor = 1.0;
    }
    if (maxP > 0.0)
    {
        x_init[0] = maxP / normFactor;
    }
    else
    {
        x_init[0] = 1.0;
    }
    x_init[1] = std::max(candidate.Redshift - zrange.GetBegin(), zrange.GetEnd() - candidate.Redshift)/2.0;
    if (verbose)
    {
        Log.LogDebug("    CPdfCandidatesZ::getCandidateSumGaussFit - init a=%e", x_init[0]);
        Log.LogDebug("    CPdfCandidatesZ::getCandidateSumGaussFit - init sigma=%e", x_init[1]);
    }

    gsl_vector_view x = gsl_vector_view_array(x_init, p);
    //    if(x.vector==0){
    //        Log.LogError( "    CPdfCandidatesZ::getCandidateSumGaussFit - Unable to
    //        allocate x"); return -1;
    //    }
    gsl_vector_view w = gsl_vector_view_array(weights, n);
    //    if(w.vector==0){
    //        Log.LogError( "    CPdfCandidatesZ::getCandidateSumGaussFit - Unable to
    //        allocate w"); return -1;
    //    }
    gsl_vector *res_f;
    double chi, chi0;

    const double xtol = 1e-8;
    const double gtol = 1e-8;
    const double ftol = 1e-8;
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

    Log.LogDebug("    CPdfCandidatesZ::getCandidateSumGaussFit - LMfit data ready");

    s = gsl_multifit_fdfsolver_alloc(T, n, p);
    if (s == 0)
    {
        Log.LogError("    CPdfCandidatesZ::getCandidateSumGaussFit - Unable to allocate the multifit solver");
        throw runtime_error("CPdfCandidatesZ::getCandidateSumGaussFit - Unable to allocate the multifit solver");
        return false;
    }

    /* initialize solver with starting point and weights */
    gsl_multifit_fdfsolver_wset(s, &f, &x.vector, &w.vector);

    /* compute initial residual norm */
    res_f = gsl_multifit_fdfsolver_residual(s);
    chi0 = gsl_blas_dnrm2(res_f);

    /* solve the system with a maximum of maxIterations iterations */
    status = gsl_multifit_fdfsolver_driver(s, maxIterations, xtol, gtol, ftol, &info);

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

        Log.LogDebug("summary from method '%s'", gsl_multifit_fdfsolver_name(s));
        Log.LogDebug("number of iterations: %zu", gsl_multifit_fdfsolver_niter(s));
        Log.LogDebug("function evaluations: %zu", f.nevalf);
        Log.LogDebug("Jacobian evaluations: %zu", f.nevaldf);
        Log.LogDebug("reason for stopping: %s", (info == 1) ? "small step size " : (info == 2) ? "small gradient" : "small change in f");
        Log.LogDebug("initial |f(x)| = %g", chi0);
        Log.LogDebug("final   |f(x)| = %g", chi);

        {
            Log.LogDebug("chisq/dof = %g", pow(chi, 2.0) / dof);

            for (Int32 k = 0; k < p; k++)
            {
                if (FIT(k) < 1e-3)
                {
                    Log.LogDebug("A %d     = %.3e +/- %.8f", k, FIT(k), c * ERR(k));
                } else
                {
                    Log.LogDebug("A %d     = %.5f +/- %.8f", k, FIT(k), c * ERR(k));
                }
            }
        }
        Log.LogDebug("status = %s (%d)", gsl_strerror(status), status);
    }

    candidate.GaussAmp = gsl_vector_get(s->x, 0) * normFactor;
    candidate.GaussAmpErr = c * ERR(0) * normFactor;
    candidate.GaussSigma = abs(gsl_vector_get(s->x, 1));
    candidate.GaussSigmaErr = c * ERR(1);

    return true;
}


