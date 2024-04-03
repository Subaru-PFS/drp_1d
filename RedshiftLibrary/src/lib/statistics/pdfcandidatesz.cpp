// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
//
// https://www.lam.fr/
//
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
//
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use,
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info".
//
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability.
//
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or
// data to be ensured and,  more generally, to use and operate it in the
// same conditions as regards security.
//
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include "RedshiftLibrary/extremum/extremum.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/operator/pdfz.h"
#include "RedshiftLibrary/statistics/deltaz.h"
#include "RedshiftLibrary/statistics/pdfcandidatesz.h"
#include "RedshiftLibrary/statistics/pdfcandidateszresult.h"

using namespace NSEpic;
using namespace std;

CPdfCandidatesZ::CPdfCandidatesZ(const TCandidateZbyID &candidates)
    : m_candidates(candidates) {}

CPdfCandidatesZ::CPdfCandidatesZ(const TFloat64List &redshifts) {
  for (Int32 i = 0; i < redshifts.size(); ++i) {
    const std::string Id = "EXT" + to_string(i);
    m_candidates[Id] = std::make_shared<TCandidateZ>();
    m_candidates[Id]->Redshift = redshifts[i];
  }
}

/**
 * 1) Fix Deltaz problems: none is passed or none could be compute --> use
 * default values 2) Check if integration windows overlap, mostly for close
 * candidates 1) if no, keep deltaz value 2) if overlapping, update the
 * half-width of the left and right sides of the integration windows Note:
 * Output range includes the multiplication by (1+z)-> Returns a list of
 * identified very close candidates, at 2*1E-4
 */
TStringList
CPdfCandidatesZ::SetIntegrationRanges(const TFloat64Range PdfZRange,
                                      TCandidateZRangebyID &ranges) {
  ranges.clear();

  for (auto &candidatePair : m_candidates) {
    const std::string &Id = candidatePair.first;
    std::shared_ptr<TCandidateZ> &cand = candidatePair.second;

    // check cases where deltaz couldnt be computed or wasnt set--> use default
    // value,
    if (cand->Deltaz == -1)
      cand->Deltaz = m_dzDefault * (1 + cand->Redshift);

    const Float64 halfWidth = integrationWindowHalfWidth(cand->Deltaz);

    // initialize range boundaries for each candidate]
    TFloat64Range range = {cand->Redshift - halfWidth,
                           cand->Redshift + halfWidth};
    ranges.emplace(Id, std::move(range));
    ranges[Id].IntersectWith(PdfZRange);
  };

  // sort candidate keys (Ids) by decreasing redshifts
  TStringList candidatesKeys(m_candidates.size());
  std::transform(m_candidates.begin(), m_candidates.end(),
                 candidatesKeys.begin(),
                 [](const pair_Id_TCandidateZ &c) { return c.first; });

  TCandidateZbyID &candidates = m_candidates;

  std::sort(candidatesKeys.rbegin(), candidatesKeys.rend(),
            [&candidates](std::string Id1, std::string Id2) {
              return candidates[Id1]->Redshift < candidates[Id2]->Redshift;
            });

  // trim overlapping ranges
  TStringList duplicates = {};
  for (auto candidateKey = candidatesKeys.begin();
       candidateKey != candidatesKeys.end() - 1; ++candidateKey) {
    std::string &candidateKeyHigh = *candidateKey;
    std::string &candidateKeyLow = *(candidateKey + 1);
    Float64 &zHigh = candidates[candidateKeyHigh]->Redshift;
    Float64 &zLow = candidates[candidateKeyLow]->Redshift;

    Float64 zSizeHighRange =
        ranges[candidateKeyHigh].GetEnd() - ranges[candidateKeyHigh].GetBegin();
    Float64 zSizeLowRange =
        ranges[candidateKeyLow].GetEnd() - ranges[candidateKeyLow].GetBegin();

    Float64 overlap =
        ranges[candidateKeyLow].GetEnd() - ranges[candidateKeyHigh].GetBegin();

    if (overlap <= 0)
      continue;

    // in the case of duplicates, split equally the overlap between both
    // candidates
    Float64 overlapHighRatio = overlap / zSizeHighRange;
    Float64 overlapLowRatio = overlap / zSizeLowRange;
    if (overlapHighRatio < OVERLAP_THRESHOLD_PDF_INTEGRATION &&
        overlapLowRatio < OVERLAP_THRESHOLD_PDF_INTEGRATION) {
      Log.LogDebug(Formatter()
                   << "    CPdfCandidatesZ::SetIntegrationRanges: integration "
                      "supports overlap for "
                   << zHigh << " and " << zLow);
      Float64 midOverlap = ranges[candidateKeyLow].GetEnd() - 0.5 * overlap;
      ranges[candidateKeyLow].SetEnd(midOverlap);
      ranges[candidateKeyHigh].SetBegin(midOverlap);
    } else {
      Log.LogInfo(Formatter()
                  << " CPdfCandidatesZ::SetIntegrationRanges: very close "
                     "candidates are identified "
                  << zHigh << " and " << zLow);
      std::string candidateKeyToRemove;
      if (candidates[candidateKeyLow]->ValProba >
          candidates[candidateKeyHigh]->ValProba) {
        candidateKeyToRemove = candidateKeyHigh;
      } else {
        candidateKeyToRemove = candidateKeyLow;
      }
      duplicates.push_back(candidateKeyToRemove);
    }
  }

  // iterate over computed ranges and check that corresponding zcandidates
  // belong to that range, otherwise throw error
  for (const auto &candidateKey : candidatesKeys) {
    // if current index belongs to the duplicates b vector, skip testing it
    // and only test the others
    if (find(duplicates.begin(), duplicates.end(), candidateKey) !=
        duplicates.end())
      continue;
    if (candidates[candidateKey]->Redshift < ranges[candidateKey].GetBegin() ||
        candidates[candidateKey]->Redshift > ranges[candidateKey].GetEnd())
      THROWG(
          INTERNAL_ERROR,
          Formatter() << "Failed to identify a range including the candidate "
                      << candidates[candidateKey]->Redshift);
  }

  return duplicates;
}

void CPdfCandidatesZ::computeCandidatesDeltaz(const TRedshiftList &PdfRedshifts,
                                              const TFloat64List &PdfProbaLog) {
  CDeltaz deltaz_op;
  for (auto &c : m_candidates)
    c.second->Deltaz =
        deltaz_op.GetDeltaz(PdfRedshifts, PdfProbaLog, c.second->Redshift);
}

inline Float64 CPdfCandidatesZ::integrationWindowHalfWidth(Float64 DeltaZ) {
  return 3 * DeltaZ;
}

/**
 * @brief CPdfCandidatesZ::Compute
 */
std::shared_ptr<PdfCandidatesZResult>
CPdfCandidatesZ::Compute(TRedshiftList const &PdfRedshifts,
                         TFloat64List const &PdfProbaLog) {
  Log.LogInfo(
      Formatter() << "    CPdfCandidatesZ::Compute pdf peaks using method= "
                  << m_optMethod
                  << " {0:direct "
                     "integration; 1:gaussian fitting}");
  computeCandidatesDeltaz(PdfRedshifts, PdfProbaLog);
  TCandidateZRangebyID zranges;
  TStringList duplicates =
      SetIntegrationRanges(TFloat64Range(PdfRedshifts), zranges);
  CDeltaz deltaz_op;
  for (auto &c : m_candidates) {
    const std::string &Id = c.first;
    std::shared_ptr<TCandidateZ> &cand = c.second;

    // common for both m_optmethod
    //  check if current candidate belongs to the identified duplicates list
    // if yes, force its pdf value to 0 and avoid calling getCandidatexxx
    if (find(duplicates.begin(), duplicates.end(), Id) != duplicates.end()) {
      cand->ValSumProba = 0.;
      continue;
    }

    if (m_optMethod == 0) {
      getCandidateSumTrapez(PdfRedshifts, PdfProbaLog, zranges[Id], cand);
    } else {
      bool GaussFitok = getCandidateRobustGaussFit(PdfRedshifts, PdfProbaLog,
                                                   zranges[Id], cand);

      cand->ValSumProba = NAN;
      if (GaussFitok)
        cand->ValSumProba = cand->GaussAmp * cand->GaussSigma * sqrt(2 * M_PI);
    }
  }

  std::shared_ptr<PdfCandidatesZResult> result =
      std::make_shared<PdfCandidatesZResult>(m_optMethod);

  SortByValSumProbaInt(result->m_ranked_candidates);

  return result;
}

void CPdfCandidatesZ::SortByValSumProbaInt(
    TCandidateZbyRank &ranked_candidates) const {
  // sort m_candidates keys (Ids) by deacreasing integ proba
  TStringList Ids(m_candidates.size());
  std::transform(m_candidates.begin(), m_candidates.end(), Ids.begin(),
                 [](const pair_Id_TCandidateZ &c) { return c.first; });
  const TCandidateZbyID &c = m_candidates;
  std::stable_sort(Ids.rbegin(), Ids.rend(),
                   [&c](std::string Id1, std::string Id2) {
                     return c.at(Id1)->ValSumProba < c.at(Id2)->ValSumProba;
                   });

  ranked_candidates.clear();
  for (const auto &Id : Ids)
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
void CPdfCandidatesZ::getCandidateSumTrapez(
    const TRedshiftList &redshifts, const TFloat64List &valprobalog,
    const TFloat64Range &zrange,
    std::shared_ptr<TCandidateZ> &candidate) const {
  // check that redshifts are sorted
  if (!std::is_sorted(redshifts.begin(), redshifts.end()))
    THROWG(INTERNAL_ERROR, Formatter() << "redshifts are not sorted");

  // find indexes kmin, kmax so that zmin and zmax are inside [
  // redshifts[kmin]:redshifts[kmax] ]
  Int32 kmin = -1;
  Int32 kmax = -1;
  zrange.getEnclosingIntervalIndices(redshifts, candidate->Redshift, kmin,
                                     kmax);
  TFloat64List ZinRange =
      TFloat64List(redshifts.begin() + kmin, redshifts.begin() + kmax + 1);
  candidate->ValSumProbaZmin = ZinRange.front();
  candidate->ValSumProbaZmax = ZinRange.back();

  TFloat64List valprobainRange =
      TFloat64List(valprobalog.begin() + kmin, valprobalog.begin() + kmax + 1);

  Float64 logSum = COperatorPdfz::logSumExpTrick(valprobainRange, ZinRange);
  candidate->ValSumProba = exp(logSum);
}

bool CPdfCandidatesZ::getCandidateRobustGaussFit(
    const TRedshiftList &redshifts, const TFloat64List &valprobalog,
    const TFloat64Range &zrange,
    std::shared_ptr<TCandidateZ> &candidate) const {
  Int32 fitSuccessful = false;
  Int32 nTry = 5;
  Int32 iTry = 0;

  TFloat64Range current_zrange = zrange;
  Float64 zwidth_max = std::max(candidate->Redshift - zrange.GetBegin(),
                                zrange.GetEnd() - candidate->Redshift);
  while (!fitSuccessful && iTry < nTry) {
    getCandidateGaussFit(redshifts, valprobalog, current_zrange, candidate);
    if (candidate->GaussSigma < zwidth_max * 2.0 &&
        std::abs(candidate->GaussSigma / candidate->GaussSigmaErr) > 1e-2) {
      fitSuccessful = true;
    } else {
      Log.LogDebug(Formatter()
                   << "    CPdfCandidatesZ::getCandidateRobustSumGaussFit - "
                      "iTry="
                   << iTry << ", for zcandidate=" << candidate->Redshift
                   << ". Found gaussAmp=" << candidate->GaussAmp
                   << " and gaussSigma=" << candidate->GaussSigma);

      Log.LogDebug("Retrying with different parameters");
    }
    zwidth_max /= 2.0;
    current_zrange.IntersectWith(TFloat64Range(
        candidate->Redshift - zwidth_max, candidate->Redshift + zwidth_max));
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

int pdfz_lmfit_f(const gsl_vector *x, void *data, gsl_vector *f) {
  size_t n = ((struct pdfz_lmfitdata *)data)->n;
  Float64 *y = ((struct pdfz_lmfitdata *)data)->y;
  Float64 *z = ((struct pdfz_lmfitdata *)data)->z;
  Float64 zcenter = ((struct pdfz_lmfitdata *)data)->zcenter;

  double a = gsl_vector_get(x, 0);
  double sigma = gsl_vector_get(x, 1);

  for (Int32 i = 0; i < n; i++) {
    Float64 t = z[i] - zcenter;
    const Float64 xsurc = t / sigma;
    Float64 Yi = a * exp(-0.5 * xsurc * xsurc);
    gsl_vector_set(f, i, Yi - y[i]);
  }

  return GSL_SUCCESS;
}

int pdfz_lmfit_df(const gsl_vector *x, void *data, gsl_matrix *J) {
  size_t n = ((struct pdfz_lmfitdata *)data)->n;
  Float64 *z = ((struct pdfz_lmfitdata *)data)->z;
  Float64 zcenter = ((struct pdfz_lmfitdata *)data)->zcenter;

  double a = gsl_vector_get(x, 0);
  double sigma = gsl_vector_get(x, 1);

  size_t i;

  for (i = 0; i < n; i++) {
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
void CPdfCandidatesZ::getCandidateGaussFit(
    const TRedshiftList &redshifts, const TFloat64List &valprobalog,
    const TFloat64Range &zrange,
    std::shared_ptr<TCandidateZ> &candidate) const {
  Log.LogDebug("    CPdfCandidatesZ::getCandidateSumGaussFit - Starting pdf "
               "peaks gaussian fitting");

  // check that redshifts are sorted
  if (!std::is_sorted(redshifts.begin(), redshifts.end()))
    THROWG(INTERNAL_ERROR, Formatter() << "redshifts are not sorted");

  // find indexes kmin, kmax so that zmin and zmax are inside [
  // redshifts[kmin]:redshifts[kmax] ]
  Int32 kmin = -1;
  Int32 kmax = -1;

  zrange.getEnclosingIntervalIndices(redshifts, candidate->Redshift, kmin,
                                     kmax);

  // initialize GSL
  const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;

  size_t n = kmax - kmin +
             1; // n samples on the support, /* number of data points to fit */
  size_t p = 2; // DOF = 1.amplitude + 2.width

  Log.LogDebug(
      Formatter() << "    CPdfCandidatesZ::getCandidateSumGaussFit - n=" << n
                  << ", p=" << p);
  if (n < p)
    THROWG(INTERNAL_ERROR, "LMfit has not enough samples on support");

  gsl_matrix *J = gsl_matrix_alloc(n, p);
  gsl_matrix *covar = gsl_matrix_alloc(p, p);

  const Float64 &zc = candidate->Redshift;

  // initialize lmfit with previously estimated values ?
  Float64 maxP = 0.0;
  std::for_each(valprobalog.begin() + kmin, valprobalog.begin() + kmin + n,
                [&maxP](const Float64 &vp) { maxP = std::max(maxP, exp(vp)); });

  TFloat64List x_init(p);
  Float64 normFactor = (maxP <= 0.0) ? 1. : maxP;
  x_init[0] = (maxP > 0.0) ? maxP / normFactor : 1.0;
  x_init[1] = std::max(candidate->Redshift - zrange.GetBegin(),
                       zrange.GetEnd() - candidate->Redshift) /
              2.0;
  Log.LogDebug(
      Formatter() << "    CPdfCandidatesZ::getCandidateSumGaussFit - init a="
                  << x_init[0] << ", sigma=" << x_init[1]);

  TFloat64List weights(n, 1.0); // no weights
  gsl_vector_view x = gsl_vector_view_array(x_init.data(), p);
  gsl_vector_view w = gsl_vector_view_array(weights.data(), n);

  // This is the data to be fitted;
  double y[n], z[n];
  for (Int32 i = 0; i < n; i++) {
    Float64 idx = i + kmin;
    y[i] = exp(valprobalog[idx]) / normFactor;
    z[i] = redshifts[idx];
  }

  struct pdfz_lmfitdata d = {n, y, z, zc};
  gsl_multifit_function_fdf f;
  f.f = &pdfz_lmfit_f;
  f.df = &pdfz_lmfit_df;
  f.n = n;
  f.p = p;
  f.params = &d;

  Log.LogDebug(
      Formatter()
      << "    CPdfCandidatesZ::getCandidateSumGaussFit - LMfit data ready");

  gsl_multifit_fdfsolver *s = gsl_multifit_fdfsolver_alloc(T, n, p);

  if (s == 0) {
    gsl_matrix_free(covar);
    gsl_matrix_free(J);
    THROWG(INTERNAL_ERROR, "Unable to allocate the multifit solver");
  }

  /* initialize solver with starting point and weights */
  gsl_multifit_fdfsolver_wset(s, &f, &x.vector, &w.vector);

  /* compute initial residual norm */
  gsl_vector *res_f = gsl_multifit_fdfsolver_residual(s);
  double chi0 = gsl_blas_dnrm2(res_f);

  /* solve the system with a maximum of maxIterations iterations */
  const double xtol = 1e-8;
  const double gtol = 1e-8;
  const double ftol = 1e-8;
  Int32 maxIterations = 500;
  Int32 info;
  Int32 status =
      gsl_multifit_fdfsolver_driver(s, maxIterations, xtol, gtol, ftol, &info);
  gsl_multifit_fdfsolver_jac(s, J);
  gsl_multifit_covar(J, 0.0, covar);

  /* compute final residual norm */
  double chi = gsl_blas_dnrm2(res_f);
  double dof = n - p;
  double c = GSL_MAX_DBL(1, chi / sqrt(dof));
#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar, i, i))

  Log.LogDebug(Formatter() << "summary from method '"
                           << gsl_multifit_fdfsolver_name(s)
                           << "': status = " << gsl_strerror(status) << " ("
                           << status << "), number of iterations: "
                           << gsl_multifit_fdfsolver_niter(s));
  Log.LogDebug(Formatter() << "function evaluations: " << f.nevalf
                           << " , Jacobian evaluations: " << f.nevaldf);
  Log.LogDebug(
      Formatter()
      << "reason for stopping: " << info
      << " where {1:small step size, 2:small gradient, small change in f }");
  Log.LogDebug(Formatter() << "initial |f(x)| = " << chi0
                           << "final   |f(x)| = " << chi);

  { Log.LogDebug(Formatter() << "chisq/dof = " << pow(chi, 2.0) / dof); }
  candidate->GaussAmp = gsl_vector_get(s->x, 0) * normFactor;
  candidate->GaussAmpErr = c * ERR(0) * normFactor;
  candidate->GaussSigma = abs(gsl_vector_get(s->x, 1));
  candidate->GaussSigmaErr = c * ERR(1);

  gsl_multifit_fdfsolver_free(s);
  gsl_matrix_free(covar);
  gsl_matrix_free(J);
}
