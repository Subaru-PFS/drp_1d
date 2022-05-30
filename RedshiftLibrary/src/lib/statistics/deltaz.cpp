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
#include "RedshiftLibrary/statistics/deltaz.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/flag.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/log/log.h"

using namespace NSEpic;
using namespace std;
#include <fstream>

#include <gsl/gsl_multifit.h>

Float64 CDeltaz::GetDeltaz(const TFloat64List &redshifts,
                           const TFloat64List &pdf, const Float64 z,
                           const Int32 gslfit) {
  Float64 dz = -1;
  if (!redshifts.size())
    THROWG(INTERNAL_ERROR, "Redshift range is empty");
  Int32 ret = -1, deltaz_i = 0, maxIter = 2;
  while (deltaz_i < maxIter) { // iterate only twice
    Int32 izmin = -1;
    Int32 iz = -1;
    Int32 izmax = -1;

    // Float64 zRangeHalf = 0.002/(deltaz_i+1);
    // Log.LogInfo("  Deltaz: Deltaz computation nb %i with zRangeHalf %f",
    // deltaz_i, zRangeHalf); TFloat64Range range = TFloat64Range(z -
    // zRangeHalf*(1+z), z + zRangeHalf*(1+z)); ret = GetRangeIndices(redshifts,
    // z, range, iz, izmin, izmax );

    Int32 half_samples_nb = 5 / (deltaz_i + 1);
    ret = GetIndices(redshifts, z, half_samples_nb, iz, izmin, izmax);

    try {
      if (gslfit)
        dz = Compute3ddl(pdf, redshifts, iz, izmin, izmax);
      else
        dz = Compute(pdf, redshifts, iz, izmin, izmax);
      break;
    } catch (GlobalException &e) {
      if (e.getErrorCode() != ErrorCode::DZ_NOT_COMPUTABLE) {
        std::string msg;
        msg = e.getMessage();
        throw GlobalException(e.getErrorCode(), msg, __FILE__, __func__,
                              __LINE__);
      } else {
        if (deltaz_i < maxIter) {
          // Log.LogWarning("  Deltaz: Deltaz computation failed for half range
          // %f", zRangeHalf);
          Flag.warning(Flag.DELTAZ_COMPUTATION_FAILED,
                       Formatter()
                           << "  CDeltaz::" << __func__
                           << ": Deltaz computation failed for half range "
                           << half_samples_nb << " samples");
          deltaz_i++;
        }
        if (deltaz_i == maxIter - 1) {
          dz = 0.001;
        }
      }
    }
  }
  return dz;
}

Int32 CDeltaz::GetIndices(const TFloat64List &redshifts, const Float64 redshift,
                          const Int32 HalfNbSamples, Int32 &iz, Int32 &izmin,
                          Int32 &izmax) {
  // find indexes: iz, izmin and izmax
  izmin = -1;
  TFloat64List::const_iterator iiz =
      std::lower_bound(redshifts.begin(), redshifts.end(), redshift);
  izmax = -1;

  iz = iiz - redshifts.begin();
  if (iiz == redshifts.end() || *iiz != redshift) {
    THROWG(INTERNAL_ERROR, Formatter()
                               << "impossible to get redshift index for z="
                               << redshift);
  }

  izmin = max(iz - HalfNbSamples, 0);
  izmax = min(iz + HalfNbSamples, Int32(redshifts.size() - 1));

  return 0;
}

Float64 CDeltaz::Compute(const TFloat64List &merits,
                         const TFloat64List &redshifts, const Int32 iz,
                         const Int32 izmin, const Int32 izmax) {
  Float64 x0 = redshifts[iz];
  Float64 y0 = merits[iz];
  Float64 xi2, yi, c0;
  Float64 sum = 0, sum2 = 0;

  Float64 sigma = -1.0;

  for (Int32 i = izmin; i < izmax + 1; i++) {
    xi2 = redshifts[i] - x0;
    xi2 *= xi2;
    yi =
        y0 -
        merits[i]; // pour re-caler les y pour que le centre soit à zero pour x0
    sum += xi2 * yi;
    sum2 += xi2 * xi2;
  }
  c0 = sum / sum2;
  if (c0 <= 0) {
    THROWG(DZ_NOT_COMPUTABLE, Formatter() << "impossible to compute sigma");
  }
  sigma = sqrt(1.0 / c0);
  return sigma;
}

// todo : merge with previous function instead of duplicating code...
Float64 CDeltaz::Compute3ddl(const TFloat64List &merits,
                             const TFloat64List &redshifts, const Int32 iz,
                             const Int32 izmin, const Int32 izmax) {
  Float64 sigma = -1.0; // default value

  // quadratic fit
  Int32 i, n;
  Float64 xi, yi, ei, chisq;
  gsl_matrix *X, *cov;
  gsl_vector *y, *w, *c;

  n = izmax - izmin + 1;
  if (n < 3) {
    THROWG(DZ_NOT_COMPUTABLE, Formatter() << "impossible to compute sigma");
  }

  X = gsl_matrix_alloc(n, 3);
  y = gsl_vector_alloc(n);
  w = gsl_vector_alloc(n);

  c = gsl_vector_alloc(3);
  cov = gsl_matrix_alloc(3, 3);

  double x0 = redshifts[iz];
  double y0 = merits[iz];
  for (i = 0; i < n; i++) {
    xi = redshifts[i + izmin];
    yi = merits[i + izmin];
    Log.LogDebug("  y = %+.5e ]\n", yi);
    Log.LogDebug("  x = %+.5e ]\n", xi);
    ei = 1.0; // todo, estimate weighting ?
    gsl_matrix_set(X, i, 0, 1.0);
    gsl_matrix_set(X, i, 1, xi - x0);
    gsl_matrix_set(X, i, 2, (xi - x0) * (xi - x0));

    gsl_vector_set(y, i, y0 - yi);
    gsl_vector_set(w, i, 1.0 / (ei * ei));
  }

  {
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n, 3);
    gsl_multifit_wlinear(X, w, y, c, cov, &chisq, work);
    gsl_multifit_linear_free(work);
  }

#define C(i) (gsl_vector_get(c, (i)))
#define COV(i, j) (gsl_matrix_get(cov, (i), (j)))

  double zcorr = x0 - C(1) / (2.0 * C(2));
  Float64 c2 = C(2);
  if (c2 > 0)
    sigma = sqrt(1.0 / c2);

  // Float64 a = (Float64)(C(0));
  // Float64 b2sur4c = (Float64)(C(1)*C(1)/((Float64)(4.0*C(2))));
  // Float64 logK = ( -(a - b2sur4c)/2.0 );
  // Float64 logarea = log(sigma) + logK + log(2.0*M_PI);
  Log.LogDebug("Center Redshift: %g", x0);
  Log.LogDebug("# best fit: Y = %g + %g X + %g X^2", C(0), C(1), C(2));
  Log.LogDebug("# covariance matrix:\n");
  Log.LogDebug("[ %+.5e, %+.5e, %+.5e  \n", COV(0, 0), COV(0, 1), COV(0, 2));
  Log.LogDebug("  %+.5e, %+.5e, %+.5e  \n", COV(1, 0), COV(1, 1), COV(1, 2));
  Log.LogDebug("  %+.5e, %+.5e, %+.5e ]\n", COV(2, 0), COV(2, 1), COV(2, 2));
  Log.LogDebug("# chisq/n = %g", chisq / n);
  Log.LogDebug("# zcorr = %g", zcorr);
  Log.LogDebug("# sigma = %g", sigma);
  // Log.LogDebug("# logarea = %g", logarea);
  Log.LogDebug("\n");

  gsl_matrix_free(X);
  gsl_vector_free(y);
  gsl_vector_free(w);
  gsl_vector_free(c);
  gsl_matrix_free(cov);

  // results.LogArea[indz] = logarea;
  // results.SigmaZ[indz] = sigma;
  // results.LogAreaCorrectedExtrema[indz] = zcorr;
  if (c2 <= 0)
    THROWG(DZ_NOT_COMPUTABLE, Formatter() << "impossible to compute sigma");
  return sigma;
}
