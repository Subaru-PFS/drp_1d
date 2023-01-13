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
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/operator/templatefittingresult.h"
#include "RedshiftLibrary/processflow/parameterstore.h"
#include "RedshiftLibrary/processflow/resultstore.h"
#include "RedshiftLibrary/statistics/deltaz.h"
#include "RedshiftLibrary/statistics/pdfcandidateszresult.h"

#include <boost/filesystem.hpp>
#include <boost/math/special_functions.hpp>
#include <fstream>
#include <iostream>
#include <istream>

#include <boost/test/unit_test.hpp>
using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(Statistics_deltaz)

Float64 precision = 1e-12;

BOOST_AUTO_TEST_CASE(GetIndices_test) {
  // simplified test
  TFloat64List Redshifts_easy = {1., 2., 3., 4., 5., 6., 7.};
  Int32 iz, izmin, izmax;
  CDeltaz deltaz;
  Int32 ret;
  Int32 half_samples_nb = 2;

  // Target outside redshifts range
  Float64 redshift_target = 8.0;
  BOOST_CHECK_THROW(deltaz.GetIndices(Redshifts_easy, redshift_target,
                                      half_samples_nb, iz, izmin, izmax),
                    GlobalException);

  // Target included in redshifts range but not defined in
  redshift_target = 4.5;
  BOOST_CHECK_THROW(deltaz.GetIndices(Redshifts_easy, redshift_target,
                                      half_samples_nb, iz, izmin, izmax),
                    GlobalException);

  // Target in redshifts range
  redshift_target = 4.0;
  ret = deltaz.GetIndices(Redshifts_easy, redshift_target, half_samples_nb, iz,
                          izmin, izmax);
  BOOST_CHECK(iz == 3 && izmin == 1 && izmax == 5);

  redshift_target = 7.0;
  ret = deltaz.GetIndices(Redshifts_easy, redshift_target, half_samples_nb, iz,
                          izmin, izmax);
  BOOST_CHECK(iz == 6 && izmin == 4 && izmax == 6);

  redshift_target = 1.0;
  ret = deltaz.GetIndices(Redshifts_easy, redshift_target, half_samples_nb, iz,
                          izmin, izmax);
  BOOST_CHECK(iz == 0 && izmin == 0 && izmax == 2);
}

BOOST_AUTO_TEST_CASE(Compute_test) {
  // simplified test
  TFloat64List Redshifts_easy = {1., 2., 3., 4., 5., 6., 7.};
  TFloat64List merits = {1., 1., 1., 8., 1., 1., 1.};
  Int32 iz, izmin, izmax;
  CDeltaz deltaz;
  Float64 sigma, deltaz_ref;

  iz = 3;
  izmin = 0;
  izmax = 6;
  deltaz_ref = 1.;
  sigma = deltaz.Compute(merits, Redshifts_easy, iz, izmin, izmax);
  BOOST_CHECK_CLOSE(sigma, deltaz_ref, precision);

  deltaz_ref = 1.7320508075688776;
  sigma = deltaz.Compute3ddl(merits, Redshifts_easy, iz, izmin, izmax);
  BOOST_CHECK_CLOSE(sigma, deltaz_ref, precision);

  // Compute : c0 <= 0
  merits = {0., 0., 0., 0., 0., 0., 0.};
  BOOST_CHECK_THROW(deltaz.Compute(merits, Redshifts_easy, iz, izmin, izmax),
                    InternalException);

  // Compute3ddl : c2 <= 0;
  BOOST_CHECK_THROW(
      deltaz.Compute3ddl(merits, Redshifts_easy, iz, izmin, izmax),
      InternalException);

  // Compute3ddl : (izmax - izmin +1) < 3
  merits = {1., 1., 1., 5., 1., 1., 1.};
  izmin = 2;
  izmax = 3;
  BOOST_CHECK_THROW(
      deltaz.Compute3ddl(merits, Redshifts_easy, iz, izmin, izmax),
      GlobalException);
}

BOOST_AUTO_TEST_CASE(GetDeltaz_test) {
  // simplified test
  TFloat64List Redshifts_easy = {1., 2., 3., 4., 5., 6., 7.};
  TFloat64List merits = {1., 1., 1., 8., 1., 1., 1.};
  Float64 redshift_target = 4.;
  CDeltaz deltaz;
  Float64 dz_1, dz_2, dz_ref;
  Int32 ret;

  // redshifts size = 0
  BOOST_CHECK_THROW(deltaz.GetDeltaz({}, merits, redshift_target),
                    GlobalException);

  // COMPUTE
  dz_ref = 1.;
  dz_1 = deltaz.GetDeltaz(Redshifts_easy, merits, redshift_target);
  BOOST_CHECK_CLOSE(dz_1, dz_ref, precision);

  // COMPUTE3DDL
  dz_ref = 1.7320508075688776;
  dz_2 = deltaz.GetDeltaz(Redshifts_easy, merits, redshift_target, true);
  BOOST_CHECK_CLOSE(dz_2, dz_ref, precision);

  // c0 <= 0 : default value
  merits = {0., 0., 0., 0., 0., 0., 0.};
  dz_2 = deltaz.GetDeltaz(Redshifts_easy, merits, redshift_target);
  BOOST_CHECK_CLOSE(dz_2, 0.001, precision);

  // c2 <= 0 : default value
  dz_2 = deltaz.GetDeltaz(Redshifts_easy, merits, redshift_target, true);
  BOOST_CHECK_CLOSE(dz_2, 0.001, precision);
}

BOOST_AUTO_TEST_SUITE_END()
