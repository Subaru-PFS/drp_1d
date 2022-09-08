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
#include "RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h"

#include <boost/test/unit_test.hpp>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(fluxcorrectioncalzetti_test)

Float64 precision = 1e-12;

BOOST_AUTO_TEST_CASE(overall_test) {
  TFloat64List lbda = {1., 100., 200., 300., 400.};
  TFloat64List flux = {0.1, 0.2, 0.5, 0.3, 0.8};
  CalzettiCorrection calzettiCorr(lbda, flux);
  CSpectrumFluxCorrectionCalzetti spcCorrCalzetti(calzettiCorr, 0., 0.1, 10);

  Float64 ebmv = spcCorrCalzetti.GetEbmvValue(1);
  BOOST_CHECK_CLOSE(ebmv, 0.1, 1e-12);

  Int32 ebmvIndex = spcCorrCalzetti.GetEbmvIndex(0.5);
  BOOST_CHECK(ebmvIndex == 5);

  Float64 dustCoeff = spcCorrCalzetti.GetDustCoeff(0, 115.);
  BOOST_CHECK_CLOSE(dustCoeff, 0.9727472237769651, 1e-12);

  Float64 lambdaMin = spcCorrCalzetti.getLambdaMin();
  BOOST_CHECK(lambdaMin == 1.);

  Float64 lambdaMax = spcCorrCalzetti.getLambdaMax();
  BOOST_CHECK(lambdaMax == 400.);

  Int32 ebmvCoeff = spcCorrCalzetti.GetNPrecomputedEbmvCoeffs();
  BOOST_CHECK(ebmvCoeff == 10.);
}

BOOST_AUTO_TEST_SUITE_END()