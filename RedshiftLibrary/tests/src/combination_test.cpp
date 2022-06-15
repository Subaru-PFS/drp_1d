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
#include "RedshiftLibrary/spectrum/combination.h"
#include "RedshiftLibrary/spectrum/spectralaxis.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <stdlib.h>
#include <time.h>

using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(Combination)

BOOST_AUTO_TEST_CASE(SpectrumCombination) {
  CSpectrumCombination combinator;
  TFloat64List lambdas = {1000, 2000, 3000};
  CSpectrumSpectralAxis spectral_axis(lambdas);

  TFloat64List lambdas_bad_1 = {1001, 2000, 3000};
  CSpectrumSpectralAxis spectral_axis_bad_1(lambdas_bad_1);

  TFloat64List lambdas_bad_2 = {1000, 2000, 3000, 4000};
  CSpectrumSpectralAxis spectral_axis_bad_2(lambdas_bad_2);

  TFloat64List flux_bad_1 = {1.2, 2.3, 3.4};
  CSpectrumFluxAxis fluxaxis_bad_1(flux_bad_1);

  TFloat64List flux_bad_2 = {1.2, 2.3, 3.4, 4.5};
  CSpectrumFluxAxis fluxaxis_bad_2(flux_bad_2);

  TFloat64List flux_a = {1.2, 2.3, 3.4};
  CSpectrumFluxAxis fluxaxis_a(flux_a);

  TFloat64List flux_b = {2, 3, 4};
  CSpectrumFluxAxis fluxaxis_b(flux_b);

  CSpectrum spectrum_out(spectral_axis, fluxaxis_a);

  std::vector<std::shared_ptr<CSpectrum>> spectrum_list;

  // not enough spectra
  BOOST_CHECK(combinator.Combine(spectrum_list, spectrum_out) == -3);

  spectrum_list.push_back(
      (std::shared_ptr<CSpectrum>)(new CSpectrum(spectral_axis, fluxaxis_a)));
  spectrum_list.push_back(
      (std::shared_ptr<CSpectrum>)(new CSpectrum(spectral_axis, fluxaxis_b)));

  // good call
  BOOST_CHECK(combinator.Combine(spectrum_list, spectrum_out) == 0);

  // different grids
  spectrum_list.push_back((std::shared_ptr<CSpectrum>)(new CSpectrum(
      spectral_axis_bad_1, fluxaxis_bad_1)));
  BOOST_CHECK(combinator.Combine(spectrum_list, spectrum_out) == -2);

  spectrum_list.pop_back();

  // different grid sizes
  spectrum_list.push_back((std::shared_ptr<CSpectrum>)(new CSpectrum(
      spectral_axis_bad_2, fluxaxis_bad_2)));
  BOOST_CHECK(combinator.Combine(spectrum_list, spectrum_out) == -2);
}

BOOST_AUTO_TEST_SUITE_END()
