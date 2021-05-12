#include <RedshiftLibrary/spectrum/combination.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/spectralaxis.h>
#include <RedshiftLibrary/spectrum/spectrum.h>

#include <time.h>
#include <iostream>
#include <stdlib.h>
#include <boost/test/unit_test.hpp>

using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(Combination)

BOOST_AUTO_TEST_CASE(SpectrumCombination)
{
  CSpectrumCombination combinator;
  TFloat64List lambdas = {1000,2000,3000};
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

  std::vector< std::shared_ptr<CSpectrum> > spectrum_list;


  // not enough spectra
  BOOST_CHECK( combinator.Combine(spectrum_list, spectrum_out) == -3 );

  spectrum_list.push_back((std::shared_ptr<CSpectrum>)(new CSpectrum(spectral_axis, fluxaxis_a)));
  spectrum_list.push_back((std::shared_ptr<CSpectrum>)(new CSpectrum(spectral_axis, fluxaxis_b)));

  // good call
  BOOST_CHECK( combinator.Combine(spectrum_list, spectrum_out) == 0 );

  // different grids
  spectrum_list.push_back((std::shared_ptr<CSpectrum>)(new CSpectrum(spectral_axis_bad_1,
								     fluxaxis_bad_1)));
  BOOST_CHECK( combinator.Combine(spectrum_list, spectrum_out) == -2 );

  spectrum_list.pop_back();

  // different grid sizes
  spectrum_list.push_back((std::shared_ptr<CSpectrum>)(new CSpectrum(spectral_axis_bad_2,
								     fluxaxis_bad_2)));
  BOOST_CHECK( combinator.Combine(spectrum_list, spectrum_out) == -2 );

}

BOOST_AUTO_TEST_SUITE_END()
