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
  Float64 lambdas[] = {1000,2000,3000};
  CSpectrumSpectralAxis spectral_axis(lambdas, 3);

  Float64 lambdas_bad_1[] = {1001, 2000, 3000};
  CSpectrumSpectralAxis spectral_axis_bad_1(lambdas_bad_1, 3);

  Float64 lambdas_bad_2[] = {1000, 2000, 3000, 4000};
  CSpectrumSpectralAxis spectral_axis_bad_2(lambdas_bad_2, 4);

  Float64 flux_bad_1[] = {1.2, 2.3, 3.4};
  CSpectrumFluxAxis fluxaxis_bad_1(flux_bad_1, 3);

  Float64 flux_bad_2[] = {1.2, 2.3, 3.4, 4.5};
  CSpectrumFluxAxis fluxaxis_bad_2(flux_bad_2, 4);

  Float64 flux_a[] = {1.2, 2.3, 3.4};
  CSpectrumFluxAxis fluxaxis_a(flux_a, 3);

  Float64 flux_b[] = {2, 3, 4};
  CSpectrumFluxAxis fluxaxis_b(flux_b, 3);

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
