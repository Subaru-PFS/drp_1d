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
#include <stdexcept>

#include <gsl/gsl_const_mksa.h>

#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/photometry/photometricband.h"

using namespace NSEpic;

CPhotometricBand::CPhotometricBand(TFloat64List trans, TFloat64List lambda)
    : m_transmission(std::move(trans)), m_lambda(std::move(lambda)) {
  // check all sizes are the same
  if (m_transmission.size() != m_lambda.size()) {
    THROWG(ErrorCode::INTERNAL_ERROR, "Transmission and "
                                      "wavelength have not the same size");
  }

  // initialize min lambda
  m_minLambda = *std::min_element(m_lambda.cbegin(), m_lambda.cend());

  // Angstrom to Herz conversion
  const Float64 c = GSL_CONST_MKSA_SPEED_OF_LIGHT * 1e10; // angstrom s-1
  m_freq.reserve(m_lambda.size());
  std::transform(m_lambda.cbegin(), m_lambda.cend(), back_inserter(m_freq),
                 [c](Float64 lambda) { return c / lambda; });

  // compute the integral of transmission in Hz
  m_transmision_sum = 0.0;
  for (auto trans = m_transmission.cbegin() + 1, E = m_transmission.cend(),
            freq = m_freq.cbegin() + 1, lbda = m_lambda.cbegin() + 1;
       trans != E; ++trans, ++freq, ++lbda) {
    // lambda factor below is to account for photon counting detector
    Float64 trapezArea = (trans[0] * lbda[0] + trans[-1] * lbda[-1]) / 2.0;
    trapezArea *= (freq[-1] - freq[0]);
    m_transmision_sum += trapezArea;
  }
}

Float64 CPhotometricBand::IntegrateFlux(const TFloat64List &inFlux) const {

  if (inFlux.size() != m_transmission.size())
    THROWG(ErrorCode::INTERNAL_ERROR, "Flux and "
                                      "transmission sizes are different");

  // convert flux from Flambda in Erg/s/cm2/Angstrom to Fnu in Erg/s/cm2/Hz
  const Float64 c = GSL_CONST_MKSA_SPEED_OF_LIGHT * 1e10; // angstrom s-1
  TAxisSampleList inFlux_nu(m_transmission.size());
  std::transform(inFlux.cbegin(), inFlux.cend(), m_lambda.cbegin(),
                 inFlux_nu.begin(), [c](Float64 Flambda, Float64 lambda) {
                   return Flambda * lambda * lambda / c;
                 });

  // multiply by transmission
  TAxisSampleList outFlux_nu(m_transmission.size());
  std::transform(inFlux_nu.cbegin(), inFlux_nu.cend(), m_transmission.cbegin(),
                 outFlux_nu.begin(), std::multiplies<Float64>{});

  Float64 sum = 0.0;
  for (auto flux = outFlux_nu.cbegin() + 1, E = outFlux_nu.cend(),
            freq = m_freq.cbegin() + 1, lbda = m_lambda.cbegin() + 1;
       flux != E; ++flux, ++freq, ++lbda) {
    // lbda factor below is to account for photon counting detector
    Float64 trapezArea = (flux[0] * lbda[0] + flux[-1] * lbda[-1]) / 2.0;
    trapezArea *= (freq[-1] - freq[0]);
    sum += trapezArea;
  }

  // convert the integrated flux to flux density per Hz
  sum /= m_transmision_sum;

  return sum;
}

TStringList CPhotBandCatalog::GetNameList() const {
  TStringList names;
  names.reserve(size());
  for (auto const &band : *this) {
    names.push_back(band.first);
  }
  return (names);
}

TStringList CPhotBandCatalog::GetNameListSortedByLambda() const {
  TStringList names = GetNameList();
  std::sort(names.begin(), names.end(),
            [this](const std::string &lhs, const std::string &rhs) {
              return at(lhs).GetMinLambda() < at(rhs).GetMinLambda();
            });
  return names;
}
