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

#include "RedshiftLibrary/photometry/photometricband.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/log/log.h"

#include <stdexcept>
#include <algorithm>

using namespace NSEpic;

CPhotometricBand::CPhotometricBand(TFloat64List trans, TFloat64List lambda)
    : m_transmission(std::move(trans)), m_lambda(std::move(lambda)) {
  // check all sizes are the same
  if (m_transmission.size() != m_lambda.size()) {
    throw GlobalException(INTERNAL_ERROR,
                          "CPhotometryBand::CPhotometryBand: transmission and "
                          "wavelength have not the same size");
  }

  // initialize min lambda
  m_minLambda = *std::min_element(m_lambda.cbegin(), m_lambda.cend());
}

Float64 CPhotometricBand::IntegrateFlux(const TFloat64List &inFlux) const {

  if (inFlux.size() != m_transmission.size())
    throw GlobalException(INTERNAL_ERROR,
                          "CPhotometryBand::IntegrateFlux: flux and "
                          "transmission sizes are different");

  TAxisSampleList outFlux(m_transmission.size());
  std::transform(inFlux.cbegin(), inFlux.cend(), m_transmission.cbegin(),
                 outFlux.begin(), std::multiplies<Float64>{});

  Float64 sum = 0.0;
  for (auto flux = outFlux.cbegin() + 1, E = outFlux.cend(),
            lambda = m_lambda.cbegin()+1;
       flux != E; ++flux, ++lambda) {
    Float64 trapezArea = (*flux + *(flux-1)) / 2.0;
    trapezArea *= (*lambda - *(lambda-1));
    sum += trapezArea;
  }

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
