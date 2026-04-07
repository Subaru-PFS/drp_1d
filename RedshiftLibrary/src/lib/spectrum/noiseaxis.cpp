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
#include <cmath>

#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/spectrum/noiseaxis.h"

using namespace NSEpic;
using namespace std;

const TBoolList CSpectrumNoiseAxis::checkNoise() const {
  TBoolList isValid(m_Samples.size(), true);
  for (std::size_t i = 0; i < m_Samples.size(); i++) {
    Float64 err2 = 1 / (m_Samples[i] * m_Samples[i]);
    if (m_Samples[i] < DBL_MIN || std::isnan(m_Samples[i]) ||
        std::isinf(m_Samples[i]) || m_Samples[i] != m_Samples[i] ||
        std::isinf(err2) || std::isnan(err2))
      isValid[i] = false;
  }
  return isValid;
}

CSpectrumNoiseAxis &
CSpectrumNoiseAxis::operator+=(CSpectrumNoiseAxis const &other) {
  if (other.GetSamplesCount() != GetSamplesCount())
    THROWG(ErrorCode::INTERNAL_ERROR,
           "Cannot sum noise axis of different sizes");
  std::transform(m_Samples.cbegin(), m_Samples.cend(), other.m_Samples.cbegin(),
                 m_Samples.begin(),
                 [](Float64 a, Float64 b) { return std::sqrt(a * a + b * b); });
  return *this;
}
