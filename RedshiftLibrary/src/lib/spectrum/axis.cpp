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
#include "RedshiftLibrary/spectrum/axis.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/log/log.h"
#include <algorithm>
#include <numeric>

using namespace NSEpic;
using namespace std;

CSpectrumAxis::CSpectrumAxis(const Float64 *samples, Int32 n) : m_Samples(n) {
  for (Int32 i = 0; i < n; i++) {
    m_Samples[i] = samples[i];
  }
}

CSpectrumAxis &CSpectrumAxis::operator*=(const Float64 op) {
  for (Int32 i = 0; i < m_Samples.size(); i++) {
    m_Samples[i] *= op;
  }
  return *this;
}

CSpectrumAxis &CSpectrumAxis::operator/=(const Float64 op) {
  for (Int32 i = 0; i < m_Samples.size(); i++) {
    m_Samples[i] /= op;
  }
  return *this;
}

void CSpectrumAxis::SetSize(Int32 s) { m_Samples.resize(s); }
void CSpectrumAxis::clear() {
  resetAxisProperties();
  m_Samples.clear();
}

/*
    maskedAxis is the output axis after applying the mask on the current object
*/
CSpectrumAxis
CSpectrumAxis::MaskAxis(const TFloat64List &mask) const // mask is 0. or 1.
{
  return CSpectrumAxis(maskVector(mask, m_Samples));
}

TFloat64List CSpectrumAxis::maskVector(const TFloat64List &mask,
                                       const TFloat64List &inputVector) {
  TFloat64List outputVector;
  if (mask.size() != inputVector.size()) {
    THROWG(INTERNAL_ERROR, "mask and vector sizes do not match");
  }
  Int32 sum = Int32(std::count(mask.begin(), mask.end(), 1));
  outputVector.clear();
  outputVector.reserve(sum);
  for (Int32 i = 0; i < mask.size(); i++) {
    if (mask[i] == 1.)
      outputVector.push_back(inputVector[i]);
  }
  return outputVector;
}
