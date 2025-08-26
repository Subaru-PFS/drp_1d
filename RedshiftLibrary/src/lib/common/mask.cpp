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
#include <cmath>

#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/common/mask.h"

using namespace NSEpic;

CMask CMask::operator&(const CMask &other) const {
  if (GetMasksCount() != other.GetMasksCount())
    THROWG(ErrorCode::INTERNAL_ERROR, Formatter()
                                          << "parameter has not the same size: "
                                          << other.GetMasksCount()
                                          << " instead of " << GetMasksCount());

  CMask result(GetMasksCount());
  std::transform(m_Mask.begin(), m_Mask.end(), other.m_Mask.begin(),
                 result.m_Mask.begin(), std::bit_and());
  return result;
}

/**
 *
 */
CMask &CMask::operator&=(const CMask &other) {
  *this = *this & other;
  return *this;
}

/**
 *
 */
Float64 CMask::ComputeOverlapFraction(const CMask &other) const {
  if (GetMasksCount() != other.GetMasksCount())
    THROWG(ErrorCode::INTERNAL_ERROR, Formatter()
                                          << "parameter has not the same size: "
                                          << other.GetMasksCount()
                                          << " instead of " << GetMasksCount());

  Float64 selfRate = GetUnMaskedSampleCount();
  Float64 otherRate = other.GetUnMaskedSampleCount();

  if (selfRate == 0.0)
    return 0;

  return otherRate / selfRate;
}

/**
 *
 */
Float64 CMask::IntersectAndComputeOverlapFraction(const CMask &other) const {
  return ComputeOverlapFraction(*this & other);
}

CMask CMask::extract(Int32 startIdx, Int32 endIdx) const {

  if (!m_Mask.size())
    return CMask();
  return CMask(
      TMaskList(m_Mask.begin() + startIdx, m_Mask.begin() + endIdx + 1));
}
