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
#ifndef _REDSHIFT_COMMON_VECTOROPERATIONS_
#define _REDSHIFT_COMMON_VECTOROPERATIONS_

#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/common/indexing.h"
#include "RedshiftLibrary/common/size.h"
#include <vector>

namespace NSEpic {

// insert source into destination with ndup overlaping elements
template <typename T>
inline void insertWithDuplicates(std::vector<T> &dest, Int32 pos,
                                 std::vector<T> src, Int32 ndup) {
  const auto src_first = src.cbegin();
  const auto src_end = src.cend();
  const auto src_last_insert = src.cend() - ndup;
  Int32 ninsert = src_last_insert - src_first;
  dest.insert(dest.begin() + pos, src_first, src_last_insert);
  std::copy(src_last_insert, src_end, dest.begin() + pos + ninsert);
}

// insert constant value into destination with ndup overlaping elements
template <typename T>
inline void insertWithDuplicates(std::vector<T> &dest, Int32 pos,
                                 std::size_t count, const T &value,
                                 Int32 ndup) {
  Int32 ninsert = count - ndup;
  dest.insert(dest.begin() + pos, ninsert, value);
  std::fill(dest.begin() + pos + ninsert, dest.begin() + pos + count,
            value); // replace existing values with the defaultV
}

template <typename T>
inline TInt32Pair find2DVectorMinIndexes(std::vector<std::vector<T>> vect) {
  Int32 iMin = 0;
  Int32 jMin = 0;
  T valMin = vect[0][0];
  for (Int32 i = 0; i < ssize(vect); i++) {
    for (Int32 j = 0; j < ssize(vect[i]); j++) {
      if (vect[i][j] < valMin) {
        iMin = i;
        jMin = j;
        valMin = vect[i][j];
      }
    }
  }
  return {iMin, jMin};
}

inline std::vector<Float64> createLinearInterpVector(Float64 v1, Float64 v2,
                                                     Float64 n) {
  // Creates a linear vector from value v1 to value v2 with n elements
  std::vector<Float64> result(n);
  if (n <= 1) {
    // Handle edge case: return a vector with a single value n1 if n <= 1
    if (n == 1)
      result.front() = v1;
    return result;
  }

  std::generate(result.begin(), result.end(),
                [i = 0, v1, v2, f0 = 1. / (n - 1)]() mutable {
                  auto f = f0 * i++;
                  return v1 * (1 - f) + v2 * f;
                });
  return result;
}

} // namespace NSEpic
#endif
