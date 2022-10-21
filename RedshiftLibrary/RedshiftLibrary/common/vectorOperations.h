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
#include <vector>
namespace NSEpic {
template <typename T>
static void insertAroundIndex(std::vector<T> &entity, Int32 idx,
                              Int32 count_smaller, Int32 count_higher,
                              T defaultVal) {
  entity.insert(entity.begin() + idx + 1, count_higher, defaultVal);
  entity.insert(entity.begin() + idx, count_smaller, defaultVal);
}

static void insertVectorAroundIndex(TFloat64List &entity, Int32 idx,
                                    const TInt32List &indices,
                                    Int32 count_smaller, Int32 count_higher,
                                    const TFloat64List &vect) {
  Float64 defaultValue = NAN;
  insertAroundIndex(entity, idx, count_smaller, count_higher,
                    defaultValue); // corresponds to a resize
  Int32 start =
      std::min(idx, indices.front()); // in case there is a common on the border
  Int32 stop = std::max(start + Int32(vect.size()) - 1, indices.back());
  for (Int32 i = start; i <= stop; i++)
    entity[i] = vect[i - start];
}

static void insertIntoRedshiftGrid(TFloat64List &entity, Int32 idx,
                                   const TInt32List &indices,
                                   Int32 &count_smaller, Int32 &count_higher,
                                   const TFloat64List &vect) {

  Int32 i = CIndexing<Float64>::getIndex(vect, entity[idx]);

  // count nb of duplicates below idx
  auto count_duplicates_low =
      std::count_if(indices.begin(), indices.end(),
                    [&](const Float64 &val) { return val < idx; });
  auto count_duplicates_high =
      std::count_if(indices.begin(), indices.end(),
                    [&](const Float64 &val) { return val > idx; });
  count_smaller =
      i - count_duplicates_low; // substract all indices smaller than i
  count_higher = vect.size() - i - count_duplicates_high - 1;
  insertVectorAroundIndex(entity, idx, indices, count_smaller, count_higher,
                          vect);
  return;
}

// insert a sub-array (src) into dest, at a specific position while avoinding
// duplicates
static void insertWithDuplicates(TFloat64List &dest, Int32 pos,
                                 const TFloat64List &src, Int32 ndup) {
  auto first = src.begin();
  auto last = src.end() - ndup;
  Int32 ninsert = last - first;
  dest.insert(dest.begin() + pos, first, last);
  std::copy(last, src.end(), dest.begin() + pos + ninsert);
  // copy into dest the missing values from src
}

// todo: check if squaching already filled values keeps things valid
template <typename T>
static void insertWithDuplicates(std::vector<T> &dest, Int32 pos,
                                 std::size_t count, const T &value,
                                 Int32 ndup) {
  Int32 ninsert = count - ndup;
  dest.insert(dest.begin() + pos, ninsert, value);
  std::fill(dest.begin() + pos + ninsert, dest.begin() + pos + count,
            value); // replace existing values with the defaultV
}

} // namespace NSEpic
#endif
