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
#ifndef _REDSHIFT_COMMON_INDEX_
#define _REDSHIFT_COMMON_INDEX_
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/log/log.h"
#include <iostream>
#include <vector>
namespace NSEpic {

/**
 * \ingroup Redshift
 * Templated INDEX manipulation class
 */
template <typename T> class CIndexing {

public:
  static Int32 getIndex(const std::vector<T> &list, const T z) {
    typename std::vector<T>::const_iterator itr =
        std::find(list.begin(), list.end(), z);
    if (itr == list.end())
      THROWG(INTERNAL_ERROR, Formatter() << "Could not find index for " << z);

    return (itr - list.begin());
  }

  // getIndex in orded_values corresponding to value:
  // value[index] can be equal or smaller than Z
  static bool getClosestLowerIndex(const std::vector<T> &ordered_values,
                                   const T &value, Int32 &i_min) {
    if (value < ordered_values.front()) {
      return false;
    }
    typename std::vector<T>::const_iterator it_min =
        std::lower_bound(ordered_values.begin(), ordered_values.end(), value);
    if (it_min == ordered_values.end())
      it_min--;
    if (*it_min > value)
      it_min = it_min - 1;

    i_min = it_min - ordered_values.begin();
    return true;
  }

  // the closest at left or right, at epsilon
  static Int32 getCloserIndex(const std::vector<T> &ordered_values,
                              const T &value) {
    typename std::vector<T>::const_iterator it =
        std::lower_bound(ordered_values.begin(), ordered_values.end(), value);

    // check if referring to the last element
    if (it == ordered_values.end())
      it = it - 1;

    else if (it != ordered_values.begin()) {
      // compare diff between value and it and it-1 --> select the it that gives
      // the minimal difference
      if (std::abs(*it - value) > std::abs(*(it - 1) - value))
        it = it - 1;
    }

    Int32 i_min = it - ordered_values.begin();
    return i_min;
  }
  // value[index] can be equal or higher than Z
  static bool getClosestUpperIndex(const std::vector<T> &ordered_values,
                                   const T &value, Int32 &i) {
    if (value > ordered_values.back()) {
      return false;
    }
    typename std::vector<T>::const_iterator it =
        std::lower_bound(ordered_values.begin(), ordered_values.end(), value);

    i = it - ordered_values.begin();
    return true;
  }
};
typedef CIndexing<Int32> TInt32Index;
typedef CIndexing<Float64> TFloat64Index;

} // namespace NSEpic

#endif
