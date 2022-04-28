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
#ifndef _REDSHIFT_COMMON_DATATYPES_
#define _REDSHIFT_COMMON_DATATYPES_

#include <memory>
#include <string>
#include <type_traits>
#include <vector>

#define __FILENAME__                                                           \
  (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)

namespace NSEpic {
#ifndef NULL
#define NULL (0)
#endif

typedef long long Int64;
typedef int Int32;
typedef short Int16;
typedef unsigned char UInt8;
typedef float Float32;
typedef double Float64;
typedef char Char;
typedef unsigned char Byte;
typedef const char *String;

typedef std::vector<Float64> TFloat64List;
typedef std::vector<Int64> TInt64List;
typedef std::vector<bool> TBoolList;
typedef std::vector<Int32> TInt32List;
typedef std::vector<std::string> TStringList;
typedef TStringList TScopeStack;

struct SPoint {
  SPoint() {
    X = 0.0;
    Y = 0.0;
  }

  SPoint(Float64 x, Float64 y) {
    X = x;
    Y = y;
  }
  Float64 X;
  Float64 Y;
};

typedef std::vector<SPoint> TPointList;

typedef UInt8 Mask;
typedef Float64 Redshift;
typedef Float64 Sample;
typedef std::vector<Mask> TMaskList;
typedef std::vector<Redshift> TRedshiftList;
typedef std::vector<Sample> TAxisSampleList;

// for C++11, defined since C++17,
template <class T>
constexpr typename std::add_const<T>::type &as_const(T &t) noexcept {
  return t;
};

template <class T> void as_const(const T &&) = delete;

typedef struct {
  Float64 x0;
  Float64 x1;
  Float64 x2;
} TPolynomCoeffs;

static const Int32 undefIdx = -1;

#include "RedshiftLibrary/common/errorcodes.i"

} // namespace NSEpic

#endif
