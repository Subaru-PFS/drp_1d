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
#ifndef _REDSHIFT_SPECTRUM_ITERATOR_
#define _REDSHIFT_SPECTRUM_ITERATOR_

#include <cstddef>  // For std::ptrdiff_t
#include <iterator> // For std::forward_iterator_tag
#include <memory>

#include "RedshiftLibrary/common/datatypes.h"

namespace NSEpic {

class CSpectraGlobalIndex {

public:
  struct Iterator {
    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = Int32;
    using pointer = std::shared_ptr<Int32>; // or also value_type*
    using reference = Int32 &;              // or also value_type&

    Iterator(pointer ptr) : m_ptr(ptr) {}

    reference operator*() const { return *m_ptr; }
    pointer operator->() { return m_ptr; }
    Iterator &operator++() {
      (*m_ptr)++;
      return *this;
    }
    Iterator operator++(int) {
      Iterator tmp = *this;
      ++(*this);
      return tmp;
    }
    friend bool operator==(const Iterator &a, const Iterator &b) {
      return *(a.m_ptr) == *(b.m_ptr);
    };
    friend bool operator!=(const Iterator &a, const Iterator &b) {
      return *(a.m_ptr) != *(b.m_ptr);
    };

  private:
    pointer m_ptr;
  };

  Iterator begin() {
    *m_curObs = 0;
    return Iterator(m_curObs);
  }
  Iterator end() { return m_endObs; }

  Int32 getCurObs() const { return *m_curObs; }

  void reset() { *m_curObs = 0; }

  CSpectraGlobalIndex(Int32 nbObs) {
    m_curObs = std::make_shared<Int32>(0);
    m_endObs = std::make_shared<Int32>(nbObs);
  }

private:
  std::shared_ptr<Int32> m_curObs;
  std::shared_ptr<Int32> m_endObs;
  // Int32 m_nbObs;
};

} // namespace NSEpic
#endif
