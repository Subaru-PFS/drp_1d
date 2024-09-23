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
#ifndef _REDSHIFT_CONTINUUM_COMPONENT_
#define _REDSHIFT_CONTINUUM_COMPONENT_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"

#include <unordered_map>

namespace continuumcomponent_test { // boost_test_suite
// all boost_auto_test_case that use private method
class tautology_test;
} // namespace continuumcomponent_test

namespace NSEpic {

class TContinuumComponent {
public:
  TContinuumComponent() : m_type(EContinuumComponent::null) {}
  TContinuumComponent(const std::string &type) { set(type); }
  ~TContinuumComponent() = default;
  TContinuumComponent(const TContinuumComponent &other) = default;
  TContinuumComponent &operator=(const TContinuumComponent &other) = default;
  TContinuumComponent(TContinuumComponent &&other) noexcept = default;
  TContinuumComponent &
  operator=(TContinuumComponent &&other) noexcept = default;

  bool operator==(const TContinuumComponent &other) const {
    return m_type == other.m_type;
  }

  bool operator!=(const TContinuumComponent &other) const {
    return m_type != other.m_type;
  }

  operator std::string() const {
    if (m_type == EContinuumComponent::null) {
      return "null";
    }
    return enumToString.at(m_type);
  }

  void set(const std::string &type) {
    auto it = stringToEnum.find(type);
    if (it != stringToEnum.end()) {
      m_type = it->second; // Assign the corresponding enum value
    } else {
      THROWG(ErrorCode::INTERNAL_ERROR,
             Formatter() << "TContinuumComponent::" << __func__
                         << ": invalid type " << type);
    }
  }

  bool isTplFit() const {
    checkInitialized();
    return m_type == EContinuumComponent::tplFit;
  }
  bool isTplFitAuto() const {
    checkInitialized();
    return m_type == EContinuumComponent::tplFitAuto;
  }
  bool isFromSpectrum() const {
    checkInitialized();
    return m_type == EContinuumComponent::fromSpectrum;
  }
  bool isNoContinuum() const {
    checkInitialized();
    return m_type == EContinuumComponent::noContinuum;
  }
  bool isPowerLaw() const {
    checkInitialized();
    return m_type == EContinuumComponent::powerLaw;
  }
  bool isTplFitxxx() const { return isTplFit() || isTplFitAuto(); }
  bool isContinuumFit() const { return isTplFitxxx() || isPowerLaw(); }

  bool isInitialized() const { return m_type != EContinuumComponent::null; }

  void checkInitialized() const {
    if (!isInitialized()) {
      THROWG(ErrorCode::INTERNAL_ERROR,
             Formatter() << "Uninitialized continuum component");
    }
  }

private:
  friend class continuumcomponent_test::tautology_test;
  enum class EContinuumComponent {
    null,
    tplFit,
    tplFitAuto,
    noContinuum,
    fromSpectrum,
    powerLaw
  };

  EContinuumComponent m_type;
  static const std::unordered_map<std::string, EContinuumComponent>
      stringToEnum;
  static const std::unordered_map<EContinuumComponent, std::string>
      enumToString;
};
} // namespace NSEpic

#endif
