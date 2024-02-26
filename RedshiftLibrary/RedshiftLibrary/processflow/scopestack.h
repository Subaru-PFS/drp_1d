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
#ifndef _SCOPE_STACK_H
#define _SCOPE_STACK_H

#include <iterator>
#include <map>
#include <ostream>
#include <string>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/defaults.h"

namespace NSEpic {

/**
 * \ingroup Redshift
 */
enum class ScopeType { UNDEFINED, SPECTRUMMODEL, STAGE, METHOD };
inline const std::map<ScopeType, std::string> scope_type_str{
    {ScopeType::UNDEFINED, "UNDEFINED"},
    {ScopeType::SPECTRUMMODEL, "SPECTRUMMODEL"},
    {ScopeType::STAGE, "STAGE"},
    {ScopeType::METHOD, "METHOD"}};

inline std::ostream &operator<<(std::ostream &os, const ScopeType &type) {
  os << scope_type_str.at(type);
  return os;
};

class CScopeStack {
public:
  CScopeStack() = default;
  CScopeStack(const TScopeStack &scope) : m_scopeStack(scope){};
  CScopeStack(TScopeStack &&scope) : m_scopeStack(std::move(scope)){};
  CScopeStack(std::initializer_list<std::string> init) : m_scopeStack(init){};

  auto size() const { return m_scopeStack.size(); };
  auto empty() const { return m_scopeStack.empty(); };

  // only const accessors
  // mutation allowed only with push/pop
  auto begin() const { return m_scopeStack.begin(); };
  auto end() const { return m_scopeStack.end(); };
  auto &front() const { return m_scopeStack.front(); };
  auto &back() const { return m_scopeStack.back(); };
  auto &operator[](std::size_t pos) const { return m_scopeStack[pos]; };
  auto &at(std::size_t pos) const { return m_scopeStack.at(pos); };

  void push_back(const std::string &value,
                 ScopeType type = ScopeType::UNDEFINED);
  void push_back(std::string &&value, ScopeType type = ScopeType::UNDEFINED);
  void pop_back();
  bool has_type(ScopeType type) const {
    return m_type_level.find(type) != m_type_level.end();
  };
  const std::string &get_type_value(ScopeType type) const {
    return m_scopeStack.at(m_type_level.at(type));
  };
  size_t get_type_level(ScopeType type) const { return m_type_level.at(type); };
  ScopeType get_current_type() const;

private:
  void set_type_level(ScopeType type);
  void reset_type_level();

  TScopeStack m_scopeStack;
  std::map<ScopeType, size_t> m_type_level;
};

} // namespace NSEpic
#endif
