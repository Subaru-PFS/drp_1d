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

#include "RedshiftLibrary/processflow/scopestack.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"

using namespace NSEpic;

void CScopeStack::push_back(const std::string &value, ScopeType type) {
  set_type_level(type);
  m_scopeStack.push_back(value);
}

void CScopeStack::push_back(std::string &&value, ScopeType type) {
  set_type_level(type);
  m_scopeStack.push_back(std::move(value));
}

void CScopeStack::pop_back() {
  m_scopeStack.pop_back();
  reset_type_level();
}

ScopeType CScopeStack::get_current_type() const {
  size_t level = size();
  auto const it = find_if(m_type_level.begin(), m_type_level.end(),
                          [level](std::pair<ScopeType, size_t> const &p) {
                            return p.second + 1 == level;
                          });
  if (it == m_type_level.end())
    return ScopeType::UNDEFINED;
  return it->first;
}

void CScopeStack::set_type_level(ScopeType type) {
  if (type == ScopeType::UNDEFINED)
    return;

  // check rules
  switch (type) {
  case ScopeType::SPECTRUMMODEL:
    // SPECTRUMMODEL should be before STAGE (ie no STAGE set)
    if (has_type(ScopeType::STAGE))
      THROWG(SCOPESTACK_ERROR,
             Formatter() << "cannot set type " << ScopeType::SPECTRUMMODEL
                         << " if " << ScopeType::STAGE << " is already set");
    break;
  case ScopeType::STAGE:
    if (has_type(ScopeType::METHOD))
      THROWG(SCOPESTACK_ERROR,
             Formatter() << "cannot set type " << ScopeType::STAGE << " if "
                         << ScopeType::METHOD << " is already set");
    break;
  case ScopeType::METHOD:
    // METHOD should follow stage
    if (!has_type(ScopeType::STAGE) ||
        m_type_level.at(ScopeType::STAGE) + 1 != m_scopeStack.size())
      THROWG(SCOPESTACK_ERROR,
             Formatter() << "cannot set type " << ScopeType::METHOD
                         << " if current scope is not " << ScopeType::STAGE);
    break;
  }

  // insert (type, rank)
  const auto &[it, success] = m_type_level.insert({type, m_scopeStack.size()});
  if (!success)
    THROWG(SCOPESTACK_ERROR, Formatter()
                                 << "scope type " << type << "already set");
}

void CScopeStack::reset_type_level() {
  for (auto it = m_type_level.begin(); it != m_type_level.end();) {
    auto const &[type, level] = *it;
    if (level >= m_scopeStack.size())
      it = m_type_level.erase(it);
    else
      ++it;
  }
}