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
#include "RedshiftLibrary/processflow/scopestore.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/log/log.h"

using namespace NSEpic;
std::string CScopeStore::GetCurrentScopeName() const {
  // TODO ugly, more elegant ways to do this
  std::string n;

  TScopeStack::const_iterator it;

  if (m_ScopeStack.size() == 0)
    return n;

  n = m_ScopeStack[0];
  it = m_ScopeStack.begin();
  ++it;

  for (; it != m_ScopeStack.end(); ++it) {
    n.append(".");
    n.append((*it));
  }

  return n;
}

std::string CScopeStore::getCurrentScopeNameAt(int depth) const {
  // TODO ugly, more elegant ways to do this
  std::string n;
  if (m_ScopeStack.size() < depth)
    THROWG(INTERNAL_ERROR, Formatter() << "Scope smaller than" << depth);

  n = m_ScopeStack[0];
  for (UInt8 i = 1; i < depth; i++) {
    n.append(".");
    n.append(m_ScopeStack[i]);
  }

  return n;
}

std::string CScopeStore::GetScopedName(const std::string &name) const {
  std::string scopedName = GetCurrentScopeName();

  if (!scopedName.empty()) {
    scopedName.append(".");
  }

  scopedName.append(name);

  return scopedName;
}

std::string CScopeStore::GetScopedNameAt(const std::string &name,
                                         int depth) const {
  std::string scopedName = getCurrentScopeNameAt(depth);

  if (!scopedName.empty()) {
    scopedName.append(".");
  }

  scopedName.append(name);

  return scopedName;
}
