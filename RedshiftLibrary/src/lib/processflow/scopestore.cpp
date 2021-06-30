#include "RedshiftLibrary/processflow/scopestore.h"

using namespace NSEpic;
std::string CScopeStore::GetCurrentScopeName() const
{
  //TODO ugly, more elegant ways to do this
  std::string n;

  TScopeStack::const_iterator it;

  if( m_ScopeStack.size() == 0 )
    return n;

  n = m_ScopeStack[0];
  it = m_ScopeStack.begin();
  it++;

  for( ; it != m_ScopeStack.end(); it++ ) {
    n.append(".");
    n.append((*it));
  }

  return n;
}

std::string CScopeStore::GetScopedName( const std::string& name ) const
{
    std::string scopedName = GetCurrentScopeName();

    if( ! scopedName.empty() ) {
        scopedName.append(".");
    }

    scopedName.append( name );

    return scopedName;

}
