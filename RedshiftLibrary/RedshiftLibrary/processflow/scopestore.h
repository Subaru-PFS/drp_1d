#ifndef _SCOPE_STORE_H
#define _SCOPE_STORE_H

#include <RedshiftLibrary/common/datatypes.h>

namespace NSEpic
{

/**
 * \ingroup Redshift
 */
// abstract class
class CScopeStore
{
 public:
 CScopeStore(const TScopeStack& scopeStack):
  m_ScopeStack(scopeStack)
    {

    }

  CScopeStore(CScopeStore const& other) = default;
  CScopeStore& operator=(CScopeStore const& other) = default;
  CScopeStore(CScopeStore&& other) = default;
  CScopeStore& operator=(CScopeStore&& other) = default;
  virtual ~CScopeStore() = default;
  
  

  
  //should be protected not public, certainly. Client code should only push and pop the scopestack to manipulate to store/get data into ScopeStores
  std::string GetCurrentScopeName() const;
  std::string GetScopedName( const std::string& name ) const ;

  //  std::string GetScope(const COperatorResult& result) const;

 protected:

  const TScopeStack& m_ScopeStack;
  
};

}
#endif
