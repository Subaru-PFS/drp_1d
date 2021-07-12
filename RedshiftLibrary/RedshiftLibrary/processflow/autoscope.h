#ifndef _AUTO_SCOPE_H
#define _AUTO_SCOPE_H

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/processflow/autoscope.h"

namespace NSEpic
{

  class CAutoScope {
  public:
    CAutoScope(TScopeStack &scopeStack,const std::string& name );
    ~CAutoScope();
  private:
    TScopeStack &m_scopeStack;

  };
}
#endif
