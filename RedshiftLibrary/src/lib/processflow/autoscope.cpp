#include "RedshiftLibrary/processflow/autoscope.h"

using namespace NSEpic;
CAutoScope::CAutoScope( TScopeStack &scopeStack, const std::string& name ):
  m_scopeStack(scopeStack)
{
  
  m_scopeStack.push_back(name);
}

CAutoScope::~CAutoScope()
{
  m_scopeStack.pop_back();
}
