#include <RedshiftLibrary/method/solveresult.h>

using namespace NSEpic;

CSolveResult::CSolveResult()
{

}

CSolveResult::~CSolveResult()
{

}

void CSolveResult::SetReliabilityLabel( std::string lbl )
{
    m_ReliabilityLabel = lbl;
}

void CSolveResult::SetTypeLabel( std::string lbl )
{
    m_TypeLabel = lbl;
}

void CSolveResult::getData(const std::string& name, std::string& v) const
{
  if (name.compare("Reliability") == 0)  v = m_ReliabilityLabel;
  else if (name.compare("Type") == 0)  v = m_TypeLabel;
  //else throw error
}
