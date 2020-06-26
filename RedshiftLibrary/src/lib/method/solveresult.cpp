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
