#include <RedshiftLibrary/processflow/result.h>

using namespace NSEpic;

COperatorResult::COperatorResult()
{

}

COperatorResult::~COperatorResult()
{

}



void COperatorResult::SetReliabilityLabel( std::string lbl )
{
    m_ReliabilityLabel = lbl;
}

void COperatorResult::SetTypeLabel( std::string lbl )
{
    m_TypeLabel = lbl;
}
