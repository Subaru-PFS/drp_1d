#include <RedshiftLibrary/processflow/result.h>

using namespace NSEpic;

COperatorResult::COperatorResult()
{

}

COperatorResult::~COperatorResult()
{

}



Void COperatorResult::SetReliabilityLabel( std::string lbl )
{
    m_ReliabilityLabel = lbl;
}

Void COperatorResult::SetTypeLabel( std::string lbl )
{
    m_TypeLabel = lbl;
}
