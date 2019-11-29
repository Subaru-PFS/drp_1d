#include <RedshiftLibrary/processflow/result.h>

using namespace NSEpic;

COperatorResult::COperatorResult()
{

}

COperatorResult::~COperatorResult()
{

}

void COperatorResult::SaveJSON(const CDataStore& store, std::ostream& stream) const
{
  // does nothing, -> no need to cast COperatorResult to LineModelResult in COperatorResultStore::SaveAllResults
}

void COperatorResult::SetReliabilityLabel( std::string lbl )
{
    m_ReliabilityLabel = lbl;
}

void COperatorResult::SetTypeLabel( std::string lbl )
{
    m_TypeLabel = lbl;
}
