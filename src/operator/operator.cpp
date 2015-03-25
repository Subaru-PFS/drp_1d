#include <epic/redshift/operator/operator.h>

using namespace NSEpic;
using namespace std;

COperator::COperator()
{

}

COperator::~COperator()
{

}

const COperator::TStatusList& COperator::GetStatus() const
{
    return m_Status;
}

const TFloat64List& COperator::GetOverlap() const
{
    return m_Overlap;
}

const TFloat64List& COperator::GetResults() const
{
    return m_Result;
}
