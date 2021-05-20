#include <RedshiftLibrary/method/reliabilityresult.h>

#include <RedshiftLibrary/processflow/context.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <RedshiftLibrary/log/log.h>


using namespace NSEpic;

/**
 * \brief Empty constructor.
 **/
CReliabilityResult::CReliabilityResult():
  CSolveResult()
{
}

void CReliabilityResult::getData(const std::string& name, std::string& v) const
{
  v = m_ReliabilityLabel;
}
