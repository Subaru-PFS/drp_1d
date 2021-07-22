#include "RedshiftLibrary/linemodel/linemodelsolution.h"

using namespace NSEpic;

CLineModelSolution::CLineModelSolution()
{
  this->m_type="CLineModelSolution";
}

void CLineModelSolution::fillRayIds()
{
  for (UInt32 j=0; j<Rays.size(); j++) rayId.emplace_back(Rays[j].GetID());
}
