#include "RedshiftLibrary/operator/tplCombinationExtremaResult.h"

#include "RedshiftLibrary/linemodel/elementlist.h"
#include "RedshiftLibrary/linemodel/modelfittingresult.h"
#include "RedshiftLibrary/operator/spectraFluxResult.h"
using namespace NSEpic;


std::shared_ptr<const COperatorResult> TplCombinationExtremaResult::getCandidate(const int& rank,const std::string& dataset) const
{ 
  if (dataset == "model_parameters" || dataset == "modeltplCombination_parameters")  return std::make_shared<const TTplCombinationResult>(this->m_ranked_candidates[rank].second);
  else if (dataset == "model")  return this->m_savedModelSpectrumResults[rank];
  else   throw GlobalException(UNKNOWN_ATTRIBUTE,"Unknown dataset");
}
    
const std::string& TplCombinationExtremaResult::getCandidateDatasetType(const std::string& dataset) const 
{
      if (dataset == "model_parameters" || dataset == "modeltplCombination_parameters")      return this->m_ranked_candidates[0].second.getType();
      else if (dataset == "model")  return this->m_savedModelSpectrumResults[0]->getType();
      else   throw GlobalException(UNKNOWN_ATTRIBUTE,"Unknown dataset");
}

bool TplCombinationExtremaResult::HasCandidateDataset(const std::string& dataset) const
{
  return (dataset == "model_parameters" || dataset == "model" || dataset == "modeltplCombination_parameters");
}
