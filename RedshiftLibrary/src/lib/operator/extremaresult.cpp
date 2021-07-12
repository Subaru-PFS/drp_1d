#include "RedshiftLibrary/operator/extremaresult.h"
#include "RedshiftLibrary/operator/spectraFluxResult.h"
#include "RedshiftLibrary/operator/modelspectrumresult.h"
using namespace NSEpic;

std::shared_ptr<const COperatorResult> ExtremaResult::getCandidate(const int& rank,const std::string& dataset) const{
      if (dataset == "model_parameters")  return std::make_shared<const TExtremaResult>(this->m_ranked_candidates[rank].second);
      else if (dataset == "model")  return this->m_savedModelSpectrumResults[rank];
      // else if (dataset == "continuum")  return this->m_savedModelContinuumSpectrumResults[rank];

      else   throw GlobalException(UNKNOWN_ATTRIBUTE,"Unknown dataset");
    }
    
const std::string& ExtremaResult::getCandidateDatasetType(const std::string& dataset) const {
      if (dataset == "model_parameters")      return this->m_ranked_candidates[0].second.getType();
      else if (dataset == "model")  return this->m_savedModelSpectrumResults[0]->getType();
      //else if (dataset == "continuum")  return this->m_savedModelContinuumSpectrumResults[0]->getType();
      else   throw GlobalException(UNKNOWN_ATTRIBUTE,"Unknown dataset");
    }

bool ExtremaResult::HasCandidateDataset(const std::string& dataset) const
{
  return (dataset == "model_parameters" || dataset == "model");
}
