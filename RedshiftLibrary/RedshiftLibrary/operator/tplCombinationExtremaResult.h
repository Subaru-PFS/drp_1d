#ifndef _REDSHIFT_TPLCOMBINATION_EXTREMARESULT_
#define _REDSHIFT_TPLCOMBINATION_EXTREMARESULT_

#include "RedshiftLibrary/operator/extremaresult.h"
#include "RedshiftLibrary/processflow/result.h"

namespace NSEpic
{
  class CModelSpectrumResult;
  class CSpectraFluxResult;

  #include "RedshiftLibrary/operator/tplCombinationExtremaResult.i"
  
  template<> class CExtremaResult<TTplCombinationResult>: public CPdfCandidateszResult<TTplCombinationResult>
  {
  public:
    std::vector<std::shared_ptr<const CModelSpectrumResult>>   m_savedModelSpectrumResults;

    CExtremaResult<TTplCombinationResult>() = default;
    
    CExtremaResult<TTplCombinationResult>(const TCandidateZbyRank& zCandidates)
    {
      this->m_type = "TplCombinationExtremaResult";

      for (std::pair<std::string,const TCandidateZ&> cand:zCandidates)
      {
	        this->m_ranked_candidates.push_back(std::make_pair<std::string,TTplCombinationResult>(std::string(cand.first),TTplCombinationResult(cand.second)));
      }
      this->Resize(zCandidates.size());
    }

    void Resize(Int32 size)
    {
      m_savedModelSpectrumResults.resize(size);
    }

    std::shared_ptr<const COperatorResult> getCandidate(const int& rank,const std::string& dataset) const;
    
    const std::string& getCandidateDatasetType(const std::string& dataset) const ;

    bool HasCandidateDataset(const std::string& dataset) const;
    
  };

  typedef CExtremaResult<TTplCombinationResult> TplCombinationExtremaResult;

}

#endif
