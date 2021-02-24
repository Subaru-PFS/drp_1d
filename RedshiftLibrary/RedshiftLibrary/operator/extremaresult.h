#ifndef _REDSHIFT_OPERATOR_EXTREMARESULT_
#define _REDSHIFT_OPERATOR_EXTREMARESULT_

#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/statistics/pdfcandidatesz.h>
#include <RedshiftLibrary/statistics/pdfcandidateszresult.h>

#include <vector>
#include <memory>

namespace NSEpic
{
class CModelSpectrumResult;
class CModelContinuumFittingResult;

#include <RedshiftLibrary/operator/extremaresult.i>

  template <class T>
  class CExtremaResult : public CPdfCandidateszResult<T>
  {
 
 
  };
    
  template<> class CExtremaResult<TExtremaResult> : public CPdfCandidateszResult<TExtremaResult>
  {
  public:
    std::vector<std::shared_ptr<const CModelSpectrumResult>  > m_savedModelSpectrumResults;
    std::vector<std::shared_ptr<const CModelContinuumFittingResult>  > m_savedModelContinuumFittingResults;

    CExtremaResult<TExtremaResult>() = default;
    CExtremaResult<TExtremaResult>(const TCandidateZbyRank& zCandidates)
    {
      this->m_type="ExtremaResult";
      this->m_ranked_candidates.resize(zCandidates.size());
      int i=0;
      for (std::pair<std::string,const TCandidateZ&> cand:zCandidates)
        {
          this->m_ranked_candidates[i].first = cand.first;
          this->m_ranked_candidates[i].second = TExtremaResult(cand.second);
          i++;
        }
      this->m_savedModelContinuumFittingResults.resize(i);
      this->m_savedModelSpectrumResults.resize(i);
    }

    std::shared_ptr<const COperatorResult> getCandidate(const int& rank,const std::string& dataset) const;
    
    const std::string& getCandidateDatasetType(const std::string& dataset) const ;

    bool HasCandidateDataset(const std::string& dataset) const;
    
   

  };
  
  typedef CExtremaResult<TExtremaResult> ExtremaResult;


}

#endif
