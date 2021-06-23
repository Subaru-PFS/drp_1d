#ifndef _REDSHIFT_LINEMODEL_LINEMODELEXTREMARESULT_
#define _REDSHIFT_LINEMODEL_LINEMODELEXTREMARESULT_

#include <RedshiftLibrary/operator/extremaresult.h>
#include <RedshiftLibrary/ray/catalog.h>
#include <RedshiftLibrary/linemodel/linemodelsolution.h>

//#include <RedshiftLibrary/operator/modelspectrumresult.h>
//#include <RedshiftLibrary/linemodel/modelfittingresult.h>
//#include <RedshiftLibrary/operator/modelcontinuumfittingresult.h>
//#include <RedshiftLibrary/linemodel/modelrulesresult.h>
//#include <RedshiftLibrary/operator/spectraFluxResult.h>
#include <RedshiftLibrary/linemodel/continuummodelsolution.h>
#include <RedshiftLibrary/processflow/result.h>



namespace NSEpic
{
class CModelSpectrumResult;
class CModelFittingResult;
class CModelContinuumFittingResult;
class CModelRulesResult;
class CSpectraFluxResult;
class CLineModelElementList;

  class CLineModelResult;
  #include <RedshiftLibrary/linemodel/linemodelextremaresult.i>
  
  template<> class CExtremaResult<TLineModelResult>: public CPdfCandidateszResult<TLineModelResult>
  {
  public:
    std::vector<TFloat64List> ExtendedRedshifts;    // z range around extrema

    std::vector<std::shared_ptr<const CModelFittingResult>  > m_savedModelFittingResults;
    std::vector<std::shared_ptr<const CModelRulesResult>  > m_savedModelRulesResults;
    std::vector<std::shared_ptr<const CSpectraFluxResult>  > m_savedModelContinuumSpectrumResults;

    //Yes there is a repetition of CExtremaResult<TExtremaResult> variables, because we inherit from CExtremaResult<TLineModelResult>, TODO find a more satisfactory solution
    std::vector<std::shared_ptr<const CModelSpectrumResult>  > m_savedModelSpectrumResults;
    std::vector<std::shared_ptr<const CModelContinuumFittingResult>  > m_savedModelContinuumFittingResults;


    CExtremaResult<TLineModelResult>() = default;
    
    CExtremaResult<TLineModelResult>(const TCandidateZbyRank& zCandidates)
    {
      this->m_type = "LineModelExtremaResult";
      int i=0;
      for (std::pair<std::string,const TCandidateZ&> cand:zCandidates)
        {
	  this->m_ranked_candidates.push_back(std::make_pair<std::string,TLineModelResult>(std::string(cand.first),TLineModelResult(cand.second)));
        }
      this->Resize(zCandidates.size());
    }

    void Resize(Int32 size)
    {
      //CExtremaResult::Resize(size);
    
      m_savedModelFittingResults.resize(size);
      m_savedModelRulesResults.resize(size);
      m_savedModelContinuumSpectrumResults.resize(size);
      m_savedModelSpectrumResults.resize(size);
      m_savedModelContinuumFittingResults.resize(size);
    }

    void setCandidateFromContinuumSolution(int rank,CContinuumModelSolution cms)
    {
      m_ranked_candidates[rank].second = TLineModelResult(cms);
    }

    std::shared_ptr<const COperatorResult> getCandidate(const int& rank,const std::string& dataset) const;
    
    const std::string& getCandidateDatasetType(const std::string& dataset) const ;

    bool HasCandidateDataset(const std::string& dataset) const;
    
  };

  typedef CExtremaResult<TLineModelResult> LineModelExtremaResult;

}

#endif
