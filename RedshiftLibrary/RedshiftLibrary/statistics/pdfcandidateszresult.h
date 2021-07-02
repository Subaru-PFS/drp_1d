#ifndef _REDSHIFT_STATISTICS_PDFCANDIDATESZRESULT_
#define _REDSHIFT_STATISTICS_PDFCANDIDATESZRESULT_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/processflow/result.h"
#include "RedshiftLibrary/statistics/pdfcandidatesz.h"

//#include <cmath>
//#include <map>
#include <ostream>
#include <string>

namespace NSEpic
{
template <class T>
class CPdfCandidateszResult : public COperatorResult
{

public:

  CPdfCandidateszResult(Int32 optMethod=0): m_optMethod(optMethod) {
    this->m_type="PdfCandidatesZResult";
  };
  


  //rule of 5 defaults
  CPdfCandidateszResult(const CPdfCandidateszResult & ) = default;
  CPdfCandidateszResult(CPdfCandidateszResult && ) = default;
  CPdfCandidateszResult & operator=(const CPdfCandidateszResult & ) = default;   
  CPdfCandidateszResult & operator=(CPdfCandidateszResult && ) = default;   
  virtual ~CPdfCandidateszResult() = default;
  
  TStringList GetIDs() const
  {
    TStringList ids;
    ids.reserve(m_ranked_candidates.size());
    for (auto c: m_ranked_candidates) ids.push_back(c.first);
    return ids;
  }
  Int32                       m_optMethod; //0: direct integration, 1:gaussian fit

  Int32 size() const{return m_ranked_candidates.size();}
  T getCandidateDataset(int rank,std::string dataset)
  {
    return m_ranked_candidates[rank].second;
  }
  
  std::vector<std::pair<std::string, T>> m_ranked_candidates;

};

  typedef CPdfCandidateszResult<TCandidateZ> PdfCandidatesZResult;
}

#endif
