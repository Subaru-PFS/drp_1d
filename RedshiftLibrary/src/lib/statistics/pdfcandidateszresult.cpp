#include <RedshiftLibrary/statistics/pdfcandidateszresult.h>
#include <RedshiftLibrary/log/log.h>

//#include <RedshiftLibrary/processflow/context.h>


using namespace NSEpic;
using namespace std;
#include <fstream>
#include <iostream>
#include <numeric>

TStringList CPdfCandidateszResult::GetIDs() const
{
    TStringList ids;
    ids.reserve(m_ranked_candidates.size());
    for (auto c: m_ranked_candidates) ids.push_back(c.first);
    return ids;
}

TFloat64List CPdfCandidateszResult::GetRedshifts() const
{
    TFloat64List redshifts;
    redshifts.reserve(m_ranked_candidates.size());
    for (auto c: m_ranked_candidates) redshifts.push_back(c.second.Redshift);
    return redshifts;
}

TFloat64List CPdfCandidateszResult::GetDeltaZs() const
{
    TFloat64List deltaZs;
    deltaZs.reserve(m_ranked_candidates.size());
    for (auto c: m_ranked_candidates) deltaZs.push_back(c.second.Deltaz);
    return deltaZs;
}

TFloat64List CPdfCandidateszResult::GetMerits() const
{
    TFloat64List merits;
    merits.reserve(m_ranked_candidates.size());
    for (auto c: m_ranked_candidates) merits.push_back(c.second.ValProba);
    return merits;
}

TFloat64List CPdfCandidateszResult::GetValSumProbas() const
{
    TFloat64List valSumProbas;
    valSumProbas.reserve(m_ranked_candidates.size());
    for  (auto c: m_ranked_candidates) valSumProbas.push_back(c.second.ValSumProba);
    return valSumProbas;
}



  void CPdfCandidateszResult::getCandidateData(const int& rank,const std::string& name, Float64& v) const
  {
    if (name.compare("Redshift") == 0) v= m_ranked_candidates[rank].second.Redshift;
    else if (name.compare("RedshiftError") == 0) v=m_ranked_candidates[rank].second.Deltaz;
    else if (name.compare("RedshiftProba") == 0) v=m_ranked_candidates[rank].second.ValSumProba;
    else if (name.compare("RedshiftProbaZmin") == 0) v=m_ranked_candidates[rank].second.ValSumProbaZmin;
    else if (name.compare("RedshiftProbaZmax") == 0) v=m_ranked_candidates[rank].second.ValSumProbaZmax;
    else if (name.compare("RedshiftProbaDensity") == 0) v=m_ranked_candidates[rank].second.ValProba;
    else {
      Log.LogError("unknown candidate data %s",name.c_str());
      throw runtime_error("CPdfCandidateszResult::getCandidateData: unknown candidate data");
    }
  }

  void CPdfCandidateszResult::getCandidateData(const int& rank,const std::string& name, Int32& v) const
  {
    if (name.compare("Rank") == 0) v=rank;
    else{
      Log.LogError("unknown candidate data %s",name.c_str());
      throw runtime_error("CPdfCandidateszResult::getCandidateData: unknown candidate data");
    }
  }

  void CPdfCandidateszResult::getCandidateData(const int& rank,const std::string& name, std::string& v) const{
    v=m_ranked_candidates[rank].first;
  }

  void CPdfCandidateszResult::getCandidateData(const int& rank,const std::string& name, double **data, int *size) const{}

  void CPdfCandidateszResult::getData(const std::string& name, Int32& v) const{
    if (name.compare("NbCandidates") == 0) v= m_ranked_candidates.size();
  
  }
  void CPdfCandidateszResult::getData(const std::string& name, Float64& v) const{}
  void CPdfCandidateszResult::getData(const std::string& name, std::string& v) const{}
  void CPdfCandidateszResult::getData(const std::string& name, double **data, int *size) const
  {

  }
