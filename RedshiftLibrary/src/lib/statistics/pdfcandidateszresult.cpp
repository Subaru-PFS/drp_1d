#include <RedshiftLibrary/statistics/pdfcandidateszresult.h>
#include <RedshiftLibrary/log/log.h>

//#include <RedshiftLibrary/processflow/context.h>


using namespace NSEpic;
using namespace std;
#include <fstream>
#include <iostream>
#include <numeric>


void CPdfCandidateszResult::Save( std::ostream& stream ) const
{
    stream  << "#fullwidth = " << "6* Deltaz" << std::endl;
    stream  << "#method = " << m_optMethod << std::endl;
    stream  << std::endl;

    //    stream  << "#" << store.GetSpectrumName() << "\t" << store.GetProcessingID() << "\t";
    stream  << std::endl;

    stream  << "#" << "rank" << "\t"  << "IDs" << "\t"<< "redshift" << "\t" << "intgProba"<< "\t" << "\t" <<"Deltaz";
    if(m_optMethod==1)
    {
        stream << "\t" << "gaussAmp" << "\t" << "gaussAmpErr" << "\t" << "gaussSigma" << "\t" << "gaussSigmaErr";
    }else{
        stream << "\t" << "gaussAmp_unused" << "\t" << "gaussAmpErr_unused" << "\t" << "gaussSigma_unused" << "\t" << "gaussSigmaErr_unused";
    }
    stream  << "\n";
    for(Int32 i=0; i< m_ranked_candidates.size(); ++i)
    {
        const TCandidateZ & cand = m_ranked_candidates[i].second;
        const std::string & Id = m_ranked_candidates[i].first;
        stream << i << "\t"; 
        stream << Id << "\t";
        stream << cand.Redshift << "\t";
        stream << cand.ValSumProba << "\t";
        stream << cand.Deltaz << "\t";
        //only for method 1, but leave columns with -1 value ste in compute()
        stream << cand.GaussAmp << "\t";
        stream << cand.GaussAmpErr << "\t";
        stream << cand.GaussSigma << "\t";
        stream << cand.GaussSigmaErr << "\t";

        stream << "\n";
    }
    stream << std::endl;
}

void CPdfCandidateszResult::SaveLine( std::ostream& stream ) const
{
  //    stream  << store.GetSpectrumName() << "\t" << store.GetProcessingID() << "\t";
    for(Int32 i=0; i< m_ranked_candidates.size(); ++i)
    {
        const TCandidateZ & cand = m_ranked_candidates[i].second;
        const std::string & Id = m_ranked_candidates[i].first;
        stream << i << "\t";
        stream << Id << "\t";
        stream << cand.Redshift << "\t";
        stream << cand.ValSumProba << "\t";
        stream << cand.Deltaz << "\t";
        stream << cand.GaussAmp << "\t";
        stream << cand.GaussSigma << "\t"; 
    }
    stream << std::endl;
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
