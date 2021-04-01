#ifndef _REDSHIFT_STATISTICS_PDFCANDIDATESZRESULT_
#define _REDSHIFT_STATISTICS_PDFCANDIDATESZRESULT_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/statistics/pdfcandidatesz.h>

//#include <cmath>
//#include <map>
#include <ostream>
#include <string>

namespace NSEpic
{

class CPdfCandidateszResult : public COperatorResult
{

public:

    CPdfCandidateszResult(Int32 optMethod=0): m_optMethod(optMethod) {};

    virtual ~CPdfCandidateszResult() = default;

    //rule of 5 defaults
    CPdfCandidateszResult(const CPdfCandidateszResult & ) = default;
    CPdfCandidateszResult(CPdfCandidateszResult && ) = default;
    CPdfCandidateszResult & operator=(const CPdfCandidateszResult & ) = default;   
    CPdfCandidateszResult & operator=(CPdfCandidateszResult && ) = default;   
  
    Int32 size() const;
 
    void Save(std::ostream& stream ) const;
    void SaveLine(std::ostream& stream ) const;
   
    Float64 getDouble(std::string name,Int32 rank) const;
    std::string getString(std::string name,Int32 rank) const;
    Int32 getInt(std::string name,Int32 rank) const;
    Int32 getNbCandidates() const;

    std::string ID(Int32 i) const;
    Float64 Redshift(Int32 i) const;
    Float64 ValProba(Int32 i) const;
    Float64 ValSumProba(Int32 i) const;
    Float64 DeltaZ(Int32 i) const;

    TStringList GetIDs() const;
    TFloat64List GetRedshifts() const;
    TFloat64List GetDeltaZs() const;
    TFloat64List GetMerits() const;
    TFloat64List GetValSumProbas() const;

    void getCandidateData(const int& rank,const std::string& name, Float64& v) const;
    void getCandidateData(const int& rank,const std::string& name, Int32& v) const;
    void getCandidateData(const int& rank,const std::string& name, std::string& v) const;
    void getCandidateData(const int& rank,const std::string& name, double **data, int *size) const;

    void getData(const std::string& name, Int32& v) const;
    void getData(const std::string& name, Float64& v) const;
    void getData(const std::string& name, std::string& v) const;
    void getData(const std::string& name, double **data, int *size) const;

    Int32                       m_optMethod; //0: direct integration, 1:gaussian fit

    TCandidateZbyRank m_ranked_candidates;

private:

};

inline Int32 CPdfCandidateszResult::size() const
{
    return m_ranked_candidates.size();
}

// for compatibility
inline std::string CPdfCandidateszResult::ID(Int32 i) const {return m_ranked_candidates[i].first;}
inline Float64 CPdfCandidateszResult::Redshift(Int32 i) const { return m_ranked_candidates[i].second.Redshift;}
inline Float64 CPdfCandidateszResult::ValProba(Int32 i) const { return m_ranked_candidates[i].second.ValProba;}
inline Float64 CPdfCandidateszResult::ValSumProba(Int32 i) const { return m_ranked_candidates[i].second.ValSumProba;}
inline Float64 CPdfCandidateszResult::DeltaZ(Int32 i) const { return m_ranked_candidates[i].second.Deltaz;}

}

#endif
