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
    ~CPdfCandidateszResult() = default;

    void Save(std::ostream& stream ) const;
    void SaveLine(std::ostream& stream ) const;
   
    Float64 getDouble(std::string name,Int32 rank) const;
    std::string getString(std::string name,Int32 rank) const;
    Int32 getInt(std::string name,Int32 rank) const;
    Int32 getNbCandidates() const;

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

}

#endif
