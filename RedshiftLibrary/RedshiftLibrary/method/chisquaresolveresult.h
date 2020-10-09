#ifndef _REDSHIFT_OPERATOR_CHISQUARESOLVERESULT_
#define _REDSHIFT_OPERATOR_CHISQUARESOLVERESULT_

#include <RedshiftLibrary/method/solveresult.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/ray/catalog.h>

#include <vector>
#include <unordered_map>

namespace NSEpic
{

class CProcessFlowContext;

/**
 * \ingroup Redshift
 */
class CChisquareSolveResult : public CSolveResult
{

public:

    enum EType
    {
             nType_raw = 1,
             nType_continuumOnly = 2,
             nType_noContinuum = 3,
             nType_all = 4,
    };

    CChisquareSolveResult(const Int32 type=nType_raw, const std::string scope="chisquaresolve");

    void Save( const CDataStore& store, std::ostream& stream ) const;
    void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    Bool GetBestRedshift(const CDataStore& store) const;
    Bool GetBestRedshiftPerTemplateString( const CDataStore& store, std::string& output ) const;
    Bool GetBestRedshiftFromPdf(const CDataStore& store) const;
    Int32 GetBestModel(const CDataStore& store, Float64 z, std::string& tplName, Int32& MeiksinIdx, Float64& DustCoeff, Float64& Amplitude) const;

    Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const;
    Bool GetRedshiftCandidates( const CDataStore& store,  std::vector<Float64>& redshiftcandidates, Int32 n_candidates, std::string outputPdfRelDir = "zPDF") const;

    Int32 m_bestRedshiftMethod = 2; //best chi2, best proba

  void preSave(const CDataStore& store);

  void getData(const std::string& name, Float64& v) const;
  void getData(const std::string& name, std::string& v) const;
  void getData(const std::string& name, Int32& v) const;
  void getCandidateData(const int& rank,const std::string& name, Float64& v) const;
  void getCandidateData(const int& rank,const std::string& name, Int32& v) const;
  void getCandidateData(const int& rank,const std::string& name, std::string& v) const;

private:

    std::unordered_map<std::string, std::string> m_scope2name = {
        {"chisquaresolve",      "ChisquareSolve"},
        {"chisquare2solve",     "Chisquare2Solve"},
        {"chisquarelogsolve",   "ChisquareLogSolve"},
        {"tplcombinationsolve", "TplcombinationSolve"}
    };

    const Int32 m_type;
    const std::string m_scope;
    std::string m_name;

  mutable Float64 m_redshift;
  mutable Float64 m_merit; 
  mutable Float64 m_evidence;
  mutable std::string m_tplName;
  //TODO: is not normal to use mutable, but kept it for test
  mutable Float64 m_amplitude;
  mutable Float64 m_amplitudeError;
  mutable Float64 m_dustCoeff;
  mutable Int32   m_meiksinIdx;

};


}

#endif
