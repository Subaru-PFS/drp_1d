#ifndef _REDSHIFT_METHOD_TEMPLATEFITTINGSOLVERESULT_
#define _REDSHIFT_METHOD_TEMPLATEFITTINGSOLVERESULT_

#include <RedshiftLibrary/method/solveresult.h>
#include <RedshiftLibrary/operator/extremaresult.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/ray/catalog.h>

#include <memory>
#include <vector>
#include <unordered_map>
#include <cmath>

namespace NSEpic
{


/**
 * \ingroup Redshift
 */
class CTemplateFittingSolveResult : public CPdfSolveResult
{

public:
    CTemplateFittingSolveResult(const std::string & scope, 
                                const std::shared_ptr<const CExtremaResult> & ExtremaResult,
                                const std::string & opt_pdfcombination,
                                Float64 evidence );
    //CTemplateFittingSolveResult(const EType type=nType_raw, const std::string scope="templatefittingsolve");

    void Save(std::ostream& stream ) const;
    void SaveLine(std::ostream& stream ) const;
/*    Bool GetBestRedshift(const CDataStore& store);
    Bool GetBestRedshiftPerTemplateString( const CDataStore& store, std::string& output ) const;
    Bool GetBestRedshiftFromPdf(const CDataStore& store);
    Int32 GetBestModel(const CDataStore& store, Float64 z);*/
    
/*  void preSave(const CDataStore& store);*/

  void getData(const std::string& name, Float64& v) const;
  void getData(const std::string& name, std::string& v) const;
  void getData(const std::string& name, Int32& v) const;
  void getCandidateData(const int& rank,const std::string& name, Float64& v) const;
  void getCandidateData(const int& rank,const std::string& name, Int32& v) const;
  void getCandidateData(const int& rank,const std::string& name, std::string& v) const;

  /*const std::string GetTemplateName();
  const Float64 GetAmplitude();
  const Float64 GetMeiksinIdx();
  const Float64 GetDustCoeff();
  const Float64 GetAmplitudeError();
  const Float64 GetMerit();
  const Float64 GetFittingSNR();*/

  //Extrema results
  std::shared_ptr<const CExtremaResult> ExtremaResult;

private:

/*    std::unordered_map<std::string, std::string> m_scope2name = {
        {"templatefittingsolve",      "TemplateFittingSolve"},
        {"templatefittinglogsolve",   "TemplateFittingLogSolve"},
        {"tplcombinationsolve", "TplcombinationSolve"}
    };*/

    //const EType m_type;
    const std::string m_scope;
    //std::string m_name;

  std::string m_tplName = "-1";
  Float64 m_amplitude = 0.0;
  Float64 m_amplitudeError = -1.0;
  Float64 m_dustCoeff = -1.0;
  Int32   m_meiksinIdx = -1.0;

  //Not sure it is necessary here
  Float64   m_fittingSNR = NAN;

};
/*
inline
const std::string CTemplateFittingSolveResult::GetTemplateName(){
  return m_tplName;
}
inline
const Float64 CTemplateFittingSolveResult::GetAmplitude(){
  return m_amplitude;
}
inline 
const Float64 CTemplateFittingSolveResult::GetMeiksinIdx(){
  return m_meiksinIdx;
}
inline
const Float64 CTemplateFittingSolveResult::GetDustCoeff(){
  return m_dustCoeff;
}
inline
const Float64 CTemplateFittingSolveResult::GetAmplitudeError(){
  return m_amplitudeError;
}
inline 
const Float64 CTemplateFittingSolveResult::GetMerit(){
  return m_merit;
}
inline
const Float64 CTemplateFittingSolveResult::GetFittingSNR(){
  return m_fittingSNR;
}
*/
}

#endif
