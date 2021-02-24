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
                                const TCandidateZ& ExtremaResult,
                                const std::string & opt_pdfcombination,
                                Float64 evidence );
    //CTemplateFittingSolveResult(const EType type=nType_raw, const std::string scope="templatefittingsolve");

/*    Bool GetBestRedshift(const CDataStore& store);
    Bool GetBestRedshiftPerTemplateString( const CDataStore& store, std::string& output ) const;
    Bool GetBestRedshiftFromPdf(const CDataStore& store);
    Int32 GetBestModel(const CDataStore& store, Float64 z);*/
    
/*  void preSave(const CDataStore& store);*/



  //Extrema results

    const std::string m_scope;
    //std::string m_name;

  std::string m_tplName = "-1";
  Float64 m_amplitude = 0.0;
  Float64 m_amplitudeError = -1.0;
  Float64 m_EbmvCoeff = -1.0;
  Int32   m_meiksinIdx = -1.0;

  //Not sure it is necessary here
  Float64   m_fittingSNR = NAN;


/*    std::unordered_map<std::string, std::string> m_scope2name = {
        {"templatefittingsolve",      "TemplateFittingSolve"},
        {"templatefittinglogsolve",   "TemplateFittingLogSolve"},
        {"tplcombinationsolve", "TplcombinationSolve"}
    };*/

    //const EType m_type;


};

}

#endif
