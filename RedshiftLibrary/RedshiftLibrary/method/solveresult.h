#ifndef _SOLVE_RESULT_
#define _SOLVE_RESULT_

#include <RedshiftLibrary/processflow/result.h>

#include <vector>

namespace NSEpic
{
  class CExtremaResult;

  class CSolveResult: public COperatorResult
  {

  public:

    CSolveResult(const std::shared_ptr<const CExtremaResult> & ExtremaResult,
                 const std::string & opt_pdfcombination,
                 Float64 evidence);
    virtual ~CSolveResult();

    void SetReliabilityLabel( std::string lbl );
    void SetTypeLabel( std::string lbl );

    //virtual void preSave(const CDataStore& store) = 0;
    Int32 m_bestRedshiftMethod = 2; //0:best chi2 or proba, 2: best marg proba

    std::string m_ReliabilityLabel="-1";
    std::string m_TypeLabel="-1";

    Float64 getMerit() {return m_merit;}
    Float64 getEvidence() {return m_evidence;}

    void getData(const std::string& name, std::string& v) const;

  protected:
    Float64 m_redshift;
    Float64 m_merit;
    Float64 m_evidence = -INFINITY;

  };
}

#endif
