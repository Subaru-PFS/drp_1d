#ifndef _SOLVE_RESULT_
#define _SOLVE_RESULT_

#include <RedshiftLibrary/processflow/result.h>

#include <vector>

namespace NSEpic
{

  class TCandidateZ;

  class CSolveResult: public COperatorResult
  {
  public:
    
    virtual ~CSolveResult()=0;
  };
  
  class CPdfSolveResult: public CSolveResult
  {

  public:

    CPdfSolveResult(const TCandidateZ& ExtremaResult,
                 const std::string & opt_pdfcombination,
                 Float64 evidence);
    CPdfSolveResult( Float64 merit, Float64 redshift,
                      const std::string & opt_pdfcombination,
                      Float64 evidence);
    virtual ~CPdfSolveResult()=default;
    CPdfSolveResult(CPdfSolveResult const& other) = default;
    CPdfSolveResult& operator=(CPdfSolveResult const& other) = default;
    CPdfSolveResult(CPdfSolveResult&& other) = default;
    CPdfSolveResult& operator=(CPdfSolveResult&& other) = default;


    //virtual void preSave(const CDataStore& store) = 0;
    Int32 m_bestRedshiftMethod = 2; //0:best chi2 or proba, 2: best marg proba


    inline Float64 getMerit() const {return m_merit;}  
    inline Float64 getEvidence() const {return m_evidence;} 


  protected:
    Float64 m_redshift;
    Float64 m_merit;
    Float64 m_evidence = -INFINITY;

  };
}

#endif
