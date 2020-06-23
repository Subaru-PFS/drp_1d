#ifndef _SOLVE_RESULT_
#define _SOLVE_RESULT_

#include <RedshiftLibrary/processflow/result.h>

#include <vector>

namespace NSEpic
{

  class CSolveResult: public COperatorResult
  {

  public:

    CSolveResult();
    virtual ~CSolveResult();

    virtual Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const = 0;
    void SetReliabilityLabel( std::string lbl );
    void SetTypeLabel( std::string lbl );

    virtual void preSave(const CDataStore& store) = 0;
    
    std::string m_ReliabilityLabel="-1";
    std::string m_TypeLabel="-1";

    Float64 getMerit() {return merit;}
  protected:
    Float64 redshift;
    Float64 merit;

  };
}

#endif
