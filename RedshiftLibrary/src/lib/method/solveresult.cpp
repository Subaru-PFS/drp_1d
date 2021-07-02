#include "RedshiftLibrary/method/solveresult.h"
#include "RedshiftLibrary/operator/extremaresult.h"

using namespace NSEpic;

CSolveResult::~CSolveResult(){} // yes pure virtual destructors needs a body !

CPdfSolveResult::CPdfSolveResult( const TCandidateZ& BestExtremumResult,
                            const std::string & opt_pdfcombination,
                            Float64 evidence):
  m_merit(opt_pdfcombination=="marg" ? BestExtremumResult.ValSumProba : BestExtremumResult.ValProba),
  m_bestRedshiftMethod(opt_pdfcombination=="marg" ? 2 :0),
  m_redshift(BestExtremumResult.Redshift),
  m_evidence(evidence)
{}
//a temporaru hack used by tplcombination
CPdfSolveResult::CPdfSolveResult( Float64 merit, Float64 redshift,
                            const std::string & opt_pdfcombination,
                            Float64 evidence):
  m_merit(merit),
  m_bestRedshiftMethod(opt_pdfcombination=="marg" ? 2 :0),
  m_redshift(redshift),
  m_evidence(evidence)
{}