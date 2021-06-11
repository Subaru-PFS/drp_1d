#include <RedshiftLibrary/method/solveresult.h>
#include <RedshiftLibrary/operator/extremaresult.h>

using namespace NSEpic;

CSolveResult::~CSolveResult(){} // yes pure virtual destructors needs a body !

CPdfSolveResult::CPdfSolveResult( const std::shared_ptr<const CExtremaResult> & ExtremaResult,
                            const std::string & opt_pdfcombination,
                            Float64 evidence):
  m_merit(opt_pdfcombination=="marg" ? ExtremaResult->ValSumProba(0) : ExtremaResult->ValProba(0)),
  m_bestRedshiftMethod(opt_pdfcombination=="marg" ? 2 :0),
  m_redshift(ExtremaResult->Redshift(0)),
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