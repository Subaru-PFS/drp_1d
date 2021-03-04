#include <RedshiftLibrary/method/solveresult.h>
#include <RedshiftLibrary/operator/extremaresult.h>

using namespace NSEpic;

CSolveResult::CSolveResult( const std::shared_ptr<const CExtremaResult> & ExtremaResult,
                            const std::string & opt_pdfcombination,
                            Float64 evidence):
  m_merit(opt_pdfcombination=="marg" ? ExtremaResult->ValSumProba(0) : ExtremaResult->ValProba(0)),
  m_bestRedshiftMethod(opt_pdfcombination=="marg" ? 2 :0),
  m_redshift(ExtremaResult->Redshift(0)),
  m_evidence(evidence)
{}

CSolveResult::~CSolveResult() {}

void CSolveResult::SetReliabilityLabel( std::string lbl )
{
    m_ReliabilityLabel = lbl;
}

void CSolveResult::SetTypeLabel( std::string lbl )
{
    m_TypeLabel = lbl;
}

void CSolveResult::getData(const std::string& name, std::string& v) const
{
  if (name.compare("Reliability") == 0)  v = m_ReliabilityLabel;
  else if (name.compare("Type") == 0)  v = m_TypeLabel;
  //else throw error
}
