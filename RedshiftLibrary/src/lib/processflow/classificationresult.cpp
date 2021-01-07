#include <RedshiftLibrary/processflow/classificationresult.h>

#include <RedshiftLibrary/processflow/context.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <RedshiftLibrary/log/log.h>


using namespace NSEpic;

/**
 * \brief Empty constructor.
 **/
CClassificationResult::CClassificationResult()
{
}

/**
 * \brief Empty destructor.
 **/
CClassificationResult::~CClassificationResult()
{

}

void CClassificationResult::SetTypeLabel( std::string lbl )
{
    m_TypeLabel = lbl;
}

void CClassificationResult::SetG(Float64 evidence, Float64 prob)
{
    m_evidence_galaxy = evidence;
    m_prob_galaxy = prob;
}

void CClassificationResult::SetS(Float64 evidence, Float64 prob)
{
    m_evidence_star = evidence;
    m_prob_star = prob;
}

void CClassificationResult::SetQ(Float64 evidence, Float64 prob)
{
    m_evidence_qso = evidence;
    m_prob_qso = prob;
}

/**
 * \brief
 **/
void CClassificationResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    stream <<  "#Type\tLogEvidenceG\tLogEvidenceS\tLogEvidenceQ\tProbaG\tProbaS\tProbaQ"<< std::endl;
    stream << m_TypeLabel << "\t"
       << m_evidence_galaxy << "\t"
       << m_evidence_star << "\t"
       << m_evidence_qso << "\t"

       << m_prob_galaxy << "\t"
       << m_prob_star << "\t"
       << m_prob_qso << "\t"
       << std::endl;
}

/**
 * \brief
 **/
void CClassificationResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    stream  << store.GetSpectrumName() << "\t"
            << store.GetProcessingID() << "\t"
            << m_TypeLabel << "\t"
            << m_evidence_galaxy << "\t"
            << m_evidence_star << "\t"
            << m_evidence_qso << "\t"
            << m_prob_galaxy << "\t"
            << m_prob_star << "\t"
            << m_prob_qso << "\t"
            << std::endl;
}

void CClassificationResult::getData(const std::string& name, std::string& v) const
{
  Log.LogDebug("CClassificationResult::getData getting type=%s",m_TypeLabel.c_str());
  v = m_TypeLabel;
}

void CClassificationResult::getData(const std::string& name, Float64& v) const
{
  if (name.compare("EvidenceGalaxy") == 0)  v = m_evidence_galaxy;
  else if (name.compare("EvidenceQSO") == 0)  v = m_evidence_qso;
  else if (name.compare("EvidenceStar") == 0)  v = m_evidence_star;
  else if (name.compare("ProbGalaxy") == 0)  v = m_prob_galaxy;
  else if (name.compare("ProbQSO") == 0)  v = m_prob_qso;
  else if (name.compare("ProbStar") == 0)  v = m_prob_star;
  //else throw error
}
