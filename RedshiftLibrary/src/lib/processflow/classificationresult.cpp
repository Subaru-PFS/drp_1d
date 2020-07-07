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

void CClassificationResult::SetG(Float64 evidence)
{
    m_evidence_galaxy = evidence;
}

void CClassificationResult::SetS(Float64 evidence)
{
    m_evidence_star = evidence;
}

void CClassificationResult::SetQ(Float64 evidence)
{
    m_evidence_qso = evidence;
}

/**
 * \brief
 **/
void CClassificationResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    Float64 sum = m_evidence_galaxy;
    if(m_evidence_star > 0 )
        sum += m_evidence_star;
    if(m_evidence_qso>0) 
        sum += m_evidence_qso;

    stream <<  "#Type\tEvidenceG\tEvidenceS\tEvidenceQ\tProbaG\tProbaS\tProbaQ"<< std::endl;
    stream << m_TypeLabel << "\t"
       << m_evidence_galaxy << "\t"
       << m_evidence_star << "\t"
       << m_evidence_qso << "\t"

       << m_evidence_galaxy/sum << "\t"
       << m_evidence_star/sum << "\t"
       << m_evidence_qso/sum << "\t"
       << std::endl;
}

/**
 * \brief
 **/
void CClassificationResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    Float64 sum = m_evidence_galaxy;
    if(m_evidence_star>0)
        sum += m_evidence_star;
    if(m_evidence_qso>0) 
        sum += m_evidence_qso;
    stream  << store.GetSpectrumName() << "\t"
            << store.GetProcessingID() << "\t"
            << m_TypeLabel << "\t"
            << m_evidence_galaxy << "\t"
            << m_evidence_star << "\t"
            << m_evidence_qso << "\t"


            << m_evidence_galaxy/sum << "\t"
            << m_evidence_star/sum << "\t"
            << m_evidence_qso/sum << "\t"
            << std::endl;
}

