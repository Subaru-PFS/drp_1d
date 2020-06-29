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
    stream <<  "#Type\tEvidenceG\tEvidenceS\tEvidenceQ"<< std::endl;
    stream << m_TypeLabel << "\t"
       << m_evidence_galaxy << "\t"
       << m_evidence_star << "\t"
       << m_evidence_qso << "\t"
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
            << std::endl;
}

