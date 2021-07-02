#include "RedshiftLibrary/method/classificationresult.h"

#include "RedshiftLibrary/processflow/context.h"
#include <stdio.h>
#include <float.h>
#include <math.h>
#include "RedshiftLibrary/log/log.h"


using namespace NSEpic;

/**
 * \brief Empty constructor.
 **/
CClassificationResult::CClassificationResult():
  CSolveResult()  
{
  this->m_type="CClassificationResult";
}

/**
 * \brief Empty destructor.
 **/

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

