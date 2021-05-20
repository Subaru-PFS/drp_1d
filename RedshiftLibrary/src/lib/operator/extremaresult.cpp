#include <RedshiftLibrary/operator/extremaresult.h>

#include <RedshiftLibrary/statistics/pdfcandidateszresult.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/common/exception.h>
#include <RedshiftLibrary/common/formatter.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision
#include <iostream>

using namespace NSEpic;

CExtremaResult::CExtremaResult(Int32 n)
{
    Resize(n);
}


void CExtremaResult::Resize(Int32 size)
{   
    m_ranked_candidates.resize(size);
    FittedTplName.resize(size);
    FittedTplAmplitude.resize(size);
    FittedTplAmplitudeError.resize(size);
    FittedTplMerit.resize(size);
    FittedTplEbmvCoeff.resize(size);
    FittedTplMeiksinIdx.resize(size);
    FittedTplDtm.resize(size);
    FittedTplMtm.resize(size);
    FittedTplLogPrior.resize(size);
    FittedTplSNR.resize(size);
    
    m_savedModelSpectrumResults.resize(size);
    m_savedModelContinuumFittingResults.resize(size);
}


void CExtremaResult::getCandidateData(const int& rank,const std::string& name, Float64& v) const
{
    if (name.compare("ContinuumIsmCoeff") == 0 || name.compare("FirstpassContinuumIsmCoeff") == 0) v = FittedTplEbmvCoeff[rank];
    else if (name.compare("FirstpassRedshift") == 0) v = m_ranked_candidates[rank].second.Redshift;
    else if (name.compare("FirstpassMerit") == 0) v = m_ranked_candidates[rank].second.ValProba;
    else if (name.compare("ContinuumAmplitude") == 0 || name.compare("FirstpassContinuumAmplitude") == 0) v = FittedTplAmplitude[rank];
    else if (name.compare("FittedTemplateDtm") == 0 || name.compare("FirstpassFittedTemplateDtm") == 0) v = FittedTplDtm[rank];
    else if (name.compare("FittedTemplateMtm") == 0 || name.compare("FirstpassFittedTemplateMtm") == 0) v = FittedTplMtm[rank];
    else if (name.compare("ContinuumAmplitudeError") == 0 || name.compare("FirstpassContinuumAmplitudeError") == 0) v = FittedTplAmplitudeError[rank];
    else if (name == "Redshift") v = Redshift(rank);
    else if (name == "RedshiftProba") v = ValSumProba(rank);
    else if (name == "RedshiftError") v = DeltaZ(rank);
    else throw GlobalException(UNKNOWN_ATTRIBUTE,Formatter() << "unknown candidate Float64 data " << name);
}

void CExtremaResult::getCandidateData(const int& rank,const std::string& name, Int32& v) const
{
  if (name.compare("ContinuumIgmIndex") == 0 || name.compare("FirstpassContinuumIgmIndex") == 0) v = FittedTplMeiksinIdx[rank];
  else throw GlobalException(UNKNOWN_ATTRIBUTE,Formatter() << "unknown candidate integer data "<<name);
}

void CExtremaResult::getCandidateData(const int& rank,const std::string& name, std::string& v) const
{
  if (name.compare("TemplateName") == 0 || name.compare("FirstpassTemplateName") == 0) v = FittedTplName[rank];
  else if (name.compare("ExtremaIds") == 0 || name.compare("FirstpassExtremaIds") == 0 ) v = m_ranked_candidates[rank].first;
  else throw GlobalException(UNKNOWN_ATTRIBUTE,Formatter() <<"unknown candidate string data "<<name);
}

void CExtremaResult::getCandidateData(const int& rank,const std::string& name, double **data, int *size) const {}
void CExtremaResult::getData(const std::string& name, Int32& v) const {}
void CExtremaResult::getData(const std::string& name, Float64& v) const {}
void CExtremaResult::getData(const std::string& name, std::string& v) const {}
void CExtremaResult::getData(const std::string& name, double **data, int *size) const {}

