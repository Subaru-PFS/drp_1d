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

void CExtremaResult::SaveJSON(std::ostream& stream) const
{
    stream << "{"<< std::endl;
    SaveJSONbody(stream);
    stream <<  "}";
}

/**
 * \brief Prints the results currently in the argument store, in the argument stream.
 * Using the argument stream to print values from the argument store:
 * Print a header as a comment.
 * Print each redshift and chisquare values.
 * Print each Extrema as a comment.
 * Print each BIC as a comment.
 * Print each POSTERIOR as a comment.
 * Print each SigmaZ as a comment.
 * Print each LogArea as a comment.
 **/
void CExtremaResult::SaveJSONbody(std::ostream& stream) const
{    
    // save extrema list, on 1 line
    SaveTFloat64List(stream,"z_extrema", GetRedshifts());
    stream << "," << std::endl;

    // save extrema list, on 1 line
    SaveStringVector(stream,"z_extremaIDs", GetIDs());
    stream << "," << std::endl;

    // save extremaPDF list, on 1 line
    SaveTFloat64List(stream,"z_extremaIntgPDF", GetValSumProbas());
    stream << "," << std::endl;
  
    // save extremaMerit list, on 1 line
    SaveTFloat64List(stream,"z_ExtremaMerit", GetMerits());
    stream << "," << std::endl;

    // save extrema Deltaz list, on 1 line
    SaveTFloat64List(stream,"z_ExtremaDeltaZ", GetDeltaZs());
    stream << "," << std::endl;

    // save FittedTplName, on 1 line
    SaveStringVector(stream,"ext_FittedTplName",FittedTplName);
    stream << "," << std::endl;
    // save FittedTplAmplitude, on 1 line
    SaveTFloat64List(stream,"ext_FittedTplAmplitude",FittedTplAmplitude);
    stream << "," << std::endl;
    // save FittedTplAmplitudeError, on 1 line
    SaveTFloat64List(stream,"ext_FittedTplAmplitudeError",FittedTplAmplitudeError);
    stream << "," << std::endl;
    // save FittedTplMerit, on 1 line
    SaveTFloat64List(stream,"ext_FittedTplMerit",FittedTplMerit);
    stream << "," << std::endl;
    // save FittedTplEbmvCoeff, on 1 line
    stream << std::setprecision(3);
    SaveTFloat64List(stream,"ext_FittedTplEbmvCoeff",FittedTplEbmvCoeff); 
    stream << "," << std::endl; 

    // save FittedTplMeiksinIdx, on 1 line
    SaveInt32Vector(stream,"ext_FittedTplMeiksinIdx",FittedTplMeiksinIdx);
    stream << "," << std::endl;
    
    // save FittedTplDtm, on 1 line
        stream << std::setprecision(8);
    SaveTFloat64List(stream,"ext_FittedTplDtm",FittedTplDtm);
    stream << "," << std::endl;
    // save FittedTplMtm, on 1 line
        stream << std::setprecision(8);
    SaveTFloat64List(stream,"ext_FittedTplMtm",FittedTplMtm);
    stream << "," << std::endl;

    // save FittedTplLogPrior, on 1 line
    stream << std::setprecision(8);
    SaveTFloat64List(stream,"ext_FittedTplLogPrior",FittedTplLogPrior);
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

