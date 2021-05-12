#include <RedshiftLibrary/linemodel/linemodelextremaresult.h>

#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/common/exception.h>
#include <RedshiftLibrary/common/formatter.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision
#include <iostream>
#include <numeric>
using namespace NSEpic;

CLineModelExtremaResult::CLineModelExtremaResult(Int32 n)
{
    Resize(n);
}

void CLineModelExtremaResult::Resize(Int32 size)
{   
    CExtremaResult::Resize(size);

    MeritContinuum.resize(size);

    mTransposeM.resize(size);
    CorrScaleMarg.resize(size);
    NDof.resize(size);
    Redshift_lmfit.resize(size);
    snrHa.resize(size);
    lfHa.resize(size);
    snrOII.resize(size);
    lfOII.resize(size);

    ExtendedRedshifts.resize(size);
    NLinesOverThreshold.resize(size);
    LogArea.resize(size);
    LogAreaCorrectedExtrema.resize(size);
    SigmaZ.resize(size);

    StrongELSNR.resize(size);
    StrongELSNRAboveCut.resize(size);
    bic.resize(size);
    ContinuumIndexes.resize(size);
    OutsideLinesMask.resize(size);
    OutsideLinesSTDFlux.resize(size);
    OutsideLinesSTDError.resize(size);

    Elv.resize(size);
    Alv.resize(size);
    GroupsELv.resize(size);
    GroupsALv.resize(size);
    for(Int32 ke=0; ke<size; ke++)
    {
        GroupsELv[ke] = std::vector<Float64>(250, -1);   //WARNING: hardcoded ngroups max
        GroupsALv[ke] = std::vector<Float64>(250, -1);   //WARNING: hardcoded ngroups max
    }

    FittedTplRedshift.resize(size);
    FittedTplpCoeffs.resize(size);

    FittedTplratioName.resize(size);
    FittedTplratioAmplitude.resize(size);
    FittedTplratioDtm.resize(size);
    FittedTplratioMtm.resize(size);
    FittedTplratioIsmCoeff.resize(size);
    
    m_savedModelFittingResults.resize(size);
    m_savedModelRulesResults.resize(size);
    m_savedModelContinuumSpectrumResults.resize(size);
}


void CLineModelExtremaResult::getCandidateData(const int& rank,const std::string& name, Float64& v) const
{
    if (name.compare("VelocityEmission") == 0 || name.compare("FirstpassVelocityEmission") == 0) v = Elv[rank];
    else if (name.compare("VelocityAbsorption") == 0 || name.compare("FirstpassVelocityAbsorption") == 0) v = Alv[rank];
    else if (name.compare("StrongEmissionLinesSNR") == 0) v = StrongELSNR[rank];
    else if (name.compare("LinesRatioIsmCoeff") == 0) v = FittedTplratioIsmCoeff[rank];
    else if (name.compare("LinesRatioAmplitude") == 0) v = FittedTplratioAmplitude[rank];
    else CExtremaResult::getCandidateData(rank, name, v);
}

void CLineModelExtremaResult::getCandidateData(const int& rank,const std::string& name, Int32& v) const
{
    CExtremaResult::getCandidateData(rank, name, v);
}

void CLineModelExtremaResult::getCandidateData(const int& rank,const std::string& name, std::string& v) const
{
    if (name.compare("LinesRatioName") == 0) v = FittedTplratioName[rank];
    else CExtremaResult::getCandidateData(rank, name, v);
}

void CLineModelExtremaResult::getCandidateData(const int& rank,const std::string& name, double **data, int *size) const
{
  if (name.compare("ContinuumIndexesColor") == 0)
    {
      *size = ContinuumIndexes[rank].size();
      if (continuumIndexesColorCopy.find(rank) != continuumIndexesColorCopy.end())
	{
	  continuumIndexesColorCopy[rank] = TFloat64List(*size);
	  for (UInt32 j=0; j<*size;j++) continuumIndexesColorCopy[rank].emplace_back(ContinuumIndexes[rank][j].Color);
	  *data = const_cast<double *>(continuumIndexesColorCopy[rank].data());
	}
    }
  else if( name.compare("ContinuumIndexesBreak") == 0)
    {
      *size = ContinuumIndexes[rank].size();
      if (continuumIndexesBreakCopy.find(rank) != continuumIndexesBreakCopy.end())
	{
	  continuumIndexesBreakCopy[rank] = TFloat64List(*size);
	  for (UInt32 j=0; j<*size;j++) continuumIndexesBreakCopy[rank].emplace_back(ContinuumIndexes[rank][j].Break);
	  *data = const_cast<double *>(continuumIndexesBreakCopy[rank].data());
	}
    }
  else throw GlobalException(UNKNOWN_ATTRIBUTE,Formatter() <<"unknown candidate string data "<<name);
}


void CLineModelExtremaResult::getData(const std::string& name, Int32& v) const{}
void CLineModelExtremaResult::getData(const std::string& name, Float64& v) const{}
void CLineModelExtremaResult::getData(const std::string& name, std::string& v) const{}
void CLineModelExtremaResult::getData(const std::string& name, double **data, int *size) const
{

}

