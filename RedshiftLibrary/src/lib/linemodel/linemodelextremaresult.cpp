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
void CLineModelExtremaResult::SaveJSON( std::ostream& stream) const
{
    stream << "{"<< std::endl;
    
    CExtremaResult::SaveJSONbody(stream);
    stream << "," << std::endl; 

    bool zeros = std::all_of(MeritContinuum.begin(), MeritContinuum.end(), [](int i) { return i==0; });

    if(!zeros){
        // save MeritContinuum list, on 1 line
        SaveTFloat64List(stream,"z_MeritContinuum",MeritContinuum);
        stream << "," << std::endl;
        // save extrema mTransposeM list, on 1 line
        SaveTFloat64List(stream,"z_mTransposeM",mTransposeM);
        stream << "," << std::endl;
        // save NDof list, on 1 line
        SaveInt32Vector(stream,"z_NDof",NDof);
        stream << "," << std::endl;
        // save CorrScaleMarg list, on 1 line
        SaveTFloat64List(stream,"z_CorrScaleMarg",CorrScaleMarg);
        stream << "," << std::endl;
        // save lmfitPass list, on 1 line
        SaveTFloat64List(stream,"z_lmfitPass", Redshift_lmfit);
        stream << "," << std::endl;
        // save snrHa list, on 1 line
        SaveTFloat64List(stream,"z_snrHa",snrHa);
        stream << "," << std::endl;
        // save lfHa list, on 1 line
        SaveTFloat64List(stream,"z_lfHa",lfHa);
        stream << "," << std::endl;

        // save snrOII list, on 1 line
        SaveTFloat64List(stream,"z_snrOII",snrOII);
        stream << "," << std::endl;
        // save lfOII list, on 1 line
        SaveTFloat64List(stream,"z_lfOII",lfOII);
        stream << "," << std::endl;
        // save bic list, on 1 line
        SaveTFloat64List(stream,"ext_BIC",bic);
        stream << "," << std::endl;
        // save posterior list, on 1 line
        SaveTFloat64List(stream,"ext_NLinesOverThreshold", NLinesOverThreshold);
        stream << "," << std::endl;
        // save ContinuumIndexes list, on 1 line
        SaveTContinuumIndexListVector(stream,"ext_ContinuumIndexes",ContinuumIndexes);
        stream << "," << std::endl;	

        // save StrongELSNR list, on 1 line
        SaveTFloat64List(stream,"ext_StrongELSNR",StrongELSNR);
        stream << "," << std::endl;

        // save StrongELSNRAboveCut list, on 1 line
        SaveStringVectorOfVector(stream,"ext_StrongELSNRAboveCut",StrongELSNRAboveCut);
        stream << "," << std::endl;
    }

    // save FittedTplRedshift, on 1 line
    stream << std::setprecision(8);
    SaveTFloat64List(stream,"ext_FittedTplRedshift",FittedTplRedshift);  
    stream << "," << std::endl;

    if(!zeros){
        // save FittedTplpCoeffs, on 1 line
        stream << std::setprecision(8);
        SaveTFloat64ListOfList(stream,"ext_FittedTplpCoeffs",FittedTplpCoeffs);
        stream << "," << std::endl;

        // save FittedTplratioName, on 1 line
        SaveStringVector(stream,"ext_FittedTplratioName",FittedTplratioName);
        stream << "," << std::endl;
        // save FittedTplratioIsmCoeff, on 1 line
        SaveTFloat64List(stream,"ext_FittedTplratioIsmCoeff",FittedTplratioIsmCoeff);
        stream << "," << std::endl;
        // save FittedTplratioAmplitude, on 1 line
        SaveTFloat64List(stream,"ext_FittedTplratioAmplitude",FittedTplratioAmplitude);
        stream << "," << std::endl;
        // save FittedTplratioDtm, on 1 line
        SaveTFloat64List(stream,"ext_FittedTplratioDtm",FittedTplratioDtm);
        stream << "," << std::endl;
        // save FittedTplratioMtm, on 1 line
        SaveTFloat64List(stream,"ext_FittedTplratioMtm",FittedTplratioMtm);
        stream << "," << std::endl;

        // save OutsideLinesSTDFlux, on 1 line
        stream << std::scientific << std::setprecision(5);
        SaveTFloat64List(stream,"ext_OutsideLinesSTDFlux",OutsideLinesSTDFlux);
        stream << "," << std::endl;
        // save OutsideLinesSTDError, on 1 line
        stream << std::scientific << std::setprecision(5);
        SaveTFloat64List(stream,"ext_OutsideLinesSTDError",OutsideLinesSTDError);
        stream << "," << std::endl;
    }
    // save Elv, on 1 line
    stream << std::fixed << std::setprecision(1);
    SaveTFloat64List(stream,"ext_Elv",Elv);
    stream << "," << std::endl;
    // save Alv, on 1 line
    stream << std::fixed << std::setprecision(1);
    SaveTFloat64List(stream,"ext_Alv",Alv);

    stream <<  "}";

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

