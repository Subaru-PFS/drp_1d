#include <RedshiftLibrary/linemodel/linemodelextremaresult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision

#include <string>

using namespace NSEpic;


/**
 * \brief Empty constructor.
 **/
CLineModelExtremaResult::CLineModelExtremaResult()
{

}

/**
 * \brief Empty destructor.
 **/
CLineModelExtremaResult::~CLineModelExtremaResult()
{

}

Void CLineModelExtremaResult::Resize(Int32 size)
{
    Extrema.resize(size);
    ExtremaMerit.resize(size);
    DeltaZ.resize(size);
    mTransposeM.resize(size);
    ExtremaLastPass.resize(size);
    lmfitPass.resize(size);
    Posterior.resize(size);
    StrongELSNR.resize(size);
    LogArea.resize(size);
    LogAreaCorrectedExtrema.resize(size);
    SigmaZ.resize(size);
    bic.resize(size);
    ContinuumIndexes.resize(size);
    OutsideLinesMask.resize(size);
    FittedTplName.resize(size);
    FittedTplAmplitude.resize(size);
    FittedTplMerit.resize(size);
    FittedTplDustCoeff.resize(size);
    FittedTplMeiksinIdx.resize(size);
    FittedTplshapeName.resize(size);
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
Void CLineModelExtremaResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    // save extrema list, on 1 line
    if(Extrema.size()>0){
        stream <<  "#Extrema for z = {";
        for ( int i=0; i<Extrema.size(); i++)
        {
            stream <<  Extrema[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save extremaMerit list, on 1 line
    if(ExtremaMerit.size()>0){
        stream <<  "#ExtremaMerit for z = {";
        for ( int i=0; i<ExtremaMerit.size(); i++)
        {
            stream <<  ExtremaMerit[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save extrema Deltaz list, on 1 line
    if(DeltaZ.size()>0){
        stream <<  "#ExtremaDeltaZ for z = {";
        for ( int i=0; i<DeltaZ.size(); i++)
        {
            stream <<  DeltaZ[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save extrema mTransposeM list, on 1 line
    if(mTransposeM.size()>0){
        stream <<  "#mTransposeM for z = {";
        for ( int i=0; i<mTransposeM.size(); i++)
        {
            stream <<  mTransposeM[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save ExtremaLastPass list, on 1 line
    if(ExtremaLastPass.size()>0){
        stream <<  "#ExtremaLastPass for z = {";
        for ( int i=0; i<ExtremaLastPass.size(); i++)
        {
            stream <<  ExtremaLastPass[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save lmfitPass list, on 1 line
    if(lmfitPass.size()>0){
        stream <<  "#lmfitPass for z = {";
        for ( int i=0; i<lmfitPass.size(); i++)
        {
            stream <<  lmfitPass[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save bic list, on 1 line
    if(bic.size()>0){
        stream <<  "#BIC for each extrema = {";
        for ( int i=0; i<bic.size(); i++)
        {
            stream <<  bic[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save posterior list, on 1 line
    if(Posterior.size()>0){
        stream <<  "#POSTERIOR for each extrema = {";
        for ( int i=0; i<Posterior.size(); i++)
        {
            stream <<  Posterior[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save SigmaZ list, on 1 line
    if(Extrema.size()>0){
        stream <<  "#SigmaZ for each extrema = {";
        for ( int i=0; i<SigmaZ.size(); i++)
        {
            stream <<  SigmaZ[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save LogArea list, on 1 line
    if(Extrema.size()>0){
        stream <<  "#LogArea for each extrema = {";
        for ( int i=0; i<LogArea.size(); i++)
        {
            stream <<  LogArea[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save ContinuumIndexes list, on 1 line
    if(Extrema.size()>0){
        stream <<  "#ContinuumIndexes Color for each extrema = {";
        for ( int i=0; i<ContinuumIndexes.size(); i++)
        {
            stream << "<";
            for(Int32 kci=0; kci<ContinuumIndexes[i].size(); kci++)
            {
                stream <<  ContinuumIndexes[i][kci].Color << "\t";
            }
            stream << ">";
        }
        stream << "}" << std::endl;
        stream <<  "#ContinuumIndexes Break for each extrema = {";
        for ( int i=0; i<ContinuumIndexes.size(); i++)
        {
            stream << "<";
            for(Int32 kci=0; kci<ContinuumIndexes[i].size(); kci++)
            {
                stream <<  ContinuumIndexes[i][kci].Break << "\t";
            }
            stream << ">";
        }
        stream << "}" << std::endl;
    }


    // save StrongELSNR list, on 1 line
    if(StrongELSNR.size()>0){
        stream <<  "#StrongELSNR for each extrema = {";
        for ( int i=0; i<StrongELSNR.size(); i++)
        {
            stream <<  StrongELSNR[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplName, on 1 line
    if(FittedTplName.size()>0){
        stream <<  "#FittedTplName for each extrema = {";
        for ( int i=0; i<FittedTplName.size(); i++)
        {
            stream <<  FittedTplName[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplAmplitude, on 1 line
    if(FittedTplAmplitude.size()>0){
        stream <<  "#FittedTplAmplitude for each extrema = {";
        for ( int i=0; i<FittedTplAmplitude.size(); i++)
        {
            stream <<  FittedTplAmplitude[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplMerit, on 1 line
    if(FittedTplMerit.size()>0){
        stream <<  "#FittedTplMerit for each extrema = {";
        for ( int i=0; i<FittedTplMerit.size(); i++)
        {
            stream <<  FittedTplMerit[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplDustCoeff, on 1 line
    if(FittedTplDustCoeff.size()>0){
        stream <<  "#FittedTplDustCoeff for each extrema = {";
        stream << std::setprecision(3);
        for ( int i=0; i<FittedTplDustCoeff.size(); i++)
        {
            stream <<  FittedTplDustCoeff[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplMeiksinIdx, on 1 line
    if(FittedTplMeiksinIdx.size()>0){
        stream <<  "#FittedTplMeiksinIdx for each extrema = {";
        for ( int i=0; i<FittedTplMeiksinIdx.size(); i++)
        {
            stream <<  FittedTplMeiksinIdx[i] << "\t";
        }
        stream << "}" << std::endl;
    }

    // save FittedTplshapeName, on 1 line
    if(FittedTplshapeName.size()>0){
        stream <<  "#FittedTplshapeName for each extrema = {";
        for ( int i=0; i<FittedTplshapeName.size(); i++)
        {
            stream <<  FittedTplshapeName[i] << "\t";
        }
        stream << "}" << std::endl;
    }

}
